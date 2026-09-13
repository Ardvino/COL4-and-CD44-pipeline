"""
MMP-8 & NAT Response - Survival Analysis

Replication of statistical methods from:
    Kesti et al. (2025). The prognostic significance of MMP-8 tissue
    immunoexpression in pancreatic ductal adenocarcinoma after neoadjuvant
    therapy. Scientific Reports.

Converted from MMP8_NAT_analysis.ipynb into a plain script. Two cohorts are
analyzed for the selected biomarker (high vs low expression):
    1. NAT cohort       - patients who received neoadjuvant therapy
                           (replicates the original paper's analyses)
    2. Upfront surgery  - patients who went straight to surgery, no NAT

All figures, CSV result tables, and a text log are written to a
per-biomarker output folder (results/<ACTIVE_BIOMARKER>/) so that runs for
different biomarkers never overwrite each other.

Endpoint: overall survival (OS) -- death from any cause, using the
pre-coded `OS` column in sorted_data.xlsx (1=died, 0=alive). This differs
from Kesti et al. (2025)'s own endpoint, disease-specific survival (DSS,
death from PDAC only) -- OS was substituted throughout by request.
"""

import re
import sys
from io import StringIO
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency, t as t_dist

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from lifelines.utils import median_survival_times

import warnings
warnings.filterwarnings('ignore')

plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.size'] = 11

# ============================================================
# Paths
# ============================================================
BASE_DIR = Path(__file__).resolve().parent
DATA_PATH = BASE_DIR / 'data' / 'sorted_data.xlsx'
SCORES_PATH = BASE_DIR / 'data' / 'patient_scores.xlsx'
RESULTS_DIR = BASE_DIR / 'results'

# Minimum number of events required before a Cox model is fit. Below this,
# HR estimates and CIs from lifelines are numerically unstable (or can fail
# to converge), so the model is skipped and flagged instead of reported.
MIN_EVENTS_FOR_COX = 10

# Valid range for each score column in patient_scores.xlsx, per biomarker.
# Four of them (plus Highest_Score, their cross-file max -- see
# score_patno.py) are raw IHC intensity scores (0 = non-detectable, 1 = low,
# 2 = moderate, 3 = high); a few readings in the source spreadsheet are
# data-entry sentinels outside that range (e.g. 7, 9) rather than real
# scores, so those are set to missing rather than treated as "very high"
# expression. CD44_SFF_percentage is scored on a different scale entirely --
# % positive tumor cells (0-100), not a 0-3 intensity grade. CD44_SFF_hscore
# and GATA6_hscore are classic H-scores (intensity x % positive cells,
# summed), 0-300.
#
# load_data() dichotomizes every column present in All_Scores_Summary, not
# just the active biomarker -- so a new score file in score_patNo/data/
# (which adds a column there automatically, see score_patno.py) needs an
# entry here (and in BIOMARKER_HIGH_CUTOFF below) before *any* biomarker can
# be analyzed, or load_data() raises a KeyError on the unregistered column.
BIOMARKER_VALID_RANGE = {
    'CD44_SFF_inflam_cells': (0, 3),
    'CD44_SFF_tumor_int':    (0, 3),
    'CD44_VFF_tumor_int':    (0, 3),
    'COL_4_stroma':          (0, 3),
    'CD44_SFF_percentage':   (0, 100),
    'CD44_SFF_hscore':       (0, 300),
    'GATA6_hscore':          (0, 300),
    'Highest_Score':         (0, 3),
}

# Fixed, pre-specified low/high cut-point per biomarker score, matching the
# convention used for MMP-8 in Kesti et al. (2025): scores 0-1 = low,
# 2-3 = high (i.e. "high" := score >= 2). This is NOT derived from this
# sample's distribution -- the four IHC scores here use the same 0-3 scale
# as MMP-8, so the same clinically-defined boundary is applied rather than a
# per-sample median split. See biomarker_cutpoint_recommendations.md for the
# analysis behind this.
# CD44_SFF_percentage's entry is only a placeholder for the 'fixed' method
# (>=10% positive, a common IHC-positivity convention) -- it has no prior
# clinical validation for this project the way score>=2 does for the IHC
# scores. Prefer dichotomization='p75' (75th-percentile split) for this
# biomarker, per Franklin et al.
# CD44_SFF_hscore and GATA6_hscore are completely novel scores with no prior
# published cut-point at all, so their entry (150, the scale's midpoint) is
# an even weaker placeholder than CD44_SFF_percentage's -- purely there so
# the 'fixed' method doesn't crash if someone selects it for these
# biomarkers. Prefer dichotomization='median' (this sample's own median
# H-score) for either instead, until/unless a clinically validated
# cut-point is established for it.
BIOMARKER_HIGH_CUTOFF = {
    'CD44_SFF_inflam_cells': 2,
    'CD44_SFF_tumor_int':    2,
    'CD44_VFF_tumor_int':    2,
    'COL_4_stroma':          2,
    'CD44_SFF_percentage':   10,
    'CD44_SFF_hscore':       150,
    'GATA6_hscore':          150,
    'Highest_Score':         2,
}

# Biomarkers with a natural continuous scale where Cox regression can also be
# fit on the raw (non-dichotomized) score, alongside the usual low/high
# binary. Franklin et al.'s own multivariable model entered CD44s
# continuously ("per 25% positive cells", their Table 2) rather than
# dichotomized, even though their KM plot used a 75th-percentile low/high
# split -- dichotomizing a continuous exposure loses statistical power, so
# this pipeline offers both for this biomarker. `scale` sets the HR's unit
# (e.g. 25 -> "per 25 percentage points").
BIOMARKER_CONTINUOUS = {
    'CD44_SFF_percentage': {'scale': 25, 'unit_label': 'per 25% positive'},
}

# ============================================================
# !! SELECT ACTIVE BIOMARKER — change this to switch markers
# ============================================================
# Options:
#   'MMP8'                  - MMP-8 (pre-coded binary in sorted_data.xlsx)
#   'CD44_SFF_inflam_cells' - CD44, SFF, inflammatory cells
#   'CD44_SFF_tumor_int'    - CD44, SFF, tumour intensity
#   'CD44_VFF_tumor_int'    - CD44, VFF, tumour intensity
#   'COL_4_stroma'          - COL-4, stroma
#   'CD44_SFF_percentage'   - CD44, SFF, % positive tumor cells (0-100 scale
#                             -- use dichotomization='p75', not 'fixed')
#   'CD44_SFF_hscore'       - CD44, SFF, H-score (0-300 scale, novel score
#                             -- use dichotomization='median', not 'fixed')
#   'GATA6_hscore'          - GATA6, H-score (0-300 scale, novel score
#                             -- use dichotomization='median', not 'fixed')
# ============================================================
ACTIVE_BIOMARKER = 'COL_4_stroma'  # change this to switch biomarkers

BIOMARKER_CONFIGS = {
    'MMP8': {
        'col':    'MMP8inCAvahvin_low_vs_high',
        'labels': {0: 'MMP-8 low', 1: 'MMP-8 high'},
        'name':   'MMP-8',
    },
    'CD44_SFF_inflam_cells': {
        'col':    'CD44_SFF_inflam_cells_binary',
        'labels': {0: 'CD44 SFF inflam low', 1: 'CD44 SFF inflam high'},
        'name':   'CD44 SFF (inflam cells)',
    },
    'CD44_SFF_tumor_int': {
        'col':    'CD44_SFF_tumor_int_binary',
        'labels': {0: 'CD44 SFF tumor low', 1: 'CD44 SFF tumor high'},
        'name':   'CD44 SFF (tumor)',
    },
    'CD44_VFF_tumor_int': {
        'col':    'CD44_VFF_tumor_int_binary',
        'labels': {0: 'CD44 VFF tumor low', 1: 'CD44 VFF tumor high'},
        'name':   'CD44 VFF (tumor)',
    },
    'COL_4_stroma': {
        'col':    'COL_4_stroma_binary',
        'labels': {0: 'COL-4 stroma low', 1: 'COL-4 stroma high'},
        'name':   'COL-4 (stroma)',
    },
    'CD44_SFF_percentage': {
        'col':    'CD44_SFF_percentage_binary',
        'labels': {0: 'CD44 SFF % positive low', 1: 'CD44 SFF % positive high'},
        'name':   'CD44 SFF (% positive tumor cells)',
    },
    'CD44_SFF_hscore': {
        'col':    'CD44_SFF_hscore_binary',
        'labels': {0: 'CD44 SFF H-score low', 1: 'CD44 SFF H-score high'},
        'name':   'CD44 SFF (H-score)',
    },
    'GATA6_hscore': {
        'col':    'GATA6_hscore_binary',
        'labels': {0: 'GATA6 H-score low', 1: 'GATA6 H-score high'},
        'name':   'GATA6 (H-score)',
    },
}

# Log file handle is set up per-cohort in main() and written to by `emit()`.
_LOG_FILE = None


def emit(text=''):
    """Print a line to the console and mirror it into the run's log file."""
    print(text)
    if _LOG_FILE is not None:
        _LOG_FILE.write(str(text) + '\n')


def emit_table(df, label=None):
    """Print a DataFrame as a readable table and mirror it into the log file."""
    if label:
        emit(f'\n{label}')
    emit(df.to_string(index=False))


# ============================================================
# Data loading
# ============================================================
# How to change the low/high split for the CD44/COL-4 scores: pick one of
# these three method keys and pass it as `dichotomization=` to load_data()
# or main() (or as the 2nd command-line argument -- see __main__ below).
# MMP-8 is unaffected either way -- it already has a pre-coded binary column
# in sorted_data.xlsx (see BIOMARKER_CONFIGS), so it never goes through this.
#   'fixed'  (default) -- BIOMARKER_HIGH_CUTOFF's fixed, pre-specified
#            cut-point per biomarker (score >= 2), matching the Kesti et al.
#            (2025) MMP-8 convention. Does NOT depend on this sample's
#            distribution -- see the module docstring for why this is the
#            default over a sample-derived split.
#   'median' -- splits at this sample's own median score for that biomarker
#            (high := score >= median).
#   'p75'    -- splits at this sample's own 75th-percentile score instead of
#            the median (high := score >= 75th percentile).
# 'median'/'p75' are sample-derived, so switching biomarkers or datasets can
# shift the cut-point itself, not just the resulting counts. The four 0-3
# IHC scores are coarse ordinal scores with heavy ties, so a meaningful
# fraction of patients can land exactly on a sample-derived cut. Despite its
# wider 0-100 range, CD44_SFF_percentage ties heavily too, since it's scored
# in coarse 5-point increments by eye rather than continuously -- in this
# sample, 26% of patients land exactly on the computed p75 value. To keep
# the resulting high-fraction close to what the quantile intends (e.g. ~25%
# for p75) despite that, load_data() picks whichever of ">="/">" the tied
# value should use rather than always defaulting to ">=" -- see the tied
# count logged below for how much it mattered for a given method/sample.
DICHOTOMIZATION_METHODS = {
    'fixed':  'Fixed clinical cut-point (score >= 2 per biomarker)',
    'median': 'Sample median split (high := score >= median)',
    'p75':    'Sample 75th-percentile split (high := score >= 75th percentile)',
}


def load_data(dichotomization='fixed'):
    """Load the raw clinical dataset and the biomarker score workbook, and
    dichotomize each CD44/COL-4 score into low/high using `dichotomization`
    (a key in DICHOTOMIZATION_METHODS -- see that dict for what each does).
    """
    if dichotomization not in DICHOTOMIZATION_METHODS:
        raise ValueError(f'Unknown dichotomization "{dichotomization}". '
                          f'Choose one of: {list(DICHOTOMIZATION_METHODS)}')

    df_raw = pd.read_excel(DATA_PATH)
    emit(f'Total rows: {len(df_raw)}')

    scores_summary = pd.read_excel(SCORES_PATH, sheet_name='All_Scores_Summary')
    emit(f'Biomarker columns: {[c for c in scores_summary.columns if c != "Patient_No"]}')

    # Align types before merging (Patient_No is string, PotNo is int)
    scores_summary['Patient_No'] = pd.to_numeric(scores_summary['Patient_No'], errors='coerce')

    # Merge biomarker scores into df_raw on patient number
    # sorted_data uses PotNo; patient_scores uses Patient_No
    df_raw = df_raw.merge(
        scores_summary,
        left_on='PotNo',
        right_on='Patient_No',
        how='left',
        suffixes=('', '_score'),
    )
    emit(f'Rows after merge: {len(df_raw)}')
    emit(f'Dichotomization method: {dichotomization} -- {DICHOTOMIZATION_METHODS[dichotomization]}')

    # The cut-point is applied to both NAT and upfront patients the same
    # way, since none of the three methods depends on either cohort's own
    # distribution (median/p75 are computed across all patients with a
    # valid score for that biomarker, not per cohort).
    biomarker_score_cols = [c for c in scores_summary.columns if c != 'Patient_No']
    for bm in biomarker_score_cols:
        # A handful of raw readings are data-entry sentinels outside the
        # valid scale for that biomarker (e.g. 7, 9 on a 0-3 IHC scale)
        # rather than real scores -- treat them as missing instead of
        # letting them masquerade as "very high".
        valid_range = BIOMARKER_VALID_RANGE[bm]
        bad = df_raw[bm].notna() & ~df_raw[bm].between(*valid_range)
        if bad.any():
            emit(f'{bm}: {bad.sum()} value(s) outside {valid_range}, '
                 f'set to NaN (PotNo: {df_raw.loc[bad, "PotNo"].tolist()})')
            df_raw.loc[bad, bm] = np.nan

        if dichotomization == 'fixed':
            # Clinically pre-specified cut-point: "high" always means
            # score >= cutoff, by definition (e.g. score >= 2).
            cutoff = BIOMARKER_HIGH_CUTOFF[bm]
            cmp_op = '>='
            high_mask = df_raw[bm] >= cutoff
        else:
            # Sample-derived cut-point. On coarse/rounded scores (the four
            # 0-3 IHC grades, or this dataset's 5-point-increment percentage
            # scores) the quantile value itself is often a heavily-repeated
            # value, so ">=" can sweep far more than the intended top
            # fraction into "high" (e.g. p75 should give ~25% high, but if
            # 26% of the sample ties at the computed cutoff, ">=" gives
            # ~48% high instead -- effectively a median split). Pick
            # whichever of ">=" / ">" lands the resulting high-fraction
            # closer to the quantile's intended fraction (1 - q).
            q = 0.5 if dichotomization == 'median' else 0.75
            cutoff = df_raw[bm].quantile(q)
            target_frac = 1 - q
            valid = df_raw.loc[df_raw[bm].notna(), bm]
            frac_ge = (valid >= cutoff).mean() if len(valid) else 0.0
            frac_gt = (valid > cutoff).mean() if len(valid) else 0.0
            # A comparison that would leave one group empty is unusable
            # (no KM/Cox model can run on a 0-patient group) regardless of
            # how numerically close its fraction happens to land to the
            # target -- e.g. on a right-skewed 0-3 scale, ">" at the cutoff
            # can leave 0 patients "high" even though 0.0 is arithmetically
            # "close" to a small target fraction. Only fall back to picking
            # by closeness when both options are actually usable.
            ge_usable = 0 < frac_ge < 1
            gt_usable = 0 < frac_gt < 1
            if ge_usable and gt_usable:
                use_ge = abs(frac_ge - target_frac) <= abs(frac_gt - target_frac)
            else:
                use_ge = ge_usable or not gt_usable
            if use_ge:
                cmp_op = '>='
                high_mask = df_raw[bm] >= cutoff
            else:
                cmp_op = '>'
                high_mask = df_raw[bm] > cutoff

        df_raw[f'{bm}_binary'] = high_mask.astype(float)
        df_raw.loc[df_raw[bm].isna(), f'{bm}_binary'] = np.nan

        if bm in BIOMARKER_CONTINUOUS:
            df_raw[f'{bm}_continuous'] = df_raw[bm] / BIOMARKER_CONTINUOUS[bm]['scale']

        n_low = (df_raw[f'{bm}_binary'] == 0).sum()
        n_high = (df_raw[f'{bm}_binary'] == 1).sum()
        tied_note = ''
        if dichotomization != 'fixed':
            n_tied = int((df_raw[bm] == cutoff).sum())
            tied_side = 'high' if cmp_op == '>=' else 'low'
            tied_note = f'  (tied at cutoff: {n_tied}, counted as {tied_side})'
        emit(f'{bm}: cutoff=score{cmp_op}{cutoff:g}  low={n_low}  high={n_high}{tied_note}')

    return df_raw


def build_column_map(biomarker_col, continuous_col=None):
    """
    Map paper concepts to dataframe column names.
    If your columns have different names, change the right-hand side here.

    `continuous_col`: the raw (non-dichotomized) column name for biomarkers
    in BIOMARKER_CONTINUOUS, or None for every other biomarker -- callers
    that want the continuous-exposure Cox tables check `col['mmp8_continuous']`.
    """
    return {
        'nat':       'NEOADJUVANTTI',              # 1=received neoadjuvant therapy, 0=upfront surgery
        'mmp8':      biomarker_col,                # set by biomarker selector above; 0=low, 1=high expression
        'mmp8_continuous': continuous_col,         # raw score, only set for biomarkers in BIOMARKER_CONTINUOUS
        'nat_resp':  'NATvaste_hyvä_012vs345',     # 1=strong(<=10% RTC), 0=weak(>=11% RTC)
        'os_event':  'os_event_binary',             # derived below: 1=died (any cause), 0=alive/censored
        'os_time':   'Survival_Months',             # follow-up time, used as the survival-analysis time axis
        'age':       'AGE_OPER',                   # age at surgery, in years (continuous)
        'sex':       'SUKUPUOLI',                  # 1=male (mies), 2=female (nainen)
        'stage':     'STAGE_8th_Binary',           # binarized AJCC 8th edition stage, 0 vs 1 (see original data prep for cut point)
        'grade':     'GRADUS',                     # histological grade, ordinal 1 (well) - 3 (poorly differentiated); raw column, used for chi-square only
        'grade_dummies': ['Grade_2', 'Grade_3'],   # dummy-coded grade for Cox models (reference = grade 1), derived in prepare_cohort_dataset()
        'logca199':  'logCA199',                   # log-transformed preoperative CA19-9 tumor marker level
        'regimen':   'GEMSITABINvsFOLFIRINOX',     # NAT regimen: 0=gemcitabine, 1=FOLFIRINOX (NAT patients only)
        'margin':    'Onkoleikkausradikaali_tunnus',  # resection margin text: 'R0'/'R1'/'muu' -- more complete than numeric RADIKALITEETTI
        'stage_full': 'STAGE_8th',                 # full AJCC 8th stage text (IA/IB/IIA/IIB/III/IV), vs the binarized 'stage' above
        'adjuvant':  'Postopsytostaatit_tunnus',   # 'Adjuvantti'/'Ei'/'Palliatiivinen' -- adjuvant vs none vs palliative treatment
    }


def prepare_cohort_dataset(df_raw, col, nat_flag):
    """
    Restrict to one treatment cohort (NAT or upfront surgery), derive the
    OS event flag, and drop rows missing biomarker/survival data.
    nat_flag: 1 to keep NAT patients, 0 to keep upfront-surgery patients.
    """
    cohort = df_raw[df_raw[col['nat']] == nat_flag].copy()
    emit(f'Patients in cohort: {len(cohort)}')

    # Drop rows missing biomarker/survival data *before* casting OS to int
    # below -- OS is a float column whenever any value is missing, and
    # .astype(int) on a NaN raises rather than letting dropna handle it.
    cohort = cohort.dropna(subset=[col['mmp8'], 'OS', col['os_time']])
    emit(f'Patients after dropping missing biomarker/survival: {len(cohort)}')

    # OS column is already coded as a binary event: 1=died (any cause), 0=alive
    cohort['os_event_binary'] = cohort['OS'].astype(int)
    emit(f'OS events (died, any cause): {cohort[col["os_event"]].sum()}')

    # Dummy-code histological grade (reference = grade 1) instead of entering
    # GRADUS directly as a continuous covariate in Cox models, which would
    # assume the 1->2 and 2->3 steps have equal, linear effects on hazard.
    grade_raw = cohort[col['grade']]
    cohort['Grade_2'] = np.where(grade_raw.isna(), np.nan, (grade_raw == 2).astype(float))
    cohort['Grade_3'] = np.where(grade_raw.isna(), np.nan, (grade_raw == 3).astype(float))

    return cohort


# ============================================================
# Chi-square: biomarker vs clinical variables
# ============================================================
def chi_square_table(df, variables, mmp8_col, output_dir, filename):
    """
    Chi-square test of each categorical variable against biomarker low/high.
    variables: list of (column_name, label) tuples.
    Returns and saves a results table as CSV.

    A significant result (p < 0.05) means the biomarker groups (low/high)
    are not evenly distributed across categories of that variable — i.e. the
    biomarker group is confounded with that clinical variable. A
    non-significant result supports comparing the biomarker groups directly,
    since the two groups are balanced with respect to that variable.
    """
    rows = []
    for c, lbl in variables:
        try:
            # Cross-tabulate counts, then test independence between the two variables
            ct = pd.crosstab(df[c], df[mmp8_col])
            chi2, p, dof, _ = chi2_contingency(ct)
            rows.append({'Variable': lbl, 'chi2': round(chi2, 3), 'df': dof, 'p': round(p, 3)})
        except Exception as e:
            rows.append({'Variable': lbl, 'chi2': np.nan, 'df': np.nan, 'p': f'ERROR: {e}'})

    result = pd.DataFrame(rows)
    emit_table(result, label='Biomarker vs clinical variables (Chi-Square)')
    result.to_csv(output_dir / filename, index=False)
    emit(f'Table saved as {output_dir / filename}')
    return result


# ============================================================
# Baseline characteristics ("Table 1")
# ============================================================
STAGE_GROUP = {'IA': 'I', 'IB': 'I', 'IIA': 'II', 'IIB': 'II', 'III': 'III', 'IV': 'IV'}


def _mean_ci(series, confidence=0.95):
    """Mean and two-sided confidence interval for a continuous series,
    dropping NaNs. Returns (mean, ci_lower, ci_upper); CI is NaN if fewer
    than 2 non-missing values (sample SD undefined for n<2).
    """
    vals = series.dropna()
    n = len(vals)
    if n < 2:
        return (vals.mean() if n else np.nan), np.nan, np.nan
    mean = vals.mean()
    sem = vals.std(ddof=1) / np.sqrt(n)
    t_crit = t_dist.ppf((1 + confidence) / 2, df=n - 1)
    return mean, mean - t_crit * sem, mean + t_crit * sem


def _add_categorical_rows(rows, series, characteristic, category_map=None, cohort_n=None):
    """Appends one {Characteristic, Category, n, %} row per distinct value in
    `series` (mapped through category_map if given -- unmapped values pass
    through unchanged), plus an explicit 'Missing' row if any value is NaN.
    % is always computed against `cohort_n` (the whole cohort, not just
    non-missing rows for this characteristic), so every characteristic's
    percentages describe the same population and a reader isn't left
    guessing what a given row's denominator was.
    """
    cohort_n = cohort_n if cohort_n is not None else len(series)
    mapped = series.map(category_map) if category_map else series
    counts = mapped.value_counts(dropna=True)
    for category, n in counts.items():
        # A float category that's a whole number (e.g. 2.0 from a grade
        # column upcast to float elsewhere due to unrelated NaNs) displays
        # as the bare integer rather than "2.0".
        if isinstance(category, float) and category.is_integer():
            category = int(category)
        rows.append({
            'Characteristic': characteristic, 'Category': str(category),
            'n': int(n), '%': round(100 * n / cohort_n, 0) if cohort_n else np.nan,
            'Value': '',
        })
    n_missing = int(series.isna().sum())
    if n_missing:
        rows.append({
            'Characteristic': characteristic, 'Category': 'Missing',
            'n': n_missing, '%': round(100 * n_missing / cohort_n, 0) if cohort_n else np.nan,
            'Value': '',
        })


def baseline_characteristics_table(df, col, output_dir, filename, is_nat=False):
    """
    Descriptive "Table 1"-style baseline characteristics for one cohort
    (already filtered via prepare_cohort_dataset()) -- sex, age, grade,
    resection margin, stage, adjuvant treatment, and follow-up time. This is
    a cohort-description table, not a comparison against the biomarker (see
    chi_square_table() for that) or a survival statistic (see
    group_summary_table() for KM-adjusted median OS).

    `is_nat`: when True, adds a NAT-regimen row (gemcitabine vs FOLFIRINOX),
    which only applies to the NAT cohort.

    Every categorical row's % is computed against the full cohort n (see
    _add_categorical_rows), with an explicit 'Missing' row rather than
    silently shrinking the denominator per characteristic.
    """
    n = len(df)
    rows = [{'Characteristic': 'N', 'Category': '', 'n': n, '%': '', 'Value': ''}]

    _add_categorical_rows(rows, df[col['sex']].map({1: 'Male', 2: 'Female'}), 'Sex', cohort_n=n)

    age_mean, age_lo, age_hi = _mean_ci(df[col['age']])
    age_value = f'{age_mean:.1f} (95% CI {age_lo:.1f}-{age_hi:.1f})' if pd.notna(age_lo) else f'{age_mean:.1f}'
    rows.append({'Characteristic': 'Age at surgery (years)', 'Category': '', 'n': '', '%': '', 'Value': age_value})

    _add_categorical_rows(rows, df[col['grade']], 'Grade', cohort_n=n)
    _add_categorical_rows(rows, df[col['margin']], 'Resection margin', cohort_n=n)
    _add_categorical_rows(rows, df[col['stage_full']].map(STAGE_GROUP), 'Stage', cohort_n=n)
    _add_categorical_rows(rows, df[col['adjuvant']].map({'Adjuvantti': 'Adjuvant', 'Ei': 'None', 'Palliatiivinen': 'Palliative'}),
                           'Adjuvant treatment', cohort_n=n)

    if is_nat:
        _add_categorical_rows(rows, df[col['regimen']].map({0: 'Gemcitabine', 1: 'FOLFIRINOX'}), 'NAT regimen', cohort_n=n)

    # Plain observed follow-up median/range -- deliberately NOT the
    # KM-adjusted median OS in group_summary_table(), which is a survival
    # statistic; this is the simple descriptive figure a paper's Table 1
    # reports for follow-up time.
    survival = df[col['os_time']].dropna()
    survival_value = f'{survival.median():.1f} (range {survival.min():.1f}-{survival.max():.1f})' if len(survival) else np.nan
    rows.append({'Characteristic': 'Survival, months', 'Category': '', 'n': '', '%': '', 'Value': survival_value})

    result = pd.DataFrame(rows)
    emit_table(result, label='Baseline characteristics')
    result.to_csv(output_dir / filename, index=False)
    emit(f'Table saved as {output_dir / filename}')
    return result


# ============================================================
# Kaplan-Meier plots (with number-at-risk table)
# ============================================================
def km_plot(df, time_col, event_col, group_col, group_labels,
            title='', ax=None, colors=('steelblue', 'tomato')):
    """
    Kaplan-Meier plot for two groups with log-rank and Wilcoxon-Breslow
    p-values, and a "number at risk" table beneath the curves.
    group_labels: dict {group_value: label_string}

    The step curves show the estimated probability of remaining event-free
    (overall survival) over time for each group; a curve that stays
    higher for longer indicates better survival for that group. The
    log-rank test compares the two full curves (not a single time point) and
    p < 0.05 means the two groups' survival experiences differ significantly.
    Log-rank weights every event equally over follow-up, so it has maximum
    power for a constant (proportional) hazard difference but can miss an
    effect concentrated early on; the Wilcoxon-Breslow test (a weighted
    log-rank variant, weighting each event by the number still at risk) is
    more sensitive to early differences and is reported alongside it for
    that reason -- the two can legitimately disagree.

    Requires exactly two groups to be present in `group_col` after dropping
    missing values (e.g. a subgroup where every remaining patient happens to
    share the same biomarker status has only one). If that's not the case,
    the plot is skipped (axis left blank with an explanatory message) rather
    than raising, and NaN p-values are returned instead.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))

    groups = sorted(df[group_col].dropna().unique())

    if len(groups) != 2:
        ax.text(0.5, 0.5, f'Skipped: found {len(groups)} group(s)\ninstead of 2 in "{group_col}"',
                ha='center', va='center', transform=ax.transAxes, fontsize=9)
        ax.set_title(title, fontsize=11)
        ax.axis('off')
        return {'logrank_p': np.nan, 'wilcoxon_p': np.nan}

    kmfs = []

    for grp, color in zip(groups, colors):
        mask = df[group_col] == grp
        kmf = KaplanMeierFitter()
        kmf.fit(
            df.loc[mask, time_col],
            event_observed=df.loc[mask, event_col],
            label=group_labels.get(grp, str(grp)),
        )
        kmf.plot_survival_function(ax=ax, ci_show=True, color=color)  # shaded band = 95% CI
        kmfs.append(kmf)

    # Number-at-risk table under the curves: how many patients are still being
    # followed (at risk), censored, or have had the event, at each time point.
    add_at_risk_counts(*kmfs, ax=ax)

    # Log-rank and Wilcoxon-Breslow tests between the two groups (each
    # compares the entire survival curves, just with different per-event
    # weighting -- see docstring above).
    g0, g1 = groups[0], groups[1]
    durations_A = df.loc[df[group_col] == g0, time_col]
    durations_B = df.loc[df[group_col] == g1, time_col]
    events_A = df.loc[df[group_col] == g0, event_col]
    events_B = df.loc[df[group_col] == g1, event_col]
    lr = logrank_test(durations_A, durations_B, event_observed_A=events_A, event_observed_B=events_B)
    wb = logrank_test(durations_A, durations_B, event_observed_A=events_A, event_observed_B=events_B,
                       weightings='wilcoxon')

    ax.set_title(f'{title}\nLog-rank p = {lr.p_value:.3f}   Wilcoxon-Breslow p = {wb.p_value:.3f}', fontsize=11)
    ax.set_xlabel('Time (months)')
    ax.set_ylabel('Overall survival')
    ax.set_ylim(0, 1.05)
    ax.legend(loc='upper right', fontsize=9)

    return {'logrank_p': lr.p_value, 'wilcoxon_p': wb.p_value}


# ============================================================
# Cox regression
# ============================================================
def univariable_cox_table(df, time_col, event_col, variables, output_dir, filename):
    """
    Run univariable Cox regression for each covariate (one model per row —
    each covariate is tested on its own, not adjusted for the others).
    variables: list of (covariate, label) tuples. `covariate`/`label` are
    normally a single column name / string, but can each be a list of the
    same length instead (e.g. dummy-coded grade) — in that case all columns
    in the list are fit together in one Cox model, and reported as separate
    rows (one contrast per dummy), since they represent one categorical
    variable rather than independent covariates.
    Returns and saves a results table as CSV.

    HR (hazard ratio) > 1 means higher values of the covariate (e.g. biomarker
    "high" coded as 1, vs "low" coded as 0) are associated with a higher risk
    of the event (worse survival); HR < 1 means lower risk (better survival).
    The 95% CI excluding 1, together with p < 0.05, indicates a statistically
    significant association.

    Models are skipped (and flagged in the table) if there are fewer than
    MIN_EVENTS_FOR_COX events, or if the fit fails to converge — both
    situations where the resulting HR/CI would be numerically unreliable.
    """
    rows = []
    for covariate, label in variables:
        cols = list(covariate) if isinstance(covariate, (list, tuple)) else [covariate]
        labels = list(label) if isinstance(label, (list, tuple)) else [label]

        # Fit one Cox model per covariate (or per dummy-coded group); dropna
        # keeps only rows with that covariate and outcome data present
        # (complete-case for this variable).
        subset = df[[time_col, event_col] + cols].dropna()
        n_events = subset[event_col].sum()

        if n_events < MIN_EVENTS_FOR_COX:
            for lbl in labels:
                rows.append({
                    'Variable': lbl, 'HR': np.nan, 'CI_lower': np.nan, 'CI_upper': np.nan,
                    'p': np.nan, 'sig': f'skipped: only {n_events} events (<{MIN_EVENTS_FOR_COX})',
                })
            continue

        try:
            cph = CoxPHFitter()
            cph.fit(subset, duration_col=time_col, event_col=event_col)
            s = cph.summary
        except Exception as e:
            for lbl in labels:
                rows.append({
                    'Variable': lbl, 'HR': np.nan, 'CI_lower': np.nan, 'CI_upper': np.nan,
                    'p': np.nan, 'sig': f'ERROR: {e}',
                })
            continue

        for c, lbl in zip(cols, labels):
            hr = np.exp(s.loc[c, 'coef'])
            ci_l = np.exp(s.loc[c, 'coef lower 95%'])
            ci_u = np.exp(s.loc[c, 'coef upper 95%'])
            p = s.loc[c, 'p']
            rows.append({
                'Variable': lbl,
                'HR': round(hr, 2),
                'CI_lower': round(ci_l, 2),
                'CI_upper': round(ci_u, 2),
                'p': round(p, 3),
                'sig': '*' if p < 0.05 else '',
            })

    result = pd.DataFrame(rows)
    emit_table(result, label='Univariable Cox regression')
    result.to_csv(output_dir / filename, index=False)
    emit(f'Table saved as {output_dir / filename}')
    return result


def multivariable_cox_table(df, time_col, event_col, covariates, output_dir, csv_filename, plot_filename, plot_title):
    """
    Run one multivariable Cox model on the given covariates (all covariates
    entered together, so each HR is adjusted for all the others in the list).
    Saves the coefficient table as CSV and a hazard-ratio forest plot as PNG.

    This uses complete-case analysis: only patients with no missing data in
    ANY of the listed covariates are included (see "Complete cases" count
    below), which is usually fewer than the full cohort size.

    The model is skipped if there are fewer than MIN_EVENTS_FOR_COX events
    or if the fit fails (e.g. convergence failure / near-perfect
    separation) — in both cases the HRs would not be reliable. If fit
    succeeds, the proportional-hazards assumption is checked via Schoenfeld
    residuals (`cph.check_assumptions`) and the result logged; this does not
    block the analysis, since a violation calls for judgement (e.g. a
    time-varying covariate) rather than an automatic pass/fail.
    """
    mv_data = df[[time_col, event_col] + covariates].dropna()
    n_events = mv_data[event_col].sum()
    emit(f'Complete cases for multivariable model: {len(mv_data)} (events: {n_events})')

    if n_events < MIN_EVENTS_FOR_COX:
        emit(f'WARNING: Only {n_events} events (<{MIN_EVENTS_FOR_COX}); '
             f'multivariable Cox model skipped as unreliable.')
        return None

    cph_mv = CoxPHFitter()
    try:
        cph_mv.fit(mv_data, duration_col=time_col, event_col=event_col)
    except Exception as e:
        emit(f'WARNING: Multivariable Cox model failed to fit ({e}); skipped.')
        return None

    # print_summary writes straight to stdout; capture it for the log file too
    buf = StringIO()
    old_stdout = sys.stdout
    sys.stdout = buf
    try:
        cph_mv.print_summary(decimals=3)
    finally:
        sys.stdout = old_stdout
    emit(buf.getvalue())

    # Proportional-hazards check (Schoenfeld residuals). Also captured from
    # stdout since check_assumptions prints its findings rather than
    # returning them.
    emit('--- Proportional hazards assumption check (Schoenfeld residuals) ---')
    buf_ph = StringIO()
    sys.stdout = buf_ph
    try:
        cph_mv.check_assumptions(mv_data, p_value_threshold=0.05, show_plots=False)
    except Exception as e:
        print(f'Could not test proportional hazards assumption: {e}')
    finally:
        sys.stdout = old_stdout
    emit(buf_ph.getvalue())

    cph_mv.summary.to_csv(output_dir / csv_filename)
    emit(f'Table saved as {output_dir / csv_filename}')

    fig, ax = plt.subplots(figsize=(8, 5))
    cph_mv.plot(hazard_ratios=True, ax=ax)
    ax.set_title(plot_title)
    ax.axvline(x=1, color='grey', linestyle='--', linewidth=0.8)
    plt.tight_layout()
    out_path = output_dir / plot_filename
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    emit(f'Figure saved as {out_path}')

    return cph_mv


# ============================================================
# Summary table
# ============================================================
def group_summary_table(groups, col, biomarker_name, output_dir, filename):
    """
    Patient counts, biomarker distribution, and Kaplan-Meier median OS per
    group.
    groups: dict {label: dataframe}

    Median OS is the Kaplan-Meier median survival time (the time at which
    the fitted survival curve crosses 0.5), not the median follow-up time
    among patients who died — that simpler statistic ignores censoring and
    is generally shorter/biased whenever a meaningful fraction of patients
    are censored before the event. 'inf' means the KM curve never drops to
    0.5 (median survival not reached within follow-up).
    """
    rows = []
    for label, df in groups.items():
        n = len(df)
        low = (df[col['mmp8']] == 0).sum()
        high = (df[col['mmp8']] == 1).sum()
        events = df[col['os_event']].sum()

        kmf = KaplanMeierFitter()
        kmf.fit(df[col['os_time']], event_observed=df[col['os_event']])
        med_surv = kmf.median_survival_time_
        med_ci = median_survival_times(kmf.confidence_interval_)
        med_lo, med_hi = med_ci.iloc[0, 0], med_ci.iloc[0, 1]

        rows.append({
            'Group': label,
            'n': n,
            f'{biomarker_name} low': low,
            f'{biomarker_name} low %': round(100 * low / n, 0) if n else np.nan,
            f'{biomarker_name} high': high,
            f'{biomarker_name} high %': round(100 * high / n, 0) if n else np.nan,
            'OS events': events,
            'Median OS, KM (mo)': round(med_surv, 1) if pd.notna(med_surv) else np.nan,
            'Median OS 95% CI lower': round(med_lo, 1) if pd.notna(med_lo) else np.nan,
            'Median OS 95% CI upper': round(med_hi, 1) if pd.notna(med_hi) else np.nan,
        })

    result = pd.DataFrame(rows)
    emit_table(result, label='Group summary')
    result.to_csv(output_dir / filename, index=False)
    emit(f'Table saved as {output_dir / filename}')
    return result


# ============================================================
# NAT cohort analysis (replicates the paper's Tables/Figures)
# ============================================================
def run_nat_cohort_analysis(df_raw, col, biomarker_name, active_biomarker, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)

    emit('\n' + '=' * 70)
    emit('NAT COHORT (patients who received neoadjuvant therapy)')
    emit('=' * 70)

    nat = prepare_cohort_dataset(df_raw, col, nat_flag=1)

    # --- Baseline characteristics ("Table 1") ---
    baseline_characteristics_table(nat, col, output_dir, filename='baseline_characteristics.csv', is_nat=True)

    # Subgroups by NAT response
    # NATvaste_hyvä_012vs345: 1=strong (scores 0-2, <=10% RTC), 0=weak (scores 3-5, >=11% RTC)
    strong = nat[nat[col['nat_resp']] == 1]
    weak = nat[nat[col['nat_resp']] == 0]
    emit(f'Strong responders (<=10% RTC): {len(strong)}')
    emit(f'Weak responders  (>=11% RTC):  {len(weak)}')

    # --- 3. Biomarker vs clinical variables (Chi-Square). Replicates Table 2. ---
    chi_square_table(
        nat,
        variables=[
            (col['sex'], 'Sex'),
            (col['stage'], 'Stage (binary)'),
            (col['grade'], 'Histological grade'),
            (col['regimen'], 'NAT regimen (gem vs FOLFIRINOX)'),
            (col['nat_resp'], 'NAT response (strong vs weak)'),
        ],
        mmp8_col=col['mmp8'],
        output_dir=output_dir,
        filename='chi_square_results.csv',
    )

    # --- 4. Kaplan-Meier survival analyses. Replicates Figure 2 (A, B, C). ---
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    fig.suptitle(f'{biomarker_name} and Overall Survival (OS) — NAT cohort', fontsize=13, y=1.02)

    labels = BIOMARKER_CONFIGS[active_biomarker]['labels']
    km_plot(nat, col['os_time'], col['os_event'], col['mmp8'], labels,
            title='(A) All NAT patients', ax=axes[0])
    km_plot(strong, col['os_time'], col['os_event'], col['mmp8'], labels,
            title='(B) Strong NAT response (<=10% RTC)', ax=axes[1])
    km_plot(weak, col['os_time'], col['os_event'], col['mmp8'], labels,
            title='(C) Weak NAT response (>=11% RTC)', ax=axes[2])

    plt.tight_layout()
    out_path = output_dir / f'KM_{active_biomarker}_OS_NAT.png'
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    emit(f'Figure saved as {out_path}')

    # --- 5. Univariable Cox regression. Replicates Table 3. ---
    emit('\n=== UNIVARIABLE COX - ALL NAT PATIENTS ===')
    univariable_cox_table(
        nat, col['os_time'], col['os_event'],
        variables=[
            (col['mmp8'], f'{biomarker_name} (low vs high)'),
            (col['nat_resp'], 'NAT response (strong vs weak)'),
            (col['age'], 'Age at surgery'),
            (col['sex'], 'Sex'),
            (col['stage'], 'Stage (binary)'),
            (col['grade_dummies'], ['Grade 2 vs 1', 'Grade 3 vs 1']),
            (col['logca199'], 'Log CA19-9 (preop)'),
        ],
        output_dir=output_dir, filename='cox_univariable_all_NAT.csv',
    )

    emit(f'\n=== UNIVARIABLE COX - STRONG RESPONDERS ONLY (n={len(strong)}) ===')
    univariable_cox_table(
        strong, col['os_time'], col['os_event'],
        variables=[(col['mmp8'], f'{biomarker_name} (low vs high)')],
        output_dir=output_dir, filename='cox_univariable_strong_responders.csv',
    )

    emit(f'\n=== UNIVARIABLE COX - WEAK RESPONDERS ONLY (n={len(weak)}) ===')
    univariable_cox_table(
        weak, col['os_time'], col['os_event'],
        variables=[(col['mmp8'], f'{biomarker_name} (low vs high)')],
        output_dir=output_dir, filename='cox_univariable_weak_responders.csv',
    )

    # --- 6. Multivariable Cox regression. Replicates Table 4 (full NAT group only). ---
    emit('\n=== MULTIVARIABLE COX - ALL NAT PATIENTS ===')
    multivariable_cox_table(
        nat, col['os_time'], col['os_event'],
        covariates=[col['mmp8'], col['age'], col['sex'], col['stage'], *col['grade_dummies'], col['logca199']],
        output_dir=output_dir,
        csv_filename='cox_multivariable_NAT.csv',
        plot_filename='ForestPlot_multivariable_NAT.png',
        plot_title='Multivariable Cox Regression - OS (all NAT patients)',
    )

    # --- 6b. Continuous-exposure Cox regression, for biomarkers in
    # BIOMARKER_CONTINUOUS only (currently just CD44_SFF_percentage).
    # Franklin et al.'s own multivariable model entered CD44s continuously
    # rather than dichotomized -- purely additive, every other biomarker
    # skips this block since col['mmp8_continuous'] is None for them.
    if col.get('mmp8_continuous'):
        unit_label = BIOMARKER_CONTINUOUS[active_biomarker]['unit_label']
        emit(f'\n=== UNIVARIABLE COX (CONTINUOUS, {unit_label.upper()}) - ALL NAT PATIENTS ===')
        univariable_cox_table(
            nat, col['os_time'], col['os_event'],
            variables=[(col['mmp8_continuous'], f'{biomarker_name} (continuous, {unit_label})')],
            output_dir=output_dir, filename='cox_univariable_continuous_NAT.csv',
        )
        emit(f'\n=== MULTIVARIABLE COX (CONTINUOUS, {unit_label.upper()}) - ALL NAT PATIENTS ===')
        multivariable_cox_table(
            nat, col['os_time'], col['os_event'],
            covariates=[col['mmp8_continuous'], col['age'], col['sex'], col['stage'], *col['grade_dummies'], col['logca199']],
            output_dir=output_dir,
            csv_filename='cox_multivariable_continuous_NAT.csv',
            plot_filename='ForestPlot_multivariable_continuous_NAT.png',
            plot_title=f'Multivariable Cox Regression, continuous exposure - OS (all NAT patients)',
        )

    # --- 7. Subgroup analysis by NAT regimen. Replicates Figure 4. ---
    emit('\n=== NAT RESPONSE AND OS BY REGIMEN ===')
    gem = nat[nat[col['regimen']] == 0]
    folf = nat[nat[col['regimen']] == 1]
    emit(f'Gemcitabine: {len(gem)}   FOLFIRINOX: {len(folf)}')

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle('NAT Response and OS by Regimen', fontsize=13)
    resp_labels = {1: 'Strong response (<=10%)', 0: 'Weak response (>=11%)'}
    km_plot(gem, col['os_time'], col['os_event'], col['nat_resp'], resp_labels,
            title='(A) Gemcitabine', ax=axes[0])
    km_plot(folf, col['os_time'], col['os_event'], col['nat_resp'], resp_labels,
            title='(B) FOLFIRINOX', ax=axes[1])
    plt.tight_layout()
    out_path = output_dir / 'KM_regimen_subgroup.png'
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    emit(f'Figure saved as {out_path}')

    # --- 8. Summary table ---
    group_summary_table(
        {
            'All NAT patients': nat,
            'Strong responders (<=10% RTC)': strong,
            'Weak responders (>=11% RTC)': weak,
        },
        col, biomarker_name, output_dir, filename='group_summary.csv',
    )


# ============================================================
# Upfront surgery cohort analysis (no NAT — biomarker high vs low only)
# ============================================================
def run_upfront_cohort_analysis(df_raw, col, biomarker_name, active_biomarker, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)

    emit('\n' + '=' * 70)
    emit('UPFRONT SURGERY COHORT (patients who did NOT receive neoadjuvant therapy)')
    emit('=' * 70)

    upfront = prepare_cohort_dataset(df_raw, col, nat_flag=0)

    # --- Baseline characteristics ("Table 1") ---
    baseline_characteristics_table(upfront, col, output_dir, filename='baseline_characteristics.csv', is_nat=False)

    # --- Biomarker vs clinical variables (Chi-Square) ---
    # NAT regimen / NAT response do not apply to this cohort, so they are omitted.
    chi_square_table(
        upfront,
        variables=[
            (col['sex'], 'Sex'),
            (col['stage'], 'Stage (binary)'),
            (col['grade'], 'Histological grade'),
        ],
        mmp8_col=col['mmp8'],
        output_dir=output_dir,
        filename='chi_square_results.csv',
    )

    # --- Kaplan-Meier survival analysis: biomarker low vs high ---
    fig, ax = plt.subplots(figsize=(7, 6))
    fig.suptitle(f'{biomarker_name} and Overall Survival (OS) — Upfront surgery', fontsize=12, y=1.02)
    labels = BIOMARKER_CONFIGS[active_biomarker]['labels']
    km_plot(upfront, col['os_time'], col['os_event'], col['mmp8'], labels,
            title='Upfront surgery patients', ax=ax)
    plt.tight_layout()
    out_path = output_dir / f'KM_{active_biomarker}_OS_upfront.png'
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    emit(f'Figure saved as {out_path}')

    # --- Univariable Cox regression ---
    emit('\n=== UNIVARIABLE COX - UPFRONT SURGERY PATIENTS ===')
    univariable_cox_table(
        upfront, col['os_time'], col['os_event'],
        variables=[
            (col['mmp8'], f'{biomarker_name} (low vs high)'),
            (col['age'], 'Age at surgery'),
            (col['sex'], 'Sex'),
            (col['stage'], 'Stage (binary)'),
            (col['grade_dummies'], ['Grade 2 vs 1', 'Grade 3 vs 1']),
            (col['logca199'], 'Log CA19-9 (preop)'),
        ],
        output_dir=output_dir, filename='cox_univariable_upfront.csv',
    )

    # --- Multivariable Cox regression ---
    emit('\n=== MULTIVARIABLE COX - UPFRONT SURGERY PATIENTS ===')
    multivariable_cox_table(
        upfront, col['os_time'], col['os_event'],
        covariates=[col['mmp8'], col['age'], col['sex'], col['stage'], *col['grade_dummies'], col['logca199']],
        output_dir=output_dir,
        csv_filename='cox_multivariable_upfront.csv',
        plot_filename='ForestPlot_multivariable_upfront.png',
        plot_title='Multivariable Cox Regression - OS (upfront surgery patients)',
    )

    # --- Continuous-exposure Cox regression (BIOMARKER_CONTINUOUS only) ---
    if col.get('mmp8_continuous'):
        unit_label = BIOMARKER_CONTINUOUS[active_biomarker]['unit_label']
        emit(f'\n=== UNIVARIABLE COX (CONTINUOUS, {unit_label.upper()}) - UPFRONT SURGERY PATIENTS ===')
        univariable_cox_table(
            upfront, col['os_time'], col['os_event'],
            variables=[(col['mmp8_continuous'], f'{biomarker_name} (continuous, {unit_label})')],
            output_dir=output_dir, filename='cox_univariable_continuous_upfront.csv',
        )
        emit(f'\n=== MULTIVARIABLE COX (CONTINUOUS, {unit_label.upper()}) - UPFRONT SURGERY PATIENTS ===')
        multivariable_cox_table(
            upfront, col['os_time'], col['os_event'],
            covariates=[col['mmp8_continuous'], col['age'], col['sex'], col['stage'], *col['grade_dummies'], col['logca199']],
            output_dir=output_dir,
            csv_filename='cox_multivariable_continuous_upfront.csv',
            plot_filename='ForestPlot_multivariable_continuous_upfront.png',
            plot_title='Multivariable Cox Regression, continuous exposure - OS (upfront surgery patients)',
        )

    # --- Summary table ---
    group_summary_table(
        {'Upfront surgery patients': upfront},
        col, biomarker_name, output_dir, filename='group_summary.csv',
    )


# ============================================================
# Patient-level data export
# ============================================================
def export_patient_data(df_raw, col, output_dir):
    """
    Build and save patient_data.csv covering all patients (NAT and upfront
    surgery) in the standard template format.

    Columns that don't apply to upfront-surgery patients (treatment_response,
    chemotherapy_regimen) are left empty for those rows.
    """
    export = df_raw.copy()

    # Keep OS event as NaN when OS is missing rather than defaulting to 0
    export['os_event_binary'] = np.where(
        export['OS'].isna(), np.nan, export['OS'].astype(float)
    )

    export['_biomarker_group']      = export[col['mmp8']].map({0: 'Low', 1: 'High'})
    export['_treatment_group']      = export[col['nat']].map({1: 'Neoadjuvant', 0: 'Upfront surgery'})
    export['_treatment_response']   = export[col['nat_resp']].map({1: 'Strong', 0: 'Weak'})
    export['_chemotherapy_regimen'] = export[col['regimen']].map({0: 'Gemcitabine', 1: 'FOLFIRINOX'})
    export['_sex']                  = export[col['sex']].map({1: 'Male', 2: 'Female'})
    export['_disease_stage']        = export[col['stage']].map({0: 'Early', 1: 'Advanced'})

    is_upfront = export[col['nat']] != 1
    export.loc[is_upfront, '_treatment_response']   = ''
    export.loc[is_upfront, '_chemotherapy_regimen'] = ''

    # Drop rows where any required field is missing
    missing_mask = (
        export[col['os_time']].isna()
        | export['os_event_binary'].isna()
        | export[col['mmp8']].isna()
    )
    if missing_mask.any():
        emit(f'export_patient_data: dropping {missing_mask.sum()} row(s) with '
             f'missing follow_up_months, event_observed, or biomarker score '
             f'(PotNo: {export.loc[missing_mask, "PotNo"].tolist()})')
        export = export[~missing_mask]

    result = export[[
        'PotNo',
        col['os_time'],
        'os_event_binary',
        '_biomarker_group',
        '_treatment_group',
        '_treatment_response',
        '_chemotherapy_regimen',
        col['age'],
        '_sex',
        '_disease_stage',
        col['grade'],
        col['logca199'],
    ]].rename(columns={
        'PotNo':                 'patient_id',
        col['os_time']:          'follow_up_months',
        'os_event_binary':       'event_observed',
        '_biomarker_group':      'biomarker_group',
        '_treatment_group':      'treatment_group',
        '_treatment_response':   'treatment_response',
        '_chemotherapy_regimen': 'chemotherapy_regimen',
        col['age']:              'patient_age',
        '_sex':                  'sex',
        '_disease_stage':        'disease_stage',
        col['grade']:            'histological_grade',
        col['logca199']:         'log_ca19_9',
    })

    out_path = output_dir / 'patient_data.csv'
    result.to_csv(out_path, index=False)
    emit(f'Patient data exported: {out_path}  ({len(result)} rows)')
    return result


# ============================================================
# Themed SVG plots (interactive display / download)
# ============================================================
# Additive only -- nothing above this section is touched, and nothing here
# is called by main() or the CLI, so `python statistics.py` keeps producing
# exactly the same PNGs it always has. These are for a UI (e.g. the Streamlit
# pipeline app) that wants a dark-themed plot on screen and a light-themed
# SVG for download, using the same validated categorical palette (blue/orange,
# slots 1-2) and chart-chrome tokens as the rest of the project's UI work.
PLOT_THEMES = {
    'dark': {
        'text':    '#ffffff',  # all text -- titles, axis/tick labels, legend, at-risk table
        'muted':   '#ffffff',
        'grid':    '#2c2c2a',
        'axis':    '#383835',  # spines
        'low':     '#d95926',  # categorical slot 2 (orange)
        'high':    '#3987e5',  # categorical slot 1 (blue)
    },
    'light': {
        'text':    '#0b0b0b',
        'muted':   '#898781',
        'grid':    '#e1e0d9',
        'axis':    '#c3c2b7',
        'low':     '#eb6834',
        'high':    '#2a78d6',
    },
}


def _style_axes(fig, theme):
    """Recolor every axes already drawn on `fig` to match `theme` (a
    PLOT_THEMES entry) -- title/tick/spine/grid/legend/text colors. Figure
    and axes backgrounds are left transparent (no fill), so the plot blends
    into whatever page/card it's placed on rather than drawing its own
    background rectangle. Applied after a figure is built, so it works on
    any figure produced by km_plot()/CoxPHFitter.plot() without changing how
    those functions draw the data itself.
    """
    fig.patch.set_alpha(0)
    for ax in fig.axes:
        ax.patch.set_alpha(0)
        if ax.title.get_text():
            ax.title.set_color(theme['text'])
        ax.xaxis.label.set_color(theme['muted'])
        ax.yaxis.label.set_color(theme['muted'])
        ax.tick_params(colors=theme['muted'])
        for spine in ax.spines.values():
            spine.set_color(theme['axis'])
        ax.grid(color=theme['grid'], linewidth=0.6, alpha=0.8)
        for text in ax.texts:  # includes the KM at-risk-count annotations
            text.set_color(theme['text'])
        legend = ax.get_legend()
        if legend is not None:
            legend.get_frame().set_alpha(0)
            legend.get_frame().set_edgecolor(theme['axis'])
            for t in legend.get_texts():
                t.set_color(theme['text'])


def fig_to_svg(fig):
    """Serializes `fig` to standalone, background-transparent SVG markup
    that scales with its container.
    """
    buf = StringIO()
    fig.savefig(buf, format='svg', bbox_inches='tight', transparent=True)
    svg = buf.getvalue()
    return re.sub(r'<svg ', '<svg style="width:100%;height:auto;display:block;" ', svg, count=1)


def themed_km_figure(panels, time_col, event_col, group_col, group_labels, theme, suptitle=None, figsize=None):
    """Builds a 1..N-panel Kaplan-Meier figure styled for `theme`, reusing
    km_plot() unmodified for the actual survival-curve fitting/log-rank test
    in each panel -- only the color scheme and chrome differ from the plain
    PNG versions run_nat_cohort_analysis()/run_upfront_cohort_analysis()
    save to disk.
    panels: list of (dataframe, panel_title) pairs, one per subplot, sharing
    the same time/event/group columns and group_labels.
    """
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=figsize or (6 * n, 5.2))
    axes = [axes] if n == 1 else list(axes)
    colors = (theme['low'], theme['high'])
    for ax, (df, title) in zip(axes, panels):
        km_plot(df, time_col, event_col, group_col, group_labels, title=title, ax=ax, colors=colors)
    if suptitle:
        fig.suptitle(suptitle, fontsize=13, y=1.02, color=theme['text'])
    plt.tight_layout()
    _style_axes(fig, theme)
    return fig


def themed_forest_figure(mv_data, time_col, event_col, covariates, title, theme):
    """Fits a fresh multivariable CoxPHFitter -- the same call
    multivariable_cox_table() makes -- and returns a themed hazard-ratio
    forest plot Figure, or None if there are too few events or the fit
    fails to converge (mirroring multivariable_cox_table()'s own checks).
    """
    subset = mv_data[[time_col, event_col] + covariates].dropna()
    if subset[event_col].sum() < MIN_EVENTS_FOR_COX:
        return None
    cph = CoxPHFitter()
    try:
        cph.fit(subset, duration_col=time_col, event_col=event_col)
    except Exception:
        return None
    fig, ax = plt.subplots(figsize=(8, max(3, 0.6 * len(covariates) + 1.5)))
    # lifelines' cph.plot() defaults its errorbar color to c="k" internally
    # (via setdefault) and draws its own HR=1 reference line in black -- pass
    # `c` (not `color`, which collides with that default) for the markers,
    # then recolor the reference line (the one LineCollection it adds to
    # ax.collections) to match the theme.
    cph.plot(hazard_ratios=True, ax=ax, c=theme['high'], markerfacecolor=theme['high'])
    ax.set_title(title, color=theme['text'])
    for coll in ax.collections:
        coll.set_color(theme['muted'])
    plt.tight_layout()
    _style_axes(fig, theme)
    return fig


# ============================================================
# Main
# ============================================================
def results_dir_for(active_biomarker, dichotomization='fixed'):
    """The output folder main() writes/reads for a given biomarker +
    dichotomization method. The 'fixed' (default) method keeps today's
    exact path (results/<biomarker>/) for CLI/back-compat; the two
    sample-derived methods get their own sub-folder so switching methods
    never overwrites another method's results for the same biomarker.
    """
    base = RESULTS_DIR / active_biomarker
    return base if dichotomization == 'fixed' else base / dichotomization


def main(active_biomarker=ACTIVE_BIOMARKER, dichotomization='fixed'):
    global _LOG_FILE

    if active_biomarker not in BIOMARKER_CONFIGS:
        raise ValueError(f'Unknown biomarker "{active_biomarker}". '
                          f'Choose one of: {list(BIOMARKER_CONFIGS)}')

    bm = BIOMARKER_CONFIGS[active_biomarker]
    biomarker_col = bm['col']
    biomarker_name = bm['name']
    continuous_col = f'{active_biomarker}_continuous' if active_biomarker in BIOMARKER_CONTINUOUS else None
    col = build_column_map(biomarker_col, continuous_col=continuous_col)

    # Output folder is named after the active biomarker (and, for a
    # sample-derived dichotomization, the method too), with one sub-folder
    # per cohort, so different runs/cohorts/methods never clobber each
    # other's figures/tables/logs.
    base_output_dir = results_dir_for(active_biomarker, dichotomization)
    nat_dir = base_output_dir / 'NAT_cohort'
    upfront_dir = base_output_dir / 'Upfront_surgery_cohort'
    nat_dir.mkdir(parents=True, exist_ok=True)
    upfront_dir.mkdir(parents=True, exist_ok=True)

    # Data is loaded once and shared between both cohort analyses.
    with open(base_output_dir / 'analysis_log.txt', 'w', encoding='utf-8') as log_file:
        _LOG_FILE = log_file

        emit(f'Active biomarker: {biomarker_name}  ->  column: {biomarker_col}')
        emit(f'Output folder: {base_output_dir}')

        df_raw = load_data(dichotomization=dichotomization)

        run_nat_cohort_analysis(df_raw, col, biomarker_name, active_biomarker, nat_dir)
        run_upfront_cohort_analysis(df_raw, col, biomarker_name, active_biomarker, upfront_dir)
        export_patient_data(df_raw, col, base_output_dir)

        emit()
        emit(f'Done. NAT cohort results:      {nat_dir}')
        emit(f'Done. Upfront surgery results: {upfront_dir}')
        emit(f'Done. Patient data export:     {base_output_dir / "patient_data.csv"}')

    _LOG_FILE = None


if __name__ == '__main__':
    # Optional: pass a biomarker key, and/or a dichotomization method key
    # (see DICHOTOMIZATION_METHODS), on the command line to override the
    # defaults, e.g.  python statistics.py MMP8            (biomarker only)
    #                 python statistics.py COL_4_stroma median
    chosen = sys.argv[1] if len(sys.argv) > 1 else ACTIVE_BIOMARKER
    method = sys.argv[2] if len(sys.argv) > 2 else 'fixed'
    main(chosen, dichotomization=method)
