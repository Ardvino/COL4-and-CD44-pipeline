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
"""

import sys
from io import StringIO
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency

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

# Valid range for the raw IHC intensity scores in patient_scores.xlsx
# (0 = non-detectable, 1 = low, 2 = moderate, 3 = high). A few readings in
# the source spreadsheet are data-entry sentinels outside this range (e.g.
# 7, 9) rather than real scores -- these are set to missing rather than
# treated as "very high" expression.
VALID_SCORE_RANGE = (0, 3)

# Fixed, pre-specified low/high cut-point per biomarker score, matching the
# convention used for MMP-8 in Kesti et al. (2025): scores 0-1 = low,
# 2-3 = high (i.e. "high" := score >= 2). This is NOT derived from this
# sample's distribution -- all four scores here use the same 0-3 IHC
# intensity scale as MMP-8, so the same clinically-defined boundary is
# applied rather than a per-sample median split. See
# biomarker_cutpoint_recommendations.md for the analysis behind this.
BIOMARKER_HIGH_CUTOFF = {
    'CD44_SFF_inflam_cells': 2,
    'CD44_SFF_tumor_int':    2,
    'CD44_VFF_tumor_int':    2,
    'COL_4_stroma':          2,
    'Highest_Score':         2,
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
# ============================================================
ACTIVE_BIOMARKER = 'CD44_VFF_tumor_int'  # change this to switch biomarkers

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
def load_data():
    """Load the raw clinical dataset and the biomarker score workbook."""
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

    # Dichotomize each score using a fixed, pre-specified cut-point (see
    # BIOMARKER_HIGH_CUTOFF) rather than a sample-derived median -- these
    # scores are coarse 0-3 IHC intensity categories with heavy ties at the
    # median, so a median split produces an arbitrary and unstable low/high
    # boundary. The same cut-point is applied to both NAT and upfront
    # patients, since it doesn't depend on either cohort's distribution.
    biomarker_score_cols = [c for c in scores_summary.columns if c != 'Patient_No']
    for bm in biomarker_score_cols:
        # A handful of raw readings are data-entry sentinels outside the
        # valid 0-3 scale (e.g. 7, 9) rather than real scores -- treat them
        # as missing instead of letting them masquerade as "very high".
        bad = df_raw[bm].notna() & ~df_raw[bm].between(*VALID_SCORE_RANGE)
        if bad.any():
            emit(f'{bm}: {bad.sum()} value(s) outside {VALID_SCORE_RANGE}, '
                 f'set to NaN (PotNo: {df_raw.loc[bad, "PotNo"].tolist()})')
            df_raw.loc[bad, bm] = np.nan

        cutoff = BIOMARKER_HIGH_CUTOFF[bm]
        df_raw[f'{bm}_binary'] = (df_raw[bm] >= cutoff).astype(float)
        df_raw.loc[df_raw[bm].isna(), f'{bm}_binary'] = np.nan
        n_low = (df_raw[f'{bm}_binary'] == 0).sum()
        n_high = (df_raw[f'{bm}_binary'] == 1).sum()
        emit(f'{bm}: cutoff=score>={cutoff}  low={n_low}  high={n_high}')

    return df_raw


def build_column_map(biomarker_col):
    """
    Map paper concepts to dataframe column names.
    If your columns have different names, change the right-hand side here.
    """
    return {
        'nat':       'NEOADJUVANTTI',              # 1=received neoadjuvant therapy, 0=upfront surgery
        'mmp8':      biomarker_col,                # set by biomarker selector above; 0=low, 1=high expression
        'nat_resp':  'NATvaste_hyvä_012vs345',     # 1=strong(<=10% RTC), 0=weak(>=11% RTC)
        'dss_event': 'dss_event_binary',           # derived below: 1=died of PDAC, 0=censored/other
        'dss_time':  'Survival_Months',            # follow-up time, used as the survival-analysis time axis
        'age':       'AGE_OPER',                   # age at surgery, in years (continuous)
        'sex':       'SUKUPUOLI',                  # 1=male (mies), 2=female (nainen)
        'stage':     'STAGE_8th_Binary',           # binarized AJCC 8th edition stage, 0 vs 1 (see original data prep for cut point)
        'grade':     'GRADUS',                     # histological grade, ordinal 1 (well) - 3 (poorly differentiated); raw column, used for chi-square only
        'grade_dummies': ['Grade_2', 'Grade_3'],   # dummy-coded grade for Cox models (reference = grade 1), derived in prepare_cohort_dataset()
        'logca199':  'logCA199',                   # log-transformed preoperative CA19-9 tumor marker level
        'regimen':   'GEMSITABINvsFOLFIRINOX',     # NAT regimen: 0=gemcitabine, 1=FOLFIRINOX (NAT patients only)
    }


def prepare_cohort_dataset(df_raw, col, nat_flag):
    """
    Restrict to one treatment cohort (NAT or upfront surgery), derive the
    DSS event flag, and drop rows missing biomarker/survival data.
    nat_flag: 1 to keep NAT patients, 0 to keep upfront-surgery patients.
    """
    cohort = df_raw[df_raw[col['nat']] == nat_flag].copy()
    emit(f'Patients in cohort: {len(cohort)}')

    # DSS column is coded 1=died of PDAC, 2=alive, 3=died of other cause -> convert to binary
    cohort['dss_event_binary'] = (cohort['DSS'] == 1).astype(int)

    # Dummy-code histological grade (reference = grade 1) instead of entering
    # GRADUS directly as a continuous covariate in Cox models, which would
    # assume the 1->2 and 2->3 steps have equal, linear effects on hazard.
    grade_raw = cohort[col['grade']]
    cohort['Grade_2'] = np.where(grade_raw.isna(), np.nan, (grade_raw == 2).astype(float))
    cohort['Grade_3'] = np.where(grade_raw.isna(), np.nan, (grade_raw == 3).astype(float))

    cohort = cohort.dropna(subset=[col['mmp8'], 'DSS', col['dss_time']])
    emit(f'Patients after dropping missing biomarker/survival: {len(cohort)}')
    emit(f'DSS events (died of PDAC): {cohort[col["dss_event"]].sum()}')

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
# Kaplan-Meier plots (with number-at-risk table)
# ============================================================
def km_plot(df, time_col, event_col, group_col, group_labels,
            title='', ax=None, colors=('steelblue', 'tomato')):
    """
    Kaplan-Meier plot for two groups with log-rank p-value and a
    "number at risk" table beneath the curves.
    group_labels: dict {group_value: label_string}

    The step curves show the estimated probability of remaining event-free
    (disease-specific survival) over time for each group; a curve that stays
    higher for longer indicates better survival for that group. The
    log-rank test compares the two full curves (not a single time point) and
    p < 0.05 means the two groups' survival experiences differ significantly.

    Requires exactly two groups to be present in `group_col` after dropping
    missing values (e.g. a subgroup where every remaining patient happens to
    share the same biomarker status has only one). If that's not the case,
    the plot is skipped (axis left blank with an explanatory message) rather
    than raising, and NaN is returned instead of a log-rank p-value.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))

    groups = sorted(df[group_col].dropna().unique())

    if len(groups) != 2:
        ax.text(0.5, 0.5, f'Skipped: found {len(groups)} group(s)\ninstead of 2 in "{group_col}"',
                ha='center', va='center', transform=ax.transAxes, fontsize=9)
        ax.set_title(title, fontsize=11)
        ax.axis('off')
        return np.nan

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

    # Log-rank test between the two groups (compares entire survival curves)
    g0, g1 = groups[0], groups[1]
    lr = logrank_test(
        df.loc[df[group_col] == g0, time_col], df.loc[df[group_col] == g1, time_col],
        event_observed_A=df.loc[df[group_col] == g0, event_col],
        event_observed_B=df.loc[df[group_col] == g1, event_col],
    )

    ax.set_title(f'{title}\nLog-rank p = {lr.p_value:.3f}', fontsize=11)
    ax.set_xlabel('Time (months)')
    ax.set_ylabel('Disease-specific survival')
    ax.set_ylim(0, 1.05)
    ax.legend(loc='upper right', fontsize=9)

    return lr.p_value


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
    Patient counts, biomarker distribution, and Kaplan-Meier median DSS per
    group.
    groups: dict {label: dataframe}

    Median DSS is the Kaplan-Meier median survival time (the time at which
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
        events = df[col['dss_event']].sum()

        kmf = KaplanMeierFitter()
        kmf.fit(df[col['dss_time']], event_observed=df[col['dss_event']])
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
            'DSS events': events,
            'Median DSS, KM (mo)': round(med_surv, 1) if pd.notna(med_surv) else np.nan,
            'Median DSS 95% CI lower': round(med_lo, 1) if pd.notna(med_lo) else np.nan,
            'Median DSS 95% CI upper': round(med_hi, 1) if pd.notna(med_hi) else np.nan,
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
    fig.suptitle(f'{biomarker_name} and Disease-Specific Survival (DSS) — NAT cohort', fontsize=13, y=1.02)

    labels = BIOMARKER_CONFIGS[active_biomarker]['labels']
    km_plot(nat, col['dss_time'], col['dss_event'], col['mmp8'], labels,
            title='(A) All NAT patients', ax=axes[0])
    km_plot(strong, col['dss_time'], col['dss_event'], col['mmp8'], labels,
            title='(B) Strong NAT response (<=10% RTC)', ax=axes[1])
    km_plot(weak, col['dss_time'], col['dss_event'], col['mmp8'], labels,
            title='(C) Weak NAT response (>=11% RTC)', ax=axes[2])

    plt.tight_layout()
    out_path = output_dir / f'KM_{active_biomarker}_DSS_NAT.png'
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    emit(f'Figure saved as {out_path}')

    # --- 5. Univariable Cox regression. Replicates Table 3. ---
    emit('\n=== UNIVARIABLE COX - ALL NAT PATIENTS ===')
    univariable_cox_table(
        nat, col['dss_time'], col['dss_event'],
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
        strong, col['dss_time'], col['dss_event'],
        variables=[(col['mmp8'], f'{biomarker_name} (low vs high)')],
        output_dir=output_dir, filename='cox_univariable_strong_responders.csv',
    )

    emit(f'\n=== UNIVARIABLE COX - WEAK RESPONDERS ONLY (n={len(weak)}) ===')
    univariable_cox_table(
        weak, col['dss_time'], col['dss_event'],
        variables=[(col['mmp8'], f'{biomarker_name} (low vs high)')],
        output_dir=output_dir, filename='cox_univariable_weak_responders.csv',
    )

    # --- 6. Multivariable Cox regression. Replicates Table 4 (full NAT group only). ---
    emit('\n=== MULTIVARIABLE COX - ALL NAT PATIENTS ===')
    multivariable_cox_table(
        nat, col['dss_time'], col['dss_event'],
        covariates=[col['mmp8'], col['age'], col['sex'], col['stage'], *col['grade_dummies'], col['logca199']],
        output_dir=output_dir,
        csv_filename='cox_multivariable_NAT.csv',
        plot_filename='ForestPlot_multivariable_NAT.png',
        plot_title='Multivariable Cox Regression - DSS (all NAT patients)',
    )

    # --- 7. Subgroup analysis by NAT regimen. Replicates Figure 4. ---
    emit('\n=== NAT RESPONSE AND DSS BY REGIMEN ===')
    gem = nat[nat[col['regimen']] == 0]
    folf = nat[nat[col['regimen']] == 1]
    emit(f'Gemcitabine: {len(gem)}   FOLFIRINOX: {len(folf)}')

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle('NAT Response and DSS by Regimen', fontsize=13)
    resp_labels = {1: 'Strong response (<=10%)', 0: 'Weak response (>=11%)'}
    km_plot(gem, col['dss_time'], col['dss_event'], col['nat_resp'], resp_labels,
            title='(A) Gemcitabine', ax=axes[0])
    km_plot(folf, col['dss_time'], col['dss_event'], col['nat_resp'], resp_labels,
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
    fig.suptitle(f'{biomarker_name} and Disease-Specific Survival (DSS) — Upfront surgery', fontsize=12, y=1.02)
    labels = BIOMARKER_CONFIGS[active_biomarker]['labels']
    km_plot(upfront, col['dss_time'], col['dss_event'], col['mmp8'], labels,
            title='Upfront surgery patients', ax=ax)
    plt.tight_layout()
    out_path = output_dir / f'KM_{active_biomarker}_DSS_upfront.png'
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    emit(f'Figure saved as {out_path}')

    # --- Univariable Cox regression ---
    emit('\n=== UNIVARIABLE COX - UPFRONT SURGERY PATIENTS ===')
    univariable_cox_table(
        upfront, col['dss_time'], col['dss_event'],
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
        upfront, col['dss_time'], col['dss_event'],
        covariates=[col['mmp8'], col['age'], col['sex'], col['stage'], *col['grade_dummies'], col['logca199']],
        output_dir=output_dir,
        csv_filename='cox_multivariable_upfront.csv',
        plot_filename='ForestPlot_multivariable_upfront.png',
        plot_title='Multivariable Cox Regression - DSS (upfront surgery patients)',
    )

    # --- Summary table ---
    group_summary_table(
        {'Upfront surgery patients': upfront},
        col, biomarker_name, output_dir, filename='group_summary.csv',
    )


# ============================================================
# Main
# ============================================================
def main(active_biomarker=ACTIVE_BIOMARKER):
    global _LOG_FILE

    if active_biomarker not in BIOMARKER_CONFIGS:
        raise ValueError(f'Unknown biomarker "{active_biomarker}". '
                          f'Choose one of: {list(BIOMARKER_CONFIGS)}')

    bm = BIOMARKER_CONFIGS[active_biomarker]
    biomarker_col = bm['col']
    biomarker_name = bm['name']
    col = build_column_map(biomarker_col)

    # Output folder is named after the active biomarker, with one
    # sub-folder per cohort, so different runs/cohorts never clobber
    # each other's figures/tables/logs.
    base_output_dir = RESULTS_DIR / active_biomarker
    nat_dir = base_output_dir / 'NAT_cohort'
    upfront_dir = base_output_dir / 'Upfront_surgery_cohort'
    nat_dir.mkdir(parents=True, exist_ok=True)
    upfront_dir.mkdir(parents=True, exist_ok=True)

    # Data is loaded once and shared between both cohort analyses.
    with open(base_output_dir / 'analysis_log.txt', 'w', encoding='utf-8') as log_file:
        _LOG_FILE = log_file

        emit(f'Active biomarker: {biomarker_name}  ->  column: {biomarker_col}')
        emit(f'Output folder: {base_output_dir}')

        df_raw = load_data()

        run_nat_cohort_analysis(df_raw, col, biomarker_name, active_biomarker, nat_dir)
        run_upfront_cohort_analysis(df_raw, col, biomarker_name, active_biomarker, upfront_dir)

        emit()
        emit(f'Done. NAT cohort results:      {nat_dir}')
        emit(f'Done. Upfront surgery results: {upfront_dir}')

    _LOG_FILE = None


if __name__ == '__main__':
    # Optional: pass a biomarker key on the command line to override
    # ACTIVE_BIOMARKER, e.g.  python MMP8_NAT_analysis.py MMP8
    chosen = sys.argv[1] if len(sys.argv) > 1 else ACTIVE_BIOMARKER
    main(chosen)
