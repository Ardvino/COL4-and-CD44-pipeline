# MMP-8 / CD44 / COL-4 Biomarker Survival Analysis

This script (`MMP8_NAT_analysis.py`) replicates the statistical analysis from:

> Kesti et al. (2025). *The prognostic significance of MMP-8 tissue
> immunoexpression in pancreatic ductal adenocarcinoma after neoadjuvant
> therapy.* Scientific Reports.

...and generalizes it so any of six immunohistochemistry biomarker scores
(MMP-8, two CD44 scores, a CD44 tumor-intensity variant, COL-4 stroma, and a
CD44 % positive tumor cells score) can be plugged in and run through the
same pipeline. It is meant to be read side-by-side with the code: each
section below matches a function or block in the script, in the order the
script executes them.

---

## 1. What question is being asked?

Does **high vs low expression of a given biomarker** in pancreatic ductal
adenocarcinoma (PDAC) tumor tissue predict **overall survival
(OS)** — i.e. time until death from any cause? (Kesti et al. (2025) used
disease-specific survival (DSS, death from PDAC only) as their endpoint;
this pipeline uses OS instead.)

The script answers this separately for two patient groups, because biology
and treatment context differ between them:

- **NAT cohort** — patients who received neoadjuvant therapy (chemo before
  surgery). This replicates the original paper's tables/figures, including
  sub-analyses by how well the tumor responded to NAT, and by chemo regimen.
- **Upfront surgery cohort** — patients who went straight to surgery with no
  NAT. This is a new addition (not in the original paper) to check whether
  the biomarker's prognostic pattern is specific to the NAT setting or holds
  in general.

---

## 2. Data inputs

| File | Role |
|---|---|
| `data/sorted_data.xlsx` | Main clinical dataset: one row per patient, treatment, staging, and survival/outcome columns. |
| `data/patient_scores.xlsx` (sheet `All_Scores_Summary`) | Biomarker immunohistochemistry scores per patient (`CD44_SFF_inflam_cells`, `CD44_SFF_tumor_int`, `CD44_VFF_tumor_int`, `COL_4_stroma`, `CD44_SFF_percentage`, `Highest_Score`). Joined onto the clinical data by patient number (`PotNo` ↔ `Patient_No`). |

`MMP8` itself is the exception: it already has a pre-computed binary column
(`MMP8inCAvahvin_low_vs_high`) in `sorted_data.xlsx`, so it does not go
through the dichotomization step described next.

### How biomarkers become "low" / "high"

Four of the scores in `patient_scores.xlsx` (plus `Highest_Score`, their
cross-file max) are a 4-tier IHC intensity grade (0 = non-detectable,
1 = low, 2 = moderate, 3 = high), the same scale used for MMP-8 in Kesti et
al. (2025). That paper dichotomizes MMP-8 as scores 0–1 = low vs. 2–3 = high
— a **fixed, pre-specified cut-point based on what the categories mean**,
not derived from the sample. `load_data()` applies that identical rule
(`BIOMARKER_HIGH_CUTOFF`, `score >= 2` for the four CD44/COL-4 IHC scores)
to those biomarkers, and the same cut-point is used for both the NAT and
upfront-surgery cohorts.

`CD44_SFF_percentage` is different: it's scored as % positive tumor cells
(0-100), not a 0-3 IHC grade, so both its valid range (`BIOMARKER_VALID_RANGE`)
and cut-point are configured separately from the other four. Its
`BIOMARKER_HIGH_CUTOFF` entry (`>=10%`) is only a placeholder for the
`'fixed'` dichotomization method — unlike `score >= 2` for the IHC scores,
it isn't a clinically pre-validated cut-point for this project. Use
`dichotomization='p75'` for this biomarker instead (75th-percentile split,
per Franklin et al.).

An earlier version of this script instead computed each score's median
*within NAT patients* and split on that. This was dropped: these are
coarse ordinal scores, not continuous measurements, and 41–64% of patients
tied exactly at the median for every score — meaning `>` vs `>=` on that
median swung some biomarkers' low/high split by up to 4x with no principled
way to choose a direction. See `biomarker_cutpoint_recommendations.md` for
the full analysis. If you want to explore a different cut-point for a
specific biomarker, edit its entry in `BIOMARKER_HIGH_CUTOFF` — just pick
it for a stated reason (e.g. a scoring protocol), not by scanning for the
best-looking split, which biases the result.

**Tie-aware `median`/`p75` cutoffs**: these two methods compute a quantile
of the sample (e.g. the 75th percentile for `p75`) and are *meant* to
produce a high-group of roughly `1 - q` of the sample (~25% for p75). On
coarse/rounded scores, though, the quantile value itself is often shared by
many patients — e.g. `CD44_SFF_percentage` is scored in 5-point increments
by eye, and in this sample 26% of patients tie exactly at the computed p75
value (80%). Always breaking ties toward "high" (`>=`) would sweep all of
them in, giving a ~48% high group — effectively a median split instead of a
quartile split. `load_data()` instead picks whichever comparison (`>=` or
`>` at the cutoff) lands the actual high-fraction closer to the quantile's
intended fraction, and logs which one it picked and how many patients were
affected (`tied at cutoff: N, counted as low/high`). This only applies to
`median`/`p75` — `fixed` always means `score >= cutoff` by definition.

**Continuous-exposure Cox regression for `CD44_SFF_percentage`**:
Franklin et al.'s own multivariable Cox model entered CD44s as a
**continuous** covariate ("per 25% positive cells", HR 1.9, p = 0.015 in
their Table 2), not the dichotomized low/high split used for their KM plot
— dichotomizing a continuous score is known to lose statistical power.
`BIOMARKER_CONTINUOUS` in `statistics.py` opts `CD44_SFF_percentage` into a
second, additive set of Cox tables fit on the raw score (scaled per 25
percentage points to match Franklin's units):
`cox_univariable_continuous_{NAT,upfront}.csv` and
`cox_multivariable_continuous_{NAT,upfront}.csv` (+ forest plot), alongside
— not replacing — the usual binarized low/high tables. No other biomarker
is affected.

**Wilcoxon–Breslow alongside log-rank**: every KM plot now reports both a
log-rank p-value and a Wilcoxon–Breslow p-value (`lifelines`'
`weightings='wilcoxon'`) in its title. Log-rank weights every event equally
over follow-up; Wilcoxon–Breslow weights earlier events more heavily, so
it's more sensitive to a survival difference concentrated early on — in
Franklin et al., osteopontin's log-rank was not significant (p = 0.0858)
but its Wilcoxon–Breslow was (p = 0.0322). The two tests can legitimately
disagree; both are shown rather than picking one.

**Data cleanup**: a few raw readings in `patient_scores.xlsx` fall outside
the valid range for their biomarker (found: PotNo 8, 56, 1345 outside 0–3 —
likely data-entry sentinels, e.g. `7`, `9`). `load_data()` now treats any
such out-of-range value as missing rather than as an extreme "high" score,
and logs which patients were affected. The valid range is per-biomarker
(`BIOMARKER_VALID_RANGE`) since `CD44_SFF_percentage` uses 0–100, not 0–3.

**Verify against the printed line**
`<biomarker>: cutoff=score>=X  low=N  high=N` — compare against a manual
count in Excel if you want a sanity check.

---

## 3. Column mapping (`build_column_map`)

All the "paper concepts" (age, sex, stage, etc.) are mapped to actual
spreadsheet column names in one place, with the coding scheme commented
next to each:

| Concept | Column | Coding |
|---|---|---|
| Received NAT | `NEOADJUVANTTI` | 1 = NAT, 0 = upfront surgery |
| Biomarker group | *(depends on `ACTIVE_BIOMARKER`)* | 0 = low, 1 = high |
| NAT response | `NATvaste_hyvä_012vs345` | 1 = strong (≤10% residual tumor cells), 0 = weak (≥11% RTC) |
| OS event | `os_event_binary` (derived) | 1 = died (any cause), 0 = alive/censored |
| Survival time | `Survival_Months` | months of follow-up |
| Age | `AGE_OPER` | years at surgery |
| Sex | `SUKUPUOLI` | 1 = male, 2 = female (verified against `Sukupuoli_tunnus` text column) |
| Stage | `STAGE_8th_Binary` | binarized AJCC 8th-edition stage (0 vs 1) |
| Grade | `GRADUS` | 1 (well) – 3 (poorly differentiated); raw column, used for the chi-square table only |
| Grade (Cox models) | `Grade_2`, `Grade_3` (derived) | dummy-coded in `prepare_cohort_dataset()`, reference = grade 1 — see §4 |
| Log CA19-9 | `logCA199` | log-transformed preop tumor marker |
| NAT regimen | `GEMSITABINvsFOLFIRINOX` | 0 = gemcitabine, 1 = FOLFIRINOX (NAT patients only) |

**If your spreadsheet's column names or codes differ, this is the only
place you need to edit them** — everything downstream just uses the `col[...]`
dictionary keys.

The raw `OS` column is already coded as a binary event (`1` = died of any
cause, `0` = alive) in `sorted_data.xlsx`, so `prepare_cohort_dataset()`
just carries it through as `os_event_binary`, used by every survival model.
**This is a key column to check**: a different OS coding scheme in a future
dataset version would silently break this. (The raw `DSS` column — `1` =
died of PDAC, `2` = alive, `3` = died of another cause — is not used by
this pipeline; it's the original paper's endpoint, not this one's.)

---

## 4. Pipeline, step by step

`main()` runs, per biomarker:

1. **`load_data()`** — reads both spreadsheets, merges them, computes the
   median-split binary column for every biomarker (see §2).
2. **`prepare_cohort_dataset(df, col, nat_flag=1 or 0)`** — filters to one
   cohort, builds `os_event_binary`, dummy-codes grade into `Grade_2` /
   `Grade_3` (reference = grade 1, so Cox models don't assume the 1→2 and
   2→3 steps have equal, linear effects on hazard — see §7), and drops any
   patient missing the biomarker score, `OS`, or survival time (prints how
   many patients were dropped at each stage — check these counts against
   your expectations).
3. **NAT cohort analysis** (`run_nat_cohort_analysis`):
   - Splits NAT patients into **strong** vs **weak** NAT responders.
   - Chi-square: biomarker group vs sex / stage / grade / regimen / response
     (checks whether the biomarker groups are balanced on these variables).
   - Kaplan-Meier: biomarker low vs high, in (A) all NAT patients, (B) strong
     responders, (C) weak responders — 3-panel figure with risk tables.
   - Univariable Cox regression: biomarker and each clinical variable, one
     at a time, in the full NAT group, then repeated within strong-only and
     weak-only responders.
   - Multivariable Cox regression: biomarker + age + sex + stage +
     grade (dummy-coded, 2 covariates) + logCA19-9, all together, in the
     full NAT group (subgroups are too small for a reliable multivariable
     model) — plus a forest plot and a proportional-hazards check.
   - Kaplan-Meier by NAT regimen (gemcitabine vs FOLFIRINOX), each showing
     strong vs weak response survival within that regimen.
   - Summary table: patient counts, biomarker distribution, OS events,
     median OS per group.
4. **Upfront surgery cohort analysis** (`run_upfront_cohort_analysis`) — the
   same logic, minus anything tied to NAT response/regimen (which don't
   apply to patients who never received NAT):
   - Chi-square: biomarker group vs sex / stage / grade.
   - Kaplan-Meier: biomarker low vs high, single panel.
   - Univariable Cox: biomarker + age + sex + stage + grade (dummy-coded) +
     logCA19-9.
   - Multivariable Cox: same covariates together, plus forest plot and a
     proportional-hazards check.
   - Summary table.

---

## 5. Output layout

```
results/<ACTIVE_BIOMARKER>/
├── analysis_log.txt                        # everything printed to console, for both cohorts
├── NAT_cohort/
│   ├── chi_square_results.csv
│   ├── KM_<biomarker>_OS_NAT.png            # 3-panel KM: all / strong / weak responders
│   ├── cox_univariable_all_NAT.csv
│   ├── cox_univariable_strong_responders.csv
│   ├── cox_univariable_weak_responders.csv
│   ├── cox_multivariable_NAT.csv
│   ├── ForestPlot_multivariable_NAT.png
│   ├── KM_regimen_subgroup.png              # gemcitabine vs FOLFIRINOX
│   └── group_summary.csv
└── Upfront_surgery_cohort/
    ├── chi_square_results.csv
    ├── KM_<biomarker>_OS_upfront.png
    ├── cox_univariable_upfront.csv
    ├── cox_multivariable_upfront.csv
    ├── ForestPlot_multivariable_upfront.png
    └── group_summary.csv
```

Every biomarker gets its own top-level folder, so re-running with a
different `ACTIVE_BIOMARKER` never overwrites a previous run's results.

## 6. Running it

Edit `ACTIVE_BIOMARKER` near the top of the script (or pass it as a command
line argument, which overrides the variable):

```bash
python MMP8_NAT_analysis.py                    # uses ACTIVE_BIOMARKER as set in the file
python MMP8_NAT_analysis.py MMP8                # override: run MMP-8 instead
python MMP8_NAT_analysis.py COL_4_stroma        # override: run COL-4 stroma instead
```

Valid keys: `MMP8`, `CD44_SFF_inflam_cells`, `CD44_SFF_tumor_int`,
`CD44_VFF_tumor_int`, `COL_4_stroma`, `CD44_SFF_percentage`.

---

## 7. How to interpret the results

### Chi-square table (`chi_square_results.csv`)
Tests whether the biomarker-low and biomarker-high groups are balanced with
respect to a clinical variable (sex, stage, grade, regimen, NAT response).

- **p ≥ 0.05**: no evidence of imbalance — the two biomarker groups look
  comparable on that variable, which is reassuring when you later attribute
  survival differences to the biomarker itself.
- **p < 0.05**: the biomarker groups differ on that variable — that
  variable may be a confounder, and should be considered in the
  multivariable Cox model (it already is, for age/sex/stage/grade/CA19-9).

### Kaplan-Meier plots (`KM_*.png`)
Each curve shows the probability of *not yet* having died (of any cause), over
time, for one biomarker group. Curves that separate and stay apart suggest
a survival difference; overlapping curves suggest none. The risk table
underneath shows, at each time point, how many patients are still being
followed ("At risk"), how many have dropped out without the event
("Censored"), and how many have had the event so far ("Events") — always
check this alongside the curve, since a tail with very few patients at risk
makes that part of the curve unreliable even if it looks visually dramatic.

- **Log-rank p** (in the panel title) tests whether the two full curves
  differ. p < 0.05 = statistically significant difference in survival
  between biomarker-low and biomarker-high groups.
- If a panel is blank with a "Skipped" message instead of curves, the
  subgroup being plotted didn't have both groups present after dropping
  missing values (e.g. every remaining patient in that subgroup happened to
  share the same biomarker status) — there's nothing to compare, so the
  plot and its log-rank test are skipped rather than crashing the run.

### Univariable Cox tables (`cox_univariable_*.csv`)
One row per variable, each from its own separate Cox model (not adjusted
for the others). Grade contributes two rows, "Grade 2 vs 1" and "Grade 3 vs
1" — both dummies are fit together in one model (grade 1 is the reference),
rather than treating `GRADUS` as a single continuous covariate, which would
force the 1→2 and 2→3 steps to have the same effect on hazard.

If a model has fewer than 10 events, or fails to converge, the row is
flagged (`sig` column explains why) instead of showing a HR — with that few
events the estimate isn't trustworthy.

- **HR (hazard ratio)**: relative risk of death (any cause) at any given
  moment, for the higher-coded group vs the reference (e.g. biomarker
  "high" [1] vs "low" [0]).
  - HR = 1: no difference.
  - HR > 1: higher risk (worse survival) in the "high"/coded-1 group.
  - HR < 1: lower risk (better survival) in the "high"/coded-1 group.
  - E.g. HR = 1.5 means roughly 50% higher instantaneous risk of death at
    any time; HR = 0.5 means roughly half the risk.
- **95% CI**: range of plausible HR values. If it excludes 1, the result is
  significant at the 5% level (this always agrees with the p-value column).
- **p**: formal significance test; `*` flags p < 0.05 in the console output.

### Multivariable Cox table (`cox_multivariable_*.csv` + forest plot)
Same HR/CI/p interpretation as above, but every covariate is entered into
*one* model together, so each HR is **adjusted for the others** — e.g. the
biomarker's HR here reflects its association with survival after
accounting for age, sex, stage, grade, and CA19-9. This is the more
rigorous test of whether the biomarker is an *independent* prognostic
factor, versus just riding on a correlated clinical variable.

- Printed alongside: **Concordance** (C-index, 0.5 = no better than chance,
  1.0 = perfect ranking of who dies first) and the **log-likelihood ratio
  test** (whether the whole model, all covariates together, explains
  survival better than no model at all).
- The forest plot shows each HR as a point with its 95% CI as a horizontal
  line; the vertical dashed line at HR = 1 is the "no effect" reference —
  a CI crossing that line means "not statistically significant."
- **Complete-case caveat**: this model only uses patients with *no* missing
  data in *any* of the covariates — check the printed "Complete cases"
  count, since it's typically well below the full cohort size (e.g. ~93 of
  115 NAT patients), which can shrink the effective sample and precision.
  If fewer than 10 events remain, the model is skipped entirely (logged as
  a warning) rather than reported.
- **Proportional-hazards check**: printed right after the coefficient table
  (`check_assumptions`, based on Schoenfeld residuals). A covariate flagged
  here has a hazard ratio that isn't constant over follow-up time, so that
  HR should be read as an average effect over the study period rather than
  a fixed risk multiplier — see the printed advice (e.g. stratifying on
  that covariate) if you need to lean on it. This is a diagnostic, not a
  gate: the model result is still reported either way.

### Summary table (`group_summary.csv`)
Descriptive: how many patients are in each group, what fraction are
biomarker-low/-high, how many OS events occurred, and median OS —
computed as the **Kaplan-Meier median survival time** (the time at which
the fitted survival curve crosses 0.5), with a 95% CI derived from the KM
curve's own confidence band. This correctly accounts for censoring, unlike
a plain median of observed times among patients who died (which ignores
censored patients entirely and is typically biased short). `inf` in the
median or CI columns means the KM curve never dropped to 0.5 within
follow-up — i.e. median survival wasn't reached.

### Practical caveats to keep in mind
- **Small subgroups**: strong NAT responders (n≈19 depending on biomarker)
  give wide, unstable confidence intervals — don't over-interpret a single
  non-significant p-value there as "no effect," it may just be underpowered.
  Cox models with fewer than 10 events (`MIN_EVENTS_FOR_COX` in the script)
  are skipped outright rather than reported, since the HR/CI would not be
  numerically reliable at that point; watch the log for "skipped"/"WARNING"
  lines.
- **No multiple-testing correction**: many tests are run across sections and
  biomarkers; a handful of p-values just under 0.05 across all of them is
  expected by chance alone. Focus on whether a result is consistent across
  cohorts/subgroups, not on any single p-value in isolation.
- **Proportional hazards assumption**: each multivariable Cox model is
  followed by a Schoenfeld-residuals check (`cph.check_assumptions()`), and
  the result is printed/logged. This is a diagnostic, not an automatic
  gate — a flagged covariate's HR should be read as an average effect over
  follow-up rather than a fixed multiplier; see the printed advice (e.g.
  stratifying on that covariate) if it matters for your conclusion.
- **Low/high cut-points** are fixed and pre-specified (`BIOMARKER_HIGH_CUTOFF`,
  §2), matching the original paper's MMP-8 convention rather than a
  per-sample median split. Group sizes can still be uneven for a given
  biomarker (e.g. `COL_4_stroma` is 15%/85% low/high in some cohorts) —
  that reflects the biology/scoring distribution, not an artifact of the
  cut-point choice; check the printed low/high counts if a result looks
  underpowered on one side.
