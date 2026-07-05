# CD44 / COL4 Project — Programming

Analysis pipeline for tissue microarray (TMA) scoring and survival statistics in pancreatic ductal adenocarcinoma (PDAC) research.

## Pipeline overview

```
formatting/          →   score_patNo/         →   statistics/
YP-number map            patient score table       KM + Cox regression
(TMA positions)          (one row per patient)     (DSS, biomarker, NAT)
```

1. **`formatting/`** — replaces anonymous YP-numbers in the raw TMA block map with real patient numbers, producing `tma_map_replaced.xlsx`.
2. **`score_patNo/`** — joins the patient-numbered map with IHC score files (CD44, COL4, …) by cell position, producing a formatted `patient_scores.xlsx`.
3. **`statistics/`** — survival analysis notebook: Kaplan-Meier curves, univariable and multivariable Cox regression, subgroup analyses by NAT regimen.

---

## Directory structure

```
programming/
│
├── README.md                        ← this file
│
├── formatting/
│   ├── formatting.py                # replaces YP-numbers with patient numbers
│   └── data/
│       ├── merged_map.xlsx          # raw TMA block map (YP-numbers)
│       ├── potno_yp_merge.xlsx      # YP → patient number lookup table
│       └── tma_map_replaced.xlsx    # output: map with patient numbers
│
├── score_patNo/
│   ├── score_patno.py               # joins map + score files by cell position
│   ├── README.md                    # detailed script documentation
│   ├── patient_scores.xlsx          # output: one sheet per score variable
│   └── data/
│       ├── tma_map_replaced.xlsx    # copy of formatting output (required)
│       ├── CD44_SFF_inflam_cells.xlsx
│       ├── CD44_SFF_tumor_int.xlsx
│       ├── CD44_VFF_tumor_int.xlsx
│       └── COL_4_stroma.xlsx
│
└── statistics/
    ├── MMP8_NAT_analysis.ipynb      # survival analysis notebook
    ├── data/
    │   └── sorted_data.xlsx         # patient-level dataset (270 rows, 122 cols)
    ├── KM_MMP8_DSS.png              # Figure: KM curves by MMP-8 expression
    ├── KM_regimen_subgroup.png      # Figure: KM curves by NAT regimen
    └── ForestPlot_multivariable.png # Figure: multivariable Cox forest plot
```

---

## Modules

### `formatting/formatting.py`

Reads the raw TMA block map (`merged_map.xlsx`) where each core position contains a YP-number (anonymised sample ID). Looks up the corresponding patient number from `potno_yp_merge.xlsx` and writes the substituted map to `tma_map_replaced.xlsx`, preserving all cell formatting.

**Run:**
```bash
cd formatting
python formatting.py
```

**Requires:** `pandas`, `openpyxl`

---

### `score_patNo/score_patno.py`

Joins the patient-numbered TMA map with any number of IHC score files placed in `data/`. The join is purely positional — sheet name + row + column — so no explicit sample ID is needed in the score files. Handles duplicate cores (same patient appearing twice on the array) by pivoting them into `Score` / `Score_2` columns and appending a `Highest_Score`. Produces a styled Excel workbook with one sheet per score variable and an `All_Scores_Summary` sheet.

See [`score_patNo/README.md`](score_patNo/README.md) for full details.

**Run:**
```bash
cd score_patNo
python score_patno.py
```

**Requires:** `pandas`, `openpyxl`

---

### `statistics/MMP8_NAT_analysis.ipynb`

Replication of statistical methods from Kesti et al. (2025) (*Scientific Reports*). Analyses disease-specific survival (DSS) in 119 NAT-treated PDAC patients stratified by a chosen biomarker and NAT response.

The notebook supports multiple biomarkers. Switch the active marker by changing `ACTIVE_BIOMARKER` in the selector cell (section 2) and re-running from there.

| `ACTIVE_BIOMARKER` | Source |
|---|---|
| `'MMP8'` | Pre-coded binary column in `sorted_data.xlsx` |
| `'CD44_SFF_inflam_cells'` | `patient_scores.xlsx` — CD44, SFF, inflammatory cells |
| `'CD44_SFF_tumor_int'` | `patient_scores.xlsx` — CD44, SFF, tumour intensity |
| `'CD44_VFF_tumor_int'` | `patient_scores.xlsx` — CD44, VFF, tumour intensity |
| `'COL_4_stroma'` | `patient_scores.xlsx` — COL-4, stroma |

CD44 and COL-4 scores are loaded from `score_patNo/patient_scores.xlsx` (`All_Scores_Summary` sheet), merged with `sorted_data.xlsx` on patient number, and dichotomised at the median within NAT patients (0 = low, 1 = high).

> **Note on cut-points:** MMP-8 uses a pathologist-defined cut-point from the original paper. The CD44/COL-4 markers use a statistical median split — if the scoring protocol defines a specific cut-point for these markers, update the dichotomisation in section 1.5 accordingly.

| Section | Content |
|---|---|
| 1.5 | Load biomarker scores from `patient_scores.xlsx`, merge, and dichotomise |
| 2 | Biomarker selector (`ACTIVE_BIOMARKER`) + column mapping |
| 3 | Chi-square: biomarker vs clinical variables |
| 4 | Kaplan-Meier curves (all NAT / strong / weak responders) |
| 5 | Univariable Cox regression |
| 6 | Multivariable Cox regression + forest plot |
| 7 | Subgroup KM by NAT regimen (gemcitabine vs FOLFIRINOX) |
| 8 | Summary table (n, biomarker distribution, DSS events) |

**Requires:** `pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`, `lifelines`, `openpyxl`, `jinja2`

```bash
pip install lifelines pandas numpy matplotlib seaborn scipy openpyxl jinja2
```
