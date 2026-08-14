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
3. **`statistics/`** — survival analysis: Kaplan-Meier curves, univariable and multivariable Cox regression, subgroup analyses by NAT regimen, with a choice of low/high dichotomization method (fixed clinical cut-point, sample median, or sample 75th percentile).

All three stages live inside **`pipeline_app/`**, along with a Streamlit app (`pipeline_app/app.py`) that runs the real pipeline scripts end-to-end and shows each stage's inputs, transformation, and outputs inline — see [`pipeline_app/README.md`](pipeline_app/README.md). Each stage also still runs standalone as a plain script (see **Modules** below).

---

## Directory structure

```
programming/
│
├── README.md                        ← this file
│
├── pipeline_app/                     Streamlit UI over the pipeline below, plus the pipeline itself
│   ├── app.py                        landing page
│   ├── pages/                        one page per stage
│   ├── requirements.txt
│   ├── README.md                     how to run the app
│   │
│   ├── formatting/
│   │   ├── formatting.py             replaces YP-numbers with patient numbers
│   │   └── data/
│   │       ├── merged_map.xlsx       raw TMA block map (YP-numbers)
│   │       ├── potno_yp_merge.xlsx   YP → patient number lookup table
│   │       └── tma_map_replaced.xlsx output: map with patient numbers
│   │
│   ├── score_patNo/
│   │   ├── score_patno.py            joins map + score files by cell position
│   │   ├── README.md                 detailed script documentation
│   │   ├── patient_scores.xlsx       output: one sheet per score variable
│   │   └── data/
│   │       ├── tma_map_replaced.xlsx copy of formatting output (required)
│   │       ├── CD44_SFF_inflam_cells.xlsx
│   │       ├── CD44_SFF_tumor_int.xlsx
│   │       ├── CD44_VFF_tumor_int.xlsx
│   │       └── COL_4_stroma.xlsx
│   │
│   └── statistics/
│       ├── statistics.py             survival analysis script
│       ├── data/
│       │   ├── sorted_data.xlsx      patient-level dataset (270 rows, 122 cols)
│       │   └── patient_scores.xlsx   copy of score_patNo output
│       └── results/                  per-biomarker, per-dichotomization-method output
│
└── tumor_marker_app/                  standalone public demo (synthetic/uploaded data — unrelated to pipeline_app)
    └── streamlit_app/
```

---

## Modules

### `pipeline_app/formatting/formatting.py`

Reads the raw TMA block map (`merged_map.xlsx`) where each core position contains a YP-number (anonymised sample ID). Looks up the corresponding patient number from `potno_yp_merge.xlsx` and writes the substituted map to `tma_map_replaced.xlsx`, preserving all cell formatting.

**Run:**
```bash
cd pipeline_app/formatting
python formatting.py
```

**Requires:** `pandas`, `openpyxl`

---

### `pipeline_app/score_patNo/score_patno.py`

Joins the patient-numbered TMA map with any number of IHC score files placed in `data/`. The join is purely positional — sheet name + row + column — so no explicit sample ID is needed in the score files. Handles duplicate cores (same patient appearing twice on the array) by pivoting them into `Score` / `Score_2` columns and appending a `Highest_Score`. Produces a styled Excel workbook with one sheet per score variable and an `All_Scores_Summary` sheet.

See [`pipeline_app/score_patNo/README.md`](pipeline_app/score_patNo/README.md) for full details.

**Run:**
```bash
cd pipeline_app/score_patNo
python score_patno.py
```

**Requires:** `pandas`, `openpyxl`

---

### `pipeline_app/statistics/statistics.py`

Replication of statistical methods from Kesti et al. (2025) (*Scientific Reports*). Analyses disease-specific survival (DSS) in NAT-treated and upfront-surgery PDAC patients stratified by a chosen biomarker and NAT response.

The script supports multiple biomarkers (`MMP8`, `CD44_SFF_inflam_cells`, `CD44_SFF_tumor_int`, `CD44_VFF_tumor_int`, `COL_4_stroma`) and three low/high dichotomization methods (`fixed`, `median`, `p75` — see `DICHOTOMIZATION_METHODS` in the script).

**Run:**
```bash
cd pipeline_app/statistics
python statistics.py                    # uses ACTIVE_BIOMARKER as set in the file, fixed cut-point
python statistics.py MMP8                # override: run MMP-8 instead
python statistics.py COL_4_stroma median # override: biomarker + dichotomization method
```

**Requires:** `pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`, `lifelines`, `openpyxl`

```bash
pip install lifelines pandas numpy matplotlib scipy openpyxl
```
