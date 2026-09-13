# score_patno

Joins TMA (Tissue Microarray) scoring files with a patient number map and produces a formatted Excel report.

## How it works

Each score file and the patient map share the same grid layout (same TMA block structure). The script joins them by cell position — sheet name + row + column — so each scored core is automatically paired with its patient number.

Patient numbers that are `x` (missing core) or `i` (identifier/control) are excluded. If a patient appears twice on the array (duplicate cores), both scores are kept as `Score` and `Score_2` on the same row.

## File structure

```
score_patNo/
├── score_patno.py        # main script
├── patient_scores.xlsx   # output (generated)
└── data/
    ├── tma_map_replaced.xlsx   # patient number map (required, fixed name)
    ├── CD44_SFF_inflam_cells.xlsx
    ├── CD44_SFF_percentage.xlsx
    ├── CD44_SFF_tumor_int.xlsx
    ├── CD44_VFF_tumor_int.xlsx
    └── COL_4_stroma.xlsx
```

### Input: `data/tma_map_replaced.xlsx`

The TMA block map with patient numbers already substituted into each cell. Produced by the `formatting/formatting.py` script upstream. Cells containing `x` or `i` are treated as non-patient entries and filtered out.

### Input: `data/*.xlsx` (score files)

Any `.xlsx` file in `data/` other than `tma_map_replaced.xlsx` is treated as a score file. The filename (without extension) becomes the column label in the output. Each file must have the **same grid layout** as the map — the join is purely positional.

To add a new scoring variable, drop its `.xlsx` file into `data/` and re-run the script.

### Output: `patient_scores.xlsx`

One sheet per score file, plus summary sheets:

| Sheet | Columns |
|---|---|
| `<score_file_name>` | `Patient_No`, `Score`, `Score_2`*, `Highest_Score` |
| `All_Scores_Summary` | `Patient_No`, `<label1>`, `<label2>`, …, `Highest_Score` |
| `H_Scores_Summary`† | `Patient_No`, `<hscore_label1>`, `<hscore_label1>_2`*, `<hscore_label1>_Highest`, `<hscore_label2>`, … |

\* `Score_2` appears only when a patient has two cores on the array. `Highest_Score` is the max across `Score` and `Score_2`.

† Only written when at least one score file's name ends in `_hscore`.

In `All_Scores_Summary`, each score column shows the **highest** value across duplicate cores for that patient. `Highest_Score` is the max across all score columns *except* any whose filename ends in `_percentage` (0-100 scale) or `_hscore` (0-300 scale) — those aren't on the shared 0-3 IHC intensity scale the other files use, so mixing them into that max would make it meaningless. They still get their own sheet/column like any other score file.

`H_Scores_Summary` gives one row per patient with every `_hscore` marker's per-core values (`<label>`, `<label>_2`, …) plus a `<label>_Highest`, so a new H-score marker's values across the whole cohort can be exported/compared without the other score types in the way. Unlike `All_Scores_Summary` (which reports only the max core per patient per marker), this sheet keeps each core's value — built from the per-marker sheets, not from `All_Scores_Summary`. Add another `_hscore`-suffixed file to `data/` and it appears here automatically, no code changes needed.

Rows are sorted by patient number (natural order, so `9` < `10` < `100`).

## Running

```bash
cd pipeline_app/score_patNo
python score_patno.py
```

Requires: `pandas`, `openpyxl`

```bash
pip install pandas openpyxl
```
