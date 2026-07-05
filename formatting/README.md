# formatting

Replaces anonymous YP-numbers in a raw TMA block map with real patient numbers, preserving all Excel cell formatting.

## How it works

The raw TMA map (`merged_map.xlsx`) uses YP-numbers (e.g. `YP1234`) as anonymised sample identifiers in each core position. A separate lookup file (`potno_yp_merge.xlsx`) maps each YP-number to a patient number (`pat_no`).

The script iterates every cell in every sheet of the map. Any cell whose value starts with `YP` (case-insensitive) is looked up and replaced in-place. All other cell content and formatting is left untouched. Unmatched YP-numbers are collected and printed as warnings without stopping the run.

The output (`tma_map_replaced.xlsx`) is used as the patient number map in the downstream `score_patNo` step.

## File structure

```
formatting/
├── formatting.py            # main script
└── data/
    ├── merged_map.xlsx      # input: raw TMA block map with YP-numbers
    ├── potno_yp_merge.xlsx  # input: YP → patient number lookup table
    └── tma_map_replaced.xlsx  # output: map with patient numbers substituted
```

### Input: `data/merged_map.xlsx`

The TMA block map exported from the pathology database. Each cell in the grid represents one core position; cells contain YP-numbers identifying the sample at that position.

### Input: `data/potno_yp_merge.xlsx`

Lookup table with at least two columns: `yp_num` and `pat_no`. Multiple sheets are supported — all are merged into one lookup dictionary. Column names are matched case-insensitively and leading/trailing whitespace is stripped.

### Output: `data/tma_map_replaced.xlsx`

A copy of `merged_map.xlsx` with every matched YP-number replaced by the corresponding patient number. Cell formatting (fonts, fills, borders) is preserved via `openpyxl`. Cells that could not be matched are left as-is; their coordinates are printed to stdout after the run.

## Running

```bash
cd formatting
python formatting.py
```

Expected output:

```
Sheet: 'Sheet1' — columns: ['yp_num', 'pat_no', ...] — rows: 250
Loaded 250 YP → pat_no mappings
Done! Replaced 480 cells → saved to data/tma_map_replaced.xlsx
```

If any YP-numbers had no match:

```
Warning: 3 YP-numbers had no match:
  Sheet 'Block_A' cell B4: YP9999
  ...
```

## Dependencies

```bash
pip install pandas openpyxl
```

| Package | Purpose |
|---|---|
| `pandas` | Reading all sheets from the lookup Excel file |
| `openpyxl` | Reading and writing the TMA map with formatting preserved |