# Pipeline app

A single Streamlit app that drives the project's real 3-stage pipeline —
[`formatting/formatting.py`](formatting/formatting.py) →
[`score_patNo/score_patno.py`](score_patNo/score_patno.py) →
[`statistics/statistics.py`](statistics/statistics.py) — and shows each
stage's inputs, transformation, and outputs inline, instead of requiring you
to run three separate scripts and dig through the resulting files by hand.
All three live inside this folder — this app is self-contained, not a UI
layer over sibling top-level project folders.

**This app is for local/internal use only.** It reads and writes real
project data files (patient numbers, IHC scores, clinical variables)
directly off disk — there is no upload flow and it is not meant to be
deployed as a public-facing web service. **The data files themselves are
not included in this repository** (see
[Bringing your own data](#bringing-your-own-data) below) — they're
gitignored patient/clinical data specific to this research project. (The
unrelated
[`tumor_marker_app/streamlit_app`](../tumor_marker_app/streamlit_app) is the
public-facing demo, built on synthetic/uploaded data — this app doesn't
touch it.)

## Running

From the repo root, with the project's dependencies installed
(`pip install -r pipeline_app/requirements.txt`):

```bash
streamlit run pipeline_app/app.py
```

## Bringing your own data

This repo ships with code only — no patient data. To run the pipeline
against your own dataset, place files in the shapes below (each stage's own
README has the full column/sheet spec):

| Stage | Where | What | Details |
|---|---|---|---|
| 1. Formatting | `formatting/data/` | `merged_map.xlsx` (raw TMA block map, YP-numbers) + `potno_yp_merge.xlsx` (YP-number → patient-number lookup) | [`formatting/README.md`](formatting/README.md) |
| 2. Score patient numbers | `score_patNo/data/` | one `.xlsx` per IHC score variable, same grid layout as the patient-numbered map | [`score_patNo/README.md`](score_patNo/README.md) |
| 3. Statistics | `statistics/data/` | `sorted_data.xlsx` (patient-level clinical dataset) + `patient_scores.xlsx` (score_patNo output) | [`statistics/README.md`](statistics/README.md) |

Each page also detects and displays a friendly message (instead of
crashing) when its required input files are missing — open stage 1 in the
sidebar for a concrete walkthrough of what's expected.

## What it does — and doesn't — change

Each page's "Run" button calls the *same* function the corresponding CLI
script calls, so results are guaranteed to match `python formatting.py` /
`python score_patno.py` / `python statistics.py <biomarker>` exactly. Running
a stage from this app **writes the same real output files** those scripts
write (`formatting/data/tma_map_replaced.xlsx`,
`score_patNo/patient_scores.xlsx`, `statistics/results/<biomarker>/...`) —
it's a UI over the real pipeline, not a sandboxed copy.

`statistics/statistics.py` itself is never modified by this app — the
Statistics page only calls its existing `main()` and reads back the CSV/PNG
files it already writes to `results/<biomarker>/`.

## Pages

| Page | Shows |
|---|---|
| `app.py` | Pipeline overview + which stages have been run and when |
| `pages/1_Formatting.py` | Lookup table & TMA map previews, replacement counts, unmatched YP-numbers, resulting map preview |
| `pages/2_Score_Patient_Numbers.py` | Formatting→here hand-off freshness check, detected score files, per-file and summary tables |
| `pages/3_Statistics.py` | Biomarker selector, chi-square/KM/Cox results per cohort, full run log, caveats |
