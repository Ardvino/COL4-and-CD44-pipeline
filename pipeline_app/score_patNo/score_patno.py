import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook, Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = BASE_DIR / "data"
DEFAULT_MAP_FILE_NAME = "tma_map_replaced.xlsx"
DEFAULT_OUTPUT_FILE = BASE_DIR / "patient_scores.xlsx"


@dataclass
class ScorePatNoResult:
    """Summary of a run_score_patno() call, for both CLI printing and app display."""
    score_labels: list
    score_files: list                 # file names used as score sources
    rows_before_filter: int
    rows_after_filter: int
    duplicate_patient_count: int      # patients with more than one core (Score_2, ...)
    per_score_tables: dict            # {label: DataFrame} — one sheet per score file
    summary_table: object             # DataFrame — All_Scores_Summary
    output_file: Path


def read_all_cells(wb):
    """Return DataFrame with columns: sheet, row, col, value."""
    rows = []
    for sheet_name in wb.sheetnames:
        ws = wb[sheet_name]
        for row in ws.iter_rows():
            for cell in row:
                if cell.value is not None:
                    rows.append((sheet_name, cell.row, cell.column, cell.value))
    return pd.DataFrame(rows, columns=["sheet", "row", "col", "value"])


def _natural_key(s):
    """Zero-pad numeric parts so string sort matches numeric order."""
    return "".join(p.zfill(10) if p.isdigit() else p.lower()
                   for p in re.split(r"(\d+)", str(s)))


def flatten_duplicates(df, patient_col, value_cols):
    """
    Pivot duplicate patients into extra columns (Score_2, Score_3, …),
    append Highest_Score, and sort by patient number.
    """
    df = df.copy()
    df["_occ"] = df.groupby(patient_col).cumcount()

    frames = []
    for col in value_cols:
        pivoted = df.pivot_table(
            index=patient_col, columns="_occ", values=col, aggfunc="first"
        )
        pivoted.columns = [
            col if i == 0 else f"{col}_{i + 1}"
            for i, _ in enumerate(pivoted.columns)
        ]
        frames.append(pivoted)

    result = pd.concat(frames, axis=1).reset_index()

    score_cols = [c for c in result.columns if c != patient_col]
    result["Highest_Score"] = result[score_cols].max(axis=1)

    result = result.sort_values(
        patient_col, key=lambda s: s.map(_natural_key)
    ).reset_index(drop=True)

    return result


# --- Styling helpers ---
HEADER_FONT = Font(name="Arial", bold=True, color="FFFFFF")
HEADER_FILL = PatternFill("solid", start_color="2F4F8F")
DATA_FONT = Font(name="Arial", size=10)
CENTER = Alignment(horizontal="center", vertical="center")
_THIN = Side(style="thin", color="CCCCCC")
BORDER = Border(left=_THIN, right=_THIN, top=_THIN, bottom=_THIN)


def write_sheet(wb, title, df):
    ws = wb.create_sheet(title=title[:31])
    headers = list(df.columns)
    for c_idx, h in enumerate(headers, 1):
        cell = ws.cell(row=1, column=c_idx, value=h)
        cell.font = HEADER_FONT
        cell.fill = HEADER_FILL
        cell.alignment = CENTER
        cell.border = BORDER
    for r_idx, row_vals in enumerate(df.itertuples(index=False), 2):
        for c_idx, val in enumerate(row_vals, 1):
            cell = ws.cell(row=r_idx, column=c_idx, value=val)
            cell.font = DATA_FONT
            cell.alignment = CENTER
            cell.border = BORDER
    for col in ws.columns:
        max_len = max((len(str(c.value)) for c in col if c.value is not None), default=10)
        ws.column_dimensions[get_column_letter(col[0].column)].width = max_len + 4
    ws.freeze_panes = "A2"


def run_score_patno(data_dir=None, map_file=None, output_file=None):
    """Join TMA scoring files with the patient-number map by cell position.

    `data_dir` holds the patient map (`map_file`, default
    "tma_map_replaced.xlsx") plus any number of score `.xlsx` files — every
    other `.xlsx` file in `data_dir` is treated as a score source, named
    after its filename. Patient numbers "x" (missing core) / "i"
    (identifier) are excluded. Duplicate cores for the same patient are
    pivoted into Score/Score_2/... plus a Highest_Score column. Writes one
    sheet per score file plus an All_Scores_Summary sheet to `output_file`.
    """
    data_dir = Path(data_dir) if data_dir else DEFAULT_DATA_DIR
    map_file_name = Path(map_file).name if map_file else DEFAULT_MAP_FILE_NAME
    map_path = data_dir / map_file_name
    output_file = Path(output_file) if output_file else DEFAULT_OUTPUT_FILE

    score_files = sorted([
        f for f in os.listdir(data_dir)
        if f.endswith(".xlsx") and not f.startswith("~$") and f != map_file_name
    ])

    # Load patient number map
    map_wb = load_workbook(map_path, data_only=True)
    map_df = read_all_cells(map_wb)
    map_df = map_df.rename(columns={"value": "Patient_No"})

    # Load scores from each file and join by position
    score_dfs = {}
    for fname in score_files:
        label = os.path.splitext(fname)[0]
        score_wb = load_workbook(data_dir / fname, data_only=True)
        sc_df = read_all_cells(score_wb).rename(columns={"value": label})
        score_dfs[label] = sc_df

    # Merge all on position key
    combined = map_df.copy()
    for label, sc_df in score_dfs.items():
        combined = combined.merge(sc_df[["sheet", "row", "col", label]],
                                  on=["sheet", "row", "col"], how="left")

    # Convert score columns to numeric
    score_labels = list(score_dfs.keys())
    for col in score_labels:
        combined[col] = pd.to_numeric(combined[col], errors="coerce")

    rows_before_filter = len(combined)

    # Filter out missing ("x") and identifier ("i") patient numbers
    EXCLUDED = {"x", "i"}
    combined = combined[
        ~combined["Patient_No"].astype(str).str.strip().str.lower().isin(EXCLUDED)
    ]
    rows_after_filter = len(combined)

    duplicate_patient_count = int((combined.groupby("Patient_No").size() > 1).sum())

    wb_out = Workbook()
    wb_out.remove(wb_out.active)

    # One sheet per score file: Patient_No | Score | Score_2 (if duplicate) | Highest_Score
    per_score_tables = {}
    for label in score_labels:
        out_df = flatten_duplicates(combined[["Patient_No", label]], "Patient_No", [label])
        rename = {label: "Score"}
        for c in out_df.columns:
            if c.startswith(f"{label}_"):
                rename[c] = f"Score_{c[len(label) + 1:]}"
        out_df = out_df.rename(columns=rename)
        per_score_tables[label] = out_df
        write_sheet(wb_out, label, out_df)

    # Summary sheet: highest score per patient per score file
    summary_df = (
        combined[["Patient_No"] + score_labels]
        .groupby("Patient_No", as_index=False)
        .max()
        .sort_values("Patient_No", key=lambda s: s.map(_natural_key))
        .reset_index(drop=True)
    )
    summary_df["Highest_Score"] = summary_df[score_labels].max(axis=1)
    write_sheet(wb_out, "All_Scores_Summary", summary_df)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    wb_out.save(output_file)

    return ScorePatNoResult(
        score_labels=score_labels,
        score_files=score_files,
        rows_before_filter=rows_before_filter,
        rows_after_filter=rows_after_filter,
        duplicate_patient_count=duplicate_patient_count,
        per_score_tables=per_score_tables,
        summary_table=summary_df,
        output_file=output_file,
    )


def _print_result(result):
    print(f"Saved: {result.output_file}")
    print(f"Sheets: {result.score_labels + ['All_Scores_Summary']}")


if __name__ == "__main__":
    _print_result(run_score_patno())
