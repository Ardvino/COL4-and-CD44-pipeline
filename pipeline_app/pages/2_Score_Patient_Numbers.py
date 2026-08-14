import shutil
import sys
from datetime import datetime
from io import BytesIO
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _pipeline_loader import load_score_patno, PIPELINE_ROOT  # noqa: E402
from _ui import inject_base_css  # noqa: E402

spn = load_score_patno()

st.set_page_config(page_title="2. Score Patient Numbers", page_icon="\U0001f9ec", layout="wide")
inject_base_css()
st.title("2. Score patient numbers")
st.write(
    "Each IHC score file shares the same grid layout as the patient-numbered TMA "
    "map, so this stage joins them purely by **cell position** (sheet + row + "
    "column) — no explicit sample ID is needed in the score files. If a patient "
    "appears twice on the array (duplicate cores), both scores are kept as "
    "`Score` / `Score_2` on the same row, and `Highest_Score` is the max of the two."
)

st.subheader("Hand-off from Formatting")
formatting_output = PIPELINE_ROOT / "formatting" / "data" / "tma_map_replaced.xlsx"
local_map = spn.DEFAULT_DATA_DIR / spn.DEFAULT_MAP_FILE_NAME

if formatting_output.exists() and local_map.exists():
    fmt_mtime = formatting_output.stat().st_mtime
    local_mtime = local_map.stat().st_mtime
    if fmt_mtime > local_mtime:
        st.warning(
            f"`score_patNo/data/{spn.DEFAULT_MAP_FILE_NAME}` is older "
            f"({datetime.fromtimestamp(local_mtime):%Y-%m-%d %H:%M}) than the Formatting "
            f"stage's output ({datetime.fromtimestamp(fmt_mtime):%Y-%m-%d %H:%M}). "
            f"This stage reads its own local copy — sync it if you've re-run Formatting since."
        )
        if st.button("Sync from Formatting output"):
            shutil.copyfile(formatting_output, local_map)
            st.success("Copied. Re-run this stage below to use the refreshed map.")
            st.rerun()
    else:
        st.success(f"`score_patNo/data/{spn.DEFAULT_MAP_FILE_NAME}` is up to date with the Formatting stage's output.")
elif not local_map.exists():
    st.error(
        f"`score_patNo/data/{spn.DEFAULT_MAP_FILE_NAME}` is missing — run the Formatting stage "
        f"first, or copy a patient-numbered map into `score_patNo/data/`."
    )
    st.stop()
else:
    st.info("Formatting stage hasn't been run yet in this session — using the existing local copy.")

st.subheader("Inputs")
score_files = sorted([
    f.name for f in spn.DEFAULT_DATA_DIR.glob("*.xlsx")
    if not f.name.startswith("~$") and f.name != spn.DEFAULT_MAP_FILE_NAME
])
st.write(f"Score files detected in `score_patNo/data/`: {', '.join(f'`{f}`' for f in score_files)}")

st.subheader("Run")
if st.button("Run score_patno step", type="primary"):
    with st.spinner("Joining score files to the patient map by position..."):
        st.session_state["score_patno_result"] = spn.run_score_patno()

result = st.session_state.get("score_patno_result")

if result is None:
    st.info("Click **Run score_patno step** to process the files above.")
    st.stop()

st.subheader("Result")

m1, m2, m3, m4 = st.columns(4)
m1.metric("Score variables found", len(result.score_labels))
m2.metric("Positions before filtering", result.rows_before_filter)
m3.metric("Patients after filtering (excl. x/i)", result.rows_after_filter)
m4.metric("Patients with duplicate cores", result.duplicate_patient_count)

st.markdown("**Output preview** — `patient_scores.xlsx`")
tabs = st.tabs(result.score_labels + ["All_Scores_Summary"])
for tab, label in zip(tabs[:-1], result.score_labels):
    with tab:
        st.dataframe(result.per_score_tables[label], width="stretch", hide_index=True)
with tabs[-1]:
    st.dataframe(result.summary_table, width="stretch", hide_index=True)

from openpyxl import Workbook  # noqa: E402

wb_out = Workbook()
wb_out.remove(wb_out.active)
for label in result.score_labels:
    spn.write_sheet(wb_out, label, result.per_score_tables[label])
spn.write_sheet(wb_out, "All_Scores_Summary", result.summary_table)
buf = BytesIO()
wb_out.save(buf)
st.download_button(
    "Download patient_scores.xlsx",
    data=buf.getvalue(),
    file_name="patient_scores.xlsx",
    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
)
st.caption(f"Also saved to `{result.output_file}` (same file the standalone script writes).")
