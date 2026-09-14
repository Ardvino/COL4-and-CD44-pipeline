import contextlib
import io
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _pipeline_loader import load_statistics, PIPELINE_ROOT  # noqa: E402
from _ui import inject_base_css, render_glossary_expander, glossary_help, render_scores_handoff  # noqa: E402

stats = load_statistics()

st.set_page_config(page_title="4. Patient Data Export", page_icon="\U0001f9ec", layout="wide")
inject_base_css()
st.title("4. Patient data export")
st.write(
    "A flat, one-row-per-patient CSV — demographics, staging, treatment, and "
    "survival fields — for use outside this app (e.g. a stats package, or a "
    "manuscript's Table 1). This is the same `patient_data.csv` the Statistics "
    "stage writes, available here on its own without needing to run the full "
    "Kaplan-Meier/Cox analysis first."
)
render_glossary_expander([
    "TNM stage", "Adjuvant treatment", "Biomarker", "NAT (neoadjuvant therapy)", "OS (overall survival)",
])

st.subheader("Hand-off from Score Patient Numbers")
render_scores_handoff(
    stats.SCORES_PATH, PIPELINE_ROOT / "score_patNo" / "patient_scores.xlsx", key="sync-page4",
)

st.subheader("Options")
col1, col2 = st.columns(2)
with col1:
    biomarker_key = st.selectbox(
        "Biomarker (used for the biomarker_group column)",
        list(stats.BIOMARKER_CONFIGS.keys()),
        format_func=lambda k: stats.BIOMARKER_CONFIGS[k]["name"],
        index=list(stats.BIOMARKER_CONFIGS.keys()).index(stats.ACTIVE_BIOMARKER),
        help=glossary_help("Biomarker", "MMP-8", "H-score"),
    )
with col2:
    dichotomization = st.selectbox(
        "Low/high split method",
        list(stats.DICHOTOMIZATION_METHODS.keys()),
        format_func=lambda k: stats.DICHOTOMIZATION_METHODS[k],
        help=glossary_help("Dichotomization / cut-point"),
    )

output_dir = stats.results_dir_for(biomarker_key, dichotomization)
csv_path = output_dir / "patient_data.csv"

if csv_path.exists():
    mtime = datetime.fromtimestamp(csv_path.stat().st_mtime)
    # export_patient_data() itself lives in statistics.py, so if that file
    # has been edited (e.g. a column added) more recently than this CSV was
    # written, the CSV predates whatever changed and may be missing columns
    # -- same "is the cached output older than the thing that produces it"
    # check used for the score_patNo/Formatting hand-offs elsewhere in this
    # app, just against the script itself rather than an upstream stage.
    script_mtime = Path(stats.__file__).stat().st_mtime
    if csv_path.stat().st_mtime < script_mtime:
        st.warning(
            f"This export was generated {mtime:%Y-%m-%d %H:%M}, before the most recent "
            f"change to `statistics.py` ({datetime.fromtimestamp(script_mtime):%Y-%m-%d %H:%M}) "
            f"— it may be missing columns added since. Click **Generate / refresh export** "
            f"below to rebuild it."
        )
    else:
        st.caption(
            f"Existing export for **{stats.BIOMARKER_CONFIGS[biomarker_key]['name']}** "
            f"({stats.DICHOTOMIZATION_METHODS[dichotomization]}) — last generated {mtime:%Y-%m-%d %H:%M}."
        )

st.subheader("Run")
if st.button("Generate / refresh export", type="primary"):
    with st.spinner("Building patient-level export..."):
        bm = stats.BIOMARKER_CONFIGS[biomarker_key]
        continuous_col = f"{biomarker_key}_continuous" if biomarker_key in stats.BIOMARKER_CONTINUOUS else None
        col = stats.build_column_map(bm["col"], continuous_col=continuous_col)
        output_dir.mkdir(parents=True, exist_ok=True)

        # export_patient_data()/load_data() call emit(), which prints to
        # stdout when no run log file is open (see statistics.py) -- capture
        # that instead of letting it only land in the terminal, since it
        # names exactly which rows got dropped and why.
        log_buf = io.StringIO()
        with contextlib.redirect_stdout(log_buf):
            df_raw = stats.load_data(dichotomization=dichotomization)
            result = stats.export_patient_data(df_raw, col, output_dir)
        st.session_state["patient_export_result"] = result
        st.session_state["patient_export_log"] = log_buf.getvalue()
    st.success("Done.")

result = st.session_state.get("patient_export_result")
if result is None and not csv_path.exists():
    st.info("Click **Generate / refresh export** to build the export for this biomarker.")
    st.stop()
if result is None:
    result = pd.read_csv(csv_path)

st.subheader("Preview")
st.dataframe(result, width="stretch", hide_index=True)

if "patient_export_log" in st.session_state:
    with st.expander("Export log (rows dropped, and why)"):
        st.text(st.session_state["patient_export_log"])

st.download_button(
    "Download patient_data.csv",
    data=result.to_csv(index=False).encode("utf-8"),
    file_name=f"patient_data_{biomarker_key}_{dichotomization}.csv",
    mime="text/csv",
)
st.caption(f"Also saved to `{csv_path}`.")
