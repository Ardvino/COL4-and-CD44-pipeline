import sys
from datetime import datetime
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _pipeline_loader import PIPELINE_ROOT  # noqa: E402
from _ui import inject_base_css  # noqa: E402

st.set_page_config(page_title="CD44/COL4 Pipeline", page_icon="\U0001f9ec", layout="wide")
inject_base_css()

st.title("CD44 / COL4 Pipeline")
st.caption("A cockpit for the real TMA scoring & survival-analysis pipeline — every stage's real inputs, transforms, and outputs, in one place.")

st.write(
    "This app runs the project's actual pipeline scripts against the real data "
    "files already in this repo (nothing synthetic, nothing uploaded) and shows "
    "what each stage read, what it did, and what it produced. Use the sidebar to "
    "open a stage."
)

st.markdown(
    """
```
formatting/          →   score_patNo/         →   statistics/
YP-number map            patient score table       KM + Cox regression
(TMA positions)          (one row per patient)     (OS, biomarker, NAT)
```
"""
)


def _status_line(path: Path, empty_hint: str) -> str:
    if not path.exists():
        return f":gray[not yet generated — {empty_hint}]"
    mtime = datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
    return f":green[ready] — last generated {mtime}"


stages = [
    (
        "1. Formatting",
        "Replaces anonymised YP-numbers in the raw TMA block map with real patient "
        "numbers, preserving cell formatting.",
        PIPELINE_ROOT / "formatting" / "data" / "tma_map_replaced.xlsx",
        "run the Formatting stage",
    ),
    (
        "2. Score patient numbers",
        "Joins the patient-numbered map with IHC score files by cell position, "
        "producing one row per patient.",
        PIPELINE_ROOT / "score_patNo" / "patient_scores.xlsx",
        "run the Score patient numbers stage",
    ),
    (
        "3. Statistics",
        "Kaplan-Meier curves, chi-square balance checks, and univariable/"
        "multivariable Cox regression, per biomarker and treatment cohort.",
        PIPELINE_ROOT / "statistics" / "results",
        "run the Statistics stage",
    ),
]

cols = st.columns(3)
for col, (title, desc, out_path, hint) in zip(cols, stages):
    with col:
        st.subheader(title)
        st.write(desc)
        if out_path.name == "results":
            biomarkers = sorted(p.name for p in out_path.glob("*") if p.is_dir()) if out_path.exists() else []
            if biomarkers:
                st.markdown(f":green[ready] — results for: {', '.join(biomarkers)}")
            else:
                st.markdown(f":gray[not yet generated — {hint}]")
        else:
            st.markdown(_status_line(out_path, hint))
