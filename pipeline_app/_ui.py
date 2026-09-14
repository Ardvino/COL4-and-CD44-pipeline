"""Shared page chrome — the small bits every page in this app injects, kept
in one place so they don't drift between pages.
"""
import shutil
from datetime import datetime

import streamlit as st

# Plain-English definitions for project-specific jargon that shows up in the
# app's own rendered text -- for a reader outside this research project, not
# just the original team. Single source of wording so pages can't drift
# apart -- see render_glossary_expander() / glossary_help() below.
GLOSSARY = {
    "YP-number": "Anonymised sample ID on the raw TMA block map, before it's matched to a real patient number.",
    "PotNo": "Abbreviation for \"patient number\" (Finnish *potilasnumero*) used in filenames/columns -- same concept as \"patient number\" elsewhere in the app.",
    "TMA (tissue microarray)": "A single slide holding many small tumor-tissue cores from different patients, arranged in a grid, so they can all be stained and scored together.",
    "IHC score": "A microscopy-based score of how strongly a protein stains in tissue, scored per TMA core.",
    "H-score": "A 0-300 IHC intensity score combining staining intensity and the percentage of positive cells -- a standard pathology scoring method.",
    "Score / Score_2 / Highest_Score": "When a patient has two TMA cores, both scores are kept as separate columns (Score, Score_2), and the higher of the two is also reported separately (Highest_Score).",
    "NAT (neoadjuvant therapy)": "Chemotherapy given *before* surgery, to shrink the tumor first.",
    "Upfront surgery": "Surgery performed as the first treatment, without prior chemotherapy.",
    "RTC (residual tumor cells)": "The percentage of viable tumor cells remaining after NAT -- a lower percentage means the chemotherapy worked better.",
    "MMP-8": "One of the biomarkers analyzed -- a protein-degrading enzyme (matrix metalloproteinase 8).",
    "Biomarker": "A measurable protein/molecular marker (here, CD44, COL-4, MMP-8, or GATA6 variants) tested for association with survival.",
    "Dichotomization / cut-point": "The rule used to split a continuous or ordinal score into \"low\" vs \"high\" groups for analysis.",
    "Chi-square (balance check)": "A statistical test checking whether the biomarker-low and biomarker-high groups are evenly distributed across a clinical variable -- an uneven split flags a possible confounder.",
    "Kaplan-Meier (KM) curve": "A step-shaped plot of the estimated probability of still being alive (or event-free) over time, accounting for patients whose follow-up ended early (censoring).",
    "Log-rank test": "A statistical test comparing two Kaplan-Meier curves over their entire follow-up period, not just at one time point.",
    "Cox regression (HR)": "A survival-analysis model estimating each variable's hazard ratio (HR) -- how much it multiplies the risk of the event. \"Univariable\" tests one variable at a time; \"multivariable\" adjusts for the others simultaneously.",
    "Proportional-hazards / Schoenfeld residuals": "Cox regression's required assumption (that a variable's effect on risk is constant over time), and Schoenfeld residuals are the standard diagnostic check for whether it holds.",
    "Concordance": "A 0.5-1.0 score measuring how well the Cox model's predicted risk ranks patients correctly -- similar in spirit to AUC.",
    "OS (overall survival)": "Time from a reference point (e.g. surgery) until death from any cause. This pipeline reports OS throughout -- the original Kesti et al. (2025) study it replicates used disease-specific survival (DSS, death from PDAC only) instead.",
    "TNM stage": "The AJCC/TNM staging system combines tumor extent (T), lymph node involvement (N), and metastasis (M) into an overall stage (e.g. IIB, III) -- exported at full sub-stage resolution, not the collapsed low/high split used elsewhere in this app.",
    "Adjuvant treatment": "Chemotherapy given *after* surgery (as opposed to NAT, given before). Reported as a category (adjuvant / none / palliative) and, where available, as free-text notes on the actual regimen given.",
}


def render_glossary_expander(terms=None, label="Glossary"):
    """Collapsed st.expander with a definition list for `terms` (a list of
    GLOSSARY keys), or the full GLOSSARY if `terms` is omitted. Pulls
    wording from the shared GLOSSARY dict so it can't drift page-to-page.
    """
    entries = GLOSSARY if terms is None else {t: GLOSSARY[t] for t in terms}
    with st.expander(label):
        for term, definition in entries.items():
            st.markdown(f"**{term}** — {definition}")


def glossary_help(*terms):
    """Joins GLOSSARY[term] for each term into a string for a widget's
    help= kwarg, e.g. help=glossary_help('Biomarker', 'MMP-8').
    """
    return "  \n".join(f"**{t}**: {GLOSSARY[t]}" for t in terms)


def render_scores_handoff(local_scores, upstream_output, *, key):
    """Hand-off check shared by every page that reads statistics/data/
    patient_scores.xlsx (via stats.load_data()) -- is that local copy in
    sync with the score_patNo stage's own output? Shows a warning + sync
    button, a success message, an info message, or (if the local copy
    doesn't exist at all) an error that also calls st.stop().

    `local_scores` / `upstream_output`: Path objects (statistics.SCORES_PATH,
    and PIPELINE_ROOT / "score_patNo" / "patient_scores.xlsx").
    `key`: unique Streamlit widget key -- required since more than one page
    renders this block, and Streamlit needs distinct keys for their buttons.
    """
    if upstream_output.exists() and local_scores.exists():
        spn_mtime = upstream_output.stat().st_mtime
        local_mtime = local_scores.stat().st_mtime
        if spn_mtime > local_mtime:
            st.warning(
                f"`statistics/data/patient_scores.xlsx` is older "
                f"({datetime.fromtimestamp(local_mtime):%Y-%m-%d %H:%M}) than the Score Patient "
                f"Numbers stage's output ({datetime.fromtimestamp(spn_mtime):%Y-%m-%d %H:%M}). "
                f"This stage reads its own local copy — sync it if you've re-run that stage since "
                f"(e.g. added a new marker's score file)."
            )
            if st.button("Sync from Score Patient Numbers output", key=key):
                shutil.copyfile(upstream_output, local_scores)
                st.success("Copied. Re-run below to use the refreshed scores.")
                st.rerun()
        else:
            st.success("`statistics/data/patient_scores.xlsx` is up to date with the Score Patient Numbers stage's output.")
    elif not local_scores.exists():
        st.error(
            f"`statistics/data/patient_scores.xlsx` is missing — run the Score Patient Numbers stage "
            f"first, or copy `patient_scores.xlsx` into `statistics/data/`."
        )
        st.stop()
    else:
        st.info("Score Patient Numbers stage hasn't been run yet in this session — using the existing local copy.")


def inject_base_css():
    """Caps the content column at 1000px, left-aligned, on every page.

    The centering isn't a margin on the block container itself (a `margin`
    override doesn't fix it, and doesn't need to) -- Streamlit's `stMain`
    section is a flex column with `align-items: center`, which centers the
    (width-capped) block container as a flex child regardless of its own
    margin. Overriding that child's `align-self` to `flex-start` is what
    actually opts it out of the parent's centering.
    """
    st.markdown(
        """
        <style>
        [data-testid="stMainBlockContainer"] {
            max-width: 1000px !important;
            align-self: flex-start !important;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
