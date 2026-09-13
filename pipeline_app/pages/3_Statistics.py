import shutil
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _pipeline_loader import load_statistics, PIPELINE_ROOT  # noqa: E402
from _ui import inject_base_css  # noqa: E402

stats = load_statistics()

st.set_page_config(page_title="3. Statistics", page_icon="\U0001f9ec", layout="wide")
inject_base_css()

st.markdown(
    """
    <style>
    [class*="st-key-plot-card-"] {
        border: 1px solid #383835;
        border-radius: 10px;
        padding: 1.1rem 1.25rem;
    }
    [class*="st-key-plot-wrap-"] [data-testid="stDownloadButton"] {
        margin-top: 0.35rem;
    }
    [class*="st-key-plot-wrap-"] [data-testid="stDownloadButton"] button {
        border: none;
        background: transparent;
        box-shadow: none;
        padding: 0.25rem;
    }
    [class*="st-key-plot-wrap-"] [data-testid="stDownloadButton"] button:hover {
        background: transparent;
        color: #3987e5;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("3. Statistics")
st.write(
    "Survival analysis, replicating the methods of Kesti et al. (2025): Kaplan-Meier "
    "curves, a biomarker-vs-clinical-variable balance check (chi-square), and "
    "univariable/multivariable Cox proportional-hazards regression — run separately "
    "for patients who received neoadjuvant therapy (NAT) and those who went straight "
    "to upfront surgery."
)

st.subheader("Hand-off from Score Patient Numbers")
score_patno_output = PIPELINE_ROOT / "score_patNo" / "patient_scores.xlsx"
local_scores = stats.SCORES_PATH

if score_patno_output.exists() and local_scores.exists():
    spn_mtime = score_patno_output.stat().st_mtime
    local_mtime = local_scores.stat().st_mtime
    if spn_mtime > local_mtime:
        st.warning(
            f"`statistics/data/patient_scores.xlsx` is older "
            f"({datetime.fromtimestamp(local_mtime):%Y-%m-%d %H:%M}) than the Score Patient "
            f"Numbers stage's output ({datetime.fromtimestamp(spn_mtime):%Y-%m-%d %H:%M}). "
            f"This stage reads its own local copy — sync it if you've re-run that stage since "
            f"(e.g. added a new marker's score file, like GATA6)."
        )
        if st.button("Sync from Score Patient Numbers output"):
            shutil.copyfile(score_patno_output, local_scores)
            st.success("Copied. Re-run the analysis below to use the refreshed scores.")
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

col1, col2 = st.columns(2)
with col1:
    biomarker_key = st.selectbox(
        "Biomarker",
        list(stats.BIOMARKER_CONFIGS.keys()),
        format_func=lambda k: stats.BIOMARKER_CONFIGS[k]["name"],
        index=list(stats.BIOMARKER_CONFIGS.keys()).index(stats.ACTIVE_BIOMARKER),
    )
with col2:
    dichotomization = st.selectbox(
        "Low/high split method",
        list(stats.DICHOTOMIZATION_METHODS.keys()),
        format_func=lambda k: stats.DICHOTOMIZATION_METHODS[k],
        help=(
            "How each CD44/COL-4 score becomes 'low' vs 'high'. Doesn't apply to MMP-8, "
            "which already has a pre-coded binary column. 'Fixed' is the default used "
            "throughout this project's write-ups; the other two are sample-derived, so "
            "switching methods can move the cut-point itself, not just the counts — see "
            "the run log below for how many patients land exactly on the cut."
        ),
    )

base_dir = stats.results_dir_for(biomarker_key, dichotomization)
nat_dir = base_dir / "NAT_cohort"
upfront_dir = base_dir / "Upfront_surgery_cohort"
log_path = base_dir / "analysis_log.txt"

# Checked via log_path (this method's own results), not base_dir -- for a
# non-'fixed' method, base_dir is a subfolder of the biomarker's root
# (results/<biomarker>/<method>/), but 'fixed' results sit directly in
# results/<biomarker>/. So the biomarker's root folder can already exist
# (as the parent of another method's subfolder) even when this exact
# method's own analysis_log.txt has never been written.
if log_path.exists():
    log_mtime = datetime.fromtimestamp(log_path.stat().st_mtime)
    st.caption(
        f"Existing results for **{stats.BIOMARKER_CONFIGS[biomarker_key]['name']}** "
        f"({stats.DICHOTOMIZATION_METHODS[dichotomization]}) — last generated {log_mtime:%Y-%m-%d %H:%M}."
    )

if st.button("Run analysis", type="primary"):
    with st.spinner(f"Running NAT + upfront-surgery cohort analysis for {stats.BIOMARKER_CONFIGS[biomarker_key]['name']}..."):
        stats.main(biomarker_key, dichotomization=dichotomization)
    st.success("Done.")

if not log_path.exists():
    st.info("Click **Run analysis** to generate results for this biomarker.")
    st.stop()

with st.expander("Full run log (includes concordance & proportional-hazards diagnostics)"):
    st.text(log_path.read_text(encoding="utf-8"))

st.info(
    "**Caveats to keep in mind:** small subgroups (e.g. strong NAT responders) give wide, "
    "unstable confidence intervals; no multiple-testing correction is applied across the "
    "many tests run per biomarker; and each multivariable Cox model's proportional-hazards "
    "assumption is checked (see the run log) but not auto-enforced — a flagged covariate's "
    "hazard ratio should be read as an average effect over follow-up."
)


@st.cache_data(show_spinner=False)
def _load_cohorts(biomarker_key, dichotomization):
    """Rebuilds the NAT/upfront cohort dataframes via the pipeline's own,
    unmodified load_data()/prepare_cohort_dataset() — the exact same calls
    main() makes internally — so the plots below are guaranteed to match the
    CSV tables above. Cached per (biomarker, dichotomization) since this
    re-reads both source spreadsheets.
    """
    bm = stats.BIOMARKER_CONFIGS[biomarker_key]
    col = stats.build_column_map(bm["col"])
    df_raw = stats.load_data(dichotomization=dichotomization)
    nat = stats.prepare_cohort_dataset(df_raw, col, nat_flag=1)
    upfront = stats.prepare_cohort_dataset(df_raw, col, nat_flag=0)
    return col, bm, nat, upfront


col, bm, nat, upfront = _load_cohorts(biomarker_key, dichotomization)
labels = bm["labels"]
strong = nat[nat[col["nat_resp"]] == 1]
weak = nat[nat[col["nat_resp"]] == 0]
gem = nat[nat[col["regimen"]] == 0]
folf = nat[nat[col["regimen"]] == 1]
mv_covariates = [col["mmp8"], col["age"], col["sex"], col["stage"], *col["grade_dummies"], col["logca199"]]


def render_themed_plot(build_fn, key_prefix, file_slug):
    """Builds the same plot in both PLOT_THEMES entries via `build_fn(theme)`,
    shows the dark version inline (matching the app's dark UI) and offers the
    light version as an SVG download (for use on a white page/print). Icon-only
    download button, matching the project's existing demo app pattern.

    Rendered via st.image (not raw st.markdown of the SVG) so Streamlit's
    built-in image toolbar gives a click-to-fullscreen preview for free.
    """
    dark_fig = build_fn(stats.PLOT_THEMES["dark"])
    if dark_fig is None:
        st.info("Not enough events in this cohort to fit a reliable model for this plot.")
        return
    light_fig = build_fn(stats.PLOT_THEMES["light"])
    dark_svg = stats.fig_to_svg(dark_fig)
    light_svg = stats.fig_to_svg(light_fig)
    plt.close(dark_fig)
    plt.close(light_fig)

    with st.container(key=f"plot-wrap-{key_prefix}"):
        with st.container(key=f"plot-card-{key_prefix}"):
            st.image(dark_svg, width="stretch")
        st.download_button(
            "", data=light_svg, file_name=f"{file_slug}.svg", mime="image/svg+xml",
            icon=":material/download:", help="Download plot (SVG, light background)",
            key=f"download-{key_prefix}",
        )


def _show_csv(dir_path, filename, caption=None):
    path = dir_path / filename
    if not path.exists():
        return
    if caption:
        st.caption(caption)
    st.dataframe(pd.read_csv(path), width="stretch", hide_index=True)


def render_nat_cohort():
    st.subheader("Baseline characteristics")
    _show_csv(nat_dir, "baseline_characteristics.csv")

    st.subheader("Biomarker vs clinical variables (chi-square)")
    _show_csv(
        nat_dir, "chi_square_results.csv",
        "p < 0.05 means the biomarker-low and biomarker-high groups are not evenly "
        "distributed across that variable — a potential confounder.",
    )

    st.subheader("Kaplan-Meier survival — all NAT patients, strong responders, weak responders")
    render_themed_plot(
        lambda theme: stats.themed_km_figure(
            [(nat, "(A) All NAT patients"),
             (strong, "(B) Strong NAT response (<=10% RTC)"),
             (weak, "(C) Weak NAT response (>=11% RTC)")],
            col["os_time"], col["os_event"], col["mmp8"], labels, theme,
            suptitle=f'{bm["name"]} and Overall Survival (OS) — NAT cohort',
        ),
        key_prefix=f"km-nat-{biomarker_key}", file_slug=f"KM_{biomarker_key}_OS_NAT",
    )

    st.subheader("Univariable Cox regression")
    st.caption("Each variable tested in its own model, one at a time, without adjusting for the others.")
    tab_all, tab_strong, tab_weak = st.tabs(["All NAT patients", "Strong responders", "Weak responders"])
    with tab_all:
        _show_csv(nat_dir, "cox_univariable_all_NAT.csv")
    with tab_strong:
        _show_csv(nat_dir, "cox_univariable_strong_responders.csv")
    with tab_weak:
        _show_csv(nat_dir, "cox_univariable_weak_responders.csv")

    st.subheader("Multivariable Cox regression (all NAT patients)")
    st.caption("All covariates entered together, so each hazard ratio is adjusted for the others.")
    _show_csv(nat_dir, "cox_multivariable_NAT.csv")
    render_themed_plot(
        lambda theme: stats.themed_forest_figure(
            nat, col["os_time"], col["os_event"], mv_covariates,
            "Multivariable Cox Regression - OS (all NAT patients)", theme,
        ),
        key_prefix=f"forest-nat-{biomarker_key}", file_slug=f"ForestPlot_{biomarker_key}_NAT",
    )

    if (nat_dir / "cox_univariable_continuous_NAT.csv").exists():
        st.subheader("Cox regression — continuous exposure")
        st.caption(
            "Fit on the raw (non-dichotomized) score rather than the low/high split above — "
            "matching how Franklin et al. modeled CD44s in their own multivariable model, "
            "since dichotomizing a continuous score can lose statistical power."
        )
        _show_csv(nat_dir, "cox_univariable_continuous_NAT.csv", "Univariable")
        _show_csv(nat_dir, "cox_multivariable_continuous_NAT.csv", "Multivariable, adjusted")

    st.subheader("NAT response by chemotherapy regimen")
    resp_labels = {1: "Strong response (<=10%)", 0: "Weak response (>=11%)"}
    render_themed_plot(
        lambda theme: stats.themed_km_figure(
            [(gem, "(A) Gemcitabine"), (folf, "(B) FOLFIRINOX")],
            col["os_time"], col["os_event"], col["nat_resp"], resp_labels, theme,
            suptitle="NAT Response and OS by Regimen",
        ),
        key_prefix=f"km-regimen-{biomarker_key}", file_slug=f"KM_{biomarker_key}_regimen_subgroup",
    )

    st.subheader("Group summary")
    st.caption("Median OS is the Kaplan-Meier median survival time, which correctly accounts for censoring.")
    _show_csv(nat_dir, "group_summary.csv")


def render_upfront_cohort():
    st.subheader("Baseline characteristics")
    _show_csv(upfront_dir, "baseline_characteristics.csv")

    st.subheader("Biomarker vs clinical variables (chi-square)")
    _show_csv(upfront_dir, "chi_square_results.csv")

    st.subheader("Kaplan-Meier survival — upfront surgery patients")
    render_themed_plot(
        lambda theme: stats.themed_km_figure(
            [(upfront, "Upfront surgery patients")],
            col["os_time"], col["os_event"], col["mmp8"], labels, theme,
            suptitle=f'{bm["name"]} and Overall Survival (OS) — Upfront surgery',
            figsize=(7, 6),
        ),
        key_prefix=f"km-upfront-{biomarker_key}", file_slug=f"KM_{biomarker_key}_OS_upfront",
    )

    st.subheader("Univariable Cox regression")
    _show_csv(upfront_dir, "cox_univariable_upfront.csv")

    st.subheader("Multivariable Cox regression")
    _show_csv(upfront_dir, "cox_multivariable_upfront.csv")
    render_themed_plot(
        lambda theme: stats.themed_forest_figure(
            upfront, col["os_time"], col["os_event"], mv_covariates,
            "Multivariable Cox Regression - OS (upfront surgery patients)", theme,
        ),
        key_prefix=f"forest-upfront-{biomarker_key}", file_slug=f"ForestPlot_{biomarker_key}_upfront",
    )

    if (upfront_dir / "cox_univariable_continuous_upfront.csv").exists():
        st.subheader("Cox regression — continuous exposure")
        st.caption(
            "Fit on the raw (non-dichotomized) score rather than the low/high split above — "
            "matching how Franklin et al. modeled CD44s in their own multivariable model, "
            "since dichotomizing a continuous score can lose statistical power."
        )
        _show_csv(upfront_dir, "cox_univariable_continuous_upfront.csv", "Univariable")
        _show_csv(upfront_dir, "cox_multivariable_continuous_upfront.csv", "Multivariable, adjusted")

    st.subheader("Group summary")
    _show_csv(upfront_dir, "group_summary.csv")


tab_nat, tab_upfront = st.tabs(["NAT cohort", "Upfront surgery cohort"])
with tab_nat:
    render_nat_cohort()
with tab_upfront:
    render_upfront_cohort()

patient_data_path = base_dir / "patient_data.csv"
if patient_data_path.exists():
    st.subheader("Exported patient-level dataset")
    st.caption(
        "One row per patient across both cohorts, in the standard template format "
        "(same schema the standalone demo app accepts as an upload)."
    )
    patient_df = pd.read_csv(patient_data_path)
    st.dataframe(patient_df, width="stretch", hide_index=True)
    st.download_button(
        "Download patient_data.csv",
        data=patient_df.to_csv(index=False),
        file_name=f"patient_data_{biomarker_key}.csv",
        mime="text/csv",
    )
