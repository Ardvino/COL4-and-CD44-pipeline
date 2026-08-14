"""Shared page chrome — the small bits every page in this app injects, kept
in one place so they don't drift between pages.
"""
import streamlit as st


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
