import sys
from io import BytesIO
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _pipeline_loader import load_formatting  # noqa: E402
from _ui import inject_base_css  # noqa: E402

fmt = load_formatting()

st.set_page_config(page_title="1. Formatting", page_icon="\U0001f9ec", layout="wide")
inject_base_css()
st.title("1. Formatting")
st.write(
    "The raw TMA block map exported from the pathology database uses anonymised "
    "**YP-numbers** (e.g. `YP1234`) as sample identifiers in each core position. "
    "This stage looks up each YP-number's real patient number in a separate "
    "lookup file and substitutes it in place, leaving every other cell and all "
    "formatting untouched."
)

st.subheader("Inputs")
col1, col2 = st.columns(2)

with col1:
    st.markdown(f"**Lookup table** — `{fmt.DEFAULT_LOOKUP_FILE.relative_to(fmt.BASE_DIR.parent)}`")
    lookup_xl = pd.ExcelFile(fmt.DEFAULT_LOOKUP_FILE)
    lookup_sheet = st.selectbox("Sheet", lookup_xl.sheet_names, key="lookup_sheet")
    st.dataframe(lookup_xl.parse(lookup_sheet).head(10), width="stretch", hide_index=True)

with col2:
    st.markdown(f"**TMA block map** — `{fmt.DEFAULT_MAP_FILE.relative_to(fmt.BASE_DIR.parent)}`")
    map_xl = pd.ExcelFile(fmt.DEFAULT_MAP_FILE)
    map_sheet = st.selectbox("Sheet", map_xl.sheet_names, key="map_sheet_in")
    st.dataframe(
        map_xl.parse(map_sheet, header=None).head(10), width="stretch",
    )

st.subheader("Run")
if st.button("Run formatting step", type="primary"):
    with st.spinner("Replacing YP-numbers with patient numbers..."):
        st.session_state["formatting_result"] = fmt.run_formatting()

result = st.session_state.get("formatting_result")

if result is None:
    st.info("Click **Run formatting step** to process the files above.")
    st.stop()

st.subheader("Result")

st.markdown("**Lookup table load**")
for sheet, columns, rows, used in result.lookup_sheets:
    if used:
        st.write(f"- Sheet `{sheet}` — columns {columns} — {rows} rows")
    else:
        st.write(f"- Sheet `{sheet}` — :orange[skipped] (missing `pat_no`/`yp_num` columns; found {columns})")

m1, m2, m3 = st.columns(3)
m1.metric("YP → patient-number mappings loaded", result.mapping_count)
m2.metric("Cells replaced", result.replaced_count)
m3.metric("Unmatched YP-numbers", len(result.unmatched))

if result.unmatched:
    with st.expander(f"Unmatched YP-numbers ({len(result.unmatched)}) — left as-is in the output"):
        st.dataframe(
            pd.DataFrame(result.unmatched, columns=["Sheet", "Cell", "YP-number"]),
            width="stretch", hide_index=True,
        )
else:
    st.success("Every YP-number in the map had a matching patient number.")

st.markdown("**Output preview** — `tma_map_replaced.xlsx`")
out_sheet = st.selectbox("Sheet", result.workbook.sheetnames, key="map_sheet_out")
ws = result.workbook[out_sheet]
preview_rows = [[c.value for c in row] for row in ws.iter_rows(max_row=10)]
st.dataframe(pd.DataFrame(preview_rows), width="stretch")

buf = BytesIO()
result.workbook.save(buf)
st.download_button(
    "Download tma_map_replaced.xlsx",
    data=buf.getvalue(),
    file_name="tma_map_replaced.xlsx",
    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
)
st.caption(f"Also saved to `{result.output_file}` (same file the standalone script writes).")
