#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# app.py
import io
import textwrap
import pandas as pd
import streamlit as st

# --- import your pipeline entrypoint from core.py ---
from ra_core.core import run_full_read_across_assessment

# Optional: if you already built a fancy Excel writer elsewhere, we try to import it.
# If not found, we fallback to a simple "dump the vertical DF" exporter below.
try:
    from ra_core.reporting_helpers import create_excel_report_bytes as _build_excel_bytes
except Exception:
    _build_excel_bytes = None


st.set_page_config(page_title="Read-Across Similarity", layout="wide")
st.title("Read-Across Assessment")

st.markdown(
    "Enter a **Target** and **Surrogate** (name + SMILES), run all modules, "
    "and download a formatted Excel report."
)

# -------- helpers --------

def _wrap_vertical_df(df: pd.DataFrame) -> pd.DataFrame:
    """Soft-wrap long text fields so the on-screen table is readable."""
    if "Target Result" in df.columns and "Surrogate Result" in df.columns:
        df = df.copy()
        for col in ["Target Result", "Surrogate Result"]:
            df[col] = df[col].apply(
                lambda v: v if isinstance(v, (int, float)) else textwrap.fill(str(v), width=80)
            )
    return df

@st.cache_data(show_spinner=False)
def _run_pair(tname, tsmi, sname, ssmi) -> pd.DataFrame:
    """Cached wrapper around your pipeline."""
    return run_full_read_across_assessment(tname, tsmi, sname, ssmi)

def _simple_excel_bytes(pairs):
    """
    Fallback Excel exporter: each comparison on its own sheet,
    writing the vertical DataFrame as-is.
    """
    bio = io.BytesIO()
    with pd.ExcelWriter(bio, engine="xlsxwriter") as writer:
        for i, (tname, tsmi, sname, ssmi) in enumerate(pairs, start=1):
            df = _run_pair(tname, tsmi, sname, ssmi).reset_index()
            sheet_name = f"{i:02d} - {tname or 'Target'} vs {sname or 'Surrogate'}"
            writer.book.add_worksheet(sheet_name[:31])  # ensure sheet exists and name <=31
            # Re-open the sheet by name to write the DF
            df.to_excel(writer, sheet_name=sheet_name[:31], index=False)
    bio.seek(0)
    return bio.getvalue()

def _build_excel(pairs):
    """Use your fancy exporter if present; otherwise fallback to the simple one."""
    if _build_excel_bytes is not None:
        return _build_excel_bytes(pairs)
    return _simple_excel_bytes(pairs)

# -------- UI --------

with st.form("inputs"):
    c1, c2 = st.columns(2)
    with c1:
        target_name = st.text_input("Target name", value="")
        target_smiles = st.text_input("Target SMILES", value="")
    with c2:
        surrogate_name = st.text_input("Surrogate name", value="")
        surrogate_smiles = st.text_input("Surrogate SMILES", value="")
    submitted = st.form_submit_button("Run assessment")

if submitted:
    if not target_smiles or not surrogate_smiles:
        st.error("Please provide both Target and Surrogate SMILES.")
    else:
        with st.spinner("Running…"):
            df = _run_pair(target_name, target_smiles, surrogate_name, surrogate_smiles)

        st.success("Done.")
        st.subheader("Report (vertical)")
        st.dataframe(_wrap_vertical_df(df), use_container_width=True)

        # Quick metrics row (if your labels match)
        try:
            total = float(df.loc["TOTAL SCORE (Sum of Modules)", "Target Result"])
            pc    = float(df.loc["P-Chem Module Score", "Target Result"])
            ms    = float(df.loc["Metabolic Similarity Module Score", "Target Result"])
            sa    = float(df.loc["Structural Alert Module Score", "Target Result"])
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Total", f"{total:.3f}")
            c2.metric("Phys Chem", f"{pc:.3f}")
            c3.metric("Metabolism", f"{ms:.3f}")
            c4.metric("Structural Alerts", f"{sa:.3f}")
        except Exception:
            # If the labels differ or a module is missing, just skip metrics.
            pass

        # Download Excel (this single pair, one sheet)
        pairs = [(target_name, target_smiles, surrogate_name, surrogate_smiles)]
        xls_bytes = _build_excel(pairs)
        st.subheader("Export")
        st.download_button(
            "Download Excel report",
            data=xls_bytes,
            file_name="read_across_report.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )

        with st.expander("Debug: raw DataFrame"):
            st.write(df)

