#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# ra_core/reporting_helpers.py
import io
import pandas as pd
import xlsxwriter

# if you implemented the fancy sectioned writer earlier, import it here.
# otherwise we'll just dump the vertical dataframe to a single sheet.

try:
    from ra_core.reporting import _df_to_sections, _write_vertical_sheet  # optional
    HAVE_SECTIONS = True
except Exception:
    HAVE_SECTIONS = False

from ra_core.core import run_full_read_across_assessment

def create_excel_report_bytes(pairs, sheet_name_fmt="{i:02d} - {tname} vs {sname}") -> bytes:
    """
    Build an in-memory .xlsx for one or more comparisons. Returns bytes you can download in Streamlit.
    """
    bio = io.BytesIO()
    wb = xlsxwriter.Workbook(bio, {'in_memory': True})

    for i, (tname, tsmi, sname, ssmi) in enumerate(pairs, start=1):
        df = run_full_read_across_assessment(tname, tsmi, sname, ssmi)
        sheet_name = (sheet_name_fmt.format(i=i, tname=(tname or "Target")[:15], sname=(sname or "Surrogate")[:15]))[:31]

        if HAVE_SECTIONS:
            sections = _df_to_sections(df)
            _write_vertical_sheet(wb, sheet_name, sections)
        else:
            # simple fallback: dump the vertical DF to a sheet
            ws = wb.add_worksheet(sheet_name)
            # write header
            ws.write(0, 0, "Parameter"); ws.write(0, 1, "Target Result"); ws.write(0, 2, "Surrogate Result")
            for r, (idx, row) in enumerate(df.reset_index().iterrows(), start=1):
                ws.write(r, 0, str(row["Parameter"]))
                ws.write(r, 1, str(row.get("Target Result","")))
                ws.write(r, 2, str(row.get("Surrogate Result","")))

    wb.close()
    bio.seek(0)
    return bio.getvalue()

