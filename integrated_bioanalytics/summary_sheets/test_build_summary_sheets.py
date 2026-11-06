#!/usr/bin/env python3
"""Regression: a sheet is stale when an image it prints is newer than the PDF, not only the chapter.

2026-11-06: the comparing-groups flowchart was redrawn (report means / medians) with the chapter
.Rmd untouched, so --check and --deploy called 11_comparing_means.pdf current and the summary sheet
kept the old flowchart. Staleness now covers the chapter plus every image the sheet prints.
Run: python3 summary_sheets/test_build_summary_sheets.py
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build_summary_sheets as b  # noqa: E402

CHAPTER = """# comparing groups {-}

```{r echo=FALSE}
knitr:::include_graphics('https://example.org/integrated_bioanalytics/images/99_demo_flowchart.png')
```

## chapter summary {-}

Some summary text.

```{r echo=FALSE}
knitr:::include_graphics('images/99_demo_table.png')
```

## exercises {-}
"""

fails = []


def check(cond, msg):
    print(("  ok    " if cond else "FAIL    ") + msg)
    if not cond:
        fails.append(msg)


with tempfile.TemporaryDirectory() as td:
    images = os.path.join(td, "images")
    os.makedirs(images)
    b.HERE, b.IMAGES = td, images
    rmd = os.path.join(td, "99_demo.Rmd")
    open(rmd, "w").write(CHAPTER)
    flow, table = os.path.join(images, "99_demo_flowchart.png"), os.path.join(images, "99_demo_table.png")
    for p in (flow, table):
        open(p, "wb").write(b"png")
    pdf = b.pdf_for(rmd)
    open(pdf, "wb").write(b"pdf")
    t = 1_700_000_000
    for p in (rmd, flow, table):
        os.utime(p, (t, t))
    os.utime(pdf, (t + 10, t + 10))

    check(set(b.sheet_inputs(rmd)) == {rmd, flow, table}, "inputs = chapter + flowchart + summary image")
    check(not b.is_stale(rmd), "PDF newer than every input: current")
    os.utime(flow, (t + 20, t + 20))
    check(b.is_stale(rmd), "flowchart redrawn after the PDF: stale")
    os.utime(flow, (t, t))
    os.utime(table, (t + 20, t + 20))
    check(b.is_stale(rmd), "summary image redrawn after the PDF: stale")

print()
if fails:
    print(f"{len(fails)} FAILURE(S)")
    sys.exit(1)
print("ALL CHECKS PASSED")
