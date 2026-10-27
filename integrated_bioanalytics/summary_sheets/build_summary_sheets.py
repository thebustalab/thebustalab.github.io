#!/usr/bin/env python3
"""Build the one-page summary PDF for each chapter from its `## chapter summary {-}` section.

The chapter is the single source of truth: this script lifts that section out of the .Rmd,
turns knitr include_graphics() chunks into plain images from ../images/, converts the markdown
with pandoc, and prints it to a one-page US-letter PDF with headless Chromium (Playwright).
Lines containing the download link itself are dropped, so the PDF does not link to itself.

    python3 summary_sheets/build_summary_sheets.py            # every chapter with a summary section
    python3 summary_sheets/build_summary_sheets.py 6b         # chapters whose file starts with 6b_
    python3 summary_sheets/build_summary_sheets.py --check    # exit 1 if any PDF is older than its chapter
    python3 summary_sheets/build_summary_sheets.py --deploy   # rebuild only stale sheets; if Playwright,
                                                              # Chromium or pandoc is missing, WARN and exit 0

--deploy is what deploy.sh runs (step 1c). A missing tool must not block a publish, so it warns
and leaves the stale PDFs in place; a sheet that overflows one page is an authoring error and
still fails.

Fails loudly if a sheet spills onto a second page: the point is one page.
"""
import glob
import html
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BOOK = os.path.dirname(HERE)
CHAPTERS = os.path.join(BOOK, "chapters")
IMAGES = os.path.join(BOOK, "images")
HEADING = "## chapter summary {-}"
LINK_MARK = "summary_sheets/"

CSS = """
@page { size: letter; margin: 0.45in; }
body { font-family: Helvetica, Arial, sans-serif; font-size: 9.5pt; line-height: 1.3; color: #222; }
h1 { font-size: 16pt; margin: 0 0 2pt 0; }
.sub { color: #666; font-size: 8.5pt; margin-bottom: 6pt; }
p { margin: 4pt 0; }
img { display: block; margin: 6pt auto; max-width: 100%; max-height: 4.3in; }
table { border-collapse: collapse; width: 100%; font-size: 8.3pt; margin-top: 6pt; }
th, td { border: 1px solid #bbb; padding: 2.5pt 4pt; text-align: left; vertical-align: top; }
th { background: #fee08b; }
code { font-family: Menlo, Consolas, monospace; font-size: 8pt; }
.cap { display: none; }
"""


class MissingTool(Exception):
    pass


def pandoc_bin():
    # Linux conda env, PATH, then the copies RStudio and Quarto bundle on a Mac.
    candidates = [os.path.expanduser("~/miniconda3/envs/phylochem/bin/pandoc"), shutil.which("pandoc")]
    candidates += glob.glob("/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/*/pandoc")
    candidates += glob.glob("/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/pandoc")
    candidates += glob.glob("/Applications/quarto/bin/tools/*/pandoc")
    for c in candidates:
        if c and os.path.exists(c):
            return c
    raise MissingTool("pandoc not found")


def launch_chromium(p):
    try:
        return p.chromium.launch()
    except Exception:
        # The pip playwright may be newer than its downloaded browsers; any already-installed
        # Chromium prints a page identically.
        found = sorted(glob.glob(os.path.expanduser("~/.cache/ms-playwright/chromium-*/chrome-linux*/chrome")))
        found += sorted(glob.glob(os.path.expanduser(
            "~/Library/Caches/ms-playwright/chromium-*/chrome-mac*/Chromium.app/Contents/MacOS/Chromium")))
        if not found:
            raise MissingTool("no Chromium for Playwright (run: python3 -m playwright install chromium)")
        return p.chromium.launch(executable_path=found[-1])


def chapter_title(rmd):
    for line in rmd.splitlines():
        m = re.match(r"#\s+(.*?)\s*\{-\}\s*$", line)
        if m and not line.startswith("##"):
            return m.group(1)
    return "chapter summary"


def summary_section(rmd):
    i = rmd.find(HEADING)
    if i < 0:
        return None
    body = rmd[i + len(HEADING):]
    m = re.search(r"^## ", body, re.M)
    return body[: m.start()] if m else body


def chunks_to_images(md):
    def repl(m):
        chunk = m.group(0)
        img = re.search(r"images/([A-Za-z0-9._-]+)", chunk)
        if not img:
            sys.exit("ERROR: an R chunk in the summary section is not an include_graphics() of images/; "
                     "the sheet builder cannot run R. Keep the summary to prose, tables and images.")
        path = os.path.join(IMAGES, img.group(1))
        if not os.path.exists(path):
            sys.exit(f"ERROR: summary image missing: {path}")
        return f"![]({path})"
    return re.sub(r"```\{r[^\n]*\}\n.*?\n```", repl, md, flags=re.S)


def build(rmd_path, pandoc):
    rmd = open(rmd_path, encoding="utf-8").read()
    section = summary_section(rmd)
    if section is None:
        return None
    section = "\n".join(l for l in section.splitlines() if LINK_MARK not in l)
    md = chunks_to_images(section)
    frag = subprocess.run([pandoc, "--from", "markdown", "--to", "html"], input=md,
                          capture_output=True, text=True, check=True).stdout
    name = os.path.splitext(os.path.basename(rmd_path))[0]
    title = html.escape(chapter_title(rmd)) + ": summary"
    page = (f"<!doctype html><html><head><meta charset='utf-8'><title>{title}</title><style>{CSS}</style></head><body>"
            f"<h1>{html.escape(chapter_title(rmd))}: summary</h1>"
            f"<div class='sub'>Integrated Bioanalytical Chemistry &middot; thebustalab.github.io/integrated_bioanalytics</div>"
            f"{frag}</body></html>")
    out = os.path.join(HERE, name + ".pdf")
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        raise MissingTool("Playwright is not installed (pip install playwright)")
    with tempfile.TemporaryDirectory() as td:
        src = os.path.join(td, "sheet.html")
        open(src, "w", encoding="utf-8").write(page)
        with sync_playwright() as p:
            b = launch_chromium(p)
            pg = b.new_page()
            pg.goto("file://" + src)
            pg.pdf(path=out, format="Letter", print_background=True, prefer_css_page_size=True)
            b.close()
    pages = len(re.findall(rb"/Type\s*/Page[^s]", open(out, "rb").read()))
    if pages != 1:
        os.remove(out)
        sys.exit(f"ERROR: {name} summary runs to {pages} pages; trim the chapter's summary section.")
    return out


def pdf_for(f):
    return os.path.join(HERE, os.path.splitext(os.path.basename(f))[0] + ".pdf")


def is_stale(f):
    return not os.path.exists(pdf_for(f)) or os.path.getmtime(pdf_for(f)) < os.path.getmtime(f)


def main():
    flags = {a for a in sys.argv[1:] if a.startswith("--")}
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    files = sorted(glob.glob(os.path.join(CHAPTERS, "[0-9]*.Rmd")))
    if args:
        files = [f for f in files if any(os.path.basename(f).startswith(a + "_") for a in args)]
    files = [f for f in files if HEADING in open(f, encoding="utf-8").read()]
    if "--check" in flags:
        stale = [f for f in files if is_stale(f)]
        for f in stale:
            print("STALE:", os.path.basename(f))
        sys.exit(1 if stale else 0)
    if "--deploy" in flags:
        files = [f for f in files if is_stale(f)]
        if not files:
            print("summary sheets: all current")
            return
    try:
        pandoc = pandoc_bin()
        for f in files:
            print("wrote", os.path.relpath(build(f, pandoc), BOOK))
    except MissingTool as e:
        if "--deploy" not in flags:
            sys.exit(f"ERROR: {e}")
        print(f"WARNING: summary sheets NOT rebuilt ({e}). Stale PDFs stay as they are:")
        for f in files:
            if is_stale(f):
                print("   STALE:", os.path.basename(f))


if __name__ == "__main__":
    main()
