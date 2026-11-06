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
    python3 summary_sheets/build_summary_sheets.py 6 --no-packet   # skip the packet rebuild

Every build also rewrites summary_packet.pdf: a contents cover, then every registered chapter's
sheet in index.Rmd reading order, for printing and handing out as one packet. --check and --deploy
treat the packet as stale whenever any registered chapter is newer than it.

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
PACKET = "summary_packet.pdf"

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
        m = re.match(r"#\s+(.*?)\s*\{(?:#[\w-]+\s+)?-\}\s*$", line)  # {-} or {#id -}
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


def sheet_body(rmd, pandoc, suffix=": summary"):
    """The sheet's inner HTML (heading, subtitle, summary), or None if the chapter has no summary.
    The packet passes suffix="" because its cover already says these are the chapter summaries."""
    section = summary_section(rmd)
    if section is None:
        return None
    section = "\n".join(l for l in section.splitlines() if LINK_MARK not in l)
    if "include_graphics" not in section:
        # The chapter-opening card holds the decision flowchart (since 2026-10-29); the sheet still
        # leads with it.
        flow = re.search(r"```\{r[^\n]*\}\n[^\n]*images/[A-Za-z0-9._-]*flowchart[A-Za-z0-9._-]*[^\n]*\n```", rmd)
        if flow:
            section = flow.group(0) + "\n\n" + section
    md = chunks_to_images(section)
    frag = subprocess.run([pandoc, "--from", "markdown", "--to", "html"], input=md,
                          capture_output=True, text=True, check=True).stdout
    return (f"<h1>{html.escape(chapter_title(rmd))}{suffix}</h1>"
            f"<div class='sub'>Integrated Bioanalytical Chemistry &middot; thebustalab.github.io/integrated_bioanalytics</div>"
            f"{frag}")


def print_pdf(page_html, out):
    """Print a full HTML page to `out` with headless Chromium; return the PDF's page count."""
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        raise MissingTool("Playwright is not installed (pip install playwright)")
    with tempfile.TemporaryDirectory() as td:
        src = os.path.join(td, "sheet.html")
        open(src, "w", encoding="utf-8").write(page_html)
        with sync_playwright() as p:
            b = launch_chromium(p)
            pg = b.new_page()
            pg.goto("file://" + src)
            pg.pdf(path=out, format="Letter", print_background=True, prefer_css_page_size=True)
            b.close()
    return len(re.findall(rb"/Type\s*/Page[^s]", open(out, "rb").read()))


def html_page(title, body):
    return (f"<!doctype html><html><head><meta charset='utf-8'><title>{title}</title><style>{CSS}</style></head>"
            f"<body>{body}</body></html>")


def build(rmd_path, pandoc):
    rmd = open(rmd_path, encoding="utf-8").read()
    body = sheet_body(rmd, pandoc)
    if body is None:
        return None
    name = os.path.splitext(os.path.basename(rmd_path))[0]
    out = os.path.join(HERE, name + ".pdf")
    pages = print_pdf(html_page(html.escape(chapter_title(rmd)) + ": summary", body), out)
    if pages != 1:
        os.remove(out)
        sys.exit(f"ERROR: {name} summary runs to {pages} pages; trim the chapter's summary section.")
    return out


def reading_order():
    """Registered chapter files with a summary section, in the order index.Rmd reads them."""
    idx = open(os.path.join(BOOK, "index.Rmd"), encoding="utf-8").read()
    files = [os.path.join(BOOK, m) for m in re.findall(r'child="(chapters/[0-9][^"]*\.Rmd)"', idx)]
    return [f for f in files if os.path.exists(f) and HEADING in open(f, encoding="utf-8").read()]


PACKET_CSS = """
.sheet { break-after: page; }
.cover h1 { font-size: 24pt; margin: 1.4in 0 6pt 0; }
.cover .sub { font-size: 11pt; margin-bottom: 0.4in; }
.cover ol { font-size: 11pt; line-height: 1.9; padding-left: 1.2em; }
.cover li span { color: #666; }
"""


def build_packet(pandoc):
    """Every registered chapter's summary, in reading order, behind a cover with a contents list:
    one PDF to print and hand out. Each sheet is the same HTML as its own one-page PDF."""
    files = reading_order()
    items, sheets = [], []
    for i, f in enumerate(files):
        rmd = open(f, encoding="utf-8").read()
        sheets.append(f"<section class='sheet'>{sheet_body(rmd, pandoc, suffix='')}</section>")
        items.append(f"<li>{html.escape(chapter_title(rmd))} <span>&middot; page {i + 2}</span></li>")
    cover = ("<section class='sheet cover'><h1>Chapter summaries</h1>"
             "<div class='sub'>Integrated Bioanalytical Chemistry &middot; thebustalab.github.io/integrated_bioanalytics</div>"
             f"<ol>{''.join(items)}</ol></section>")
    out = os.path.join(HERE, PACKET)
    page = html_page("Chapter summaries", cover + "".join(sheets)).replace("</style>", PACKET_CSS + "</style>", 1)
    pages = print_pdf(page, out)
    if pages != len(files) + 1:
        os.remove(out)
        sys.exit(f"ERROR: the packet came out at {pages} pages, expected {len(files) + 1} "
                 "(cover plus one per chapter); a sheet is spilling over.")
    return out


def sheet_inputs(f):
    """The chapter plus every image its sheet prints: the summary's own include_graphics() images and the
    chapter-opening flowchart the sheet leads with. A redrawn flowchart with an unchanged chapter must
    still rebuild the sheet (test_build_summary_sheets.py)."""
    rmd = open(f, encoding="utf-8").read()
    names = set(re.findall(r"images/([A-Za-z0-9._-]+)", summary_section(rmd) or ""))
    flow = re.search(r"images/([A-Za-z0-9._-]*flowchart[A-Za-z0-9._-]*)", rmd)
    if flow:
        names.add(flow.group(1))
    return [f] + [p for p in (os.path.join(IMAGES, n) for n in sorted(names)) if os.path.exists(p)]


def newest_input(f):
    return max(os.path.getmtime(p) for p in sheet_inputs(f))


def packet_stale():
    out = os.path.join(HERE, PACKET)
    return not os.path.exists(out) or any(os.path.getmtime(out) < newest_input(f) for f in reading_order())


def pdf_for(f):
    return os.path.join(HERE, os.path.splitext(os.path.basename(f))[0] + ".pdf")


def is_stale(f):
    return not os.path.exists(pdf_for(f)) or os.path.getmtime(pdf_for(f)) < newest_input(f)


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
        if packet_stale():
            print("STALE:", PACKET)
        sys.exit(1 if stale or packet_stale() else 0)
    if "--deploy" in flags:
        files = [f for f in files if is_stale(f)]
        if not files and not packet_stale():
            print("summary sheets: all current")
            return
    try:
        pandoc = pandoc_bin()
        for f in files:
            print("wrote", os.path.relpath(build(f, pandoc), BOOK))
        # The packet follows every sheet build, so it can never lag the sheets it is made of.
        if "--no-packet" not in flags:
            print("wrote", os.path.relpath(build_packet(pandoc), BOOK))
    except MissingTool as e:
        if "--deploy" not in flags:
            sys.exit(f"ERROR: {e}")
        print(f"WARNING: summary sheets NOT rebuilt ({e}). Stale PDFs stay as they are:")
        for f in files:
            if is_stale(f):
                print("   STALE:", os.path.basename(f))


if __name__ == "__main__":
    main()
