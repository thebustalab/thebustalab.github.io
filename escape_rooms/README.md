# WebR Escape Rooms

Narrative, escape-room-style R exercises that run entirely in the student's
browser (no server, no install) via [WebR](https://docs.r-wasm.org/webr/latest/).
A student works through a two-screen mystery, does real R analysis in an
embedded console, answers a short series of multiple-choice questions with
feedback, and receives an alphanumeric **submission code** to paste into Canvas.
A small R decoder turns the codes back into per-student answer paths and points.

Prototype scenario: **`alaska/`** — "Signal in the Cold" (the Noatak lake
mystery, from the CHEM 5725 data-visualization exercise).

## Layout

```
escape_rooms/
  shared/
    webr-console.js    WebR boot + package install + data load + run
    escape-engine.js   two-screen flow, MC steps, feedback, code generation
    codec.js           submission-code encoder (mirrors the R decoder)
  alaska/
    index.html         the room shell
    style.css
    scenario.js        ALL the Alaska content, as data
  decoder/
    decode_codes.R     decode + grade submitted codes (run on your machine)
  README.md
```

To add a mystery, copy `alaska/` to a new folder, give the scenario a unique
`id` (1–15) in its `scenario.js`, and rewrite the story, datasets, packages,
and steps. The engine and codec are shared and unchanged.

## Running / hosting

WebR needs to be served over http/https (ES modules + `fetch` don't work from a
`file://` path). To try it locally from the site root:

```
cd websites/thebustalab.github.io
python3 -m http.server 8000
# then open http://localhost:8000/escape_rooms/alaska/
```

First load takes ~20–40 s while WebR downloads and installs ggplot2/dplyr/readr;
after that it's cached. On the live site it's just static files under
`thebustalab.github.io/escape_rooms/alaska/` — link to it from the book.

## The submission code

`codec.js` (browser) and `decode_codes.R` (your machine) implement the **same**
scheme; if you change one, change the other and re-run the R self-test
(`Rscript decode_codes.R`).

What a code encodes, per run:

- a header byte: version + scenario id
- one byte per step: chosen answer index (5 bits) + attempts taken (3 bits)
- a checksum byte

The payload is XOR-scrambled with a keystream derived from a shared secret
**and the student's x500**, then Crockford-base32 encoded. Consequences:

- Two students with the same answer path get different-looking codes.
- A code decoded with the **wrong** x500 fails its checksum — so a shared code
  is detectable (it won't validate against the friend's identity).
- Typos are caught by the checksum.

Grading, after downloading a Canvas assignment export with x500 + code columns:

```r
source("decode_codes.R")
roster <- readr::read_csv("canvas_export.csv")
scored <- grade_submissions(roster, ALASKA_KEY, id_col = "x500", code_col = "code")
readr::write_csv(scored, "graded.csv")
```

`ALASKA_KEY$score_step` is where you set points-per-step and the attempts
penalty — edit it to taste.

### Honest security note

The secret lives in `escape-engine.js`, which ships to the browser. This is
**obfuscation, not security**: a determined student who reads the JavaScript
could forge a code, and anyone can read an AI's answer on a second screen and
type it in. It's a speed bump appropriate for low-stakes practice, not an exam
control. The paste-warning banner is likewise a nudge, not a block. The real
defence is that answers are the *product of running the analysis*, and the code
records attempts — so shortcutting still shows.

Change the secret per course by editing `SECRET` in **both**
`escape-engine.js` and `decode_codes.R`.

## WebR compatibility (what runs, what doesn't)

Verified against the WebR wasm package repo for the CHEM 5725 techniques:

- **Runs cleanly:** data visualization, data wrangling, PCA, flat clustering,
  comparing means, modeling (tidymodels + ranger/randomForest). All required
  packages have wasm builds.
- **Needs a tweak:** hierarchical clustering draws its tree with **ggtree**,
  which has no wasm build. The `hclust` maths is base R; only the drawing needs
  swapping to `ape`/`ggdendro`, **or** building ggtree+treeio for wasm (both are
  pure-R Bioconductor packages, so this is feasible — see project notes).
- **Does not fit this format:** the embeddings exercises call live APIs with
  secret keys (`searchPubMed`, `embedText`) and can't run safely client-side.

You cannot `source()` the full `phylochemistry.R` in WebR (it loads shiny/httpuv
and a Bioconductor stack). Each scenario instead loads only the specific
packages and datasets it needs, and any custom phylochemistry functions should
be extracted into a slim per-scenario bundle.

## Status

The Alaska prototype and the code loop are built; the R codec is verified
(encode/decode/grade round-trip, wrong-id rejection, agreement with a
JS-faithful reference port for both short and long codes). The **live WebR run**
— booting R and rendering the ggplot in a real browser — has not yet been
executed and is the next thing to confirm by opening the page.
