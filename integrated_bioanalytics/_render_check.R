# _render_check.R — render SOME chapters locally and report every figure's real dimensions.
#
# The full `Rscript build.R` renders all 19 chapters and hits live APIs; that is the right tool for a
# release canary and the wrong one for "did my edit to chapter 7 still knit, and do the figures look
# the right shape". This renders only the chapters you name, in one session, with the chunk cache off,
# and then measures every PNG knitr produced.
#
#   Rscript _render_check.R 7                 # one chapter
#   Rscript _render_check.R 7 8 9 10          # several
#   Rscript _render_check.R --figs-only 7     # skip the render, just re-measure the last one
#
# Use the conda R, which has the whole stack:
#   ~/miniconda3/envs/phylochem/bin/Rscript _render_check.R 7
#
# Output: book_test/check.html plus a figure table on stdout. book_test/ is gitignored and separate
# from the published book/, so this can never overwrite a real render.
#
# It picks up _local_harness.R through _common.R, so it tests the LOCAL phylochemistry.R. That is the
# whole point: without the harness a local render silently tests the published library instead.

args <- commandArgs(trailingOnly = TRUE)
figs_only <- "--figs-only" %in% args
nums <- setdiff(args, "--figs-only")

here <- normalizePath(".")
outdir <- file.path(here, "book_test")
dir.create(outdir, showWarnings = FALSE)

## ── locate the chapter files from the numbers given ───────────────────────────────────────────────
all_ch <- list.files("chapters", pattern = "^[0-9]+_.*\\.Rmd$", full.names = FALSE)
pick <- function(n) {
  hit <- all_ch[grepl(paste0("^", n, "_"), all_ch)]
  if (length(hit) != 1) stop("no single chapter file for number ", n,
                             " (found: ", paste(hit, collapse = ", "), ")")
  hit
}

if (!figs_only) {
  if (!length(nums)) stop("give at least one chapter number, e.g. Rscript _render_check.R 7 8")
  files <- vapply(nums, pick, character(1))
  cat("rendering:", paste(files, collapse = ", "), "\n")

  ## conda pandoc 3.x — system pandoc 2.9 wants the removed pandoc-citeproc filter and fails assembly
  conda_pandoc <- "/home/bustalab/miniconda3/envs/phylochem/bin"
  if (dir.exists(conda_pandoc)) Sys.setenv(RSTUDIO_PANDOC = conda_pandoc)
  options(repos = c(CRAN = "https://cloud.r-project.org"))

  ## Build a merged doc from the chapters requested. Chunk labels repeat across chapters, and the
  ## cache is unsafe across a child merge (a hit replays output without re-materialising objects),
  ## so both are dealt with up front.
  rmd <- c(
    "---",
    "title: \"render check\"",
    "output:",
    "  html_document:",
    "    self_contained: false",
    "---",
    "",
    "```{r setup-check, include=FALSE}",
    "options(knitr.duplicate.label = 'allow')",
    "knitr::opts_chunk$set(echo = TRUE, prompt = FALSE, eval = TRUE,",
    "                      warning = FALSE, comment = '##', cache = FALSE,",
    "                      fig.width = 4, fig.height = 3,",
    "                      collapse = TRUE, results = 'markup', max.print = 6,",
    "                      out.width = '100%', fig.align = 'center')",
    "knitr::opts_hooks$set(cache = function(o) { o$cache <- FALSE; o })",
    "source('_common.R')",
    "```",
    ""
  )
  for (f in files) {
    rmd <- c(rmd,
             paste0("```{r child='chapters/", f, "'}"),
             "```",
             "")
  }
  writeLines(rmd, file.path(outdir, "check.Rmd"))

  ## Clear figures from any PREVIOUS run. Chunk labels are positional across the merged document, so
  ## rendering a different set of chapters renumbers them and leaves the old PNGs behind under new
  ## names -- which then show up in the table below as phantom duplicates of the right size. Caught
  ## doing exactly that on 2026-09-22.
  unlink(file.path(outdir, "check_files"), recursive = TRUE)
  unlink(file.path(outdir, "check_cache"), recursive = TRUE)

  t0 <- Sys.time()
  rmarkdown::render(file.path(outdir, "check.Rmd"),
                    knit_root_dir = here,
                    output_dir = outdir,
                    envir = new.env(),
                    quiet = FALSE)
  cat("\n===== RENDER_OK in", round(difftime(Sys.time(), t0, units = "secs")), "s =====\n")
}

## ── measure every figure knitr produced ───────────────────────────────────────────────────────────
figdir <- file.path(outdir, "check_files")
pngs <- list.files(figdir, pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
if (!length(pngs)) {
  cat("\nno figures found under", figdir, "\n")
} else {
  png_dim <- function(f) {
    con <- file(f, "rb"); on.exit(close(con))
    raw <- readBin(con, "raw", 33)
    if (length(raw) < 33) return(c(NA, NA))
    # IHDR width/height are big-endian 4-byte ints at offsets 16 and 20
    be <- function(b) sum(as.integer(b) * 256^((length(b) - 1):0))
    c(be(raw[17:20]), be(raw[21:24]))
  }
  info <- do.call(rbind, lapply(pngs, function(f) {
    d <- png_dim(f)
    data.frame(figure = basename(f), px_w = d[1], px_h = d[2],
               aspect = round(d[1] / d[2], 2),
               in_w = round(d[1] / 192, 2), in_h = round(d[2] / 192, 2),
               stringsAsFactors = FALSE)
  }))
  info <- info[order(info$figure), ]
  cat("\n===== FIGURES (", nrow(info), ") =====\n", sep = "")
  print(info, row.names = FALSE)
  cat("\nAspect ratio is width/height. Under ~0.7 is a tall figure, over ~2.5 a wide one -- check\n",
      "those render legibly at page width. in_w/in_h are back-computed at 192 dpi (knitr renders at\n",
      "2x here), so they should match the fig.width/fig.height set on the chunk. A mismatch means\n",
      "the chunk is falling back to the 4x3 global default -- i.e. you forgot to size it.\n", sep = "")
}
