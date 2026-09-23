# example R options set globally
options(width = 60)

# Ensure a CRAN mirror is always set. Each chapter renders in its own
# non-interactive session (new_session: true), where repos defaults to
# "@CRAN@". phylochemistry.R calls install.packages() for any missing
# package, which errors ("trying to use CRAN without setting a mirror")
# unless a real mirror is set here first.
if (is.null(getOption("repos")) || any(getOption("repos") == "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# Bootstrap phylochemistry (packages, custom functions, datasets) at the start of
# every chapter session. This used to happen as a side effect of reading order:
# the installation chapter's live source() chunk sat in "GETTING STARTED", ahead
# of every data chapter. That chapter moved to the appendix (2026-08-27), so the
# bootstrap is explicit here instead of implicit in chapter position. The guard
# keeps it a no-op when the library is already attached.
# LOCAL BUILD HARNESS (optional, gitignored, absent on the Mac). When present it sources the LOCAL
# phylochemistry instead of the published URL and translates the chapters' hardcoded Mac paths onto
# this machine, so a local render actually tests the working copy. It bootstraps the library itself,
# hence the exists() guard below still short-circuits afterwards. No-op if the file is not there.
if (file.exists("_local_harness.R")) try(source("_local_harness.R"), silent = TRUE)

# WHICH phylochemistry.R? The LOCAL one in this same working tree, whenever it is there; the
# published URL only as a fallback. This matters because the library and the chapters are edited
# together and deploy.sh compiles the book BEFORE it pushes the site: sourcing the URL means a
# compile always tests the LIVE (i.e. previous) library, so any chapter written against a
# same-commit library change cannot compile. That is exactly how ch7 broke on 2026-09-22 —
# `runMatrixAnalysis(analysis = "dist", output_format = "long")` returned a bare `dist` from the
# published library and ggplot() could not fortify it. Sourcing the tree removes the whole class.
.pc_local <- file.path("..", "phylochemistry", "phylochemistry.R")

# NOT folded into cache.extra on purpose. Keying the chunk cache on this file's checksum would be
# more correct, but it makes every library edit a COLD rebuild of the whole book, which re-runs the
# deliberately-live API chunks in ch13/ch14 and can fail a deploy on a spent class budget. If a
# library change needs to show up in already-cached chapters, force it by hand:
#   rm -rf _bookdown_files      (then Rscript build.R)
# The Linux canary harness does set cache.extra, because there a full cold render is the point.

# ANY source() OF A phylochemistry URL IS REDIRECTED TO THE LOCAL WORKING COPY (2026-09-22).
#
# WHY. _common.R loading the local library is not enough, because CHAPTERS re-source the published
# one mid-render and the whole book knits in ONE session (index.Rmd + children). ch3 line ~327 is an
# `include = FALSE` chunk -- hidden output, but it still EVALUATES -- holding a bare
# source("https://thebustalab.github.io/phylochemistry/phylochemistry.R"). ch3 runs before ch7, so it
# quietly replaced the freshly-loaded working copy with the PREVIOUS RELEASE, and ch7 then called a
# runMatrixAnalysis() with no long-format `dist` branch and died inside ggplot(), ~12000 lines from
# the cause. ch2 (appendix) and ch13/ch14 (the module) do the same thing later in the book.
# This cost the 2026-09-22 deploy several hours. The Linux canary harness had solved it with exactly
# this override -- but _local_harness.R is gitignored and Mac-absent, so the Mac deploy, the one that
# actually publishes, was the only machine without the fix. Hence it lives HERE now.
#
# Applies to modules/ and bustalabfunctions/ too, so a local library can no longer be half-shadowed
# by published modules. No-op wherever the local file is absent. Render-only: _common.R is the book's
# before_chapter_script and is not part of what students source.
.pc_root <- file.path("..", "phylochemistry")
if (dir.exists(.pc_root) && !isTRUE(attr(get0("source", envir = globalenv()), "pc_redirect"))) {
  .pc_url <- "^https://thebustalab\\.github\\.io/phylochemistry/"
  .pc_redirect <- function(file, ...) {
    if (is.character(file) && length(file) == 1L && grepl(.pc_url, file)) {
      .local <- file.path(.pc_root, sub(.pc_url, "", file))
      if (file.exists(.local)) file <- .local
    }
    base::source(file, ...)
  }
  attr(.pc_redirect, "pc_redirect") <- TRUE
  assign("source", .pc_redirect, envir = globalenv())
}

# A HALF-LOADED LIBRARY MUST NOT PASS SILENTLY (2026-09-22).
# The local source() used to be try(..., silent = TRUE) with `exists("algae_data")` as the
# did-it-work test. Both are wrong together: phylochemistry.R defines algae_data at ~line 311 (the
# datasets module) but runMatrixAnalysis at ~line 14104, so ANY error in between leaves algae_data
# present, the fallback skipped, and the library half-loaded -- with the modules/ sourced near the
# top still in scope. The render then dies thousands of lines from the cause, in a chapter chunk,
# with an error that says nothing about the library. Fail loudly at the point of failure instead.
if (!exists("algae_data")) {
  bustalab <- TRUE
  if (file.exists(.pc_local)) {
    .pc_err <- tryCatch({ suppressMessages(source(.pc_local)); NULL },
                        error = function(e) conditionMessage(e))
    if (!is.null(.pc_err)) {
      stop(
        "phylochemistry.R failed to load from ", normalizePath(.pc_local, mustWork = FALSE),
        "\n  ", .pc_err,
        "\nThe book compiles against the LOCAL library on purpose (see AGENTS.md). Fix the library;",
        " do NOT fall back to the published one, which is a different version."
      )
    }
  } else {
    # No local working copy (a machine that only has the book). Published library, as before.
    try(
      suppressMessages(
        source("https://thebustalab.github.io/phylochemistry/phylochemistry.R")
      ),
      silent = TRUE
    )
  }
}

# Post-load assertion. Cheap, and it names the real problem in one line rather than letting a stale
# or partial library surface as a ggplot()/fortify() error twelve thousand lines away.
if (!exists("runMatrixAnalysis") || !is.function(runMatrixAnalysis)) {
  stop("phylochemistry loaded but runMatrixAnalysis() is missing -- the library did not finish loading.")
}
# Assert on BEHAVIOUR, not on the formals. Every stale fork of runMatrixAnalysis() still HAS an
# `output_format` argument -- it just ignores it for analysis = "dist" and hands back a bare `dist`
# object, which is precisely the ch7 failure. So look for the branch that actually builds the long
# pair table (2026-09-22: a formals-only check passed happily while the wrong function was loaded).
if (!any(grepl('"sample_1", "sample_2", "distance"',
               deparse(body(runMatrixAnalysis)), fixed = TRUE))) {
  stop(
    "The runMatrixAnalysis() in scope is a STALE copy: it has no long-format branch for\n",
    "  analysis = \"dist\", so it returns a bare `dist` object and chapter 7 cannot plot it.\n",
    "The working copy at ", normalizePath(.pc_local, mustWork = FALSE), " does have that branch,\n",
    "so what loaded is either the PUBLISHED library or a fork left inside modules/ by a\n",
    "phylochemistry.R load that failed part-way through."
  )
}

# example chunk options set globally
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE
  )

try({
  if (requireNamespace("bookdown", quietly = TRUE)) {
    ns <- asNamespace("bookdown")
    if (exists("tweak_part_screwup", envir = ns, inherits = FALSE)) {
      unlockBinding("tweak_part_screwup", ns)
      assign(
        "tweak_part_screwup",
        function(html) {
          sidebar <- xml2::xml_find_first(html, "//div[contains(@class, 'sidebar-chapter')]")
          if (inherits(sidebar, "xml_missing")) return()
          parent <- xml2::xml_parent(sidebar)
          if (inherits(parent, "xml_missing") || is.null(parent)) return()
          parent_class <- xml2::xml_attr(parent, "class")
          if (!is.na(parent_class) && parent_class == "row") return()
          main <- xml2::xml_find_first(html, "//main")
          if (inherits(main, "xml_missing") || is.null(main)) return()
          xml2::xml_add_sibling(main, sidebar)
          xml2::xml_remove(sidebar)
        },
        envir = ns
      )
      lockBinding("tweak_part_screwup", ns)
    }
  }
}, silent = TRUE)

# BUILD STAMP (temporary, 2026-09-22). Printed to stderr so it lands in the deploy log. Its only job
# is to answer "is the render actually using THIS file, on THIS machine?" -- a question that cost
# several rounds of fixing things that were never being executed. Remove once ch7 compiles.
cat(
  "[_common.R] stamp 2026-09-22-B | file: ", normalizePath("_common.R", mustWork = FALSE), "\n",
  "[_common.R] local library: ", normalizePath(.pc_local, mustWork = FALSE),
  " exists=", file.exists(.pc_local), "\n",
  "[_common.R] source() redirect installed: ",
  isTRUE(attr(get0("source", envir = globalenv()), "pc_redirect")), "\n",
  "[_common.R] runMatrixAnalysis has long-dist branch: ",
  any(grepl('"sample_1", "sample_2", "distance"',
            deparse(body(runMatrixAnalysis)), fixed = TRUE)), "\n",
  sep = "", file = stderr()
)
