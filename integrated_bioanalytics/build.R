# Build entry point for the Integrated Bioanalytics bookdown.
#
# Running the from-scratch Linux harness or the core-scope validation? Read
# _agent_reference/build_harness.md first — method, conda env, and the traps that misreport success.
#
# WHY THIS EXISTS: _common.R monkeypatches bookdown's tweak_part_screwup() to
# guard against a missing "class" attribute. That patch must be in scope in the
# MAIN render_book() process (that is where the bs4_book HTML-tweak phase runs),
# but _common.R is only the `before_chapter_script`, so it is otherwise sourced
# only in the per-chapter knit sessions. Calling render_book() directly therefore
# crashes in the tweak phase with "missing value where TRUE/FALSE needed".
# Sourcing _common.R here first puts the patch in the main process.
# Full failure mode + rationale: AGENTS.md -> "Build process".
#
# Run from this directory:  Rscript build.R   (or:  source("build.R") in an R session)

source("_common.R")

# PRUNE THE CACHE'S PACKAGE LIST BEFORE RENDERING.
# knitr keeps a `__packages` file inside the chunk cache dir listing every package that was attached
# while the cache was written, and re-attaches ALL of them at the first cache hit of the next render
# — before any chapter runs. So one `library(ggtern)` in an appendix chapter contaminates the WHOLE
# book on every later build, and keeps doing so after that line is removed, because the list is only
# ever added to. ggtern replaces ggplot2::ggplot(), which is how chapter 7 came to fail inside
# `ggtern:::ggplot.default` on 2026-09-22. Drop the known maskers, and drop anything not installed
# on this machine (a stale entry from another box otherwise halts the render outright). Packages
# genuinely attached during the render get recorded again.
local({
  maskers <- c("ggtern", "plyr")
  files <- list.files("_bookdown_files", pattern = "^__packages$",
                      recursive = TRUE, full.names = TRUE)
  for (f in files) {
    pkgs <- readLines(f, warn = FALSE)
    keep <- setdiff(pkgs, maskers)
    keep <- keep[nzchar(keep) & nzchar(vapply(keep, function(p) system.file(package = p), character(1)))]
    if (!identical(keep, pkgs)) {
      writeLines(keep, f)
      message("[build] pruned from ", f, ": ", paste(setdiff(pkgs, keep), collapse = ", "))
    }
  }
})

# STALE-CACHE GUARD (2026-09-22). The chunk cache is ON (`cache = TRUE` in index.Rmd's setup chunk)
# and is deliberately NOT keyed on the library's checksum — see AGENTS.md for why (a cold rebuild
# re-runs the live ch13/ch14 API chunks and can fail a deploy on a spent class budget). The cost of
# that trade-off is this failure mode, and it burned most of an evening on 2026-09-22:
#
#   The library gains a feature. A chapter is written against it. The chunk that CALLS the new feature
#   is a cache HIT, so it replays the OLD object; the chunk that PLOTS it errored last time, was never
#   cached, and re-runs. The render then fails identically on every attempt, no matter what you fix,
#   because the fixed code is never executed. Nothing in the output says "cache".
#
# So: if the library is newer than the newest cache entry, stop and say so. One line to clear, and the
# message names it. Override with BOOK_ALLOW_STALE_CACHE=1 when you know the cache predates a
# library change that cannot affect it.
# Both `phylochemistry.R` AND `_common.R` count as inputs: the 2026-09-22 cache was poisoned not by
# a library edit but by a render in which ch3 re-sourced the PUBLISHED library over the local one, so
# the cache entry was NEWER than the library and a library-only mtime check would have waved it
# through. `_common.R` is where the loading logic lives, so touching it invalidates the cache too.
#
# SELF-DISABLING (2026-09-23). `_common.R` now keys each chunk's cache on the fingerprints of the
# library functions THAT chunk calls (transitive closure), so a library edit invalidates exactly the
# affected chunks and leaves the rest cached. Where that is live, this blunt mtime guard would only
# force needless cold rebuilds -- so it stands down. If the hook ever fails to install,
# `.pc_fingerprint` is absent and the blunt guard comes back automatically. Fail safe, not silent.
local({
  fp <- get0(".pc_fingerprint", envir = globalenv())
  if (!is.null(fp) && length(fp) > 0) {
    message("[build] per-chunk library fingerprints active (", length(fp),
            " functions) — blunt stale-cache guard stood down")
    return(invisible(NULL))
  }
  message("[build] WARNING: per-chunk library fingerprints NOT active — falling back to the ",
          "mtime guard. Check _common.R.")
  lib <- c(file.path("..", "phylochemistry", "phylochemistry.R"), "_common.R")
  lib <- lib[file.exists(lib)]
  cache_files <- list.files("_bookdown_files", recursive = TRUE, full.names = TRUE,
                            pattern = "\\.(rdb|rdx|RData)$")
  if (length(lib) == 0 || length(cache_files) == 0) return(invisible(NULL))
  if (identical(Sys.getenv("BOOK_ALLOW_STALE_CACHE"), "1")) {
    message("[build] stale-cache guard overridden by BOOK_ALLOW_STALE_CACHE=1")
    return(invisible(NULL))
  }
  lib_time   <- max(file.info(lib)$mtime)
  newest     <- lib[which.max(file.info(lib)$mtime)]
  cache_time <- max(file.info(cache_files)$mtime)
  if (lib_time > cache_time) {
    stop(
      "\n\n  STALE CHUNK CACHE.\n",
      "  ", basename(newest), " was modified ", format(lib_time),   "\n",
      "  newest cache entry is from    ", format(cache_time), "\n\n",
      "  Cached chunks will replay objects built by the OLD library, so any chapter written\n",
      "  against a library change cannot compile and will fail the same way every time.\n\n",
      "  FIX (from anywhere — this is what --cold-book is for):\n\n",
      "      ", normalizePath("..", mustWork = FALSE), "/deploy.sh --cold-book\n\n",
      "  It clears the cache and compiles in one step, so there is no directory to get wrong.\n",
      "  Doing it by hand needs the BOOK dir, not the site root — `rm -rf _bookdown_files` from the\n",
      "  site root silently removes nothing and you land right back here:\n\n",
      "      rm -rf ", normalizePath(".", mustWork = FALSE), "/_bookdown_files\n\n",
      "  (A cold render re-runs the live API chunks in ch13/ch14 — expected.)\n",
      "  To clear only from chapter 7 onward and keep the early chapters cached:\n\n",
      "      for n in $(seq 237 634); do rm -f ", normalizePath(".", mustWork = FALSE),
      "/_bookdown_files/index_cache/html/unnamed-chunk-${n}_*; done\n\n",
      "  Or set BOOK_ALLOW_STALE_CACHE=1 if you are certain the library change cannot reach the cache.\n",
      call. = FALSE
    )
  }
})

bookdown::render_book("index.Rmd")
