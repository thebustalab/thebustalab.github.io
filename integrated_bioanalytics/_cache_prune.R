# Prune knitr's `__packages` lists before a render. Sourced by build.R; tested by
# tests/cache_prune_test.R.
#
# knitr re-attaches every package in `__packages` at the first cache hit, before chapter 1 runs, and
# the list is only ever added to. A package nothing in the book attaches any more can therefore still
# mask a function the book uses: after the 2026-11 package switch the list still carried the old
# phylochemistry.R's Hmisc, whose summarize() hid dplyr's and broke ch6. So keep only what the book's
# own sources attach (plus what those packages attach in turn), minus the known maskers and anything
# not installed on this machine.

book_attached_packages <- function(source_files) {
  txt <- unlist(lapply(source_files, readLines, warn = FALSE))
  m <- regmatches(txt, gregexpr("\\b(?:library|require)\\(\\s*[\"']?([A-Za-z][A-Za-z0-9.]*)", txt, perl = TRUE))
  pkgs <- unique(sub("^(?:library|require)\\(\\s*[\"']?", "", unlist(m), perl = TRUE))
  # what attaching those also attaches: Depends (recursively), and tidyverse's core set
  seen <- character(0)
  queue <- pkgs
  while (length(queue)) {
    p <- queue[1]; queue <- queue[-1]
    if (p %in% seen) next
    seen <- c(seen, p)
    dep <- suppressWarnings(utils::packageDescription(p, fields = "Depends"))
    if (is.character(dep) && !is.na(dep)) {
      d <- trimws(sub("\\(.*", "", strsplit(dep, ",")[[1]]))
      queue <- c(queue, setdiff(d[nzchar(d) & d != "R"], seen))
    }
    if (p == "tidyverse") {
      queue <- c(queue, c("ggplot2", "dplyr", "tidyr", "readr", "purrr", "tibble",
                          "stringr", "forcats", "lubridate"))
    }
  }
  seen
}

prune_cache_packages <- function(pkgs, allowed, maskers = c("ggtern", "plyr")) {
  keep <- pkgs[pkgs %in% allowed]
  keep <- setdiff(keep, maskers)
  keep[nzchar(keep) & nzchar(vapply(keep, function(p) system.file(package = p), character(1)))]
}

prune_cache_dir <- function(book_dir = ".") {
  sources <- c(file.path(book_dir, c("index.Rmd", "_common.R")),
               list.files(file.path(book_dir, "chapters"), pattern = "\\.Rmd$", full.names = TRUE))
  allowed <- book_attached_packages(sources[file.exists(sources)])
  files <- list.files(file.path(book_dir, "_bookdown_files"), pattern = "^__packages$",
                      recursive = TRUE, full.names = TRUE)
  for (f in files) {
    pkgs <- readLines(f, warn = FALSE)
    keep <- prune_cache_packages(pkgs, allowed)
    if (!identical(keep, pkgs)) {
      writeLines(keep, f)
      message("[build] pruned from ", f, ": ", paste(setdiff(pkgs, keep), collapse = ", "))
    }
  }
}
