# example R options set globally
options(width = 60)

# Ensure a CRAN mirror is always set. The book renders in a non-interactive session
# (new_session: true), where repos defaults to "@CRAN@", and anything that installs or checks a
# package then errors ("trying to use CRAN without setting a mirror") unless a real mirror is set.
if (is.null(getOption("repos")) || any(getOption("repos") == "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# ---------------------------------------------------------------------------------------------
# LOAD THE PACKAGES: library(tidyverse), then phylochemistry -- the same two lines the book teaches.
#
# WHICH phylochemistry? The package SOURCE in this same synced tree whenever it is there, loaded
# with pkgload::load_all(); the INSTALLED package (r-universe) only as a fallback. The chapters and
# the package are edited together and deploy.sh compiles the book BEFORE anything is published, so
# compiling against the installed release would test the PREVIOUS version and a chapter written
# against a same-commit package change could not compile. (Same rule, same reason, as when this file
# sourced ../phylochemistry/phylochemistry.R: that is how ch7 broke on 2026-09-22.)
#
# The source lives at Tools/public_tools/phylochemistry/phylochemistry, i.e. three levels up from
# this folder on both the Mac and host2. Overrides, for a test or an unusual checkout:
#   PHYLOCHEMISTRY_PKG_DIR=<path to the package source>
#   PHYLOCHEMISTRY_BOOK_LOAD=installed   (use library(phylochemistry) even when the source is there)
#
# There is no bustalab <- TRUE any more: no chapter needs the lab package, and the book must render
# with exactly what a student installs.
.pc_pkg_dir <- Sys.getenv("PHYLOCHEMISTRY_PKG_DIR",
                          file.path("..", "..", "..", "public_tools", "phylochemistry", "phylochemistry"))
.pc_use_source <- !identical(Sys.getenv("PHYLOCHEMISTRY_BOOK_LOAD"), "installed") &&
  file.exists(file.path(.pc_pkg_dir, "DESCRIPTION"))

if (!"package:tidyverse" %in% search()) suppressPackageStartupMessages(library(tidyverse))

# A HALF-LOADED LIBRARY MUST NOT PASS SILENTLY (2026-09-22): stop with the real load error here,
# rather than letting it surface in a chapter chunk with an error that says nothing about the library.
if (!"package:phylochemistry" %in% search()) {
  .pc_err <- tryCatch({
    if (.pc_use_source) {
      if (!requireNamespace("pkgload", quietly = TRUE)) {
        stop("the book loads phylochemistry from its source folder and needs pkgload: ",
             "install.packages(\"pkgload\")")
      }
      suppressMessages(pkgload::load_all(.pc_pkg_dir, export_all = FALSE, helpers = FALSE,
                                         attach_testthat = FALSE, quiet = TRUE))
    } else {
      suppressPackageStartupMessages(library(phylochemistry))
    }
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(.pc_err)) {
    stop(
      "phylochemistry failed to load from ",
      if (.pc_use_source) normalizePath(.pc_pkg_dir, mustWork = FALSE) else "the installed library",
      "\n  ", .pc_err,
      "\nThe book compiles against the package SOURCE on purpose (see AGENTS.md). Fix the package;",
      " do NOT fall back to the published release, which is a different version."
    )
  }
}

# LOCAL BUILD HARNESS (optional, gitignored, absent on the Mac). Translates the Mac paths and env vars
# the chapters expect onto this machine. Sourced AFTER the packages load, because it wraps package
# functions. No-op if the file is not there.
if (file.exists("_local_harness.R")) try(source("_local_harness.R"), silent = TRUE)

# ANY source() OF A thebustalab.github.io/phylochemistry/... URL IS REDIRECTED TO THE LOCAL WORKING
# COPY when one exists (2026-09-22). No chapter sources the toolkit any more, but the appendix still
# shows `source(".../modules/gcms.R")` (the standalone GC-MS app), and a render must never pull a
# published file over the local tree. No-op wherever the local file is absent. Render-only.
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

# The two large datasets the book reads with getDataset() are downloaded on first use, and that first
# use prints a message naming the local cache folder -- which would land in the rendered HTML (a
# private path; the pre-push scan blocks it). Fetch them quietly here so the chapter chunks always
# hit the cache.
try(local({
  cat_ <- phylochemistry::getDataset()
  for (nm in intersect(c("ny_trees", "lake_superior_shoreline"), cat_$name[!cat_$cached])) {
    suppressMessages(phylochemistry::getDataset(nm))
  }
}), silent = TRUE)

# Post-load assertion, on BEHAVIOUR (2026-09-22): the runMatrixAnalysis() in scope must be the
# package's, and must have the long-format `dist` branch chapter 7 depends on. A formals-only check
# is useless (every stale fork still HAS an output_format argument), and a stray copy in the global
# environment (someone sourcing the old phylochemistry.R) would silently shadow the package.
.pc_rma <- get0("runMatrixAnalysis", mode = "function")
if (is.null(.pc_rma)) {
  stop("phylochemistry loaded but runMatrixAnalysis() is missing -- the package did not finish loading.")
}
if (!identical(environmentName(environment(.pc_rma)), "phylochemistry")) {
  stop("The runMatrixAnalysis() in scope does not come from the phylochemistry package (found in ",
       environmentName(environment(.pc_rma)), "). Something sourced an old copy over it.")
}
if (!any(grepl('"sample_1", "sample_2", "distance"', deparse(body(.pc_rma)), fixed = TRUE))) {
  stop(
    "The runMatrixAnalysis() in scope has no long-format branch for analysis = \"dist\", so it\n",
    "  returns a bare `dist` object and chapter 7 cannot plot it. Check which phylochemistry loaded."
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

# ---------------------------------------------------------------------------------------------
# PER-CHUNK CACHE KEYING ON THE PACKAGE FUNCTIONS AND DATASETS A CHUNK ACTUALLY USES (2026-09-23;
# moved from the sourced script to the package namespace with the package switch)
#
# THE PROBLEM. `cache = TRUE` is global (index.Rmd setup). knitr keys a chunk on its own CODE only,
# so a package change that alters what a function RETURNS leaves every chunk a cache hit and the new
# behaviour never runs. The blunt fix -- key every chunk on a checksum of the whole package -- makes
# any package edit a COLD rebuild of the whole book, which re-runs the deliberately-live ch13/ch14
# API chunks. Most of the compile time is in chunks that never touch the package.
#
# WHAT THIS DOES. Fingerprint every function in the phylochemistry NAMESPACE (exported and internal),
# then give each chunk a `cache.extra` built from just the fingerprints of the exported functions and
# lazy datasets THAT chunk names. Change `runMatrixAnalysis()` and only the chunks calling it re-run.
#
# INDIRECT CALLS ARE COVERED. A function's fingerprint is the transitive closure: its own body plus the
# bodies of every package function it calls (internal helpers included), and everything they call.
#
# DATASETS ARE COVERED TOO (new with the package): each lazy dataset (`algae_data`, ...) is hashed,
# so a chunk naming one re-runs when its data/ file changes. getDataset() downloads are NOT covered.
#
# Names are taken from the chunk's PARSED code (`all.names`), not a regex, so a mention in a string
# or a comment does not count and a call does. Source references are stripped before hashing, so the
# key does not depend on file paths or on load_all() versus library().
#
# LIMITS: package code that runs at LOAD time (.onLoad), and getDataset() downloads, are invisible to
# this. `../deploy.sh --cold-book` remains the belt-and-braces answer before a semester starts.
local({
  ns <- tryCatch(asNamespace("phylochemistry"), error = function(e) NULL)
  if (is.null(ns)) return(invisible(NULL))
  fn_names <- Filter(function(n) {
    obj <- get0(n, envir = ns, inherits = FALSE)
    is.function(obj) && !is.primitive(obj)
  }, ls(ns, all.names = TRUE))
  if (length(fn_names) == 0) return(invisible(NULL))

  hash_of <- if (requireNamespace("digest", quietly = TRUE)) {
    function(x) digest::digest(x)
  } else {
    function(x) if (is.character(x)) paste(x, collapse = "\n") else NA_character_
  }

  own <- vapply(fn_names, function(n) {
    f <- utils::removeSource(get(n, envir = ns))
    hash_of(paste(c(deparse(formals(f)), deparse(body(f))), collapse = "\n"))
  }, character(1))

  ## direct package-callees of each function, from the language object (no re-parsing)
  deps <- lapply(fn_names, function(n) {
    b <- body(get(n, envir = ns))
    if (is.null(b)) character(0) else intersect(all.names(b), fn_names)
  })
  names(deps) <- fn_names

  ## transitive closure, iterated to a fixed point (cheap: a few passes over a few hundred names)
  repeat {
    grown <- lapply(deps, function(d) unique(c(d, unlist(deps[d], use.names = FALSE))))
    if (identical(lengths(grown), lengths(deps))) { deps <- grown; break }
    deps <- grown
  }

  fingerprint <- vapply(fn_names, function(n) {
    hash_of(paste(own[sort(unique(c(n, deps[[n]])))], collapse = "|"))
  }, character(1))
  fingerprint <- fingerprint[names(fingerprint) %in% getNamespaceExports(ns)]

  ## lazy datasets: hash each object (only when digest is available; raw text would be huge)
  lazy <- tryCatch(getNamespaceInfo(ns, "lazydata"), error = function(e) NULL)
  if (is.environment(lazy) && requireNamespace("digest", quietly = TRUE)) {
    ds <- ls(lazy, all.names = TRUE)
    if (length(ds)) {
      fingerprint <- c(fingerprint, vapply(ds, function(d) digest::digest(get(d, envir = lazy)),
                                           character(1)))
    }
  }

  assign(".pc_fingerprint", fingerprint, envir = globalenv())

  knitr::opts_hooks$set(cache = function(options) {
    if (is.null(options$cache) || identical(options$cache, FALSE) ||
        identical(options$cache, 0)) return(options)
    ids <- tryCatch(all.names(parse(text = options$code), unique = TRUE),
                    error = function(e) character(0))
    used <- intersect(ids, names(.pc_fingerprint))
    if (length(used)) {
      options$cache.extra <- c(options$cache.extra,
                               unname(.pc_fingerprint[sort(used)]))
    }
    options
  })
})

# BUILD STAMP -- one line, kept deliberately (2026-09-22). Written to stderr so it lands in the deploy
# log: which phylochemistry loaded (source folder or installed, and its version), whether the URL
# redirect is armed, and whether the function ch7 depends on is the current one. If it ever reads
# long_dist=FALSE, stop and fix that before touching a chapter.
cat(sprintf(
  "[_common.R] phylochemistry=%s %s | redirect=%s | long_dist=%s | fingerprints=%d\n",
  if (.pc_use_source) paste0("source:", normalizePath(.pc_pkg_dir, mustWork = FALSE)) else "installed",
  as.character(utils::packageVersion("phylochemistry")),
  isTRUE(attr(get0("source", envir = globalenv()), "pc_redirect")),
  any(grepl('"sample_1", "sample_2", "distance"', deparse(body(.pc_rma)), fixed = TRUE)),
  length(get0(".pc_fingerprint", envir = globalenv(), ifnotfound = character(0)))
), file = stderr())
