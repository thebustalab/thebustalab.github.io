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
if (!exists("algae_data")) {
  bustalab <- TRUE
  try(
    suppressMessages(
      source("https://thebustalab.github.io/phylochemistry/phylochemistry.R")
    ),
    silent = TRUE
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
