# example R options set globally
options(width = 60)

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
