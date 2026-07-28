# Build entry point for the Integrated Bioanalytics bookdown.
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
bookdown::render_book("index.Rmd")
