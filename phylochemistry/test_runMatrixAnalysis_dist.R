## Regression test for runMatrixAnalysis(analysis = "dist", output_format = "long")
##
## FAILURE MODE (found in WebR by Lucas, 2026-09-23; fixed same day)
## ----------------------------------------------------------------
## A call with exactly ONE `columns_w_sample_ID_info` died inside
##   data_wide[, match(c(columns_w_sample_ID_info, "sample_unique_ID"), colnames(data_wide))]
## while the SAME call with two ID columns worked. The error named a subsetting
## expression and neither the argument at fault nor the reason.
##
## WHY: the pre-processing step branches on the number of ID columns. With two or more it
## keeps them as they are and PASTES a new `sample_unique_ID` alongside. With exactly one it
## RENAMES that column to `sample_unique_ID` -- so the name the caller passed no longer
## exists, `match()` returns NA, and the subset blows up.
##
## THE FIX: the long-format `dist` branch selects only the ID columns that actually survived
## pre-processing (`intersect(columns_w_sample_ID_info, colnames(data_wide))` plus
## `sample_unique_ID`). In the one-column case `sample_unique_ID` IS the caller's column, so
## nothing is lost.
##
## THE SIGNATURE THIS TEST KEYS ON: a one-ID-column call returns a data frame whose first two
## columns are sample_unique_ID_sample_1 / sample_unique_ID_sample_2 -- NOT an error, and NOT
## a bare `dist` object (that second failure is the 2026-09-22 ch7 one; see
## ../integrated_bioanalytics/_agent_reference/render_failure_modes.md).
##
## Run:  ~/miniconda3/envs/phylochem/bin/Rscript test_runMatrixAnalysis_dist.R
##       (from the phylochemistry/ directory; needs dplyr + tibble, nothing bioinformatic --
##        it extracts the function from phylochemistry.R instead of sourcing the toolkit)

suppressMessages({library(dplyr); library(tibble)})

## --- load only the function under test -------------------------------------
exprs <- parse("phylochemistry.R")
fn <- NULL
for (e in exprs) {
    if (is.call(e) && length(e) >= 2 && identical(e[[1]], as.name("<-")) &&
        identical(e[[2]], as.name("runMatrixAnalysis"))) fn <- e
}
if (is.null(fn)) stop("runMatrixAnalysis definition not found in phylochemistry.R")
eval(fn)

## --- a small, fully deterministic fixture ----------------------------------
d <- tibble(
  lake = c("Alpha", "Beta", "Gamma", "Delta"),
  park = c("P1", "P1", "P2", "P2"),
  a    = c(1, 2, 3, 4),
  b    = c(4, 3, 2, 1),
  c    = c(0, 1, 0, 1)
)
vals <- c("a", "b", "c")

## --- test 1: THE BUG -- one ID column must work ----------------------------
one <- runMatrixAnalysis(
  data = d, analysis = "dist", output_format = "long", scale_variance = TRUE,
  columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = c("lake")
)
stopifnot(
  is.data.frame(one),
  length(dim(one)) == 2L,
  identical(colnames(one)[1:3],
            c("sample_unique_ID_sample_1", "sample_unique_ID_sample_2", "distance")),
  nrow(one) == 12,                                  # 4 samples -> 4*3 ordered pairs
  all(sort(unique(one$sample_unique_ID_sample_1)) == sort(d$lake))
)

## --- test 2: two ID columns still behave, and carry their columns through ---
two <- runMatrixAnalysis(
  data = d, analysis = "dist", output_format = "long", scale_variance = TRUE,
  columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = c("lake", "park")
)
stopifnot(
  identical(colnames(two)[1:3],
            c("sample_unique_ID_sample_1", "sample_unique_ID_sample_2", "distance")),
  all(c("lake_sample_1", "park_sample_1", "lake_sample_2", "park_sample_2") %in% colnames(two)),
  nrow(two) == 12
)

## --- test 3: the two calls agree on the NUMBERS ----------------------------
## The ID columns differ, the distances must not. Pin one against the other so a future
## change to either branch cannot silently move the values in just one of them.
key1 <- paste(one$sample_unique_ID_sample_1, one$sample_unique_ID_sample_2)
key2 <- paste(two$lake_sample_1, two$lake_sample_2)
stopifnot(all.equal(one$distance[order(key1)], two$distance[order(key2)]))

## --- test 4: wide is still a `dist`, i.e. the default did not shift ---------
wide <- runMatrixAnalysis(
  data = d, analysis = "dist", scale_variance = TRUE,
  columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = c("lake")
)
stopifnot(inherits(wide, "dist"))

cat("runMatrixAnalysis dist/long: OK (one ID column, two ID columns, agreement, wide default)\n")
