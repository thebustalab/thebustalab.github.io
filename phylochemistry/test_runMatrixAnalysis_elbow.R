## Regression test for runMatrixAnalysis(analysis = "kmeans", parameters = "elbow")
##
## WHAT THIS GUARDS (added 2026-09-30)
## -----------------------------------
## Book ch.10 teaches choosing k by the elbow, but the only thing that ever DREW the elbow
## was the Shiny picker in findClusterParameters(), which needs desktop R. Every WebR
## surface -- the escape rooms, the book's cells, the sandbox -- therefore taught the idea
## and withheld the instrument. `parameters = "elbow"` returns the table the chapter plots:
## one row per candidate k, columns `k` and `within_cluster_variance`.
##
## The flat-clustering escape room (rooms/flat_clustering/waterfalls, station2) is a
## click-the-elbow puzzle GRADED at k = 3, so the numbers below are not decoration -- they
## are the puzzle's keyed answer. If this test starts failing, a room's answer moved.
##
## Run:  Rscript test_runMatrixAnalysis_elbow.R
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

## --- fixture: three well-separated families of four, in two dimensions ------
## Deterministic, and structured so the elbow is unmistakably at k = 3: the drop from
## 2 to 3 is large, everything after 3 is small.
d <- tibble(
  specimen = paste0("s", 1:12),
  x = c(0, 1, 0, 1,   20, 21, 20, 21,   0,  1,  0,  1),
  y = c(0, 0, 1, 1,    0,  0,  1,  1,  20, 20, 21, 21)
)
vals <- c("x", "y")

elbow_for <- function(parameters) {
    suppressMessages(runMatrixAnalysis(
        data = d, analysis = "kmeans", parameters = parameters,
        columns_w_values_for_single_analyte = vals,
        columns_w_sample_ID_info = c("specimen")
    ))
}

## --- test 1: shape ----------------------------------------------------------
## The column must be called `k`: the escape-room picker tags clickable marks by the
## column named in pick.idColumn, and the flat-clustering room names "k".
e <- elbow_for(c("elbow", 6))
stopifnot(
  is.data.frame(e),
  identical(colnames(e), c("k", "within_cluster_variance")),
  nrow(e) == 6L,
  identical(as.integer(e$k), 1:6),
  is.numeric(e$within_cluster_variance)
)

## --- test 2: it is an elbow, and it bends where it should -------------------
## Monotonically non-increasing (adding a cluster can never increase total spread),
## and the k=2 -> k=3 drop dwarfs every drop after it.
drops <- -diff(e$within_cluster_variance)
stopifnot(
  all(drops >= -1e-8),
  drops[2] > 10 * max(drops[3:length(drops)])
)

## --- test 3: default range, capped at nrow - 1 ------------------------------
## 12 samples, default maximum 10 -> 10 rows (the cap does not bite here).
stopifnot(nrow(elbow_for("elbow")) == 10L)
## A maximum above nrow - 1 is capped: k = n is a degenerate zero-spread fit.
stopifnot(nrow(elbow_for(c("elbow", 50))) == 11L)

## --- test 4: it is case-insensitive and accepts a bare string ---------------
stopifnot(identical(elbow_for("Elbow")$k, elbow_for("elbow")$k))

## --- test 5: the guards ------------------------------------------------------
## "elbow" with any other analysis must STOP, not silently return that analysis.
bad <- try(suppressMessages(runMatrixAnalysis(
    data = d, analysis = "pca", parameters = "elbow",
    columns_w_values_for_single_analyte = vals,
    columns_w_sample_ID_info = c("specimen")
)), silent = TRUE)
stopifnot(
  inherits(bad, "try-error"),
  grepl("only makes sense with analysis", conditionMessage(attr(bad, "condition")), fixed = TRUE)
)

## A nonsense maximum must STOP rather than produce an empty or one-row table.
bad2 <- try(elbow_for(c("elbow", "banana")), silent = TRUE)
stopifnot(inherits(bad2, "try-error"))
bad3 <- try(elbow_for(c("elbow", 1)), silent = TRUE)
stopifnot(inherits(bad3, "try-error"))

## --- test 6: the ordinary k-means path is untouched -------------------------
## parameters = c(3) must still return a clustering, one row per sample.
km <- suppressMessages(runMatrixAnalysis(
    data = d, analysis = "kmeans", parameters = c(3),
    columns_w_values_for_single_analyte = vals,
    columns_w_sample_ID_info = c("specimen")
))
stopifnot(
  is.data.frame(km),
  nrow(km) == 12L,
  "cluster" %in% colnames(km),
  length(unique(km$cluster)) == 3L,
  all(table(km$cluster) == 4L)
)

cat("test_runMatrixAnalysis_elbow.R: all checks passed\n")
