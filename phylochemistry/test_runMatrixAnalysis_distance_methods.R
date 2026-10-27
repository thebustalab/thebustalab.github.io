## Regression test for runMatrixAnalysis(distance_method = ...)  (added 2026-10-27)
##
## WHAT IT GUARDS
## --------------
## distance_method used to offer only euclidean, manhattan and gower. For the distances-and-scaling
## chapter it gained correlation, cosine, bray_curtis, jaccard and aitchison, all computed in base R
## by ONE helper (.rma_distance) that the dist output, neighbour joining and the bootstrapped linkage
## tree all call. The ways this can quietly go wrong, each pinned below:
##   - a metric computes something other than its textbook definition (checked against vegan where
##     vegan has it, and against hand-worked numbers otherwise);
##   - the tree branches fall back to stats::dist(), which knows none of the new names (Bclust's own
##     default FUN does exactly that), so asking for a tree with bray_curtis errors or differs;
##   - bray_curtis / jaccard / aitchison are fed SCALED data (the default), which has negative values
##     and turns them into nonsense; they must stop and say "scale_variance = FALSE";
##   - the long output drops every pair at distance 0 instead of only self-pairs (it did until this
##     change), which deletes genuinely identical pairs under jaccard / correlation / cosine;
##   - an unknown method name reaches stats::dist and silently computes some other metric.
##
## Run:  ~/miniconda3/envs/phylochem/bin/Rscript test_runMatrixAnalysis_distance_methods.R
##       (from phylochemistry/; extracts the function from phylochemistry.R rather than sourcing the
##        toolkit; needs dplyr, tibble, cluster, ape, shipunov, and vegan for the cross-checks)

suppressMessages({library(dplyr); library(tibble)})

exprs <- parse("phylochemistry.R")
fn <- NULL
for (e in exprs) {
    if (is.call(e) && length(e) >= 2 && identical(e[[1]], as.name("<-")) &&
        identical(e[[2]], as.name("runMatrixAnalysis"))) fn <- e
}
if (is.null(fn)) stop("runMatrixAnalysis definition not found in phylochemistry.R")
eval(fn)
Bclust <- shipunov::Bclust

## Fixture: abundances. S2 is exactly 2 x S1 (same pattern and composition, different size);
## S4 has the same presence/absence as S1 but different amounts.
d <- tibble(
  sample = c("S1", "S2", "S3", "S4"),
  a = c(1, 2, 0, 5),
  b = c(2, 4, 3, 1),
  c = c(3, 6, 1, 2),
  e = c(0, 0, 4, 0)
)
vals <- c("a", "b", "c", "e")
x <- as.matrix(d[, vals]); rownames(x) <- d$sample

rma <- function(method, scale = FALSE, ...) {
  out <- runMatrixAnalysis(
    data = d, analysis = "dist", distance_method = method, scale_variance = scale,
    columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = "sample", ...
  )
  if (inherits(out, "dist")) as.matrix(out)[d$sample, d$sample] else out
}
same <- function(a, b) isTRUE(all.equal(unname(as.matrix(a)), unname(as.matrix(b)), tolerance = 1e-10))
fails <- 0
check <- function(label, ok) {
  if (isTRUE(ok)) cat("ok   ", label, "\n") else { cat("FAIL ", label, "\n"); fails <<- fails + 1 }
}

## 1. euclidean unchanged (scaled default, as before the change)
check("euclidean, scaled = stats::dist(scale(x))",
      same(rma("euclidean", scale = TRUE), as.matrix(dist(scale(x)))))

## 2. each new metric against its definition
check("manhattan = stats::dist manhattan", same(rma("manhattan"), as.matrix(dist(x, "manhattan"))))
check("correlation = 1 - cor(t(x))",       same(rma("correlation"), 1 - cor(t(x))))
check("correlation: S1 vs S2 (proportional) is 0", abs(rma("correlation")["S1", "S2"]) < 1e-12)
check("cosine: S1 vs S2 (proportional) is 0",      abs(rma("cosine")["S1", "S2"]) < 1e-12)
check("cosine S1 vs S3 by hand",
      abs(rma("cosine")["S1", "S3"] - (1 - sum(x[1, ] * x[3, ]) / sqrt(sum(x[1, ]^2) * sum(x[3, ]^2)))) < 1e-12)
check("bray_curtis = vegan::vegdist bray",  same(rma("bray_curtis"), as.matrix(vegan::vegdist(x, "bray"))))
check("bray_curtis S1 vs S3 by hand (= 8/14)", abs(rma("bray_curtis")["S1", "S3"] - 8 / 14) < 1e-12)
check("jaccard = vegan::vegdist jaccard, binary",
      same(rma("jaccard"), as.matrix(vegan::vegdist(x, "jaccard", binary = TRUE))))
check("jaccard: S1 vs S4 (same presences) is 0", rma("jaccard")["S1", "S4"] == 0)

## aitchison: scale-invariant, and equals euclidean on clr of the zero-replaced matrix
xa <- x; xa[xa == 0] <- min(x[x > 0]) / 2
clr <- log(xa) - rowMeans(log(xa))
ait <- suppressMessages(rma("aitchison"))
check("aitchison = euclidean on clr", same(ait, as.matrix(dist(clr))))
## scale invariance holds only without zeros: the pseudocount is the same absolute value in every
## sample, so after replacement S2 is no longer exactly 2 x S1. Check on the zero-free columns.
ait_abc <- as.matrix(runMatrixAnalysis(
  data = d, analysis = "dist", distance_method = "aitchison", scale_variance = FALSE,
  columns_w_values_for_single_analyte = c("b", "c"), columns_w_sample_ID_info = "sample"
))
check("aitchison: S1 vs S2 (proportional, no zeros) is 0", abs(ait_abc["S1", "S2"]) < 1e-12)
msg <- tryCatch({ rma("aitchison"); "" }, message = function(m) conditionMessage(m))
check("aitchison announces its zero replacement", grepl("zeros were replaced", msg))

## 3. raw-value metrics refuse scaled data, and say how to fix it
for (m in c("bray_curtis", "jaccard", "aitchison")) {
  err <- tryCatch({ rma(m, scale = TRUE); "" }, error = function(e) conditionMessage(e))
  check(paste0(m, " with scaling stops and names scale_variance"), grepl("scale_variance = FALSE", err, fixed = TRUE))
}

## 4. unknown names are refused rather than passed to stats::dist
err <- tryCatch({ rma("maximum"); "" }, error = function(e) conditionMessage(e))
check("unknown method refused", grepl("is not available", err))

## 5. long output keeps distinct pairs at distance 0, drops only self-pairs
lj <- rma("jaccard", output_format = "long")
check("long: 12 rows (4 samples, no self-pairs)", nrow(lj) == 12)
check("long: S1-S4 jaccard 0 pair kept",
      any(lj$sample_unique_ID_sample_1 == "S1" & lj$sample_unique_ID_sample_2 == "S4" & lj$distance == 0))

## 6. the tree branches use the same helper
nj <- runMatrixAnalysis(
  data = d, analysis = "hclust_phylo", tree_method = "neighbor_joining",
  distance_method = "bray_curtis", scale_variance = FALSE,
  columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = "sample"
)
check("neighbour joining with bray_curtis = ape::nj(vegdist)",
      same(ape::cophenetic.phylo(nj)[d$sample, d$sample],
           ape::cophenetic.phylo(ape::nj(vegan::vegdist(x, "bray")))[d$sample, d$sample]))

lk <- suppressMessages(runMatrixAnalysis(
  data = d, analysis = "hclust_phylo", distance_method = "bray_curtis", scale_variance = FALSE,
  agglomeration_method = "average",
  columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = "sample"
))
ref <- ape::as.phylo(hclust(vegan::vegdist(x, "bray"), method = "average"))
check("linkage tree with bray_curtis = hclust(vegdist, average)",
      same(ape::cophenetic.phylo(lk)[d$sample, d$sample], ape::cophenetic.phylo(ref)[d$sample, d$sample]))

## 7. gower takes MIXED columns (character + logical + numeric) with scaling at its default.
##    Until 2026-10-27 daisy() rejected character columns ("invalid type character") and scale()
##    could not take them either, so the one metric meant for mixed data could not be given any.
dm <- tibble(
  item  = c("A", "B", "C"),
  kind  = c("x", "x", "y"),
  flag  = c(TRUE, FALSE, TRUE),
  value = c(1, 2, 5)
)
gm <- as.matrix(runMatrixAnalysis(
  data = dm, analysis = "dist", distance_method = "gower",
  columns_w_values_for_single_analyte = c("kind", "flag", "value"), columns_w_sample_ID_info = "item"
))
ref_df <- data.frame(kind = factor(dm$kind), flag = factor(dm$flag), value = dm$value, row.names = dm$item)
check("gower on mixed columns = daisy on factors",
      same(gm[dm$item, dm$item], as.matrix(cluster::daisy(ref_df, metric = "gower"))))
## A vs B: kind same (0), flag differs (1), value |1-2|/4 = 0.25 -> mean 1.25/3
check("gower A vs B by hand (= 1.25/3)", abs(gm["A", "B"] - 1.25 / 3) < 1e-12)

if (fails > 0) stop(fails, " check(s) failed")
cat("test_runMatrixAnalysis_distance_methods.R: all checks passed\n")
