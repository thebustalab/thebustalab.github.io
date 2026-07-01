## Regression test for codonAlignFromAA() (phylochemistry.R)
##
## codonAlignFromAA() is the pure-R replacement for the former
## orthologr::codon_aln() / PAL2NAL dependency in alignSequences(mode="codon_align").
## It back-translates a gapped amino-acid alignment onto its codons.
##
## Run:  Rscript test_codonAlignFromAA.R   (from the phylochemistry/ directory)
## Needs only base R -- it extracts the function definition from phylochemistry.R
## without sourcing the whole toolkit, so it runs on a machine with no bioinformatics
## packages installed.

## --- load only the function under test -------------------------------------
exprs <- parse("phylochemistry.R")
fn <- NULL
for (e in exprs) {
    if (is.call(e) && length(e) >= 2 && identical(e[[1]], as.name("<-")) &&
        identical(e[[2]], as.name("codonAlignFromAA"))) fn <- e
}
if (is.null(fn)) stop("codonAlignFromAA definition not found in phylochemistry.R")
eval(fn)

## --- test 1: exact-translation case with a gap -----------------------------
## This is the canonical PAL2NAL example; the reference tool produces the same
## output (NUC2 -> ATGGCT---GAT), so we pin to it.
aa <- c(NUC1 = "MACD", NUC2 = "MA-D")
nt <- c(NUC1 = "ATGGCTTGTGAT", NUC2 = "ATGGCTGAT")
res <- codonAlignFromAA(aa, nt)
stopifnot(
    res[["NUC1"]] == "ATGGCTTGTGAT",
    res[["NUC2"]] == "ATGGCT---GAT"
)

## --- test 2: FAILURE MODE -- peptide/nucleotide length mismatch -------------
## The whole point of the guard: if the peptide is not an exact in-frame
## translation of its nucleotides (duplicate accession, off-by-one ORF, stop
## codon on one side only), fail loudly instead of silently misaligning.
e2 <- tryCatch(codonAlignFromAA(aa, c(NUC1 = "ATGGCTTGTGAT", NUC2 = "ATGGCTGA")),
               error = function(e) conditionMessage(e))
stopifnot(grepl("inconsistency between the peptide and nucleotide", e2))

## --- test 3: FAILURE MODE -- no nucleotide sequence for a peptide name ------
e3 <- tryCatch(codonAlignFromAA(aa, c(NUC1 = "ATGGCTTGTGAT")),
               error = function(e) conditionMessage(e))
stopifnot(grepl("no nucleotide sequence named", e3))

cat("codonAlignFromAA: all regression tests PASS\n")
