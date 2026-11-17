## permanovaTest() regression test (ANOVA chapter's PERMANOVA section, 2026-11-17).
## Run:  ~/miniconda3/envs/phylochem/bin/Rscript test_permanovaTest.R   (from phylochemistry/)
## Extracts anovaTest + permanovaTest alone, as the WebR shim does, and pins them against vegan::adonis2,
## anovaTest(), long-vs-wide input, the chapter's two worked examples, and the input guards.
## Break check: PERMANOVA_TEST_BREAK=1 nudges one plant and must fail.
suppressMessages({ library(dplyr); library(vegan) })
src <- readLines("phylochemistry.R")
grab <- function(nm) {
  i <- grep(paste0("^\\s*", nm, " <- function"), src); j <- i; depth <- 0
  repeat {
    l <- gsub('"[^"]*"', "", src[j])
    depth <- depth + lengths(regmatches(l, gregexpr("\\{", l))) - lengths(regmatches(l, gregexpr("\\}", l)))
    if (depth == 0 && j > i) break
    j <- j + 1
  }
  eval(parse(text = src[i:j]), envir = globalenv())
}
grab("anovaTest"); grab("permanovaTest")
ft <- read.csv("../integrated_bioanalytics/data/fertilizer_trial.csv")
if (Sys.getenv("PERMANOVA_TEST_BREAK") == "1") ft$limonene[1] <- ft$limonene[1] + 5
ok <- function(cond, what) { if (!isTRUE(cond)) stop("FAIL: ", what, call. = FALSE); cat("ok  ", what, "\n") }

p2 <- ft %>% permanovaTest(cbind(limonene, pinene) ~ species * fertilizer)
ok(identical(p2$term, c("species", "fertilizer", "species:fertilizer", "noise", "total")), "two-way rows")
ok(isTRUE(all.equal(p2$share, c(0.747, 0.089, 0.040, 0.124, 1))), "two-way shares (chapter: 75/9/4%)")
ok(isTRUE(all.equal(p2$F[1:3], c(72.15, 17.12, 3.85))), "two-way F")
ref <- adonis2(ft[c("limonene", "pinene")] ~ species * fertilizer, data = ft, method = "euclidean", by = "terms", permutations = 0)
ok(isTRUE(all.equal(p2$sum_of_squares, signif(ref$SumOfSqs, 4))), "sums of squares = adonis2")

a <- ft %>% anovaTest(height_cm ~ species * fertilizer)
p1 <- ft %>% permanovaTest(height_cm ~ species * fertilizer)
cols <- c("term", "sum_of_squares", "free_values", "variance", "F", "share")
ok(isTRUE(all.equal(a[cols], p1[cols])), "one response, euclidean = anovaTest()")

long <- rbind(data.frame(plant = ft$plant, species = ft$species, fertilizer = ft$fertilizer, compound = "limonene", amount = ft$limonene),
              data.frame(plant = ft$plant, species = ft$species, fertilizer = ft$fertilizer, compound = "pinene", amount = ft$pinene))
set.seed(3); long <- long[sample(nrow(long)), ]
pl <- long %>% group_by(species) %>% permanovaTest(amount ~ species * fertilizer, profile = "compound", sample = "plant")
ok(isTRUE(all.equal(pl[c("sum_of_squares", "F", "share")], p2[c("sum_of_squares", "F", "share")])), "long (shuffled, grouped) = wide")

h <- read.csv("sample_data/hawaii_aquifers.csv") %>% filter(!analyte %in% c("dissolved_solids", "HCO3"))
hb <- h %>% permanovaTest(abundance ~ aquifer_code, profile = "analyte", sample = c("well_name", "longitude", "latitude"), distance = "bray")
he <- h %>% permanovaTest(abundance ~ aquifer_code, profile = "analyte", sample = c("well_name", "longitude", "latitude"))
ok(hb$share[1] == 0.473 && he$share[1] == 0.343 && hb$free_values[3] == 105, "hawaii: bray 47%, euclidean 34%, 106 wells (chapter)")

err <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))
ok(grepl("more than one value", err(h %>% permanovaTest(abundance ~ aquifer_code, profile = "analyte", sample = "well_name"))), "shared well names are caught")
ok(grepl("column names", err(ft %>% permanovaTest(log(limonene) ~ species))), "a transformed response is refused, not ignored")
f2 <- ft; f2$species[1] <- NA
ok(grepl("missing values", err(f2 %>% permanovaTest(cbind(limonene, pinene) ~ species))), "NA in a factor is named")
ok(grepl("full profile", err(long[-1, ] %>% permanovaTest(amount ~ species, profile = "compound", sample = "plant"))), "incomplete profile is caught")
cat("ALL OK\n")
