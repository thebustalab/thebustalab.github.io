# datasets.R

message("Loading datasets module...")

## Sample datasets

    sample_datasets <- as.data.frame(rbind(
        c("algae_data", "https://thebustalab.github.io/phylochemistry/sample_data/algae_data.csv"),
        c("alaska_lake_data", "https://thebustalab.github.io/phylochemistry/sample_data/alaska_lake_data.csv"),
        c("solvents", "https://thebustalab.github.io/phylochemistry/sample_data/solvents.csv"),
        c("periodic_table", "https://thebustalab.github.io/phylochemistry/sample_data/per_table.csv"),
        c("fadb_sample", "https://thebustalab.github.io/phylochemistry/sample_data/fadb_sample.csv"),
        c("periodic_table_subset", "https://thebustalab.github.io/phylochemistry/sample_data/per_table_small.csv"),
        c("ny_trees", "https://thebustalab.github.io/phylochemistry/sample_data/ny_trees.csv"),
        c("metabolomics_data", "https://thebustalab.github.io/phylochemistry/sample_data/metabolomics_data.csv"),
        c("wine_grape_data", "https://thebustalab.github.io/phylochemistry/sample_data/wine_grape_data.csv"),
        c("hawaii_aquifers", "https://thebustalab.github.io/phylochemistry/sample_data/hawaii_aquifers.csv"),
        c("beer_components", "https://thebustalab.github.io/phylochemistry/sample_data/beer_components.csv"),
        c("wood_smoke", "https://thebustalab.github.io/phylochemistry/sample_data/wood_smoke_data.csv"),
        c("unknown_smoke", "https://thebustalab.github.io/phylochemistry/sample_data/unknown_smoke.csv"),
        c("hops_components", "https://thebustalab.github.io/phylochemistry/sample_data/hops_components.csv"),
        c("tequila_chemistry", "https://thebustalab.github.io/phylochemistry/sample_data/tequila_chemistry.csv"),
        c("chemical_blooms", "https://thebustalab.github.io/phylochemistry/sample_data/chemical_blooms.csv"),
        c("metabolomics_unknown", "https://thebustalab.github.io/phylochemistry/sample_data/metabolomics_unknown.csv"),
        c("wine_quality", "https://thebustalab.github.io/phylochemistry/sample_data/wine_quality.csv"),
        c("lake_superior_shoreline", "https://thebustalab.github.io/phylochemistry/sample_data/lake_superior_shoreline.csv"),
        c("rice_proteins", "https://thebustalab.github.io/phylochemistry/sample_data/rice_proteins.csv"),
        c("mushrooms", "https://thebustalab.github.io/phylochemistry/sample_data/mushrooms.csv"),
        c("unknown_mushroom", "https://thebustalab.github.io/phylochemistry/sample_data/unknown_mushroom.csv")
    ))

    # pb <- progress::progress_bar$new(total = dim(sample_datasets)[1])
    for (i in 1:dim(sample_datasets)[1]) {
        temp_obj <- readr::read_csv(sample_datasets[i,2], show_col_types = FALSE)
        gdata::mv("temp_obj", as.character(sample_datasets[i,1]))
        # pb$tick()
    }

    tmpfile <- tempfile(fileext = ".fasta")
    utils::download.file(
      "https://raw.githubusercontent.com/thebustalab/thebustalab.github.io/refs/heads/master/phylochemistry/sample_data/OSCs.fasta",
      tmpfile,
      quiet = TRUE
    )
    OSC_sequences <- Biostrings::readAAStringSet(tmpfile)

## Busta lab specific datasets

    if (exists("bustalab")) {
        
        if (bustalab == TRUE) {

            busta_spectral_library <- read_csv("https://thebustalab.github.io/phylochemistry/sample_data/busta_spectral_library_v1.csv", col_types = c(Compound_common_name = "c"))
            plant_phylogeny <- read.tree("https://thebustalab.github.io/data/plant_phylogeny.newick")
            plant_species <- readMonolist("https://thebustalab.github.io/data/plant_species.csv")
        
        }
    }


## Load color schemes

    cont_1 <- c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")
    cont_2 <- c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6")
    discrete_palette <- c(
        "dodgerblue2", "#E31A1C", # red
        "green4",
        "#6A3D9A", # purple
        "#FF7F00", # orange
        "black", "gold1",
        "skyblue2", "#FB9A99", # lt pink
        "palegreen2",
        "#CAB2D6", # lt purple
        "#FDBF6F", # lt orange
        "gray70", "khaki2",
        "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
        "darkturquoise", "green1", "yellow4", "yellow3",
        "darkorange4", "brown"
    )