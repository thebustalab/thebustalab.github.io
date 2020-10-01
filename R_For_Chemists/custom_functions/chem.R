#### libraries

    library(tidyverse)
    library(ape)
    library(ggtree)

#### datasets
    
    algae_data <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/algae_data.csv")
    alaska_lake_data <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/alaska_lake_data.csv")
    solvents <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/solvents.csv")
    periodic_table <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/per_table.csv")
    periodic_table_small <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/per_table_small.csv")
    NY_trees <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/NY_trees.csv")
    ckd_data <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/ckd_metabolomics.csv")
    wine_grape_data <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/wine_grape_data.csv")

#### readCSV

    #'

    readCSV <- function() { return(readr::read_csv(file.choose())) }


#### runMatrixAnalysis

    #' Runs a matrix analysis (clustering, kmeans, pca).
    #'
    #' @param data The data.frame or tibble to use.
    #' @param analysis  
    #' @param column_w_names_of_multiple_analytes
    #' @param column_w_values_for_multiple_analytes
    #' @param columns_w_values_for_single_analyte
    #' @param columns_w_additional_analyte_info
    #' @param columns_for_sample_unique_ID
    #' @param columns_w_sample_annotation_info
    #' @examples
    #' @export
    #' runMatrixAnalysis

        runMatrixAnalysis <-    function(
                                
                                    data,

                                    analysis = c("hclust", "pca", "pca-ord", "pca-dim"),
                                    
                                    column_w_names_of_multiple_analytes,
                                    column_w_values_for_multiple_analytes,

                                    columns_w_values_for_single_analyte,

                                    columns_w_additional_analyte_info,

                                    columns_w_sample_ID_info

                                ) {

            # Check that column names are spelled correctly

                if( any(
                    !c(
                        column_w_names_of_multiple_analytes,
                        column_w_values_for_multiple_analytes,
                        columns_w_values_for_single_analyte,
                        columns_w_additional_analyte_info,
                        columns_w_sample_ID_info
                    ) %in% colnames(data)
                    ) == TRUE
                ) {
                    stop("There is a mismatch in the column names delivered to the command and the column names in your data. Please double check the spelling of your column names you gave to the command.")
                }

            # Add analyte_unique_ID_column if necessary

                # if( length(columns_for_analyte_unique_ID) > 1 ) {
                    #add analyte_unique_ID column if necessary
                # }

            # Remove columns that are not included in input column lists

                if (length(
                        which(!colnames(data) %in% 
                            c(
                                column_w_names_of_multiple_analytes,
                                column_w_values_for_multiple_analytes,

                                columns_w_values_for_single_analyte,
                                columns_w_additional_analyte_info,

                                columns_w_sample_ID_info
                            )
                        )
                    ) > 0 
                ) {
                    data <- data[,-which(!colnames(data) %in% 
                        c(
                            column_w_names_of_multiple_analytes,
                            column_w_values_for_multiple_analytes,

                            columns_w_values_for_single_analyte,
                            columns_w_additional_analyte_info,

                            columns_w_sample_ID_info
                        )
                    )]
                }

            # Check for duplicate analyte names

                if( length(columns_w_values_for_single_analyte) > 0 ) {
                    if( any(duplicated(columns_w_values_for_single_analyte)) ) {
                        stop("There are duplicate analyte names.")
                    }
                }

                if( length(column_w_names_of_multiple_analytes) > 0 ) {
                    x <- table(select(data, column_w_names_of_multiple_analytes))
                    if( all(range(x) / mean(x) != c(1,1)) ) {
                        stop("There are duplicate analyte names.")
                    }
                } 

            # Remove analyte annotation columns before pivoting

                if( length(columns_w_additional_analyte_info) > 0 ) {
                    analyte_annotation_free_data <- data[,-match(columns_w_additional_analyte_info, colnames(data))]
                } else {
                    analyte_annotation_free_data <- data
                }

            # If no pivot required, skip pivoting

                if( length(column_w_names_of_multiple_analytes) == 0 & length(columns_w_values_for_single_analyte) >= 1 ) {
                    data_wide <- analyte_annotation_free_data
                    analyte_columns <- columns_w_values_for_single_analyte
                }

            # If pivoting required, pivot_wider any long-style data

                if( length(column_w_names_of_multiple_analytes) == 1 ) {
                    data_wide <- pivot_wider(
                        analyte_annotation_free_data,
                        names_from = all_of(column_w_names_of_multiple_analytes),
                        values_from = all_of(column_w_values_for_multiple_analytes)
                    )
                    analyte_columns <- unlist(unique(analyte_annotation_free_data[,colnames(analyte_annotation_free_data) == column_w_names_of_multiple_analytes]))
                    analyte_columns <- c(columns_w_values_for_single_analyte, analyte_columns)
                }

            # Convert analyte columns to numeric

                which_analyte_columns <- which(colnames(data_wide) %in% analyte_columns)
                for( i in which_analyte_columns ) {
                  data_wide[[i]] <- as.numeric(data_wide[[i]])
                }

            # Add sample_unique_ID_column if necessary, or just change column name of existing sample_unique_ID column

                if( length(columns_w_sample_ID_info) > 1 ) {
                    sample_unique_IDs <- apply(
                        data_wide[,match(columns_w_sample_ID_info, colnames(data_wide))],
                        1, paste, collapse = "_"
                    )
                    if( any(duplicated(sample_unique_IDs)) ) {stop("columns_w_sample_ID_info specified do not lead to unique sample IDs")}
                    data_wide$sample_unique_ID <- sample_unique_IDs
                } else {
                    colnames(data_wide)[colnames(data_wide) == columns_w_sample_ID_info] <- "sample_unique_ID"
                    if( any(duplicated(data_wide$sample_unique_ID)) ) {stop("columns_w_sample_ID_info specified do not lead to unique sample IDs")}
                }

                # Make sure "sample_unique_ID" is character
                    data_wide$sample_unique_ID <- as.character(data_wide$sample_unique_ID)

            # Clustering analysis

                # Prepare the matrix

                    matrix <- as.data.frame(data_wide[,match(analyte_columns, colnames(data_wide))])
                    rownames(matrix) <- data_wide$sample_unique_ID

                # Run hclust, if requested

                    if( analysis == "hclust" ) {
                        clustering <- ggtree::fortify(ape::as.phylo(stats::hclust(stats::dist(matrix))))
                        clustering$sample_unique_ID <- clustering$label 
                    }

                # Run PCA, if requested

                    as.data.frame(matrix)

                    if( analysis == "pca" ) {
                        coords <- FactoMineR::PCA(matrix, graph = FALSE)$ind$coord[,c(1:2)]
                        clustering <- as_tibble(coords)
                        clustering$sample_unique_ID <- rownames(coords)
                    }

                # Run PCA and return ordination plot coordinates, if requested

                    if( analysis == "pca-ord" ) {
                        coords <- FactoMineR::PCA(matrix, graph = FALSE)$var$coord[,c(1,2)]
                        clustering <- as_tibble(coords)
                        clustering$analyte <- rownames(coords)
                        clustering <- select(clustering, analyte, Dim.1, Dim.2)
                        return(clustering)
                        stop("Returning ordination plot coordinates.")
                    }

                # Run PCA and return eigenvalues, if requested

                    if( analysis == "pca-dim" ) {
                        coords <- FactoMineR::PCA(matrix, graph = FALSE)$eig[,2]
                        clustering <- tibble::enframe(coords, name = NULL)
                        clustering$principal_component <- names(coords)
                        clustering$principal_component <- as.numeric(gsub("comp ", "", clustering$principal_component))
                        colnames(clustering)[colnames(clustering) == "value"] <- "percent_variance_explained"
                        clustering <- select(clustering, principal_component, percent_variance_explained)
                        return(clustering)
                        stop("Returning eigenvalues.")
                    }
                    
                # Add back annotations to the output

                    clustering <- full_join(
                        data_wide[,match(
                            if( length(columns_w_sample_ID_info) == 1 ) {
                                "sample_unique_ID"
                            } else {
                                c(columns_w_sample_ID_info, "sample_unique_ID")
                            }, colnames(data_wide))],
                        clustering,
                        by = "sample_unique_ID"
                    )

                    rownames_matrix <- tibble::enframe(rownames(matrix), name = NULL)
                    colnames(rownames_matrix)[1] <- "sample_unique_ID"

                    clustering <- full_join(
                        clustering,
                        as_tibble(cbind(rownames_matrix, as_tibble(matrix))),
                        by = "sample_unique_ID"
                    )
                    clustering

            # Return results

                return( clustering )

        }