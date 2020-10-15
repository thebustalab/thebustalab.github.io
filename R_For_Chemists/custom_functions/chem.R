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
    # NY_trees <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/NY_trees.csv")
    ckd_data <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/ckd_metabolomics.csv")
    wine_grape_data <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/wine_grape_data.csv")
    data <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/housing.csv")

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
                                    
                                    column_w_names_of_multiple_analytes = NULL,
                                    column_w_values_for_multiple_analytes = NULL,

                                    columns_w_values_for_single_analyte = NULL,

                                    columns_w_additional_analyte_info = NULL,

                                    columns_w_sample_ID_info = NULL,

                                    transpose = FALSE,

                                    kmeans = c("none", "auto", "elbow", "1", "2", "3", "etc."),

                                    na_replacement = c("none", "mean", "zero", "drop"),

                                    ...

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
                    data_wide <- unique(data_wide)
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

                # Replace NAs with colmeans

                    if( na_replacement[1] == "none") {
                        cat("Not replacing any NAs in your data set \n")
                    }
                    if( na_replacement[1] == "drop" ) {
                        cat("Dropping any variables in your dataset that have NA as a value.\nVariables dropped:\n")
                        cat(names(which(apply(is.na(matrix), 2, any))))
                        cat("\n")
                        matrix <- matrix[,!apply(is.na(matrix), 2, any)]
                    }
                    if( na_replacement[1] %in% c("zero", "mean") ) {

                        if( any(is.na(matrix)) ) {
                            
                            cat(paste0("Replacing NAs in your data with ", na_replacement), "\n")

                                for( column in 1:dim(matrix)[2]) {
                                    
                                    if( any(is.na(matrix[,column])) ) {

                                        if( na_replacement == "mean" ) {
                                            replacement <- mean(matrix[,column], na.rm = TRUE)
                                        }
                                        if( na_replacement == "zero" ) {
                                            replacement <- 0
                                        }
                                        if( !any(na_replacement %in% c("mean", "zero")) ) {
                                            stop("Your data contains NAs. Please specify how to deal with them using na_replacement. \n")
                                        }
                                        
                                        matrix[,column][is.na(matrix[,column])] <- replacement

                                    } else {}
                                }
                        }
                    }

                    if( transpose == TRUE ) {

                        matrix <- t(matrix)

                    }

                # Run hclust, if requested

                    if( analysis == "hclust" ) {
                        phylo <- ape::as.phylo(stats::hclust(stats::dist(matrix)))
                        clustering <- ggtree::fortify(phylo)
                        clustering$sample_unique_ID <- clustering$label
                    }

                # Run PCA, if requested

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

                # Run kMeans, if requested

                    if( kmeans[1] == "none" ) {}

                    if( kmeans[1] != "none" ) {

                        ## Check for NAs

                            if( any(is.na(matrix)) == TRUE ) {
                                stop("kmeans cannot handle NA. Please choose an option for na_replacement.")
                            }

                        ## Run k-means
                            kmeans_results <- list()
                            for( i in 1:(dim(matrix)[1]-1) ) {
                                kmeans_results[[i]] <- sum(stats::kmeans(x = matrix, centers = i, nstart = 25, iter.max = 1000)$withinss)
                            }
                            kmeans_results <- do.call(rbind, kmeans_results)
                            kmeans_results

                        ## If auto, determine sharpest angle and return clusters

                            if( kmeans[1] == "auto" ) {
                                angles <- list()
                                for( i in 1:(length(kmeans_results)-2) ) {
                                    slope_1 <- kmeans_results[i+1] - kmeans_results[i]
                                    slope_2 <- kmeans_results[i+2] - kmeans_results[i+1]
                                    angles[[i]] <- atan( (slope_1 - slope_2) / (1 + slope_1*slope_2) )
                                }
                                angles <- do.call(rbind, angles)

                                cluster_number <- which(angles == min(angles)) + 1

                                kmeans_clusters <- stats::kmeans(x = matrix, centers = cluster_number, nstart = 25, iter.max = 1000)$cluster
                            }

                        ## If scree, return raw results

                            if( kmeans[1] == "elbow" ) {
                                results <- as_tibble(data.frame(
                                    cluster_number = seq(1, dim(kmeans_results)[1], 1),
                                    variance_within_cluster = kmeans_results
                                ))
                                return( results )
                                stop("Returning elbow plot.")
                            }

                        ## If a number, return that many clusters

                            if( !kmeans[1] %in% c("auto", "elbow") ) {
                                if( is.na(as.numeric(kmeans[1])) ) {stop("For kmeans, please use 'none', 'auto', 'elbow', or a number.")}
                                kmeans_clusters <- stats::kmeans(x = matrix, centers = as.numeric(kmeans[1]), nstart = 25, iter.max = 1000)$cluster
                            }

                        ## Bind kmeans cluster numbers to output

                            kmeans_clusters <- as_tibble(data.frame(sample_unique_ID = names(kmeans_clusters), kmeans_cluster = paste0("cluster_", kmeans_clusters)))
                            clustering$kmeans_cluster <- kmeans_clusters$kmeans_cluster[match(clustering$sample_unique_ID, kmeans_clusters$sample_unique_ID)]

                    }

                # Return without annotations if transpose = TRUE

                    if( transpose == TRUE ) {
                        return(clustering)
                        stop("Returning transposed cluster output. Make sure all your variables have the same units!")
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

                # Annotate internal branches if tree output

                    if( analysis == "hclust" ) {
                        for( node in filter(clustering, isTip == FALSE)$node ) {
                            for (sample_property in colnames(clustering)[colnames(clustering) %in% columns_w_sample_ID_info] ) {
                                descends <- clustering[clustering$node %in% descendants(phylo, node),]
                                if ( table(duplicated(descends[,colnames(descends) == sample_property]))[1] == 1 ) {
                                    clustering[
                                        which(clustering$node == node),
                                        which(colnames(clustering) == sample_property)
                                    ] <- descends[,colnames(descends) == sample_property][1,1]
                                }
                            }
                        }
                    }

            # Return results

                return( clustering )

        }

##### descendants

    descendants <- function (phy, node, type = "t", ignore.tip = TRUE, labels = FALSE) {
            
            if (inherits(phy, "phylo")) {
                edge <- phy$edge
            }
            else {
                if (!is.matrix(phy)) {
                    stop("'phy' must be of classes 'phylo' or 'matrix'")
                }
                else {
                    edge <- phy
                    labels <- FALSE
                }
            }
            if (length(node) > 1) 
                stop("'node' must be vector of length 1")
            type <- match.arg(type, c("all", "daughter", "internal", 
                "terminal"))
            tips <- setdiff(edge[, 2], edge[, 1])
            if (node <= max(tips)) {
                if (ignore.tip) {
                    x <- node
                }
                else {
                    stop("node ", node, " is not an internal node")
                }
            }
            else {
                x <- edge[edge[, 1] == node, 2]
                if (type %in% c("internal", "terminal", "all")) {
                    repeat {
                        xx <- x
                        x <- sort(unique(c(x, edge[, 2][edge[, 1] %in% 
                          x])))
                        if (identical(x, xx)) 
                          break
                    }
                    if (type == "internal") 
                        x <- setdiff(x, tips)
                }
            }
            if (type == "terminal") {
                x <- intersect(x, tips)
                if (labels) {
                    x <- phy$tip.label[x]
                }
            }
            x
        }
