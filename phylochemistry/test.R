###### Libraries

    if ( !exists("packages") ) {

        ## Define necessary libraries
            
            CRAN_packages <- c(
                "stats",
                "dplyr",
                "shiny",
                "data.table",
                "ape",
                "rhandsontable",
                "ggplot2",
                "utils",
                "tidyr",
                "tibble",
                "readr",
                "ips",
                "FactoMineR",
                "remotes"
            )

            Bioconductor_packages <- c(
                "ggtree",
                "xcms"
            )

            Github_packages <- c(
                # "HajkD/orthologr"
            )

            packages_needed <- c(CRAN_packages, Bioconductor_packages, Github_packages)[!c(CRAN_packages, Bioconductor_packages, gsub(".*/", "", Github_packages)) %in% rownames(installed.packages())]

        ## Determine if anything needs to be installed
            
            if (length(packages_needed) > 0) {

                message <- paste0(
                    "You need to install the following packages before proceeding: ",
                    paste(packages_needed, collapse = ", "),
                    " Is it okay if phylochemistry installs them for you?"
                )

                if (.Platform$OS.type == "unix"){ response <- askYesNo( message ) }

                if (.Platform$OS.type == "windows") { response <- askYesNo("yesnocancel", msg = message) }

                if(response) {
                    
                    if (length(CRAN_packages[CRAN_packages %in% packages_needed]) > 0) {
                        install.packages(CRAN_packages[CRAN_packages %in% packages_needed], dependencies = TRUE)
                    }

                    if (length(Bioconductor_packages[Bioconductor_packages %in% packages_needed]) > 0) {
                        BiocManager::install(Bioconductor_packages[Bioconductor_packages %in% packages_needed], dependencies = TRUE)
                    }

                    if (length(Github_packages[Github_packages %in% packages_needed]) > 0) {
                        remotes::install_github(Github_packages[Github_packages %in% packages_needed], dependencies = TRUE)
                    }

                } else {
                    stop("Cannot load phylochemistry without the required packages. Exiting.")
                }
            }

            message("Loading packages...")
            invisible(suppressMessages(suppressWarnings(lapply(c(CRAN_packages, Bioconductor_packages), library, character.only = TRUE))))
    
    } else {
    
        message("Object 'packages' exists, not loading phylochemistry packages...")
    
    }

    ## Set up prioriy functions

        mutate <- dplyr::mutate
        summarize <- dplyr::summarize
        filter <- dplyr::filter
        select <- dplyr::select

    ## Set up options

        options(readr.show_progress = FALSE)
        options(dplyr.summarise.inform = FALSE)

#### Datasets

    message("Loading MS Library...")

    busta_spectral_library <- read_csv("https://thebustalab.github.io/R_For_Chemists_2/sample_data/busta_spectral_library_v1.csv", col_types = c(Compound_common_name = "c"))


#### writeMonolist

    #' Writes a monolist
    #'
    #' @param monolist The monolist to write out
    #' @param monolist_out_path The path to where the monolist should be written
    #' @examples
    #' @export
    #' writeMonolist

    writeMonolist <- function( monolist, monolist_out_path, ... ) {
            write.table(monolist, file = monolist_out_path, sep = ",", row.names = FALSE, ...)
    }

#### writeCSV

    #' Interactive selection of a CSV file to write
    #'
    #' @param monolist The monolist to write out
    #' @examples
    #' @export
    #' writeCSV

    writeCSV <- function(monolist) {
        writeMonolist(monolist = monolist, monolist_out_path = file.choose(new = TRUE))
    }

#### runMatrixAnalysis

    #' Runs a matrix analysis (clustering, kmeans, pca).
    #'
    #' @param data The data frame or tibble to analyze
    #' @param analysis The type of analysis to run. Can be one of: "hclust" (heirarchical clustering), "pca" (principal components analysis), "pca-ord" (principal components analysis ordination plot), or "pca-dim" (principal components scree plot)
    #' @param column_w_names_of_multiple_analytes  
    #' @param column_w_values_for_multiple_analytes
    #' @param columns_w_additional_analyte_info
    #' @param columns_w_sample_ID_info
    #' @param kmeans
    #' @param na_replacement
    #' @param output_format
    #' @examples
    #' @export
    #' runMatrixAnalysis

        runMatrixAnalysis <-    function(
                                    data,
                                    analysis = c(
                                        "hclust", "hclust_phylo",
                                        "pca", "pca_ord", "pca_dim",
                                        "mca", "mca_ord", "mca_dim",
                                        "mds", "mds_ord", "mds_dim",
                                        "tSNE", "DBSCAN"
                                    ),
                                    column_w_names_of_multiple_analytes = NULL,
                                    column_w_values_for_multiple_analytes = NULL,
                                    columns_w_values_for_single_analyte = NULL,
                                    columns_w_additional_analyte_info = NULL,
                                    columns_w_sample_ID_info = NULL,
                                    transpose = FALSE,
                                    distance_method = c("euclidean", "coeff_unlike"),
                                    unknown_sample_ID_info = NULL,
                                    scale_variance = TRUE,
                                    kmeans = c("none", "auto", "elbow", "1", "2", "3", "etc.", "pca"),
                                    na_replacement = c("none", "mean", "zero", "drop"),
                                    output_format = c("wide", "long"),
                                    ...
                                ) {

            # Pre-process data

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
                            stop("There are duplicate analyte names in columns_w_values_for_single_analyte.")
                        }
                    }

                    if( length(column_w_names_of_multiple_analytes) == 1 ) {
                        n_unique_samples <- table(duplicated(data[,colnames(data) %in% columns_w_sample_ID_info]))[1]

                        analyte_breakdown <- table(duplicated(data[,colnames(data) == column_w_names_of_multiple_analytes]))
                        n_unique_groups_of_analytes <- sum(analyte_breakdown)/analyte_breakdown[1]

                        if (n_unique_samples != n_unique_groups_of_analytes) {
                            stop("Either:\n1. The variables you have specified in columns_w_sample_ID_info do not define groups of samples with unique sets of analytes and you probably need to add more variables to columns_w_sample_ID_info.\n\nOR\n\n2. Not all your samples have the same number of analytes in column_w_names_of_multiple_analytes/column_w_values_of_multiple_analytes. Please fix that before continuing.")
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

                # Check to see if analyte columns are numeric

                    which_analyte_columns <- which(colnames(data_wide) %in% analyte_columns)
                    are_they_numeric <- list()
                    for( i in which_analyte_columns ) {
                      are_they_numeric <- c(are_they_numeric, is.numeric(data_wide[[i]]))
                    }

                    ## Convert all to numeric?
                        # for( i in which_analyte_columns ) {
                        #   data_wide[[i]] <- as.numeric(data_wide[[i]])
                        # }

                    # Should selected analysis proceed?

                        if ( all(unlist(are_they_numeric)) ) {
                            if ( analysis %in% c("pca", "pca_dim", "pca_ord", "hclust", "hclust_phylo") ) {
                                # cat("Analytes are all numeric and compatible with the analysis selected.\n")
                            }
                            if ( analysis %in% c("mca", "mca_ord", "mca_dim") ) {
                                stop("Analytes are all numeric, but the analysis selected is for categorical variables. Please choose a different analysis method.\n")
                            }
                        }

                        if ( !all(unlist(are_they_numeric)) ) {
                            if (analysis %in% c("mca", "mca_ord", "mca_dim")) {
                                # cat("Analytes are all categorical and compatible with the analysis selected.\n")
                            }
                            if ( analysis %in% c("pca", "pca_dim", "pca_ord", "hclust", "hclust_phylo") ) {
                                stop("Analytes are all categorical, but the analysis selected is for numeric variables. Please choose a different analysis method.\n")
                            }

                        }

                        # if (    all( 
                        #             c( !all(unlist(are_they_numeric)), all(unlist(are_they_numeric)) )
                        #         )
                        # ) {
                        #     if (analysis %in% c("mca")) {
                        #         cat("Analytes are mixed and ...")
                        #     }
                        # }

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

                # Prepare the Matrix for analysis - Clustering analysis

                    # Prepare the matrix

                        matrix <- as.data.frame(data_wide[,match(analyte_columns, colnames(data_wide))])
                        rownames(matrix) <- data_wide$sample_unique_ID

                    # Replace NAs with colmeans

                        if( na_replacement[1] == "none") {
                        }
                        if( na_replacement[1] == "drop" ) {
                            cat("Dropping any variables in your dataset that have NA as a value.\nVariables dropped:\n")
                            if (length(names(which(apply(is.na(matrix), 2, any)))) > 0) {
                                cat(names(which(apply(is.na(matrix), 2, any))))    
                            } else {
                                cat("none")
                            }
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

                    # Transpose matrix, if requested
                        if( transpose == TRUE ) { matrix <- t(matrix) }

                # Run unknown, if requested

                    if( length(unknown_sample_ID_info) > 0 ) {
                    
                        ## Use na_replacement != "drop"

                            if( na_replacement != "drop" ) {
                              stop("It is highly recommended that you use na_replacement = \"drop\" when matching an unknown, since anything with an NA value will be sent to the bottom of the matching list.")
                            }

                        ## Separate unknown and knowns matrix

                            sample_ID_of_unknown <- unknown_sample_ID_info
                            index_of_unknown <- which(rownames(matrix) == sample_ID_of_unknown)
                            unknown <- matrix[index_of_unknown,]
                            matrix_minus_unknown <- matrix[-c(index_of_unknown),]

                        ## Identify how many knowns there are to test

                            indices_to_test <- seq(1, dim(matrix_minus_unknown)[1], 1)

                        ## Loop and select 20 nearest neighbors until there are just 100 left

                            while( length(indices_to_test) > 100 ) { ## Has been optimized for 50, 100, and 1000 indeces bites. 100 and 1000 both take ~5s.

                                ## Select 100 entries at random, subset matrix_minus_unknown for that, append the unknowns
                                
                                    indices_being_tested <- sample(indices_to_test, 100, replace = FALSE)

                                    matrix_to_be_tested <- matrix_minus_unknown[indices_being_tested,]
                                    matrix_to_be_tested <- rbind(unknown, matrix_to_be_tested)
                                    dist <- Rfast::Dist(matrix_to_be_tested)
                                    results <- data.frame(
                                      distance = dist[2:dim(dist)[1], 1],
                                      index = indices_being_tested
                                    )

                                ## Keep closest 20 indices, modify indices_to_test

                                    indices_to_keep <- results[order(results$distance),]$index[1:20]
                                    indices_to_toss <- indices_being_tested[!indices_being_tested %in% indices_to_keep]
                                    indices_to_test <- indices_to_test[!indices_to_test %in% indices_to_toss]

                            }

                            ## Test the last 100 neighbors and return the closest 10

                                final_matrix_to_be_tested <- rbind(unknown, matrix_minus_unknown[indices_to_test,])
                                dist <- Rfast::Dist(final_matrix_to_be_tested)
                                results <- data.frame(
                                  distance = dist[2:dim(dist)[1], 1],
                                  index = indices_to_test
                            )

                            if( length(indices_to_test) < 10 ) {
                                number_for_final_tree <- length(indices_to_test)
                            } else {
                                number_for_final_tree <- 10
                            }

                            indices_to_keep <- results[order(results$distance),]$index[1:number_for_final_tree]
                            indices_to_keep

                        ## Subset the matrix

                            matrix <- rbind(unknown, matrix_minus_unknown[indices_to_keep,])

                            ### THINK CAREFULLY HERE! THINGS SOMEWHAT CLOSE TO THE UNKNOWN 
                            ### MIGHT NOT CLUSTER DIRECTLY WITH IT
                            ### DEPENDING ON HOW RELATED THEY ARE TO OTHER THINGS

                }

            # Run the matrix analysis selected

                ## Scale data if requested

                    if( scale_variance == TRUE ) {
                        matrix <- scale(matrix)
                    }

                ## Distance-based methods. First, get distance matrix:

                    if( distance_method[1] == "euclidean" ) {
                        dist_matrix <- stats::dist(matrix, method = "euclidean")
                    } else {
                        stop("Please specify distance method")
                    }

                    ## HCLUST, HCLUST_PHYLO ##

                        if( analysis == "hclust" ) {
                            phylo <- ape::as.phylo(stats::hclust(dist_matrix))
                            clustering <- ggtree::fortify(phylo)
                            clustering$sample_unique_ID <- clustering$label
                        }

                        if( analysis == "hclust_phylo" ) {                        
                            phylo <- ape::as.phylo(stats::hclust(dist_matrix))
                            return(phylo)
                            stop("Returning hclust_phylo.")
                        }

                    ## MDS

                        if( analysis == "mds" ) {
                            coords <- stats::cmdscale(dist_matrix)
                            colnames(coords) <- c("Dim.1", "Dim.2")
                            clustering <- as_tibble(coords)
                            clustering$sample_unique_ID <- rownames(coords)
                        }

                ## Non distance-based methods
                
                    ## tSNE

                        # library(Rtsne)
                        # out <- Rtsne(matrix, theta = 0.0, perplexity = 2)
                        # out <- out$Y
                        # colnames(out) <- c("x","y")
                        # out <- as.data.frame(out)
                        # ggplot(out) +
                        #     geom_point(aes(x = x, y = y), shape = 21, size= 4)

                    ## MCA, MCA_ORD, MCA_DIM ##

                        if( analysis == "mca" ) {
                            cat("Running Multiple Correspondence Analysis, extracting sample coordinates...\n")
                            coords <- FactoMineR::MCA(matrix, graph = FALSE)$ind$coord[,c(1:2)]
                            clustering <- as_tibble(coords)
                            clustering$sample_unique_ID <- rownames(coords)
                            colnames(clustering) <- c("Dim_1", "Dim_2", "sample_unique_ID")
                            cat("Done!\n")
                        }

                        if( analysis == "mca_ord" ) {
                            cat("Running Multiple Correspondence Analysis, extracting ordination plot...\n")
                            coords <- FactoMineR::MCA(matrix, graph = FALSE)$var$coord[,c(1,2)]
                            clustering <- as_tibble(coords)
                            clustering$analyte <- rownames(coords)
                            colnames(clustering) <- c("Dim_1", "Dim_2", "analyte")
                            clustering <- select(clustering, analyte, Dim_1, Dim_2)
                            return(clustering)
                            stop("Returning ordination plot coordinates. \nDone!")
                        }

                        if( analysis == "mca_dim" ) {
                            cat("Running Multiple Correspondence Analysis, extracting dimensional contributions...\n")
                            coords <- FactoMineR::MCA(matrix, graph = FALSE)$eig[,2]
                            clustering <- tibble::enframe(coords, name = NULL)
                            clustering$principal_component <- names(coords)
                            clustering$principal_component <- as.numeric(gsub("dim ", "", clustering$principal_component))
                            colnames(clustering)[colnames(clustering) == "value"] <- "percent_variance_explained"
                            clustering <- select(clustering, principal_component, percent_variance_explained)
                            return(clustering)
                            stop("Returning eigenvalues. \nDone!")
                        }

                ## PCA, PCA_ORD, PCA_DIM ## 

                    if( analysis == "pca" ) {
                        coords <- FactoMineR::PCA(matrix, graph = FALSE, scale.unit = FALSE)$ind$coord[,c(1:2)]
                        clustering <- as_tibble(coords)
                        clustering$sample_unique_ID <- rownames(coords)
                    }

                    if( analysis == "pca_ord" ) {
                        coords <- FactoMineR::PCA(matrix, graph = FALSE, scale.unit = FALSE)$var$coord[,c(1,2)]
                        clustering <- as_tibble(coords)
                        clustering$analyte <- rownames(coords)
                        clustering <- select(clustering, analyte, Dim.1, Dim.2)
                        return(clustering)
                        stop("Returning ordination plot coordinates.")
                    }

                    if( analysis == "pca_dim" ) {
                        coords <- FactoMineR::PCA(matrix, graph = FALSE, scale.unit = FALSE)$eig[,2]
                        clustering <- tibble::enframe(coords, name = NULL)
                        clustering$principal_component <- names(coords)
                        clustering$principal_component <- as.numeric(gsub("comp ", "", clustering$principal_component))
                        colnames(clustering)[colnames(clustering) == "value"] <- "percent_variance_explained"
                        clustering <- select(clustering, principal_component, percent_variance_explained)
                        return(clustering)
                        stop("Returning eigenvalues.")
                    }

                
            # Clustering

                ## Density-based clustering = DBSCAN and OPTICS

                    # out <- dbscan::optics(matrix, minPts = 5)
                    # out <- data.frame(
                    #     order = out$order,
                    #     reach_dist = out$reachdist,
                    #     name = rownames(matrix)
                    # )
                    # out$name <- factor(out$name, levels = rev(rownames(matrix)[out$order]))
                    # ggplot(out[2:19,]) +
                    #     geom_col(aes(x = name, y = reach_dist))

                ## K-MEANS ##

                    if( kmeans[1] %in% c("none", FALSE) ) {}

                    if( !kmeans[1] %in% c("none", FALSE) ) {

                        ## Check for NAs

                            if( any(is.na(matrix)) == TRUE ) {
                                stop("kmeans cannot handle NA. Please choose an option for na_replacement.")
                            }

                        ## If kmeans = "pca", substitute the matrix kmeans sees with the first two dimensions of the PCA

                            if ( kmeans[1] == "pca") {
                                matrix <- as.matrix(clustering[,1:2])
                                rownames(matrix) <- as.character(as.data.frame(clustering)[,3])
                            }

                        ## Run k-means

                            kmeans_results <- list()
                            for( i in 1:(dim(matrix)[1]-1) ) {
                                kmeans_results[[i]] <- sum(stats::kmeans(x = matrix, centers = i, nstart = 25, iter.max = 1000)$withinss)
                            }
                            kmeans_results <- do.call(rbind, kmeans_results)
                            kmeans_results

                        ## If auto, determine sharpest angle and return clusters

                            if( kmeans[1] == "auto" | kmeans[1] == "pca" ) {
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

                            if( !kmeans[1] %in% c("auto", "elbow", "pca") ) {
                                if( is.na(as.numeric(kmeans[1])) ) {stop("For kmeans, please use 'none', 'auto', 'elbow', 'pca' or a number.")}
                                kmeans_clusters <- stats::kmeans(x = matrix, centers = as.numeric(kmeans[1]), nstart = 25, iter.max = 1000)$cluster
                            }

                        ## Bind kmeans cluster numbers to output

                            kmeans_clusters <- as_tibble(data.frame(sample_unique_ID = names(kmeans_clusters), kmeans_cluster = paste0("cluster_", kmeans_clusters)))
                            clustering$kmeans_cluster <- kmeans_clusters$kmeans_cluster[match(clustering$sample_unique_ID, kmeans_clusters$sample_unique_ID)]

                    }

            # If transpose = TRUE, then return without annotations

                if( transpose == TRUE ) {
                    return(clustering)
                    stop("Returning transposed cluster output. Make sure all your variables have the same units!")
                }
                
            # Post processing and return.

                if( !analysis %in% c("hclust_phylo")) {

                    ## Add back annotations to the output

                        if( length(columns_w_sample_ID_info) == 1 ) {
                        } else {
                        clustering <-   right_join(
                                            data_wide[,match(
                                                c(columns_w_sample_ID_info, "sample_unique_ID"),
                                                colnames(data_wide))
                                            ], clustering, by = "sample_unique_ID"
                                        )
                        }

                        rownames_matrix <- tibble::enframe(rownames(matrix), name = NULL)
                        colnames(rownames_matrix)[1] <- "sample_unique_ID"

                        # if (analysis != "pca") { ## Don't do this for pca for some reason?? I don't understand why...
                            clustering <- full_join(
                                clustering,
                                as_tibble(cbind(rownames_matrix, as_tibble(matrix))),
                                by = "sample_unique_ID"
                            )
                        # }
                        clustering

                        ## Order the returned matrix so that the sample_unique_ID comes first

                            clustering <- select(clustering, sample_unique_ID, everything())

                    # Annotate internal branches if tree output

                        if( analysis == "hclust" ) {
                            for( node in dplyr::filter(clustering, isTip == FALSE)$node ) {
                                for (sample_property in colnames(clustering)[colnames(clustering) %in% columns_w_sample_ID_info] ) {
                                    descends <- clustering[clustering$node %in% ips::descendants(phylo, node),]
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

                        if( output_format[1] == "long" ) {
                            clustering <- pivot_longer(
                                clustering,
                                cols = c(which(colnames(clustering) == analyte_columns[1]): dim(clustering)[2]),
                                names_to = "analyte_name", 
                                values_to = "value"
                            )
                        }

                        return( clustering )
                }
        }

#### readMonolist

    #' Reads a monolist
    #'
    #' @param monolist_in_path The path to the monolist (in .csv format) to be read
    #' @examples
    #' @export
    #' readMonolist

    readMonolist <- function( monolist_in_path ) {
        monolist <- as.data.frame(data.table::fread(file = monolist_in_path))
        return( monolist )
    }

#### normalize

    #' Normalizes a vector of numbers to a range of zero to one.
    #'
    #' @param x The vector to normalize
    #' @param old_min The minimum of the old range
    #' @param old_max The maximum of the old range
    #' @param new_min The minimum of the new range, defaults to 0
    #' @param new_max The maximum of the new range, defaults to 1
    #' @examples
    #' @export
    #' normalize

    normalize <- function( x, old_min = NULL, old_max = NULL, new_min = 0, new_max = 1 ) {

        if ( length(old_min) == 0 & length(old_max) == 0 ) {
            
            (x - min(x)) * (new_max - new_min) / (max(x) - min(x)) + new_min                
        
        } else {

            (x - old_min) * (new_max - new_min) / (old_max - old_min) + new_min

        }

    }

#### dropNA

    #' drops NAs from a vector
    #'
    #' @param x The vector to drop NAs from 
    #' @examples
    #' @export
    #' normalize

    dropNA <- function( x ) {

        return(x[!is.na(x)])

    }

#### integrationAppLite2

    #' A Shiny app to integrate GC-FID and GC-MS data
    #'
    #' @param chromatograms A data frame containing columns: "rt", "tic", and "path_to_cdf_csv", which contain retention time, total ion chromatogram intensities, and paths to CDF.csv files generated by the convertCDFstoCSVs function.
    #' @param x_axis_start A numeric value for the lower x-axis bounds on the plot generated by the app. Defaults to full length.
    #' @param x_axis_end A numeric value for the upper x-axis bounds on the plot generated by the app. Defaults to full length.
    #' @param samples_monolist_path A path to a .csv file containing metadata for the samples you wish to analyze. Requied columns are: "rt_offset", "baseline_window", and "path_to_cdf_csv", which are for aligning chromatograms, adjusting baseline determination, and defining the path to the CDF.csv files for each sample, respectively.
    #' @param samples_monolist_subset Optional, a numeric vector (for example, "c(1:10)"), defining a subset of samples to be loaded.
    #' @param peaks_monolist_path A path to a .csv file containing metadata for all peaks in the sample set. Required columns are: peak_start", "peak_end", "path_to_cdf_csv", "area", "peak_number_within_sample", "rt_offset", "peak_start_rt_offset", "peak_end_rt_offset". This file is automatically generated by the app.
    #' @param zoom_and_scroll_rate Defines intervals of zooming and scrolling movement while running the app
    #' @examples
    #' @export
    #' integrationAppLite2

    integrationAppLite2 <- function(
            CDF_directory_path = getwd(),
            zoom_and_scroll_rate = 100,
            baseline_window = 400,
            x_axis_start_default = NULL,
            x_axis_end_default = NULL,
            path_to_reference_library = busta_spectral_library,
            samples_monolist_subset = NULL,
            ions = "tic"
        ) {

            setwd(CDF_directory_path)
            paths_to_cdfs <- dir()[grep("*.CDF$", dir())]
            paths_to_cdf_csvs <- paste0(paths_to_cdfs, ".csv")

            ## PREPARE DATA: Check for CDF to CSV conversion, check chromatograms

                if ( length(paths_to_cdfs) == 0 ) {

                    stop("The directory specified does not contain any .CDF files.")

                } else {

                    if ( !file.exists("chromatograms.csv") ) {

                        chromatograms <- list()

                    } else {

                        chromatograms <- readMonolist("chromatograms.csv")

                    }

                    chromatograms_to_add <- list()

                    for (file in 1:length(paths_to_cdfs)) {

                        ## If the cdf.csv doesn't exist for this cdf, create it.

                            if ( !file.exists(paths_to_cdf_csvs[file]) ) {

                                cat(paste("CDF to CSV conversion. Reading data file ", paths_to_cdfs[file], "\n", sep = ""))
                                    rawDataFile <- xcms::loadRaw(xcms::xcmsSource(paths_to_cdfs[file]))

                                cat("   Framing data file ... \n")
                                    rt <- rawDataFile$rt
                                    scanindex <- rawDataFile$scanindex

                                    filteredRawDataFile <- list()
                                    for ( i in 1:(length(rt)-1) ) {
                                        filteredRawDataFile[[i]] <- data.frame(
                                            mz = rawDataFile$mz[(scanindex[i]+1):(scanindex[i+1])],
                                            intensity = rawDataFile$intensity[(scanindex[i]+1):(scanindex[i+1])],
                                            rt = rt[i]
                                        )
                                    }
                                    framedDataFile <- do.call(rbind, filteredRawDataFile)
                                    framedDataFile$mz <- round(framedDataFile$mz, digits = 1)

                                cat("   Merging duplicate rows ...\n")
                                    if ( dim(table(duplicated(paste(framedDataFile$mz, framedDataFile$rt, sep = "_")))) > 1 ) {
                                        framedDataFile %>% group_by(mz,rt) %>% summarize(intensity = sum(intensity), .groups = "drop") -> framedDataFile
                                        framedDataFile <- as.data.frame(framedDataFile)
                                    }

                                cat("   Writing out data file as CSV... \n")
                                    data.table::fwrite(framedDataFile, file = paste(paths_to_cdfs[file], ".csv", sep = ""), col.names = TRUE, row.names = FALSE)

                            }

                        ## If any chromatograms (tic and ion) are not present for this csv, extract them

                            if ( file.exists("chromatograms.csv") ) {

                                ions_for_this_cdf_csv <- unique(filter(readMonolist("chromatograms.csv"), path_to_cdf_csv == paths_to_cdf_csvs[file])$ion)
                                missing_ions <- as.numeric(as.character(ions[!ions %in% ions_for_this_cdf_csv]))
                                missing_ions <- dropNA(missing_ions)

                            } else {

                                missing_ions <- ions

                            }

                            if (length(missing_ions) > 0) {
                                
                                cat(paste("Chromatogram extraction. Reading data file ", paths_to_cdf_csvs[file], "\n", sep = ""))    
                                    framedDataFile <- as.data.frame(data.table::fread(paths_to_cdf_csvs[file]))
                                
                                cat("   Extracting chromatograms...\n")
                                    
                                    if ("tic" %in% ions) {
                                        framedDataFile$row_number <- seq(1,dim(framedDataFile)[1],1)
                                        framedDataFile %>% 
                                            group_by(rt) %>% summarize(
                                            abundance = sum(intensity),
                                            ion = "tic",
                                            rt_first_row_in_raw = min(row_number),
                                            rt_last_row_in_raw = max(row_number)
                                        ) -> chromatogram
                                        chromatogram <- as.data.frame(chromatogram)
                                        chromatogram$rt <- as.numeric(chromatogram$rt)
                                        chromatogram$path_to_cdf_csv <- paste(paths_to_cdfs[file], ".csv", sep = "")
                                        chromatograms_to_add <- rbind(chromatograms_to_add, chromatogram)
                                    }

                                    if ( length(ions[ions != "tic"]) > 0 ) {

                                        numeric_ions <- as.numeric(as.character(ions[ions != "tic"]))
                                        for ( ion in 1:length(numeric_ions) ){
                                            framedDataFile$row_number <- seq(1,dim(framedDataFile)[1],1)
                                            framedDataFile %>% 
                                                group_by(rt) %>% 
                                                filter(mz > (numeric_ions[ion] - 0.6)) %>%
                                                filter(mz < (numeric_ions[ion] + 0.6)) %>%
                                                summarize(
                                                    abundance = sum(intensity),
                                                    ion = numeric_ions[ion],
                                                    rt_first_row_in_raw = min(row_number),
                                                    rt_last_row_in_raw = max(row_number)
                                                ) -> chromatogram
                                            chromatogram <- as.data.frame(chromatogram)
                                            chromatogram$rt <- as.numeric(chromatogram$rt)
                                            chromatogram$path_to_cdf_csv <- paste(paths_to_cdfs[file], ".csv", sep = "")
                                            chromatograms_to_add <- rbind(chromatograms_to_add, chromatogram)   
                                        }
                                    }
                            }

                        ## If the chromatograms file already exists, append to it and re-write out, else create it

                            if ( file.exists("chromatograms.csv") ) {
                                
                                print("writing it")

                                writeMonolist(
                                    monolist = rbind( chromatograms, chromatograms_to_add ),
                                    monolist_out_path = "chromatograms.csv"
                                )

                            } else {

                                writeMonolist(chromatograms_to_add, "chromatograms.csv")

                            }
                    }

                    print("done")
                }

                ## Read in chromatograms

                    chromatograms <- readMonolist("chromatograms.csv")
                            
                ## Set up new samples monolist

                    ## If it doesn't exist, create it
                    
                        if ( !file.exists("samples_monolist.csv") ) {

                            samples_monolist <- data.frame(
                                Sample_ID = unique(chromatograms$path_to_cdf_csv),
                                rt_offset = 0,
                                baseline_window = baseline_window,
                                path_to_cdf_csv = unique(chromatograms$path_to_cdf_csv)
                            )

                            write.table(
                                x = samples_monolist,
                                file = "samples_monolist.csv",
                                row.names = FALSE,
                                sep = ","
                            )

                    # If it exists, check to see if all cdfs in this folder are in it, if not, add them

                        } else {

                            samples_monolist <- readMonolist("samples_monolist.csv")

                            missing_from_samples_monolist <- paths_to_cdf_csvs[!paths_to_cdf_csvs %in% unique(readMonolist("samples_monolist.csv")$path_to_cdf_csv)]

                            if ( length(missing_from_samples_monolist) > 0 ) {

                                samples_monolist_additions <- data.frame(
                                    Sample_ID = gsub("\\..*$", "", gsub(".*/", "", missing_from_samples_monolist)),
                                    rt_offset = 0,
                                    baseline_window = baseline_window,
                                    path_to_cdf_csv = missing_from_samples_monolist
                                )

                                write.table(
                                    x = samples_monolist_additions,
                                    file = "samples_monolist.csv",
                                    row.names = FALSE,
                                    col.names = FALSE,
                                    sep = ",",
                                    append = TRUE
                                )
                            }
                        }

                ## Filter chromatograms so only the CDFs in this folder are included

                    chromatograms <- chromatograms[chromatograms$path_to_cdf_csv %in% dir()[grep(".CDF.csv", dir())],]

                ## Set up several variables, plot_height, and x_axis limits if not specified in function call
                    
                    peak_data <- NULL
                    peak_points <- NULL
                    plot_height <- 200 + 100*length(unique(chromatograms$path_to_cdf_csv))
                    
                ## Set up new peak monolist if it doesn't exist
                
                if ( !file.exists("peaks_monolist.csv") ) {
                    
                    peak_data <- data.frame(
                      peak_start = 0,
                      peak_end = 0,
                      peak_ID = "unknown",
                      path_to_cdf_csv = "a",
                      area = 0
                    )

                    write.table(
                      x = peak_data[-1,],
                      file = "peaks_monolist.csv",
                      append = FALSE,
                      row.names = FALSE,
                      col.names = TRUE,
                      sep = ","
                    )

                }

            ## SET UP USER INTERFACE

                ui <- fluidPage(

                    tags$script('
                    $(document).on("keypress", function (e) {
                       Shiny.onInputChange("keypress", e.which);
                    });
                    '), 

                    tabsetPanel(type = "tabs",

                        tabPanel("Main",

                            verticalLayout(

                                plotOutput(
                                    outputId = "massSpectra_1",
                                    brush = brushOpts(
                                        id = "massSpectra_1_brush"
                                    ),
                                    height = 150
                                ),

                                plotOutput(
                                    outputId = "massSpectra_2",
                                    brush = brushOpts(
                                        id = "massSpectra_2_brush"
                                    ),
                                    height = 150
                                ),

                                plotOutput(
                                    outputId = "chromatograms",
                                    brush = brushOpts(
                                        id = "chromatogram_brush"
                                    ),
                                  click = "chromatogram_click", 
                                  dblclick = "chromatogram_double_click",
                                  height = plot_height
                                ),

                                verbatimTextOutput("key", placeholder = TRUE),

                                rhandsontable::rHandsontableOutput("peak_table")
                            )
                        ),

                        tabPanel("MS Library",

                            verticalLayout(
                                
                                plotOutput(
                                    outputId = "massSpectrumLookup",
                                    height = 1200
                                )

                            )
                        )
                    )
                )

            ## SET UP SERVER

                server <- function(input, output, session) {

                    ## Check keystoke value
                        output$key <- renderPrint({
                            input$keypress
                        })

                    ## Keys to move chromatogram view - zoom in and out, move L and R
                        observeEvent(input$keypress, {
                            if( input$keypress == 70 ) { x_axis_start_default <<- x_axis_start_default + zoom_and_scroll_rate } # Forward on "F"
                            if( input$keypress == 70 ) { x_axis_end_default <<- x_axis_end_default + zoom_and_scroll_rate } # Forward on "F"
                            if( input$keypress == 68 ) { x_axis_start_default <<- x_axis_start_default - zoom_and_scroll_rate } # Backward on "D"
                            if( input$keypress == 68 ) { x_axis_end_default <<- x_axis_end_default - zoom_and_scroll_rate } # Backward on "D"
                            if( input$keypress == 86 ) { x_axis_start_default <<- x_axis_start_default - zoom_and_scroll_rate } # Wider on "V"
                            if( input$keypress == 86 ) { x_axis_end_default <<- x_axis_end_default + zoom_and_scroll_rate } # Wider on "V"
                            if( input$keypress == 67 ) { x_axis_start_default <<- x_axis_start_default + zoom_and_scroll_rate } # Closer on "C"
                            if( input$keypress == 67 ) { x_axis_end_default <<- x_axis_end_default - zoom_and_scroll_rate } # Closer on "C"
                        })

                    ## Save manual changes to table on "Z" (90) keystroke

                        observeEvent(input$keypress, {

                            if (input$keypress == 90 ) {

                                ## Write out any modifications to peak table (i.e. sample IDs)
                                    
                                    hot = isolate(input$peak_table)
                                    if (!is.null(hot)) {
                                        writeMonolist(rhandsontable::hot_to_r(input$peak_table), "peaks_monolist.csv")
                                        cat("Peak list saved!\n")
                                    }

                            }
                        
                        })

                    ## Update chromatogram on "Q" (81) keystroke
                        
                        observeEvent(input$keypress, {      
                            
                            if( input$keypress == 81 ) { # Update on "Q"

                                ## Read in samples monolist and put chromatograms into chromatograms_updated
                                    
                                    samples_monolist <- read.csv("samples_monolist.csv")
                                    if ( length(samples_monolist_subset) > 0 ) {
                                        samples_monolist <- samples_monolist[samples_monolist_subset,]    
                                    }
                                    chromatograms_updated <- dplyr::filter(chromatograms, path_to_cdf_csv %in% samples_monolist$path_to_cdf_csv)

                                ## Calculate baseline for each sample

                                    baselined_chromatograms <- list()

                                    for ( chrom in 1:length(unique(chromatograms_updated$path_to_cdf_csv)) ) {
                              
                                        chromatogram <- dplyr::filter(chromatograms_updated, path_to_cdf_csv == unique(chromatograms_updated$path_to_cdf_csv)[chrom])
                                        tic <- filter(chromatogram, ion == "tic")

                                        prelim_baseline_window <- samples_monolist$baseline_window[match(chromatogram$path_to_cdf_csv[1], samples_monolist$path_to_cdf_csv)]

                                        n_prelim_baseline_windows <- floor(length(tic$rt)/prelim_baseline_window)
                                        prelim_baseline <- list()
                                        for ( i in 1:n_prelim_baseline_windows ) {
                                            abundances_in_window <- tic$abundance[((prelim_baseline_window*(i-1))+1):(prelim_baseline_window*i)]
                                            prelim_baseline[[i]] <- data.frame(
                                                rt = tic$rt[(which.min(abundances_in_window)+((i-1)*prelim_baseline_window))],
                                                min = min(abundances_in_window)
                                            )
                                        }
                                        prelim_baseline <- do.call(rbind, prelim_baseline)
                                        tic$in_prelim_baseline <- FALSE
                                        tic$in_prelim_baseline[tic$rt %in% prelim_baseline$rt] <- TRUE

                                        y = prelim_baseline$min
                                        x = prelim_baseline$rt

                                        baseline2 <- data.frame(
                                            rt = tic$rt,
                                            y = approx(x, y, xout = tic$rt)$y
                                        )
                                        baseline2 <- baseline2[!is.na(baseline2$y),]
                                        tic <- tic[tic$rt %in% baseline2$rt,]
                                        tic$baseline <- baseline2$y

                                        baselined_chromatograms[[chrom]] <- data.frame(
                                            rt = tic$rt,
                                            abundance = tic$baseline,
                                            ion = "baseline",
                                            path_to_cdf_csv = tic$path_to_cdf_csv,
                                            rt_first_row_in_raw = tic$rt_first_row_in_raw,
                                            rt_last_row_in_raw = tic$rt_last_row_in_raw
                                        )
                                    }

                                    baselined_chromatograms <- do.call(rbind, baselined_chromatograms)
                                    chromatograms_updated <- rbind(chromatograms, baselined_chromatograms)

                                ## Add rt offset information for all chromatograms

                                    chromatograms_updated$rt_offset <- samples_monolist$rt_offset[match(chromatograms_updated$path_to_cdf_csv, samples_monolist$path_to_cdf_csv)]
                                    chromatograms_updated$rt_rt_offset <- chromatograms_updated$rt + chromatograms_updated$rt_offset
                                    chromatograms_updated <<- chromatograms_updated

                                ## Subset x_axis according to selection in chromatogram

                                    ## If null from initial start up, assign extreme values

                                        if ( is.null(x_axis_start_default) ) {
                                            x_axis_start_default <<- min(chromatograms$rt)
                                            cat(paste("x_axis_start_default is ", x_axis_start_default, "\n"))
                                        }

                                        if ( is.null(x_axis_end_default) ) {
                                            x_axis_end_default <<- max(chromatograms$rt)
                                            cat(paste("x_axis_end_default is ", x_axis_end_default, "\n"))
                                        }

                                    ## If brush is null, assign default values to start and end

                                        if ( is.null(input$chromatogram_brush) ) {
                                            x_axis_start <<- x_axis_start_default
                                            x_axis_end <<- x_axis_end_default
                                        }

                                    ## If brush is not null, assign brush values to start and end

                                        if ( !is.null(input$chromatogram_brush) ) {
                                            peak_points <<- isolate(brushedPoints(chromatograms_updated, input$chromatogram_brush))
                                            x_axis_start <<- min(peak_points$rt)
                                            x_axis_end <<- max(peak_points$rt)
                                        }
                                    
                                    ## Filter chromatogram
                                        
                                        chromatograms_updated <- dplyr::filter(
                                            chromatograms_updated, rt_rt_offset > x_axis_start & rt_rt_offset < x_axis_end
                                        )

                                ## Plot
                                        
                                    facet_labels <- gsub(".CDF.csv", "", gsub(".*/", "", chromatograms_updated$path_to_cdf_csv))
                                    names(facet_labels) <- chromatograms_updated$path_to_cdf_csv

                                    chromatogram_plot <- ggplot() +
                                        geom_line(
                                            data = filter(chromatograms_updated, ion == "baseline"),
                                            mapping = aes(x = rt_rt_offset, y = abundance), color = "grey"
                                        ) +
                                        # geom_line(
                                        #     data = filter(chromatograms_updated, ion != "baseline"),
                                        #     mapping = aes(x = rt_rt_offset, y = abundance, color = ion),
                                        #     alpha = 0.8
                                        # ) +
                                        scale_x_continuous(limits = c(x_axis_start, x_axis_end), name = "Retention (Scan number)") +
                                        scale_y_continuous(name = "Abundance (counts)") +
                                        facet_grid(path_to_cdf_csv~., scales = "free_y", labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                        theme_classic() +
                                        guides(fill = "none") +
                                        scale_fill_continuous(type = "viridis") +
                                        scale_color_brewer(palette = "Paired")

                                ## Add peaks, if any
                                    
                                    peak_table <- read.csv("peaks_monolist.csv")
                            
                                    if (dim(peak_table)[1] > 0) {

                                        ## Filter out duplicate peaks and NA peaks
                                            
                                            peak_table <- peak_table %>% group_by(path_to_cdf_csv) %>% mutate(duplicated = duplicated(peak_start))
                                            peak_table <- as.data.frame(peak_table)
                                            peak_table <- dplyr::filter(peak_table, duplicated == FALSE)
                                            peak_table <- peak_table %>% group_by(path_to_cdf_csv) %>% mutate(duplicated = duplicated(peak_end))
                                            peak_table <- as.data.frame(peak_table)
                                            peak_table <- dplyr::filter(peak_table, duplicated == FALSE)
                                            peak_table <- peak_table[,!colnames(peak_table) == "duplicated"]
                                            peak_table <- peak_table[!is.na(peak_table$peak_start),]

                                        ## Update with peak_number_within_sample
                                            
                                            peak_table_updated <- list()
                                            for (sample_number in 1:length(unique(peak_table$path_to_cdf_csv))) {
                                              peaks_in_this_sample <- peak_table[peak_table$path_to_cdf_csv == unique(peak_table$path_to_cdf_csv)[sample_number],]
                                              peaks_in_this_sample <- peaks_in_this_sample[order(peaks_in_this_sample$peak_start),]
                                              peaks_in_this_sample$peak_number_within_sample <- seq(1,length(peaks_in_this_sample$path_to_cdf_csv),1)
                                              peak_table_updated[[sample_number]] <- peaks_in_this_sample
                                            }
                                            peak_table_updated <- do.call(rbind, peak_table_updated)
                                            peak_table <- peak_table_updated

                                        ## Modify peaks with RT offset
                                            
                                            # for (sample_number in 1:length(unique(samples_monolist$path_to_cdf_csv))) {
                                              
                                            #   peaks_in_this_sample <- peak_table[peak_table$path_to_cdf_csv == samples_monolist$path_to_cdf_csv[sample_number],]
                                              
                                            #   rt_offsets <- samples_monolist$rt_offset[match(peaks_in_this_sample$path_to_cdf_csv, samples_monolist$path_to_cdf_csv)]
                                            #   peak_start_rt_offsets <- peak_table$peak_start + peak_table$rt_offset
                                            #   peak_end_rt_offsets <- peak_table$peak_end + peak_table$rt_offset
                                              
                                            #   peak_table$rt_offset[peak_table$path_to_cdf_csv == as.character(samples_monolist$path_to_cdf_csv[sample_number])] <- rt_offsets
                                            #   peak_table$peak_start_rt_offset[peak_table$path_to_cdf_csv == as.character(samples_monolist$path_to_cdf_csv[sample_number])] <- peak_start_rt_offsets
                                            #   peak_table$peak_end_rt_offset[peak_table$path_to_cdf_csv == as.character(samples_monolist$path_to_cdf_csv[sample_number])] <- peak_end_rt_offsets

                                            # }

                                            peak_table$rt_offset <- samples_monolist$rt_offset[match(peak_table$path_to_cdf_csv, samples_monolist$path_to_cdf_csv)]
                                            peak_table$peak_start_rt_offset <- peak_table$peak_start + peak_table$rt_offset
                                            peak_table$peak_end_rt_offset <- peak_table$peak_end + peak_table$rt_offset
                                            peak_table$path_to_cdf_csv <- as.character(peak_table$path_to_cdf_csv)
                                
                                        ## Update all peak areas in case baseline was adjusted
                                            
                                            for (sample_number in 1:length(unique(samples_monolist$path_to_cdf_csv))) {
                                              
                                              peaks_in_this_sample <- peak_table[peak_table$path_to_cdf_csv == samples_monolist$path_to_cdf_csv[sample_number],]
                                              
                                              areas <- vector()
                                              for (peak in 1:length(peaks_in_this_sample$peak_number_within_sample)) {
                                                areas <- append(areas, 
                                                  sum(dplyr::filter(
                                                    chromatograms_updated[chromatograms_updated$path_to_cdf_csv == as.character(peaks_in_this_sample$path_to_cdf_csv[peak]),], 
                                                    rt >= peaks_in_this_sample$peak_start[peak] & rt <= peaks_in_this_sample$peak_end[peak],
                                                    ion == "tic")$abundance
                                                  ) - 
                                                  sum(dplyr::filter(
                                                    chromatograms_updated[chromatograms_updated$path_to_cdf_csv == as.character(peaks_in_this_sample$path_to_cdf_csv[peak]),], 
                                                    rt >= peaks_in_this_sample$peak_start[peak] & rt <= peaks_in_this_sample$peak_end[peak])$baseline
                                                  )
                                                )
                                              }

                                              peak_table$area[peak_table$path_to_cdf_csv == as.character(samples_monolist$path_to_cdf_csv[sample_number])] <- areas

                                              # peaks_in_this_sample$area <- areas
                                              # peak_table_updated[[sample_number]] <- peaks_in_this_sample
                                            }
                                            # peak_table_updated <- do.call(rbind, peak_table_updated)
                                            # peak_table <- peak_table_updated
                              
                                        ## Write out peaks now with assigned peak_number_within_sample and RT offset, update the peak_table in ui
                                
                                            write.table(peak_table, file = "peaks_monolist.csv", col.names = TRUE, sep = ",", row.names = FALSE)

                                            output$peak_table <- rhandsontable::renderRHandsontable(rhandsontable::rhandsontable({
                                                peak_table2 <- read.csv("peaks_monolist.csv")
                                                peak_table2
                                            }))

                                        ## Add peaks

                                            print(x_axis_start)
                                            print(x_axis_end)
                                            if (length(x_axis_start) == 0) {x_axis_start <<- min(chromatograms$rt)}
                                            if (length(x_axis_end) == 0) {x_axis_end <<- max(chromatograms$rt)}
                                            print(x_axis_start)
                                            print(x_axis_end)

                                            peak_table <- dplyr::filter(peak_table, peak_start_rt_offset > x_axis_start & peak_end_rt_offset < x_axis_end)
                                            print("filter passed")

                                            for (peak in 1:dim(peak_table)[1]) {
                                                
                                                signal_for_this_peak <- dplyr::filter(
                                                    chromatograms_updated[chromatograms_updated$path_to_cdf_csv == peak_table[peak,]$path_to_cdf_csv,], 
                                                    rt_rt_offset > peak_table[peak,]$peak_start_rt_offset, 
                                                    rt_rt_offset < peak_table[peak,]$peak_end_rt_offset
                                                )

                                                if (dim(signal_for_this_peak)[1] > 0) {

                                                    signal_for_this_peak$peak_number_within_sample <- peak_table$peak_number_within_sample[peak]
                                                    
                                                    ribbon <- filter(signal_for_this_peak, ion == "tic")
                                                    ribbon$baseline <- filter(signal_for_this_peak, ion == "baseline")$abundance

                                                    chromatogram_plot <- chromatogram_plot +
                                                        geom_vline(data = signal_for_this_peak[1,], mapping = aes(xintercept = rt_rt_offset), alpha = 0.3) +
                                                        geom_ribbon(
                                                            data = ribbon,
                                                            mapping = aes(x = rt_rt_offset, ymax = abundance, ymin = baseline, fill = peak_number_within_sample),
                                                            alpha = 0.8
                                                        ) +
                                                        geom_text(
                                                            data = filter(signal_for_this_peak, ion == "tic"),
                                                            mapping = aes(label = peak_number_within_sample, x = median(rt_rt_offset), y = max(abundance))
                                                        )
                                                }
                                            }
                                    }

                                ## Draw the plot and communicate

                                    chromatogram_plot <- chromatogram_plot +
                                        # geom_line(
                                        #     data = filter(chromatograms_updated, ion == "baseline"),
                                        #     mapping = aes(x = rt_rt_offset, y = abundance), color = "grey"
                                        # ) +
                                        geom_line(
                                            data = filter(chromatograms_updated, ion != "baseline"),
                                            mapping = aes(x = rt_rt_offset, y = abundance, color = ion),
                                            alpha = 0.8
                                        )
                                    #     scale_x_continuous(limits = c(x_axis_start, x_axis_end), name = "Retention (Scan number)") +
                                    #     scale_y_continuous(name = "Abundance (counts)") +
                                    #     facet_grid(path_to_cdf_csv~., scales = "free_y", labeller = labeller(path_to_cdf_csv = facet_labels))
                                    
                                    output$chromatograms <- renderPlot({chromatogram_plot})

                                    cat("Chromatogram updated.\n")
                            }
                        })

                    ## Transfer chromatogram_brush info to selected_peak table
                        
                        output$selected_peak <- DT::renderDataTable(DT::datatable({

                            if ( !is.null(input$chromatogram_brush )) {
                                peak_points <- brushedPoints(chromatograms_updated, input$chromatogram_brush)
                                peak_data <-  data.frame(
                                    peak_start = min(peak_points$rt),
                                    peak_end = max(peak_points$rt),
                                    peak_ID = "unknown",
                                    path_to_cdf_csv = peak_points$path_to_cdf_csv[1],
                                    area = sum(peak_points$abundance[peak_points$ion == "tic"])
                                )
                                peak_data
                            } else {
                                NULL
                            }

                        }))

                    ## Remove selected peaks with "R" (82) keystroke
                        
                        observeEvent(input$keypress, {
                            if( input$keypress == 82 ) { # Update on "R"

                                if ( !is.null(input$chromatogram_brush )) {

                                    peak_points <- brushedPoints(chromatograms_updated, input$chromatogram_brush)
                                    selection_start = min(peak_points$rt)
                                    selection_end = max(peak_points$rt)
                                    path_to_cdf_csv = peak_points$path_to_cdf_csv[1]
                                    peak_table <- read.csv("peaks_monolist.csv")

                                    peak_table <- peak_table[
                                        !apply(cbind(
                                            peak_table$peak_start > selection_start,
                                            peak_table$peak_end < selection_end,
                                            peak_table$path_to_cdf_csv == as.character(peak_points$path_to_cdf_csv[1])
                                        ), 1, all)
                                    ,]

                                    write.table(
                                        x = peak_table,
                                        file = "peaks_monolist.csv",
                                        append = FALSE,
                                        row.names = FALSE,
                                        col.names = TRUE,
                                        sep = ","
                                    )
                                    cat("Removed peaks.\n")
                                }
                            }
                        })

                    ## Append single peak with "A" (65) keystroke 
                        
                        observeEvent(input$keypress, {

                            # Do nothing if no selection
                                if(is.null(input$chromatogram_brush)) {
                                    return()
                                }

                            # If selection and "a" is pressed, add the selection to the peak table
                                if( input$keypress == 65 ) {
                                
                                    write.table(
                                        x = data.frame(
                                                peak_start = min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
                                                peak_end = max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
                                                peak_ID = "unknown",
                                                path_to_cdf_csv = brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1],
                                                area = sum(brushedPoints(chromatograms_updated, input$chromatogram_brush)$tic) - sum(brushedPoints(chromatograms_updated, input$chromatogram_brush)$baseline)
                                            ),
                                        file = "peaks_monolist.csv",
                                        append = TRUE,
                                        row.names = FALSE,
                                        col.names = FALSE,
                                        sep = ","
                                    )

                                    output$peak_table <- rhandsontable::renderRHandsontable(rhandsontable::rhandsontable({
                                        peak_table2 <- read.csv("peaks_monolist.csv")
                                        peak_table2
                                    }))
                                    cat("Added peak.\n")
                                }
                        })

                    ## Global append peak with "G" (71) keystroke

                        observeEvent(input$keypress, {

                            # Do nothing if no selection
                                if(is.null(input$chromatogram_brush)) {
                                    return()
                                }

                            # If selection and "G" is pressed, add the selection to the peak table
                                if( input$keypress == 71 ) {
                                
                                    x_peaks <-  data.frame(
                                                    peak_start = min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt_rt_offset),
                                                    peak_end = max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt_rt_offset),
                                                    peak_ID = "unknown",
                                                    path_to_cdf_csv = unique(chromatograms_updated$path_to_cdf_csv),
                                                    area = sum(brushedPoints(chromatograms_updated, input$chromatogram_brush)$tic) - sum(brushedPoints(chromatograms_updated, input$chromatogram_brush)$baseline)
                                                )

                                    x_peaks$peak_start <- x_peaks$peak_start - chromatograms_updated$rt_offset[match(x_peaks$path_to_cdf_csv, chromatograms_updated$path_to_cdf_csv)]
                                    x_peaks$peak_end <- x_peaks$peak_end - chromatograms_updated$rt_offset[match(x_peaks$path_to_cdf_csv, chromatograms_updated$path_to_cdf_csv)]

                                    write.table(
                                        x = x_peaks,
                                        file = "peaks_monolist.csv",
                                        append = TRUE,
                                        row.names = FALSE,
                                        col.names = FALSE,
                                        sep = ","
                                    )

                                    output$peak_table <- DT::renderDataTable(DT::datatable({
                                        peak_table <- read.csv("peaks_monolist.csv")
                                        peak_table
                                    }))
                                    cat("Added global peak.\n")
                                }
                        })

                    ## [MS1 extract ("shift+1"), update ("shift+2"), subtract ("shift+3"), library search ("shift+4"), save ("shift+5")]
                        
                        observeEvent(input$keypress, {

                            ## If "shift+1", MS from chromatogram brush -> MS_out_1
                                
                                if( input$keypress == 33 ) {

                                    ret_start_MS <- min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
                                    ret_end_MS <- max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
                                    sample_name_MS <- as.character(brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1])

                                    chromatogram_updated_MS <- filter(chromatograms_updated, path_to_cdf_csv == sample_name_MS)

                                    MS_ret_start_line <<- chromatogram_updated_MS$rt_first_row_in_raw[which.min(abs(
                                        chromatogram_updated_MS$rt - ret_start_MS
                                    ))] + 1

                                    MS_ret_end_line <<- chromatogram_updated_MS$rt_last_row_in_raw[which.min(abs(
                                        chromatogram_updated_MS$rt - ret_end_MS
                                    ))] + 1

                                    if (.Platform$OS.type == "unix") {
                                        
                                        system(paste0("head -1"," ",CDF_directory_path,"/",sample_name_MS," > ",CDF_directory_path,"/temp_MS.csv"))
                                        system(paste0("sed -n ",MS_ret_start_line,",",MS_ret_end_line,"p ",sample_name_MS," >> ",CDF_directory_path,"/temp_MS.csv"))    
                                        framedDataFile <- readMonolist(paste0(CDF_directory_path, "/temp_MS.csv"))
                                    
                                    }

                                    if (.Platform$OS.type == "windows") {
                                        
                                        # system(paste0("head -1"," ",CDF_directory_path,"/",sample_name_MS," > ",CDF_directory_path,"/temp_MS.csv"))
                                        # system(paste0("sed -n ",MS_ret_start_line,",",MS_ret_end_line,"p ",sample_name_MS," >> ",CDF_directory_path,"/temp_MS.csv"))    
                                        framedDataFile <-   isolate(
                                            as.data.frame(
                                                data.table::fread(sample_name_MS)
                                            )
                                        )
                                    
                                    }

                                    framedDataFile <- isolate(dplyr::filter(framedDataFile, rt > ret_start_MS, rt < ret_end_MS))
                                    framedDataFile$mz <- round(framedDataFile$mz, 1)
                                    framedDataFile <- framedDataFile %>% group_by(mz) %>% summarize(intensity = sum(intensity))
                                    MS_out_1 <<- as.data.frame(framedDataFile)
                                }

                            ## If "shift+3" subtract chromatogram brush from MS_out_1 -> MS_out_1

                                if ( input$keypress == 35 ) {

                                    if (!exists("MS_out_1")) {
                                        cat("No mass spectrum extracted yet.\n")
                                        return()
                                    } else {
                                    
                                        framedDataFile_to_subtract <- isolate(as.data.frame(
                                                            data.table::fread(as.character(
                                                                brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1]
                                                            ))
                                        ))
                                        framedDataFile_to_subtract <- isolate(dplyr::filter(
                                                                framedDataFile_to_subtract, 
                                                                rt > min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
                                                                rt < max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
                                                            ))
                                        framedDataFile_to_subtract$mz <- round(framedDataFile_to_subtract$mz, 1)
                                        framedDataFile_to_subtract <- framedDataFile_to_subtract %>% group_by(mz) %>% summarize(intensity = sum(intensity))
                                        framedDataFile_to_subtract <<- as.data.frame(framedDataFile_to_subtract)
                                        joined <- left_join(MS_out_1, framedDataFile_to_subtract, by = "mz")
                                        joined$intensity.y[is.na(joined$intensity.y)] <- 0
                                        MS_out_1$intensity <- joined$intensity.x - joined$intensity.y
                                        MS_out_1$intensity[MS_out_1$intensity < 0] <- 0
                                        MS_out_1 <<- MS_out_1
                                    }
                                }

                            ## If "shift+1", or "shift+2", or "shift+3, make plot based on brush if any

                                if ( input$keypress == 33 | input$keypress == 64 | input$keypress == 35 ) {

                                    if (!exists("MS_out_1")) {
                                        cat("No mass spectrum extracted yet.\n")
                                        return()
                                    } else {

                                        # Normalize to max abu 100
                                            MS_out_1$intensity <- MS_out_1$intensity*100/max(MS_out_1$intensity)

                                        # Set ranges on mass spec (allow zooming in with selection), Sometimes brush returns Inf or -Inf if the range gets too small, so if that happens, just go up to the main view
                                    
                                            if (isolate(is.null(input$massSpectra_1_brush))) {
                                                MS1_low_x_limit <- 0; MS1_high_x_limit <- 800; MS1_high_y_limit <- 110
                                            } else {
                                                MS1_low_x_limit <- isolate(min(brushedPoints(MS_out_1, input$massSpectra_1_brush)$mz))
                                                MS1_high_x_limit <- isolate(max(brushedPoints(MS_out_1, input$massSpectra_1_brush)$mz))
                                                MS1_high_y_limit <- max(dplyr::filter(MS_out_1, mz > MS1_low_x_limit & mz < MS1_high_x_limit)$intensity) + 8
                                            }
                                            if (MS1_low_x_limit %in% c(Inf, -Inf) | MS1_high_x_limit %in% c(Inf, -Inf) | MS1_high_y_limit %in% c(Inf, -Inf)) {
                                                MS1_low_x_limit <- 0; MS1_high_x_limit <- 800; MS1_high_y_limit <- 110
                                            }

                                        # Make plot

                                            output$massSpectra_1 <- renderPlot({

                                                ggplot() + 
                                                    geom_bar(
                                                        data = MS_out_1,
                                                        mapping = aes(x = mz, y = intensity),
                                                        stat = "identity", width = 0.1,
                                                        color = "black", fill = "grey"
                                                    ) +
                                                    theme_classic() +
                                                    scale_x_continuous(expand = c(0,0)) +
                                                    scale_y_continuous(expand = c(0,0)) +
                                                    coord_cartesian(xlim = c(MS1_low_x_limit, MS1_high_x_limit), ylim = c(0, MS1_high_y_limit)) +
                                                    geom_text(
                                                        data = # If it's less than 10 bars, label all. Otherwise, label 10 biggest ones
                                                            if (MS1_high_x_limit - MS1_low_x_limit >= 1) {
                                                                dplyr::filter(MS_out_1, mz > MS1_low_x_limit & mz < MS1_high_x_limit)[
                                                                    order(
                                                                        dplyr::filter(MS_out_1, mz > MS1_low_x_limit & mz < MS1_high_x_limit)$intensity,
                                                                        decreasing = TRUE
                                                                    )[1:10]
                                                                ,]
                                                            } else {
                                                                MS_out_1
                                                            },
                                                        mapping = aes(x = mz, y = intensity + 5, label = mz)
                                                    )
                                            })
                                    }
                                }

                            ## If "shift+4" library lookup

                                if (input$keypress == 36) {
                                        
                                    message("Searching reference library for unknown spectrum...\n")

                                        ## Round to nominal mass spectrum
                                            MS_out_1$mz <- floor(MS_out_1$mz)
                                            MS_out_1 %>%
                                                group_by(mz) %>%
                                                dplyr::summarize(intensity = sum(intensity)) -> MS_out_1
                                            MS_out_1$intensity <- MS_out_1$intensity/max(MS_out_1$intensity)*100

                                        ## Add zeros for mz values missing from unknown spectrum if necessary
                                            mz_missing <- seq(1, 800, 1)[!seq(1, 800, 1) %in% MS_out_1$mz]
                                            if (length(mz_missing) > 0) {
                                                MS_out_1 <- rbind(MS_out_1, data.frame(mz = mz_missing, intensity = 0))
                                                MS_out_1 <- MS_out_1[order(MS_out_1$mz),]
                                            }

                                        ## Bind it with metadata
                                            unknown <- cbind(
                                                data.frame(
                                                    Accession_number = "unknown",
                                                    Compound_systematic_name = "unknown",
                                                    Compound_common_name = "unknown",
                                                    SMILES = NA,
                                                    Source = "unknown"
                                                ),
                                                t(MS_out_1$intensity)
                                            )
                                            colnames(unknown)[6:805] <- paste("mz_", seq(1,800, 1), sep = "")

                                        ## Bind it to the library and run the lookup
                                            lookup_data <- rbind(busta_spectral_library, unknown)
                                            
                                            hits <- runMatrixAnalysis(
                                                data = lookup_data,
                                                analysis = c("hclust"),
                                                column_w_names_of_multiple_analytes = NULL,
                                                column_w_values_for_multiple_analytes = NULL,
                                                columns_w_values_for_single_analyte = colnames(lookup_data)[6:805],
                                                columns_w_additional_analyte_info = NULL,
                                                columns_w_sample_ID_info = c("Accession_number", "Compound_systematic_name"),
                                                transpose = FALSE,
                                                unknown_sample_ID_info = c("unknown_unknown"),
                                                scale_variance = FALSE,
                                                kmeans = "none",
                                                na_replacement = "drop",
                                                output_format = "long"
                                            )

                                        ## Find distance between tips

                                            phylo <- runMatrixAnalysis(
                                                data = hits[!is.na(hits$sample_unique_ID),],
                                                analysis = c("hclust_phylo"),
                                                column_w_names_of_multiple_analytes = "analyte_name",
                                                column_w_values_for_multiple_analytes = "value",
                                                columns_w_values_for_single_analyte = NULL,
                                                columns_w_additional_analyte_info = NULL,
                                                columns_w_sample_ID_info = c("Accession_number", "Compound_systematic_name"),
                                                transpose = FALSE,
                                                unknown_sample_ID_info = NULL,
                                                scale_variance = FALSE,
                                                kmeans = "none",
                                                na_replacement = "drop",
                                                output_format = "long"
                                            )

                                            distances_1 <- as.data.frame(cophenetic.phylo(phylo)[,colnames(cophenetic.phylo(phylo)) == "unknown_unknown"])
                                            distances_2 <- data.frame(
                                                sample_unique_ID = rownames(distances_1),
                                                distance = distances_1[,1]
                                            )
                                            distances_2$distance <- normalize(distances_2$distance, old_min = min(distances_2$distance), old_max = max(distances_2$distance), new_min = 100, new_max = 0)

                                        ## Make the bar data

                                            hits[!is.na(hits$sample_unique_ID),] %>%
                                                select(analyte_name, value, sample_unique_ID) %>%
                                                unique() -> bars

                                            bars$analyte_name <- gsub(".*_", "", bars$analyte_name)

                                            bars %>%
                                                group_by(sample_unique_ID) %>%
                                                arrange(desc(value)) -> bar_labels

                                            bar_labels <- bar_labels[1:110,]

                                        ## Order bar data

                                            bars$sample_unique_ID <- factor(
                                                bars$sample_unique_ID, levels = distances_2$sample_unique_ID[order(distances_2$distance, decreasing = TRUE)]
                                            )

                                            distances_2$sample_unique_ID <- factor(
                                                distances_2$sample_unique_ID, levels = distances_2$sample_unique_ID[order(distances_2$distance, decreasing = TRUE)]
                                            )

                                            bar_labels$sample_unique_ID <- factor(
                                                bar_labels$sample_unique_ID, levels = distances_2$sample_unique_ID[order(distances_2$distance, decreasing = TRUE)]
                                            )

                                        ## Make the bar plot
                                            
                                            plot <- ggplot() +
                                                geom_col(data = bars, aes(x = as.numeric(as.character(analyte_name)), y = value)) +
                                                geom_text(data = unique(select(bars, sample_unique_ID)), aes(label = sample_unique_ID, x = 400, y = 90), hjust = 0.5, size = 4) +
                                                facet_grid(sample_unique_ID~.) +
                                                theme_bw() +
                                                scale_x_continuous(name = "m/z") +
                                                scale_y_continuous(name = "Relative intensity (%)") +
                                                geom_text(
                                                        data = bar_labels,
                                                        mapping = aes(
                                                            x = as.numeric(as.character(analyte_name)),
                                                            y = value + 5, label = analyte_name
                                                        )
                                                    ) +
                                                geom_text(
                                                        data = distances_2,
                                                        mapping = aes(
                                                            x = 800,
                                                            y = 75,
                                                            label = paste0(
                                                                "Relative similarity to unknown: ",
                                                                round(distance, 2), "%"
                                                            ), hjust = 1
                                                        )
                                                    )

                                            output$massSpectrumLookup <- renderPlot({plot})

                                    message("Done.\n")

                                }
                        
                            ## If "shift+5" save mass spectrum

                                if( input$keypress == 37 ) {

                                    # Do nothing if no MS extracted
                                        if (!exists("MS_out_1")) {
                                            cat("No mass spectrum extracted yet.\n")
                                            return()
                                        } else {
                                    
                                        MS_out_1_to_write <- MS_out_1
                                        MS_out_1_to_write$mz <- round(MS_out_1_to_write$mz)
                                        MS_out_1_to_write %>% 
                                            group_by(mz) %>%
                                            dplyr::summarize(intensity = sum(intensity)) -> MS_out_1_to_write
                                        MS_out_1_to_write$intensity <- MS_out_1_to_write$intensity*100/max(MS_out_1_to_write$intensity)
                                        MS_out_1_to_write <- data.frame(
                                            Compound_common_name = NA,
                                            Compound_systematic_name = NA,
                                            SMILES = NA,
                                            Source = "Busta",
                                            mz = MS_out_1_to_write$mz,
                                            abu = MS_out_1_to_write$intensity
                                        )
                                        writeCSV(MS_out_1_to_write)
                                    }
                                }
                        })

                }

            ## Call the app
                
                shinyApp(ui = ui, server = server)

        }