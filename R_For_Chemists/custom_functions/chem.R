#### libraries

    library(tidyverse)
    library(ape)
    library(ggtree)
    library(rstatix)
    library(multcompView)
    library(imager)

#### datasets
    
    algae_data <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/algae_data.csv")
    alaska_lake_data <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/alaska_lake_data.csv")
    solvents <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/solvents.csv")
    periodic_table <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/per_table.csv")
    periodic_table_small <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/per_table_small.csv")
    # NY_trees <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/NY_trees.csv")
    ckd_data <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/ckd_metabolomics.csv")
    wine_grape_data <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/wine_grape_data.csv")
    # data <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/housing.csv")
    hawaii_aquifers <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/hawaii_aquifer_data.csv")
    beer_components <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/beer_components.csv")
    hops_components <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/hops_components.csv")

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
                                    analysis = c("hclust", "pca", "pca-ord", "pca-dim"),
                                    column_w_names_of_multiple_analytes = NULL,
                                    column_w_values_for_multiple_analytes = NULL,
                                    columns_w_values_for_single_analyte = NULL,
                                    columns_w_additional_analyte_info = NULL,
                                    columns_w_sample_ID_info = NULL,
                                    transpose = FALSE,
                                    unknown_sample_ID_info = NULL,
                                    scale_variance = TRUE,
                                    kmeans = c("none", "auto", "elbow", "1", "2", "3", "etc."),
                                    na_replacement = c("none", "mean", "zero", "drop"),
                                    output_format = c("wide", "long"),
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

                # if( length(column_w_names_of_multiple_analytes) > 0 ) {
                #     x <- table(select(data, column_w_names_of_multiple_analytes))
                #     if( all(range(x) / mean(x) != c(1,1)) ) {
                #         stop("There are duplicate analyte names.")
                #     }
                # } 

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

                # Run hclust, if requested

                    if( analysis == "hclust" ) {                        
                        phylo <- ape::as.phylo(stats::hclust(stats::dist(matrix)))
                        clustering <- ggtree::fortify(phylo)
                        clustering$sample_unique_ID <- clustering$label
                    }

                # Run PCA, if requested

                    if( analysis == "pca" ) {
                        if( scale_variance == TRUE ) {
                            coords <- FactoMineR::PCA(matrix, graph = FALSE)$ind$coord[,c(1:2)]    
                        } else {
                            coords <- FactoMineR::PCA(matrix, graph = FALSE, scale.unit = FALSE)$ind$coord[,c(1:2)]
                        }
                        clustering <- as_tibble(coords)
                        clustering$sample_unique_ID <- rownames(coords)
                    }

                # Run PCA and return ordination plot coordinates, if requested

                    if( analysis == "pca-ord" ) {
                        if( scale_variance == TRUE ) {
                            coords <- FactoMineR::PCA(matrix, graph = FALSE)$var$coord[,c(1,2)]
                        } else {
                            coords <- FactoMineR::PCA(matrix, graph = FALSE, scale.unit = FALSE)$var$coord[,c(1,2)]
                        }
                        clustering <- as_tibble(coords)
                        clustering$analyte <- rownames(coords)
                        clustering <- select(clustering, analyte, Dim.1, Dim.2)
                        return(clustering)
                        stop("Returning ordination plot coordinates.")
                    }

                # Run PCA and return eigenvalues, if requested

                        if( analysis == "pca-dim" ) {
                            if( scale_variance == TRUE ) {
                                coords <- FactoMineR::PCA(matrix, graph = FALSE)$eig[,2]
                            } else {
                                coords <- FactoMineR::PCA(matrix, graph = FALSE, scale.unit = FALSE)$eig[,2]
                            }
                            clustering <- tibble::enframe(coords, name = NULL)
                            clustering$principal_component <- names(coords)
                            clustering$principal_component <- as.numeric(gsub("comp ", "", clustering$principal_component))
                            colnames(clustering)[colnames(clustering) == "value"] <- "percent_variance_explained"
                            clustering <- select(clustering, principal_component, percent_variance_explained)
                            return(clustering)
                            stop("Returning eigenvalues.")
                        }

                # Run kMeans, if requested

                    if( kmeans[1] %in% c("none", FALSE) ) {}

                    if( !kmeans[1] %in% c("none", FALSE) ) {

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

                    clustering <- right_join(
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
                        for( node in dplyr::filter(clustering, isTip == FALSE)$node ) {
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

#### descendants

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

#### p_groups

    p_groups <- function(data) {
        
        p <- data$p.adj
        names(p) <- paste(data$group1, data$group2, sep = "-")

        output <- data.frame(
            treatment = names(multcompLetters(p)$Letters),
            group = multcompLetters(p)$Letters,
            spaced_group = multcompLetters(p)$monospacedLetters
        )
        return(output)
    }

#### HSD.test

    HSD.test <- function (y, trt, DFerror, MSerror, alpha = 0.05, group = TRUE, 
            main = NULL, unbalanced = FALSE, console = FALSE) 
        {
            name.y <- paste(deparse(substitute(y)))
            name.t <- paste(deparse(substitute(trt)))
            if (is.null(main)) 
                main <- paste(name.y, "~", name.t)
            clase <- c("aov", "lm")
            if ("aov" %in% class(y) | "lm" %in% class(y)) {
                if (is.null(main)) 
                    main <- y$call
                A <- y$model
                DFerror <- df.residual(y)
                MSerror <- deviance(y)/DFerror
                y <- A[, 1]
                ipch <- pmatch(trt, names(A))
                nipch <- length(ipch)
                for (i in 1:nipch) {
                    if (is.na(ipch[i])) 
                        return(if (console) cat("Name: ", trt, "\n", 
                          names(A)[-1], "\n"))
                }
                name.t <- names(A)[ipch][1]
                trt <- A[, ipch]
                if (nipch > 1) {
                    trt <- A[, ipch[1]]
                    for (i in 2:nipch) {
                        name.t <- paste(name.t, names(A)[ipch][i], sep = ":")
                        trt <- paste(trt, A[, ipch[i]], sep = ":")
                    }
                }
                name.y <- names(A)[1]
            }
            junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
            Mean <- mean(junto[, 1])
            CV <- sqrt(MSerror) * 100/Mean
            medians <- tapply.stat(junto[, 1], junto[, 2], stat = "median")
            for (i in c(1, 5, 2:4)) {
                x <- tapply.stat(junto[, 1], junto[, 2], function(x) quantile(x)[i])
                medians <- cbind(medians, x[, 2])
            }
            medians <- medians[, 3:7]
            names(medians) <- c("Min", "Max", "Q25", "Q50", "Q75")
            means <- tapply.stat(junto[, 1], junto[, 2], stat = "mean")
            sds <- tapply.stat(junto[, 1], junto[, 2], stat = "sd")
            nn <- tapply.stat(junto[, 1], junto[, 2], stat = "length")
            means <- data.frame(means, std = sds[, 2], r = nn[, 2], medians)
            names(means)[1:2] <- c(name.t, name.y)
            ntr <- nrow(means)
            Tprob <- qtukey(1 - alpha, ntr, DFerror)
            nr <- unique(nn[, 2])
            nr1 <- 1/mean(1/nn[, 2])
            if (console) {
                cat("\nStudy:", main)
                cat("\n\nHSD Test for", name.y, "\n")
                cat("\nMean Square Error: ", MSerror, "\n\n")
                cat(paste(name.t, ",", sep = ""), " means\n\n")
                print(data.frame(row.names = means[, 1], means[, 2:6]))
                cat("\nAlpha:", alpha, "; DF Error:", DFerror, "\n")
                cat("Critical Value of Studentized Range:", Tprob, "\n")
            }
            HSD <- Tprob * sqrt(MSerror/nr)
            statistics <- data.frame(MSerror = MSerror, Df = DFerror, 
                Mean = Mean, CV = CV, MSD = HSD)
            if (group & length(nr) == 1 & console) 
                cat("\nMinimun Significant Difference:", HSD, "\n")
            if (group & length(nr) != 1 & console) 
                cat("\nGroups according to probability of means differences and alpha level(", 
                    alpha, ")\n")
            if (length(nr) != 1) 
                statistics <- data.frame(MSerror = MSerror, Df = DFerror, 
                    Mean = Mean, CV = CV)
            comb <- utils::combn(ntr, 2)
            nn <- ncol(comb)
            dif <- rep(0, nn)
            sig <- NULL
            LCL <- dif
            UCL <- dif
            pvalue <- rep(0, nn)
            for (k in 1:nn) {
                i <- comb[1, k]
                j <- comb[2, k]
                dif[k] <- means[i, 2] - means[j, 2]
                sdtdif <- sqrt(MSerror * 0.5 * (1/means[i, 4] + 1/means[j, 
                    4]))
                if (unbalanced) 
                    sdtdif <- sqrt(MSerror/nr1)
                pvalue[k] <- round(1 - ptukey(abs(dif[k])/sdtdif, ntr, 
                    DFerror), 4)
                LCL[k] <- dif[k] - Tprob * sdtdif
                UCL[k] <- dif[k] + Tprob * sdtdif
                sig[k] <- " "
                if (pvalue[k] <= 0.001) 
                    sig[k] <- "***"
                else if (pvalue[k] <= 0.01) 
                    sig[k] <- "**"
                else if (pvalue[k] <= 0.05) 
                    sig[k] <- "*"
                else if (pvalue[k] <= 0.1) 
                    sig[k] <- "."
            }
            if (!group) {
                tr.i <- means[comb[1, ], 1]
                tr.j <- means[comb[2, ], 1]
                comparison <- data.frame(difference = dif, pvalue = pvalue, 
                    signif. = sig, LCL, UCL)
                rownames(comparison) <- paste(tr.i, tr.j, sep = " - ")
                if (console) {
                    cat("\nComparison between treatments means\n\n")
                    print(comparison)
                }
                groups = NULL
            }
            if (group) {
                comparison = NULL
                Q <- matrix(1, ncol = ntr, nrow = ntr)
                p <- pvalue
                k <- 0
                for (i in 1:(ntr - 1)) {
                    for (j in (i + 1):ntr) {
                        k <- k + 1
                        Q[i, j] <- p[k]
                        Q[j, i] <- p[k]
                    }
        }
            groups <- orderPvalue(means[, 1], means[, 2], alpha, 
                Q, console)
            names(groups)[1] <- name.y
            if (console) {
                cat("\nTreatments with the same letter are not significantly different.\n\n")
                print(groups)
            }
        }
        parameters <- data.frame(test = "Tukey", name.t = name.t, 
            ntr = ntr, StudentizedRange = Tprob, alpha = alpha)
        rownames(parameters) <- " "
        rownames(statistics) <- " "
        rownames(means) <- means[, 1]
        means <- means[, -1]
        output <- list(statistics = statistics, parameters = parameters, 
            means = means, comparison = comparison, groups = groups)
        class(output) <- "group"
        invisible(output)
    }

#### orderPvalue

    orderPvalue <- function (treatment, means, alpha, pvalue, console) 
        {
            n <- length(means)
            z <- data.frame(treatment, means)
            letras <- c(letters[1:26], LETTERS[1:26], 1:9, c(".", "+", 
                "-", "*", "/", "#", "$", "%", "&", "^", "[", "]", ":", 
                "@", ";", "_", "?", "!", "=", "#", rep(" ", 2000)))
            w <- z[order(z[, 2], decreasing = TRUE), ]
            M <- rep("", n)
            k <- 1
            k1 <- 0
            j <- 1
            i <- 1
            cambio <- n
            cambio1 <- 0
            chequeo = 0
            M[1] <- letras[k]
            q <- as.numeric(rownames(w))
            while (j < n) {
                chequeo <- chequeo + 1
                if (chequeo > n) 
                    break
                for (i in j:n) {
                    s <- pvalue[q[i], q[j]] > alpha
                    if (s) {
                        if (lastC(M[i]) != letras[k]) 
                          M[i] <- paste(M[i], letras[k], sep = "")
                    }
                    else {
                        k <- k + 1
                        cambio <- i
                        cambio1 <- 0
                        ja <- j
                        for (jj in cambio:n) M[jj] <- paste(M[jj], "", 
                          sep = "")
                        M[cambio] <- paste(M[cambio], letras[k], sep = "")
                        for (v in ja:cambio) {
                          if (pvalue[q[v], q[cambio]] <= alpha) {
                            j <- j + 1
                            cambio1 <- 1
                          }
                          else break
                        }
                        break
                    }
                }
                if (cambio1 == 0) 
                    j <- j + 1
            }
            w <- data.frame(w, stat = M)
            trt <- as.character(w$treatment)
            means <- as.numeric(w$means)
            output <- data.frame(means, groups = M)
            rownames(output) <- trt
            if (k > 81) 
                cat("\n", k, "groups are estimated.The number of groups exceeded the maximum of 81 labels. change to group=FALSE.\n")
            invisible(output)
        }

#### tapply.stat

    tapply.stat <- function (y, x, stat = "mean") 
        {
            k <- 0
            numerico <- NULL
            if (is.null(ncol(x))) {
                if (is.numeric(x)) {
                    k <- 1
                    numerico[1] <- 1
                }
            }
            else {
                ncolx <- ncol(x)
                for (i in 1:ncolx) {
                    if (is.numeric(x[, i])) {
                        k <- k + 1
                        numerico[k] <- i
                    }
                }
            }
            cx <- deparse(substitute(x))
            cy <- deparse(substitute(y))
            x <- data.frame(c1 = 1, x)
            y <- data.frame(v1 = 1, y)
            nx <- ncol(x)
            ny <- ncol(y)
            namex <- names(x)
            namey <- names(y)
            if (nx == 2) 
                namex <- c("c1", cx)
            if (ny == 2) 
                namey <- c("v1", cy)
            namexy <- c(namex, namey)
            for (i in 1:nx) {
                x[, i] <- as.character(x[, i])
            }
            z <- NULL
            for (i in 1:nx) {
                z <- paste(z, x[, i], sep = "&")
            }
            w <- NULL
            for (i in 1:ny) {
                m <- tapply(y[, i], z, stat)
                m <- as.matrix(m)
                w <- cbind(w, m)
            }
            nw <- nrow(w)
            c <- rownames(w)
            v <- rep("", nw * nx)
            dim(v) <- c(nw, nx)
            for (i in 1:nw) {
                for (j in 1:nx) {
                    v[i, j] <- strsplit(c[i], "&")[[1]][j + 1]
                }
            }
            rownames(w) <- NULL
            junto <- data.frame(v[, -1], w)
            junto <- junto[, -nx]
            names(junto) <- namexy[c(-1, -(nx + 1))]
            if (k == 1 & nx == 2) {
                junto[, numerico[1]] <- as.character(junto[, numerico[1]])
                junto[, numerico[1]] <- as.numeric(junto[, numerico[1]])
                junto <- junto[order(junto[, 1]), ]
            }
            if (k > 0 & nx > 2) {
                for (i in 1:k) {
                    junto[, numerico[i]] <- as.character(junto[, numerico[i]])
                    junto[, numerico[i]] <- as.numeric(junto[, numerico[i]])
                }
                junto <- junto[do.call("order", c(junto[, 1:(nx - 1)])), 
                    ]
            }
            rownames(junto) <- 1:(nrow(junto))
            return(junto)
        }

#### readPolylist

    #' Reads a polylist
    #'
    #' Reads a spreadsheet in polylist format (wide format, multiple column and row headers), and converts it into tidy format
    #' @param polylist_in_path The path to the polylist (in .csv format) to be read
    #' @param centroid The characters in the cell that defines the boundaries of the multiple row and column headers
    #' @param table_value_unit The name to be given to the column that will contain the values in the polylist table
    #' @examples
    #' @export
    #' readPolylist

    readPolylist <- function(   
                        polylist_in_path,
                        centroid = "~~~",
                        table_value_unit = "abundance"
                    ) {

        # Check for centroid
            if ( length(centroid) != 1 ) {
                stop("Please provide a centroid")
            }

        # Import polylist
            polylist <- as.data.frame(data.table::fread(polylist_in_path, header = FALSE))

        # Identify location of centroid
            center_column <- unlist(apply(polylist, 1, function(x) grep(centroid, x)))
            center_row <- grep(centroid, polylist[,center_column])

        # Use centroid to extract vertical_monolist
            vertical_monolist <- polylist[(center_row+1):dim(polylist)[1], 1:(center_column-1)]
            colnames(vertical_monolist) <- as.character(unlist(polylist[center_row, 1:(center_column-1)]))
            vertical_monolist$URI_URI_URI <- apply(vertical_monolist, 1, function(x) paste(x, collapse = ""))

        # Use centroid to extract horizontal_monolist
            horizontal_monolist <- as.data.frame(t(polylist[-c(center_row), (center_column+1):dim(polylist)[2]]))
            rownames(horizontal_monolist) <- NULL
            colnames(horizontal_monolist) <-    c(
                                                as.character(unlist(polylist[c(1:(center_row-1), (center_row+1):dim(polylist)[1]), center_column]))[1:(center_row-1)],
                                                as.character(vertical_monolist$URI_URI_URI)
                                            )
            horizontal_monolist <- tidyr::gather(horizontal_monolist, URI_URI_URI, table_value_unit, (center_row):dim(horizontal_monolist)[2])
            colnames(horizontal_monolist)[colnames(horizontal_monolist) == "table_value_unit"] <- table_value_unit

        # Bind the two monolists, drop the URI column
            polylist <- cbind(horizontal_monolist,vertical_monolist[match(horizontal_monolist$URI_URI_URI, vertical_monolist[,colnames(vertical_monolist) == "URI_URI_URI"]),])
            polylist <- polylist[,colnames(polylist) != "URI_URI_URI"]

        # Make the table_value_unit column numeric
            polylist[,colnames(polylist) == table_value_unit] <- as.numeric(as.character(polylist[,colnames(polylist) == table_value_unit]))

        # Add Genus_species column if not present
            # if ( any(colnames(polylist) == "Genus_species") == FALSE) {
            #     print("Adding Genus_species column")
            #     polylist$Genus_species <- paste(polylist$Genus, polylist$species, sep="_")
            # }

        # Reset row numbers, return the polylist
            rownames(polylist) <- NULL
            return(as_tibble(polylist))
    }

#### drawMolecules

    drawMolecules <- function(path_to_csv) {

        data <- read_csv(path_to_csv)

        data$bond_start_x <- data$x[match(data$bond_start_atom, data$atom_number)]
        data$bond_start_y <- data$y[match(data$bond_start_atom, data$atom_number)]
        data$bond_end_x <- data$x[match(data$bond_end_atom, data$atom_number)]
        data$bond_end_y <- data$y[match(data$bond_end_atom, data$atom_number)]
        data$bond_start_element_or_group <- data$element_or_group[match(data$bond_start_atom, data$atom_number)]
        data$bond_end_element_or_group <- data$element_or_group[match(data$bond_end_atom, data$atom_number)]

        # keep <- !apply(cbind(
        #   data$element_or_group == "H",
        #   data$bond_start_element_or_group == "H",
        #   data$bond_end_element_or_group == "H"
        # ), 1, any)
        # keep[is.na(keep)] <- TRUE
        # plot_data <- data[keep,]
        plot_data <- data

        atom_colors_pre <- as.data.frame(rbind(
            c("C", "black"),
            c("CH2", "grey"),
            c("CH2OH", "#4daf4a"), # green
            c("CH3", "#ff7f00"), # orange
            c("COOH", "#ffff33"), # yellow
            c("H", "white"),
            c("O", "#e41a1c"), # red
            c("OH", "#377eb8") # blue
        ))
        atom_colors <- atom_colors_pre[,2]
        names(atom_colors) <- atom_colors_pre[,1]

        #e41a1c red
        #377eb8 blue
        #4daf4a green
        #984ea3 purple
        #ff7f00 orange
        #ffff33 yellow
        #a65628 brown

        ## Plotting

          plot <- ggplot() +

            ## Add achiral single bonds
              geom_link(
                data = filter(
                  plot_data, 
                  molecule_component == "bond" &  
                  bond_type == "single" & 
                  bond_direction == "flat"
                ),
                aes(
                  x = bond_start_x, y = bond_start_y,
                  xend = bond_end_x, yend = bond_end_y
                ), size = 2, color = "grey30"
              ) +

            ## Add down chiral single bonds
              geom_link(
                data = filter(
                  plot_data, 
                  molecule_component == "bond" &
                  bond_type == "single" &
                  bond_direction == "down"
                ),
                aes(
                  x = bond_start_x, y = bond_start_y,
                  xend = bond_end_x, yend = bond_end_y,
                  size = stat(index)
                ), color = "grey60"
              ) +
            
            ## Add up chiral single bonds
              geom_link(
                data = filter(
                  plot_data, 
                  molecule_component == "bond" &
                  bond_type == "single" &
                  bond_direction == "up"
                ),
                aes(
                  x = bond_start_x, y = bond_start_y,
                  xend = bond_end_x, yend = bond_end_y,
                  size = stat(index)
                ), color = "grey0"
              ) +

            ## Add nearly horizontal double bonds
              geom_segment(
                data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x != bond_end_x),
                aes(x = bond_start_x, y = bond_start_y-0.2, xend = bond_end_x, yend = bond_end_y-0.2),
                size = 1.5, alpha = 0.7
              ) +
              geom_segment(
                data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x != bond_end_x),
                aes(x = bond_start_x, y = bond_start_y+0.2, xend = bond_end_x, yend = bond_end_y+0.2),
                size = 1.5, alpha = 0.7
              ) +

            ## Add vertical double bonds
              geom_segment(
                data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x == bond_end_x),
                aes(x = bond_start_x-0.18, y = bond_start_y, xend = bond_end_x-0.18, yend = bond_end_y),
                size = 1.5, alpha = 0.7
              ) +
              geom_segment(
                data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x == bond_end_x),
                aes(x = bond_start_x+0.18, y = bond_start_y, xend = bond_end_x+0.18, yend = bond_end_y),
                size = 1.5, alpha = 0.7
              ) +

            ## Add atom number labels
              geom_point(
                data = filter(plot_data, molecule_component == "atom"),
                aes(x = x, y = y, fill = element_or_group),
                shape = 21, size = 4
              ) +
              geom_text(
                data = filter(plot_data, molecule_component == "atom", element_or_group %in% c("H", "COOH")),
                aes(x = x, y = y, label = atom_number),
                size = 2, color = "black"
              ) +
              geom_text(
                data = filter(plot_data, molecule_component == "atom", element_or_group != "H" & element_or_group != "COOH"),
                aes(x = x, y = y, label = atom_number),
                size = 2, color = "white"
              ) +

            ## Add molecule name as label
              geom_text(
                data = drop_na(unique(select(plot_data, molecule_name, skeleton))),
                aes(x = 1, y = 13, label = molecule_name),
                size = 4, color = "black", hjust = 0
              ) +

            ## Scales and theme
              scale_color_gradient(low = "black", high = "black") +
              scale_size_continuous(range = c(0,4)) +
              scale_x_continuous(breaks = seq(0,20,1)) +
              # scale_linetype_manual(values = bond_line_types) +
              scale_y_continuous(breaks = seq(0,20,1)) +
              scale_fill_manual(values = atom_colors, name = "") +
              facet_wrap(.~molecule_name, ncol = 2) +
              theme_void() +
              guides(size = "none", alpha = "none") +
              coord_fixed() +
              theme(
                strip.text = element_blank()
              )

            # print(plot)
              return(plot)
    }

#### SSexp

    SSexp <- structure(function (input, A, rc) {
      .expr2 <- exp(rc * input)
      .value <- A * .expr2
      .actualArgs <- as.list(match.call()[c("A", "rc")])
      if (all(unlist(lapply(.actualArgs, is.name)))) {
        .grad <- array(0, c(length(.value), 2), list(NULL, c("A", "rc")))
        .grad[, "A"] <- .expr2
        .grad[, "rc"] <- A * (.expr2 * input)
        dimnames(.grad) <- list(NULL, .actualArgs)
        attr(.value, "gradient") <- .grad
      }
      .value
    }
    , initial = function (mCall, data, LHS) {
      xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
      if (nrow(xy) < 3)
        stop("Too few distinct input values to fit an exponential")
      xy$logy <- log(xy$y)
      ## Keep only finite cases (if there are y <= 0)
      xyfinite <- xy[is.finite(xy$logy), ]
      if (nrow(xyfinite) < 2)
        stop("Too few distinct LHS values > 0 to fit an exponential")
      res <- lsfit(xyfinite$x, xyfinite$logy)$coef
      value <- c(exp(res[1]), res[2])
      setNames(value, mCall[c("A", "rc")])
    }
    , pnames = c("A", "rc"), class = "selfStart")

#### analyzeMassSpectralImages

    # library(imager)
    # Read a png with left-most pixel low mass, right-most pixel high mass
    # Filename: reference~compoundName_lowmz-rightmz

    analyzeMassSpectralImages <- function(image_directory_path_in) {

        images_to_analyze <- dir(image_directory_path_in)[grep(".png", dir(image_directory_path_in))]
        
        ## Skip comparison images, if necessary
            if (length(grep("comparison.png", images_to_analyze)) > 0) {
                images_to_analyze <- images_to_analyze[-c(grep("comparison.png", images_to_analyze))]
            }

        merged_spectral_output <- list()
        for (image_to_analyze in 1:length(images_to_analyze)) {

            ## Load and rotate image, filter for just one color channel

                cat(paste0("\n", "Analyzing ", images_to_analyze[image_to_analyze], "\n", "\n"))

                left_mz <- as.numeric(gsub("-.*$", "", gsub(".*_", "", images_to_analyze[image_to_analyze])))
                right_mz <- as.numeric(gsub("\\..*$", "", gsub(".*-", "", images_to_analyze[image_to_analyze])))
                reference <- gsub("~.*$", "", images_to_analyze[image_to_analyze])

                image <- imager::load.image(paste0(image_directory_path_in, "/", images_to_analyze[image_to_analyze]))
                image <- as_tibble(as.data.frame(as.cimg(image)))
                image <- image[image$cc == 1,]
                image$y <- -as.numeric(image$y)
                image$y <- image$y + -min(image$y)
                
                # ggplot(image, aes(x = x, y = y, fill = value)) + geom_tile()

            ## Set to monochrome, normalize coordinates, sort output
                
                filtered_image <- filter(image, value < 0.8)
                filtered_image$value <- 1
                filtered_image$x <- normalize(
                    filtered_image$x, 
                    old_min = min(filtered_image$x),
                    old_max = max(filtered_image$x),
                    new_min = left_mz,
                    new_max = right_mz
                )
                filtered_image <- filtered_image[order(filtered_image$x),]

                # ggplot(filtered_image, aes(x = x, y = y)) + geom_point(size = 0.1) + coord_cartesian(xlim = c(330,350))

            ## Evaluate continuity and create processed image, do some clean up with median and negative numbers

                processed_image_y <- list()
                processed_image_x <- list()
                pb <- progress::progress_bar$new(total = length(unique(filtered_image$x)))
                for( i in 1:length(unique(filtered_image$x)) ) {
                    
                    data <- filtered_image[filtered_image$x == unique(filtered_image$x)[i],]
                    
                    # # One way to find gaps
                    # if (any(!diff(data$y) == -1)) {
                    #   start <- data$y[sum(head(rle(diff(data$y))$lengths, -1))+1]
                    # } else {
                    #   start <- max(data$y)
                    # }

                    # Another way.. both are slow
                    start <- min(data$y)
                    while((start + 1) %in% data$y) {
                        start <- start + 1
                    }

                    processed_image_y[[i]] <- start
                    processed_image_x[[i]] <- unique(filtered_image$x)[i]
                    pb$tick()
                }

                processed_image <- data.frame(
                    x = do.call(rbind, processed_image_x),
                    y = do.call(rbind, processed_image_y)
                )

                processed_image$y <- processed_image$y - median(processed_image$y)
                processed_image$y[processed_image$y < 0] <- 0

                # ggplot(processed_image, aes(x = x, y = y)) + geom_col()

            ## Bin x coordinates to nominal resolution, set negatives to zero

                processed_image$x_round <- round(processed_image$x, 0)
                processed_image <- group_by(processed_image, x_round)
                processed_image <- summarize(processed_image, y = sum(y))
                processed_image$y <- processed_image$y / max(processed_image$y) * 100
                colnames(processed_image)[colnames(processed_image) == "x_round"] <- "x"
                

                # ggplot(processed_image, aes(x = x, y = y)) + geom_col()

            ## Deal with duplicate peaks by resolving them in the way that creates the highest number of odd peaks

                badness <- 0
                for( i in 1:(length(processed_image$x)-1) ) {
                    if (processed_image$y[i] > 0 ) {
                        if( processed_image$y[i]/processed_image$y[i+1] > 0.9 & processed_image$y[i]/processed_image$y[i+1] < 1.1 ) {
                            if (processed_image$y[i] > 25) {
                                badness <- badness + 1
                            }
                        }
                    }
                }
                # warning(paste0("Badness: ", badness, "\n"))

            ## Return

                colnames(processed_image) <- c("mz", "relative_abundance")
                write_csv(processed_image, paste0(image_directory_path_in, "/", gsub(".png", "", images_to_analyze[image_to_analyze]), ".csv"))

            ## Make plot of image versus processed data

                plot_output <- ggplot() +
                    geom_col(data = processed_image, aes(x = mz, y = -relative_abundance), width = 0.8, color = "black") +
                    geom_text(data = 
                        processed_image[apply(cbind(
                            processed_image$relative_abundance > 2,
                            processed_image$mz > max(processed_image$mz)*0.999,
                            processed_image$mz < min(processed_image$mz)*1.001
                        ), 1, any),],
                        aes(x = mz, y = -relative_abundance - 5, label = mz)
                    ) +
                    geom_point(data = filtered_image, aes(x = x, y = y/max(y)*100), size = 0.1) +
                    theme_bw()
                print(plot_output)
                ggsave(filename = paste0(image_directory_path_in, "/", gsub(".png", "", images_to_analyze[image_to_analyze]), "_comparison.png"), plot = plot_output)

            ## Add to merged_spectral_output

                merged_spectral_output[[image_to_analyze]] <- data.frame(
                    reference = reference,
                    compound_name = gsub("_.*$", "", gsub(".*~", "", gsub(".png", "", images_to_analyze[image_to_analyze]))),
                    mz = processed_image$mz,
                    relative_abundance = processed_image$relative_abundance
                )
        }

            ## Write out merged_spectral_output

                merged_spectral_output <- do.call(rbind, merged_spectral_output)
                write_csv(merged_spectral_output, paste0(image_directory_path_in, "/_merged_spectral_output.csv"))
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