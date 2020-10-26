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
    hawaii_aquifer_data <- read_csv("https://thebustalab.github.io/R_For_Chemists/sample_data/hawaii_aquifer_data.csv")

#### readCSV

    #'asdf

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

##### tukey_groups

    tukey_groups <- function(data, formula) {
        
        # unique_groups <- unique(c(unique(data$group1), unique(data$group2)))

        # output <- data.frame( group = unique_groups, tukey_group = "-" )
        # output$tukey_group <- as.character(output$tukey_group)
        
        # for( i in 1:length(unique_groups) ) {
            
        #     groups_in_this_tukey_group <- c(
        #         unique_groups[i],
        #         c(
        #             dplyr::filter(data, group1 == unique_groups[i] & p.adj.signif == "ns")$group2,
        #             dplyr::filter(data, group2 == unique_groups[i] & p.adj.signif == "ns")$group1
        #         )
        #     )
            
        #     output$tukey_group[output$group %in% groups_in_this_tukey_group] <- paste0(
        #         output$tukey_group[output$group %in% groups_in_this_tukey_group],
        #         LETTERS[i]
        #     )

        # }

        # output$tukey_group <- gsub("-", "", output$tukey_group)
        # # output$tukey_group <- LETTERS[as.numeric(as.factor(output$tukey_group))]

        # return(output)

        grouping_output <- eval(
            parse(
                text = paste0(
                    "agricolae::HSD.test(aov(",
                    # gsub(" ", "", gsub("~.*$", "", formula)),
                    as.character(formula)[2],
                    "~",
                    # gsub(" ", "", gsub(".*~", "", formula)),
                    as.character(formula)[3],
                    ", data = data), trt = '",
                    # gsub(" ", "", gsub(".*~", "", formula)),
                    as.character(formula)[3],
                    "')$groups"
                )
            )
        )

        output <- data.frame(
            treatment = rownames(grouping_output),
            group = grouping_output$groups
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

##### orderPvalue

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