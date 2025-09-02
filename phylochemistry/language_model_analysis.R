################################
## LANGUAGE MODEL ANALYSIS v1 ##
################################

###### Libraries

    Bioconductor_packages <- c(
        "Biostrings", "ggtree", "ips", "rhdf5"
    )
    CRAN_packages <- c(
        "BiocManager", "ggrepel", "ggplot2", "rentrez",
        "XML", "lubridate", "tibble", "httr",
        "jsonlite", "dplyr", "tidyr", "FactoMineR", "fpc",
        "Rtsne", "Rfast", "umap", "cluster", "ape",
        "bootstrap", "shipunov"
    )
    packages_needed <- c(CRAN_packages, Bioconductor_packages)[!c(CRAN_packages, Bioconductor_packages) %in% rownames(installed.packages())]

    ## Determine if anything needs to be installed
        
        if (length(packages_needed) > 0) {

            message <- paste0(
                "You need to install the following packages before proceeding: ",
                paste(packages_needed, collapse = ", "),
                " Is it okay if this script installs them for you?"
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

            } else {
                stop("Cannot load this script without the required packages. Exiting.")
            }
        }
        
        message("Loading language model packages...")

        invisible(suppressMessages(suppressWarnings(lapply(c(CRAN_packages, Bioconductor_packages), library, character.only = TRUE))))

        message("Loading language model functions...")

###### Functions

    #### searchNCBI

        searchNCBI <- function(search_term, retmax = 5) {

            # Search, message if no results
                search_results <- rentrez::entrez_search(db = "protein", term = search_term, retmax = retmax)
                if (length(search_results$ids) == 0) {
                    message("No proteins found for the search term.")
                    return(NULL)
                }

            # Write results to temp file then read in as stringset
                temp_fasta <- tempfile(fileext = ".fasta")
                write(
                    rentrez::entrez_fetch(db = "protein", id = search_results$ids, rettype = "fasta"),
                    file = temp_fasta
                )
                    return(readAAStringSet(temp_fasta))
        }

    #### searchPubMed

        searchPubMed <- function(search_terms, pubmed_api_key, sort = c("date", "relevance"), retmax_per_term = 20) {

              pm_entries <- character()
              term_vector <- character()

              for (i in 1:length(search_terms)) { 
                  search_output <- rentrez::entrez_search(
                    db = "pubmed", term = as.character(search_terms[i]), 
                    retmax = retmax_per_term, use_history = TRUE, sort = sort[1]
                  )
                  
                  # Initialize variables for retry mechanism
                  success <- FALSE
                  attempts <- 0
                  max_attempts <- 3  # Maximum number of retry attempts
                  
                  while (!success && attempts < max_attempts) {
                      attempts <- attempts + 1
                      
                      # Attempt to fetch data
                      query_output <- try(rentrez::entrez_fetch(
                        db = "pubmed", web_history = search_output$web_history, 
                        rettype = "xml", retmax = retmax_per_term, 
                        api_key = pubmed_api_key, timeout = 60), silent = TRUE)
                      
                      # Check if the attempt was successful
                      if (inherits(query_output, "try-error")) {
                          message("Error encountered. Attempt ", attempts, " of ", max_attempts, ". Retrying in 5 seconds...")
                          Sys.sleep(5)  # Wait before retrying
                      } else {
                          # Parse and store the entries if successful
                          current_pm_entries <- XML::xmlToList(XML::xmlParse(query_output))
                          pm_entries <- c(pm_entries, current_pm_entries)
                          term_vector <- c(term_vector, rep(as.character(search_terms[i]), length(current_pm_entries)))
                          success <- TRUE
                      }
                  }
                  
                  if (!success) {
                      message("Failed to fetch data for term: ", as.character(search_terms[i]), " after ", max_attempts, " attempts.")
                  }
                  
                  Sys.sleep(4)  # Wait between different search terms to respect API rate limits
              }

              unique_indices <- !duplicated(pm_entries)
              pm_entries <- pm_entries[unique_indices]
              term_vector <- term_vector[unique_indices]

              pm_results <- list()
              for (i in 1:length(pm_entries)) { # i=1

                  if (length(pm_entries[[i]]) == 1) { next }
                  if (is.null(pm_entries[[i]]$MedlineCitation$Article$ELocationID$text)) { next }

                  options <- which(names(pm_entries[[i]]$MedlineCitation$Article) == "ELocationID")
                  for (option in options) { # option = 4
                      if (grepl("10\\.", pm_entries[[i]]$MedlineCitation$Article[[option]]$text)) {
                          doi <<- pm_entries[[i]]$MedlineCitation$Article[[option]]$text
                          break
                      } else {next}
                  }

                  pm_results[[i]] <- data.frame(
                      entry_number = as.numeric(i),
                      term = term_vector[[i]],
                      date = lubridate::as_date(paste(
                          pm_entries[[i]]$MedlineCitation$DateRevised$Year,
                          pm_entries[[i]]$MedlineCitation$DateRevised$Month,
                          pm_entries[[i]]$MedlineCitation$DateRevised$Day,
                      sep = "-"
                      )),
                      journal = pm_entries[[i]]$MedlineCitation$Article$Journal$Title,
                      title = paste0(pm_entries[[i]]$MedlineCitation$Article$ArticleTitle, collapse = ""),
                      doi = doi,
                      abstract = paste0(pm_entries[[i]]$MedlineCitation$Article$Abstract$AbstractText, collapse = "")
                  )
              }
              return(as_tibble(do.call(rbind, pm_results)))
        }

    #### embedText

        embedText <- function(
          df,
          column_name,
          hf_api_key = Sys.getenv("HF_TOKEN"),
          path_to_glove_file = "glove.6B.50d.txt",
          local = FALSE,
          model_id = "BAAI/bge-small-en-v1.5",
          batch_size = 16,
          max_retries = 5,
          timeout_sec = 60
        ) {
          # deps
          suppressPackageStartupMessages({
            library(jsonlite)
            library(httr)
            library(dplyr)
            library(tibble)
            library(data.table)
            library(stringr)
          })

          # --- prep input ----
          df <- as_tibble(df)
          if (!column_name %in% colnames(df)) stop("`column_name` not found in df.")
          text_vector <- as.character(df[[column_name]])
          text_vector[is.na(text_vector)] <- ""  # avoid NAs to API

          if (local) {
            # ---------- LOCAL (GloVe) ----------
            if (!file.exists(path_to_glove_file)) {
              stop("GloVe file not found at: ", path_to_glove_file)
            }
            # load GloVe
            glove <- fread(path_to_glove_file, data.table = FALSE, quote = "")
            rownames(glove) <- glove[, 1]
            glove <- glove[, -1, drop = FALSE]

            # simple tokenizer
            word_tokenizer <- function(x) {
              # split on non-letters/numbers, lowercase
              toks <- str_split(tolower(x), "[^a-z0-9_]+", simplify = FALSE)[[1]]
              toks[nzchar(toks)]
            }

            pb <- txtProgressBar(min = 0, max = length(text_vector), style = 3)
            embeddings_list <- vector("list", length(text_vector))
            for (i in seq_along(text_vector)) {
              tokens <- word_tokenizer(text_vector[i])
              if (length(tokens) == 0) {
                # empty -> zeros
                embeddings_list[[i]] <- rep(0, ncol(glove))
              } else {
                present <- tokens[tokens %in% rownames(glove)]
                if (length(present) == 0) {
                  embeddings_list[[i]] <- rep(0, ncol(glove))
                } else {
                  embeddings_list[[i]] <- colMeans(glove[present, , drop = FALSE], na.rm = TRUE)
                }
              }
              setTxtProgressBar(pb, i)
            }
            close(pb)

            embeddings_df <- as.data.frame(do.call(rbind, embeddings_list))
            colnames(embeddings_df) <- paste0("embedding_", seq_len(ncol(embeddings_df)))
            out <- bind_cols(df, embeddings_df)
            return(out)
          }

          # ---------- REMOTE (Hugging Face Router) ----------
          if (is.null(hf_api_key) || hf_api_key == "") {
            stop("Please provide an HF API key via `hf_api_key` or set HF_TOKEN env var.")
          }

          base_url <- sprintf(
            "https://router.huggingface.co/hf-inference/models/%s/pipeline/feature-extraction",
            URLencode(model_id, reserved = TRUE)
          )

          # helper: POST with retries/backoff
          post_with_retries <- function(payload_json) {
            delay <- 1
            for (attempt in seq_len(max_retries)) {
              resp <- httr::POST(
                url = base_url,
                httr::add_headers(
                  Authorization = paste0("Bearer ", hf_api_key),
                  `Content-Type` = "application/json"
                ),
                body = payload_json,
                encode = "raw",
                timeout(timeout_sec)
              )

              code <- resp$status_code
              if (code >= 200 && code < 300) return(resp)

              # handle 429 / 5xx with backoff
              if (code == 429 || (code >= 500 && code < 600)) {
                retry_after <- as.numeric(httr::headers(resp)[["retry-after"]])
                if (is.finite(retry_after)) {
                  Sys.sleep(retry_after)
                } else {
                  Sys.sleep(delay)
                  delay <- min(delay * 2, 30)
                }
                next
              }

              # Other errors: stop with message body
              msg <- tryCatch(httr::content(resp, as = "text", encoding = "UTF-8"), error = function(e) "")
              stop(sprintf("HF request failed (HTTP %s): %s", code, msg))
            }
            stop("HF request failed after retries.")
          }

          # batching
          n <- length(text_vector)
          idx <- split(seq_len(n), ceiling(seq_len(n) / batch_size))

          pb <- txtProgressBar(min = 0, max = length(idx), style = 3)
          all_rows <- vector("list", length(idx))
          out_dim <- NULL

          for (b in seq_along(idx)) {
            batch_idx <- idx[[b]]
            batch_text <- as.list(text_vector[batch_idx])

            # Body: either single string or list of strings is accepted.
            # We always send a list to keep parsing simple.
            body <- list(inputs = batch_text)
            payload <- jsonlite::toJSON(body, auto_unbox = TRUE)

            resp <- post_with_retries(payload)
            txt <- httr::content(resp, as = "text", encoding = "UTF-8")
            parsed <- jsonlite::fromJSON(txt, simplifyVector = FALSE)

            # Expected: list(list(numeric ...), list(numeric ...), ...)
            # Some models may return matrix-like nested lists; handle both.
            # Ensure we end up with a matrix rows = items, cols = dims.
            # If single item came back as a single vector, coerce to list-of-one.
            if (!is.list(parsed[[1]]) && is.numeric(unlist(parsed))) {
              # single vector -> wrap
              parsed <- list(parsed)
            }

            # Convert each element to numeric vector
            mat <- do.call(rbind, lapply(parsed, function(v) as.numeric(unlist(v, recursive = TRUE))))
            if (is.null(out_dim)) out_dim <- ncol(mat)
            if (ncol(mat) != out_dim) {
              stop("Inconsistent embedding dimensions returned by the model.")
            }

            all_rows[[b]] <- mat
            setTxtProgressBar(pb, b)
          }
          close(pb)

          embeddings_mat <- do.call(rbind, all_rows)
          colnames(embeddings_mat) <- paste0("embedding_", seq_len(ncol(embeddings_mat)))
          embeddings_df <- as.data.frame(embeddings_mat, stringsAsFactors = FALSE)

          # preserve original row order
          out <- bind_cols(df, embeddings_df)
          return(out)
        }



        # embedText <- function(df, column_name, hf_api_key, path_to_glove_file = "glove.6B.50d.txt", local = FALSE) {

        #     ## Prep input
        #     embeddings_list <- list()
        #     df <- as_tibble(df)
        #     text_vector <- unlist(df[,which(colnames(df) == column_name)])

        #     if (local == TRUE) {

        #     }
          
        #     ## Run HF embeddings, process and return output
        #     for (i in 1:length(text_vector)) { # i=1
        #         response <- httr::POST(
        #             url = "https://api-inference.huggingface.co/models/BAAI/bge-small-en-v1.5",
        #             httr::add_headers(Authorization = paste0("Bearer ", hf_api_key)),
        #             body = toJSON(list(inputs = text_vector[i])),
        #             encode = "json"
        #         )
        #         if (response$status_code == 429) {stop("Warning: you have (probably) exceeded your HuggingFace rate limit.")}
        #         embeddings_list[[i]] <- as.numeric(as.character(jsonlite::fromJSON(httr::content(response, as = "text", encoding = "UTF-8"))))
        #     }

        #     embeddings_df <- as.data.frame(do.call(rbind, embeddings_list))
        #     colnames(embeddings_df) <- paste0("embedding_", seq_len(ncol(embeddings_df)))
        #     df <- bind_cols(df, embeddings_df)
        #     return(df)
        # }

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
                                            "pca", "pca_ord", "pca_dim",
                                            "mca", "mca_ord", "mca_dim",
                                            "mds", "mds_ord", "mds_dim",
                                            "tsne", "dbscan", "kmeans",
                                            "hclust", "hclust_phylo", "hclust_cat", "dist"
                                        ),
                                        parameters = NULL,
                                        column_w_names_of_multiple_analytes = NULL,
                                        column_w_values_for_multiple_analytes = NULL,
                                        columns_w_values_for_single_analyte = NULL,
                                        columns_w_additional_analyte_info = NULL,
                                        columns_w_sample_ID_info = NULL,
                                        transpose = FALSE,
                                        distance_method = c("euclidean", "manhattan", "gower"),
                                        agglomeration_method = c(
                                            "ward.D2", "ward.D", "single", "complete",
                                            "average", # (= UPGMA)
                                            "mcquitty", # (= WPGMA)
                                            "median", # (= WPGMC)
                                            "centroid" # (= UPGMC)
                                        ),
                                        unknown_sample_ID_info = NULL,
                                        components_to_return = 2,
                                        scale_variance = NULL, ## default = TRUE, except for hclust, then default = FALSE
                                        na_replacement = c("mean", "none", "zero", "drop"),
                                        output_format = c("wide", "long"),
                                        ...
                                    ) {

                # Check that argument names are spelled correctly

                    passed_args <- names(c(as.list(environment()), list(...)))

                    if (!all(passed_args %in%
                        c(
                            "data",
                            "analysis",
                            "parameters",
                            "column_w_names_of_multiple_analytes",
                            "column_w_values_for_multiple_analytes",
                            "columns_w_values_for_single_analyte",
                            "columns_w_additional_analyte_info",
                            "columns_w_sample_ID_info",
                            "transpose",
                            "distance_method",
                            "agglomeration_method",
                            "unknown_sample_ID_info",
                            "components_to_return",
                            "scale_variance",
                            "kmeans",
                            "na_replacement",
                            "output_format",
                            "..."
                        ))
                    ) {stop("One of your argument names is misspelled, please double check spelling.")}

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

                    if (all(c(
                        is.null(column_w_values_for_multiple_analytes),
                        is.null(columns_w_values_for_single_analyte)
                    ))) { stop("You need to specify at least one column with values for analytes.")}

                # Pre-process data

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

                    # Check to see if analyte columns are numeric and compatible with analysis

                        which_analyte_columns <- which(colnames(data_wide) %in% analyte_columns)
                        are_they_numeric <- list()
                        for( i in which_analyte_columns ) {
                          are_they_numeric <- c(are_they_numeric, is.numeric(data_wide[[i]]))
                        }

                        # Should selected analysis proceed?

                            if ( all(unlist(are_they_numeric)) ) {
                                if ( analysis %in% c("pca", "pca_dim", "pca_ord") ) {
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
                                if ( analysis %in% c("pca", "pca_dim", "pca_ord") ) {
                                    stop("Analytes are all categorical, but the analysis selected is for numeric variables. Please choose a different analysis method.\n")
                                }

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

                    # Prepare the matrix

                        matrix <- as.data.frame(data_wide[,match(analyte_columns, colnames(data_wide))])
                        rownames(matrix) <- data_wide$sample_unique_ID

                    # Handle NAs

                        if( na_replacement[1] == "none") {
                        }

                        if( na_replacement[1] == "drop" ) {
                            message("Dropping any variables in your dataset that have NA as a value.\nVariables dropped:\n")
                            if (length(names(which(apply(is.na(matrix), 2, any)))) > 0) {
                                message(names(which(apply(is.na(matrix), 2, any))))    
                            } else {
                                message("none")
                            }
                            cat("\n")
                            matrix <- matrix[,!apply(is.na(matrix), 2, any)]
                        }
                        if( na_replacement[1] %in% c("zero", "mean") ) {

                            if( any(is.na(matrix)) ) {
                                
                                message(paste0("Replacing NAs in your data with ", na_replacement[1]), "\n")

                                    for( column in 1:dim(matrix)[2]) {
                                        
                                        if( any(is.na(matrix[,column])) ) {

                                            if( na_replacement[1] == "mean" ) {
                                                replacement <- mean(matrix[,column], na.rm = TRUE)
                                            }
                                            if( na_replacement[1] == "zero" ) {
                                                replacement <- 0
                                            }
                                            if( !any(na_replacement %in% c("mean", "zero")) ) {
                                                stop("Your data contains NAs. Please specify how to deal with them using na_replacement. \n")
                                            }
                                            
                                            matrix[,column][is.na(matrix[,column])] <- as.numeric(replacement)

                                        } else {}
                                    }
                            }
                        }

                    # Transpose matrix, if requested
                        if( transpose == TRUE ) { matrix <- t(matrix) }

                    # Run unknown, if requested

                        if( length(unknown_sample_ID_info) > 0 ) {
                        
                            ## Use na_replacement != "drop"

                                if( na_replacement[1] != "drop" ) {
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

                    # Scale data, unless not requested

                        if ( is.null(scale_variance) ) {
                            if (analysis != "hclust") {scale_variance <- TRUE} else {scale_variance <- FALSE}
                        }

                        if( scale_variance == TRUE & !analysis %in% c("mca", "mca_ord", "mca_dim")) {
                            
                            scaled_matrix <- scale(matrix)

                            if( any(is.na(scaled_matrix)) ) {
                                message("Some analytes have zero variance and will be assigned a value of zero in the scaled matrix.")
                                scaled_matrix[is.na(scaled_matrix)] <- 0
                            }

                        }

                        if( scale_variance == FALSE ) {
                            scaled_matrix <- matrix
                        }

                    ## HCLUST, HCLUST_PHYLO ##

                        if( analysis == "hclust" | analysis == "hclust_phylo") {
                            
                            ## BClust approach to bootstrapped hclust

                                bclust <- Bclust(
                                    scaled_matrix, method.d = distance_method[1],
                                    method.c = agglomeration_method[1],
                                    monitor = FALSE
                                )
                                # print(bclust$value)
                                # plot(bclust)
                                phylo <- ape::as.phylo(bclust$hclust)

                            if( analysis == "hclust_phylo" ) {
                                return(phylo)
                                stop("Returning hclust_phylo.")
                            }
                            clustering <- ggtree::fortify(phylo)
                            clustering$sample_unique_ID <- clustering$label
                            clustering$bootstrap <- NA

                            ## Add bootstrap values starting from the furthest node to the highest node
                                bs_vals <- data.frame(
                                    xval = clustering$x[clustering$isTip != TRUE],
                                    bs_val = NA
                                )
                                for (i in 1:length(bclust$value)) { # i=1
                                    bs_vals$bs_val[
                                        order(bs_vals$xval, decreasing = TRUE)[i]
                                    ] <- bclust$values[i]
                                }

                            clustering$bootstrap[clustering$isTip != TRUE] <- bs_vals$bs_val
                        }

                        if (analysis == "hclust_cat") {

                            ## "bootstrap" approach
                                temp <- as.data.frame(lapply(scaled_matrix, function(x) if(is.character(x)) factor(x) else x))
                                rownames(temp) <- rownames(scaled_matrix)
                                scaled_matrix <- temp
                                createHclustObject <- function(x)hclust(cluster::daisy(x, metric = distance_method[1]), method = agglomeration_method[1])
                                b <- bootstrap(scaled_matrix, fun = createHclustObject, n = 100L)
                                phylo <- ape::as.phylo(createHclustObject(scaled_matrix))

                                clustering <- ggtree::fortify(phylo)
                                clustering$sample_unique_ID <- clustering$label
                                clustering$bootstrap <- NA

                            ## Add bootstrap values starting from the furthest node to the highest node
                                bs_vals <- data.frame( xval = clustering$x[clustering$isTip != TRUE], bs_val = NA )
                                for (i in 1:length(b)) { bs_vals$bs_val[order(bs_vals$xval, decreasing = TRUE)[i]] <- b[i] }

                            clustering$bootstrap[clustering$isTip != TRUE] <- bs_vals$bs_val

                        }

                        # if( analysis == "hclust" | analysis == "hclust_phylo" ) {
                        #     scaled_matrix <- as.data.frame(lapply(scaled_matrix, function(x) if(is.character(x)) factor(x) else x))
                        #     createHclustObject <- function(x)hclust(cluster::daisy(x, metric = distance_method[1]), method = agglomeration_method[1])
                        #     b <- bootstrap(scaled_matrix, fun = createHclustObject, n = 100L)
                        #     phylo <- ape::as.phylo(createHclustObject(scaled_matrix))
                        #     if( analysis == "hclust_phylo" ) {
                        #         return(phylo)
                        #         stop("Returning hclust_phylo.")
                        #     }
                        #     clustering <- ggtree::fortify(phylo)
                        #     clustering$sample_unique_ID <- clustering$label
                        #     clustering$bootstrap <- NA

                        #     ## Add bootstrap values starting from the furthest node to the highest node
                        #         bs_vals <- data.frame(
                        #             xval = clustering$x[clustering$isTip != TRUE],
                        #             bs_val = NA
                        #         )
                        #         for (i in 1:length(b)) { # i=1
                        #             bs_vals$bs_val[
                        #                 order(bs_vals$xval, decreasing = TRUE)[i]
                        #             ] <- b[i]
                        #         }

                        #     clustering$bootstrap[clustering$isTip != TRUE] <- bs_vals$bs_val
                        # }

                    # Generate distance matrix

                        if(  !analysis %in% c("mca", "mca_ord", "mca_dim") ) {

                            dist_matrix <- stats::dist(scaled_matrix, method = distance_method[1])
                            
                            if( analysis == "dist") {
                                return(dist_matrix)
                                stop()
                            }

                        }

                        if( analysis %in% c("mca", "mca_ord", "mca_dim") ) { scaled_matrix <- matrix }

                    ## Dimensionality reduction

                        ## MDS

                            if( analysis == "mds" ) {
                                coords <- stats::cmdscale(dist_matrix)
                                colnames(coords) <- c("Dim_1", "Dim_2")
                                clustering <- as_tibble(coords)
                                clustering$sample_unique_ID <- rownames(coords)
                            }
                    
                        ## tSNE

                            if( analysis == "tsne" ) {
                                clustering <- Rtsne(scaled_matrix, theta = 0.0, perplexity = 2)
                                clustering <- as.data.frame(clustering$Y)
                                clustering$sample_unique_ID <- rownames(scaled_matrix)
                                colnames(clustering) <- c("Dim_1", "Dim_2", "sample_unique_ID")
                                clustering <- select(clustering, sample_unique_ID, Dim_1, Dim_2)
                                rownames(clustering) <- NULL
                            }

                        ## umap

                            if( analysis == "umap" ) {
                                clustering <- as.data.frame(umap(scaled_matrix)$layout)
                                clustering$sample_unique_ID <- rownames(clustering)
                                colnames(clustering) <- c("Dim_1", "Dim_2", "sample_unique_ID")
                                clustering <- select(clustering, sample_unique_ID, Dim_1, Dim_2)
                                rownames(clustering) <- NULL
                            }

                        ## MCA, MCA_ORD, MCA_DIM ##

                            if( analysis == "mca" ) {
                                message("Running Multiple Correspondence Analysis, extracting sample coordinates...\n")
                                coords <- FactoMineR::MCA(matrix, graph = FALSE)$ind$coord[,c(1:components_to_return)]
                                clustering <- as_tibble(coords)
                                clustering$sample_unique_ID <- rownames(coords)
                                colnames(clustering) <- c("Dim_1", "Dim_2", "sample_unique_ID")
                                message("Done!\n")
                            }

                            if( analysis == "mca_ord" ) {
                                message("Running Multiple Correspondence Analysis, extracting ordination plot...\n")
                                coords <- FactoMineR::MCA(matrix, graph = FALSE)$var$eta2[,c(1,components_to_return)]
                                clustering <- as_tibble(coords)
                                clustering$analyte <- rownames(coords)
                                colnames(clustering) <- c("Dim_1", "Dim_2", "analyte")
                                clustering <- select(clustering, analyte, Dim_1, Dim_2)
                                return(clustering)
                                stop("Returning ordination plot coordinates. \nDone!")
                            }

                            if( analysis == "mca_dim" ) {
                                message("Running Multiple Correspondence Analysis, extracting dimensional contributions...\n")
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
                                coords <- FactoMineR::PCA(scaled_matrix, graph = FALSE, scale.unit = FALSE)$ind$coord[,c(1:components_to_return)]
                                clustering <- as_tibble(coords)
                                clustering$sample_unique_ID <- rownames(coords)
                                # colnames(clustering) <- c("Dim_1", "Dim_2", "sample_unique_ID")
                            }

                            if( analysis == "pca_ord" ) {
                                coords <- FactoMineR::PCA(scaled_matrix, graph = FALSE, scale.unit = FALSE)$var$coord[,c(1:components_to_return)]
                                clustering <- as_tibble(coords)
                                clustering$analyte <- rownames(coords)
                                clustering <- select(clustering, analyte, c(paste0("Dim.", seq(1,components_to_return,1))))
                                # colnames(clustering) <- c("analyte", "Dim_1", "Dim_2")
                                return(clustering)
                                stop("Returning ordination plot coordinates.")
                            }

                            if( analysis == "pca_dim" ) {
                                coords <- FactoMineR::PCA(scaled_matrix, graph = FALSE, scale.unit = FALSE)$eig[,2]
                                clustering <- tibble::enframe(coords, name = NULL)
                                clustering$principal_component <- names(coords)
                                clustering$principal_component <- as.numeric(gsub("comp ", "", clustering$principal_component))
                                colnames(clustering)[colnames(clustering) == "value"] <- "percent_variance_explained"
                                clustering <- select(clustering, principal_component, percent_variance_explained)
                                return(clustering)
                                stop("Returning eigenvalues.")
                            }
                
                    ## Clustering

                        if(  !analysis %in% c("mca", "mca_ord", "mca_dim") ) {

                            if( any(is.na(scaled_matrix)) == TRUE ) {
                                stop("clustering cannot handle NA. Please choose an option for na_replacement.")
                            }

                        }

                        ## DBSCAN

                            if( analysis == "dbscan" ) {

                                if ( length(parameters) > 0 ) {
                                    cluster_k <- parameters[1]
                                    cluster_threshold <- parameters[2]
                                }

                                if ( length(parameters) == 0 ) {
                                    findClusterParameters(dist_matrix = dist_matrix, matrix = matrix, analysis = "dbscan")
                                }

                                message("Using", cluster_k, "as a value for k.\n")
                                message("Using", cluster_threshold, "as a value for threshold.\n")
                                clustering <- as_tibble(data.frame(
                                    sample_unique_ID = colnames(as.matrix(dist_matrix)),
                                    cluster = paste0("cluster_", fpc::dbscan(dist_matrix, eps = as.numeric(cluster_threshold), MinPts = as.numeric(cluster_k), scale = FALSE, method = "dist")[[1]])
                                ))
                                clustering$cluster[clustering$cluster == "cluster_0"] <- NA

                            }

                        ## k-means

                            if (analysis == "kmeans") {

                                if ( length(parameters) > 0 ) {
                                    n_clusters <- parameters[1]
                                }

                                if ( length(parameters) == 0 ) {
                                    findClusterParameters(dist_matrix = dist_matrix, matrix = matrix, analysis = "kmeans")
                                }

                                message("Using", n_clusters, "as a value for cluster_number.\n")
                                clustering <- as_tibble(data.frame(
                                    sample_unique_ID = colnames(as.matrix(dist_matrix)),
                                    cluster = stats::kmeans(x = matrix, centers = as.numeric(n_clusters), nstart = 25, iter.max = 1000)$cluster
                                ))

                            }

                        ## OPTICS

                            # out <- dbscan::optics(scaled_matrix, minPts = 5)
                            # out <- data.frame(
                            #     order = out$order,
                            #     reach_dist = out$reachdist,
                            #     name = rownames(scaled_matrix)
                            # )
                            # out$name <- factor(out$name, levels = rev(rownames(scaled_matrix)[out$order]))
                            # ggplot(out[2:19,]) +
                            #     geom_col(aes(x = name, y = reach_dist))                                

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

                            rownames_matrix <- tibble::enframe(rownames(scaled_matrix), name = NULL)
                            colnames(rownames_matrix)[1] <- "sample_unique_ID"

                            # if (analysis != "pca") { ## Don't do this for pca for some reason?? I don't understand why...
                                clustering <- full_join(
                                    clustering,
                                    as_tibble(cbind(rownames_matrix, as_tibble(matrix))),
                                    by = "sample_unique_ID"
                                )
                            # }
                            # clustering

                            ## Order the returned matrix so that the sample_unique_ID comes first

                                clustering <- select(clustering, sample_unique_ID, everything())

                        # Annotate internal nodes in tree output if all its descendants share a property

                            if( analysis == "hclust" ) {
                                for( node in dplyr::filter(clustering, isTip == FALSE)$node ) {
                                    for (sample_property in colnames(clustering)[colnames(clustering) %in% columns_w_sample_ID_info] ) {
                                        descends <- clustering[clustering$node %in% ips::descendants(phylo, node),]
                                        if (length( unlist(unique(descends[,colnames(descends) == sample_property])) ) == 1 ) {
                                            clustering[
                                                which(clustering$node == node),
                                                which(colnames(clustering) == sample_property)
                                            ] <- unlist(descends[,colnames(descends) == sample_property])[1]
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
                                analyte_annotation_frame <- unique(select(ungroup(data), all_of(c(column_w_names_of_multiple_analytes, columns_w_additional_analyte_info))))
                                clustering <- left_join(clustering, analyte_annotation_frame, by = c("analyte_name" = column_w_names_of_multiple_analytes))
                            }

                            return( clustering )
                    }
            }

    #### embedAminoAcids

        embedAminoAcids <- function(
                amino_acid_stringset = NULL,
                dataframe = NULL,
                biolm_api_key = NULL,
                nvidia_api_key = NULL,
                platform = c("nvidia", "biolm", "local"),
                input_type = c("amino_acid_stringset", "dataframe"),
                model_name = c("esm2-650m", "esm2-8m", "esm2-35m", "esm2-150m", "esm2_t30_150M_UR50D"),
                seq_column = NULL
            ) {
              # Select first option for parameters
              platform <- platform[1]
              input_type <- input_type[1]
              model_name <- model_name[1]
              
              # Validate that the provided model_name is allowed for the chosen platform
              if (platform == "biolm") {
                allowed_models <- c("esm2-650m", "esm2-8m", "esm2-35m", "esm2-150m")
                if (!model_name %in% allowed_models) {
                  stop("For platform 'biolm', model_name must be one of: ", paste(allowed_models, collapse = ", "))
                }
              }
              if (platform == "nvidia") {
                if (model_name != "esm2-650m") {
                  stop("For platform 'nvidia', only model_name 'esm2-650m' is allowed.")
                }
              }
              if (platform == "local") {
                if (model_name != "esm2_t30_150M_UR50D") {
                  stop("For platform 'local', only model_name 'esm2_t30_150M_UR50D' is allowed.")
                }
              }
              
              # Validate input type
              if (input_type == "amino_acid_stringset") { 
                if (!inherits(amino_acid_stringset, "XStringSet")) {
                  stop("For input_type 'amino_acid_stringset', amino_acid_stringset must be of class 'XStringSet'.")
                }
                # Convert the XStringSet to a character vector and preserve names
                sequences <- as.character(amino_acid_stringset)
                seq_names <- names(amino_acid_stringset)
              }
              
              # BIOLM    
              if (platform == "biolm") {
                # Build the items string manually by iterating over the sequences.
                # If needed, escape any internal quotes in the sequence.
                items_str <- paste0(
                  sapply(sequences, function(seq) {
                    seq_escaped <- gsub("\"", "\\\\\"", seq)
                    sprintf("{\"sequence\": \"%s\"}", seq_escaped)
                  }),
                  collapse = ","
                )
                
                # Create the final JSON payload string.
                # The formatting (line breaks and indentations) here matches the manual example.
                payload <- sprintf("{\n    \"params\": {\n        \"include\": [\n            \"mean\",\n            \"contacts\",\n            \"logits\",\n            \"attentions\"\n        ]\n    },\n    \"items\": [\n        %s\n    ]\n}", items_str)
                
                headers <- c(
                  'Authorization' = paste('Token', biolm_api_key),
                  'Content-Type'  = 'application/json',
                  'Accept'        = 'application/json'
                )
                
                # Construct the API endpoint URL dynamically based on model_name
                biolm_url <- sprintf("https://biolm.ai/api/v3/%s/encode/", model_name)
                response_content <- postForm(biolm_url, .opts = list(postfields = payload, httpheader = headers, followlocation = TRUE), style = "httppost")
                
                # Parse the JSON response
                result_data <- jsonlite::fromJSON(response_content)
                
                # Check for any error messages in the API response
                if (!is.null(result_data$error)) {
                  stop("API Error: ", result_data$error)
                }
                
                embeddings <- list()
                for (i in 1:length(result_data$results$embeddings)) { # i=1
                  embeddings[[i]] <- result_data$results$embeddings[[i]]$embedding[[1]]
                }
                embeddings <- as_tibble(as.data.frame(do.call(rbind, embeddings)))
                colnames(embeddings) <- paste0("embedding_", seq(1, dim(embeddings)[2], 1))
                embeddings <- as_tibble(as.data.frame(cbind(data.frame(name = seq_names), embeddings)))
                return(embeddings)
              }
              
              if (platform == "nvidia") {
                # NVIDIA API accepts a list of sequences
                payload <- list(
                  sequences = as.list(as.data.frame(amino_acid_stringset)[[1]]),
                  format = "h5"
                )
                response <- httr::POST(
                  url = "https://health.api.nvidia.com/v1/biology/meta/esm2-650m",
                  httr::add_headers(
                    `Content-Type` = "application/json",
                    Authorization = paste("Bearer", nvidia_api_key)
                  ),
                  body = jsonlite::toJSON(payload, auto_unbox = TRUE),
                  encode = "json"
                )
                # Save response to a temporary h5 file
                file_path <- tempfile(fileext = ".h5")
                writeBin(httr::content(response, as = "raw"), file_path)
                temp_dir <- tempdir()
                unzipped_files <- unzip(file_path, exdir = temp_dir)
                embeddings <- as_tibble(t(rhdf5::h5read(unzipped_files, "embeddings")))
                colnames(embeddings) <- paste0("embedding_", seq(1, dim(embeddings)[2], 1))
                embeddings <- as_tibble(as.data.frame(cbind(data.frame(name = seq_names), embeddings)))
                return(embeddings)
              }
              # 
              # if (platform == "local") {
              #   # Create a temporary FASTA file with all sequences
              #   temp_dir <- tempdir()
              #   fasta_file <- file.path(temp_dir, "temp_sequences.fasta")
              #   output_file <- file.path(temp_dir, "temp_embeddings.csv")
              #   Biostrings::writeXStringSet(amino_acid_stringset, fasta_file)
              #   
              #   # Run the Python script to compute embeddings
              #   python_path <- "/usr/bin/python3"  # Adjust path if necessary
              #   python_script <- "/home/bustalab/Documents/protein_lm/encode_proteins.py"  # Adjust script path as needed
              #   command <- sprintf("%s %s %s %s %s", python_path, python_script, fasta_file, output_file, model_name)
              #   system(command, intern = TRUE)
              #   
              #   if (!file.exists(output_file)) {
              #     stop("Python script did not produce an output file.")
              #   }
              #   embeddings <- read.csv(output_file)
              #   unlink(c(fasta_file, output_file))
              # }
              # 
              # # Combine embeddings with original sequence names
              # if (platform %in% c("biolm", "nvidia")) {
              #   
              #   colnames(embeddings_df) <- paste0("embedding_", seq_len(ncol(embeddings_df)))
              #   embeddings_df <- cbind(name = seq_names, embeddings_df)
              #   return(as_tibble(embeddings_df))
              # } else if (platform == "local") {
              #   embeddings_df <- as.data.frame(embeddings)
              #   embeddings_df <- cbind(name = seq_names, embeddings_df)
              #   return(as_tibble(embeddings_df))
              # }
            } 

message("Done with language model loading!")