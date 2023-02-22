########################
## PHYLOCHEMISTRY 1.0 ##
########################

###### Libraries

    ## Load bustalab-specific libraries and functions

        if (exists("bustalab")) {
            
            if (bustalab) {

                message("bustalab is true")
                
                CRAN_packages <- c(
                    # "imager",
                    "minpack.lm"
                )

                Bioconductor_packages <- c(
                    "xcms",
                    "DESeq2",
                    "msa",
                    "rtracklayer",
                    "Biostrings",
                    "GenomicRanges",
                    "GenomicFeatures",
                    "Rsamtools"
                )

                source("https://thebustalab.github.io/phylochemistry/genomescope.R")
            }
        }

    ## Load packages in general

        if ( !exists("packages") ) {

            ## Define necessary libraries
                
                if (!exists("CRAN_packages")) {CRAN_packages <- vector()}
                CRAN_packages <- c(CRAN_packages, 
                    "gridExtra",
                    "ape",
                    "multcompView",
                    # "imager",
                    "shiny",
                    "DT",
                    "RColorBrewer",
                    "data.table",
                    "rhandsontable",
                    "ips",
                    "eulerr",
                    "phangorn",
                    "seqinr",
                    "Rfast",
                    "picante",
                    "BiocManager",
                    "googlesheets4",
                    "Hmisc",
                    "ggforce",
                    "network",
                    "pracma",
                    "ggnetwork",
                    "patchwork",
                    "FactoMineR",
                    "dplyr",
                    "stringr", 
                    "progress",
                    "tidyverse",
                    "ggrepel",
                    "cowplot",
                    "rstatix",
                    "agricolae",
                    "ggpmisc",
                    "exifr",
                    "lubridate",
                    "bio3d",
                    "remotes",
                    "gdata",
                    "treemapify",
                    "viridis",
                    "umap",
                    "ggside",
                    "fpc",
                    "dbscan",
                    "Rtsne",
                    "readxl",
                    "httpuv",
                    "ggdist"
                )

                if (!exists("Bioconductor_packages")) {Bioconductor_packages <- vector()}
                Bioconductor_packages <- c(
                    Bioconductor_packages,
                    "ggtree",
                    "ggtreeExtra"
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
        slice <- dplyr::slice

    ## Set up options

        options(readr.show_progress = FALSE)
        options(dplyr.summarise.inform = FALSE)

    ## Set up better warning messages for some functions

        shapiroTest <- function( data, ... ) {

            if ( any(summarize(data, size = n())$size < 3) ) {
                stop("One of the groups defined has fewer than 3 members. A Shapiro test cannot be run on such a group. Please filter your data or choose new groups.")
            } else {
                rstatix::shapiro_test(data = data, ...)
            }

        }

    ## Rename some stats functions so they use camelCase

        leveneTest <- function( data, ... ) {rstatix::levene_test( data = data, ... )}
        tTest <- function( data, ... ) {rstatix::t_test( data = data, ... )}
        wilcoxTest <- function( data, ... ) {rstatix::wilcox_test( data = data, ... )}
        anovaTest <- function( data, ... ) {rstatix::anova_test( data = data, ... )}
        tukeyTest <- function( data, ... ) {rstatix::tukey_hsd( data = data, ... )}
        kruskalTest <- function( data, ... ) {rstatix::kruskal_test( data = data, ... )}
        dunnTest <- function( data, ... ) {rstatix::dunn_test( data = data, ... )}
        pairwiseWilcoxTest <- function( data, ... ) {rstatix::pairwise_wilcox_test( data = data, ... )}
        pairwiseTTest <- function( data, ... ) {rstatix::pairwise_t_test( data = data, ... )}
        # tukeyHSD <- function( data, ... ) {rstatix::tukey_hsd( data = data, ... )}

    ## Set up lookups

        ## Phred_ascii_33 lookup
            phred33_lookup <- data.frame(rbind(
                c("!","0"),c("\"","1"),c("#","2"),c("$","3"),c("%","4"),c("&","5"),c("'","6"),c("(","7"),c(")","8"),c("*","9"),c("+","10"),c(",","11"),c("-","12"),c(".","13"),c("/","14"),c("0","15"),c("1","16"),c("2","17"),c("3","18"),c("4","19"),c("5","20"),c("6","21"),c("7","22"),c("8","23"),c("9","24"),c(":","25"),c(";","26"),c("<","27"),c("=","28"),c(">","29"),c("?","30"),c("@","31"),c("A","32"),c("B","33"),c("C","34"),c("D","35"),c("E","36"),c("F","37"),c("G","38"),c("H","39"),c("I","40"),c("J","41"),c("K","42"),c("L","43"),c("M","44"),c("N","45"),c("O","46"),c("P","47"),c("Q","48"),c("R","49"),c("S","50"),c("T","51"),c("U","52"),c("V","53"),c("W","54"),c("X","55"),c("Y","56"),c("Z","57"),c("[","58"),c("\\","59"),c("]","60"),c("^","61"),c("_","62"),c("`","63"),c("a","64"),c("b","65"),c("c","66"),c("d","67"),c("e","68"),c("f","69"),c("g","70"),c("h","71"),c("i","72"),c("j","73"),c("k","74"),c("l","75"),c("m","76"),c("n","77"),c("o","78"),c("p","79"),c("q","80"),c("r","81"),c("s","82"),c("t","83"),c("u","84"),c("v","85"),c("w","86"),c("x","87"),c("y","88"),c("z","89"),c("{","90"),c("|","91"),c("}","92"),c("~","93")
            ))

###### Functions

    message("Loading functions...")

    ##### Read and write functions

        #### readCSV

            #' Interactive selection of a CSV file to read
            #'
            #' @param sep The delimiter to use when writing the file out. Default is a comma, i.e. CSV.
            #' @examples
            #' @export
            #' readCSV

            readCSV <- function(sep = c(",", " ")) {
                if (sep[1] == ",") {
                    return(readr::read_csv(file.choose()))
                }  
                if (sep[1] == " ") {
                    return(utils::read.delim(file.choose(), delim = " "))
                }
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

        #### readMonolist

            #' Reads a monolist
            #'
            #' @param monolist_in_path The path to the monolist (in .csv format) to be read. URLs are also accepted.
            #' @examples
            #' @export
            #' readMonolist

            readMonolist <- function( monolist_in_path ) {
                if (length(grep("http", monolist_in_path)) > 0) {
                    monolist <- read_csv(monolist_in_path, show_col_types = FALSE)        
                } else {
                    monolist <- as.data.frame(data.table::fread(file = monolist_in_path))
                }
                return( monolist )
            }

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
                    center_row <- as.numeric(grep(centroid, polylist[,center_column]))

                # Use centroid to extract vertical_monolist
                    vertical_monolist <- as.data.frame(polylist[(center_row+1):dim(polylist)[1], 1:(center_column-1)])
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

                # Add Genus_species column if not present
                    # if ( any(colnames(polylist) == "Genus_species") == FALSE) {
                    #     print("Adding Genus_species column")
                    #     polylist$Genus_species <- paste(polylist$Genus, polylist$species, sep="_")
                    # }

                # Convert to numeric
                    for ( i in 1:dim(polylist)[2] ) {
                        if ( numbersOnly(polylist[,i]) ) {
                            polylist[,i] <- as.numeric(polylist[,i])
                        }
                    }

                # Reset row numbers, return the polylist
                    rownames(polylist) <- NULL
                    return(polylist)
            }

        #### writePolylist

            #' Writes a polylist
            #'
            #' @param polylist_out_path The path to where the polylist should be written
            #' @examples
            #' @export
            #' writePolylist

            writePolylist <- function( polylist, polylist_out_path ) {
                write.table( polylist, file = polylist_out_path, sep = ",", row.names = FALSE, col.names = FALSE)
            }

        #### readAlignment

            #' Reads an alignment into a tidy dataframe
            #'
            #' @param alignment_in_path The path to the alignment
            #' @examples
            #' @export
            #' readAlignment

            readAlignment <- function( alignment_in_path, type = c("DNA", "AA")) {

                ## DNA v AA

                    if (length(type) > 1) {
                        stop("Please specify a single type of alignment: DNA or AA")
                    }

                ## Read and return the alignment
                    
                    if ( type == "DNA" ) {
                        alignment <- phangorn::read.phyDat(
                            file = alignment_in_path,
                            format = "fasta", type = "DNA"
                        )
                    }
                    if ( type == "AA" ) {
                        alignment <- phangorn::read.phyDat(
                            file = alignment_in_path,
                            format = "fasta", type = "AA"
                        )   
                    }
                    alignment <- t(as.matrix(as.data.frame(alignment)))
                    alignment <- as_tibble(data.frame(cbind(rownames(alignment), alignment)))
                    alignment <- pivot_longer(alignment, cols = c(2:dim(alignment)[2]), names_to = "position", values_to = "state")
                    alignment$position <- as.numeric(as.character(gsub("X", "", alignment$position)))
                    colnames(alignment)[1] <- "name"
                    
                    return(alignment)
            }

        #### readTree

            #' Reads a phylogenetic tree
            #'
            #' @param tree_in_path The path to the phylogenetic tree
            #' @importFrom ape read.tree
            #' @examples
            #' @export
            #' readTree

            readTree <- function( tree_in_path ) {

                ## Read and return the tree
                    return( ape::read.tree(file = tree_in_path) )
            }

        #### writeTree

            #' Writes a phylogenetic tree
            #'
            #' @param tree_out_path The path to where the tree should be written
            #' @importFrom ape write.tree
            #' @examples
            #' @export
            #' writeTree

            writeTree <- function( tree, tree_out_path ) {
                ape::write.tree( phy = tree, file = tree_out_path )
            }

        #### downloadFasta

            #' Returns contents of one fasta file as a StringSet object.
            #'
            #' @param share_link The path to the fasta to read. Can be local or a google drive share link.
            #' @importFrom Biostrings readBStringSet readDNAStringSet readRNAStringSet readAAStringSet
            #' @examples
            #' @export
            #' downloadFasta

            downloadFasta <- function(share_link) {
                cat("Downloading... this may take a while...")
                suppressMessages(googledrive::drive_download(
                    file = googledrive::as_id(share_link),
                    path = "TEMP___fasta",
                    overwrite = TRUE
                ))
            }

        #### readFasta

        	#' Returns contents of one fasta file as a StringSet object.
            #'
            #' @param fasta_in_path The path to the fasta to read. Can be local or a google drive share link.
            #' @param fasta_type The type of sequence data contained in the fasta file. Options: "nonspecific", "DNA", "RNA", or "AA"
            #' @importFrom Biostrings readBStringSet readDNAStringSet readRNAStringSet readAAStringSet
            #' @examples
            #' @export
            #' readFasta

            readFasta <- function( fasta_in_path, fasta_type = "nonspecific" ) {

                ## If its on google drive, download it and make a temp file
                    if (length(grep("http", fasta_in_path)) > 0) {
                        
                        downloadFasta(share_link = fasta_in_path)

                        temp_fasta <- Biostrings::readBStringSet("TEMP___fasta")
                        
                        file.remove("TEMP___fasta")

                        return(temp_fasta)
                    }

                ## Read the file in
                	if ( fasta_type == "nonspecific" ) {
                		return(Biostrings::readBStringSet( filepath = fasta_in_path ))
                	}

                	if ( fasta_type == "DNA" ) {
                		return(Biostrings::readDNAStringSet( filepath = fasta_in_path ))
                	}

                	if ( fasta_type == "RNA" ) {
                		return(Biostrings::readRNAStringSet( filepath = fasta_in_path ))
                	}

                	if ( fasta_type == "AA" ) {
                		return(Biostrings::readAAStringSet( filepath = fasta_in_path ))
                	}

            }

        #### readManyFastas

        	#' Returns the contents of many fasta files as a StringSet object.
            #'
            #' @param fasta_in_paths The paths to the fastas to read
            #' @param fasta_type The type of sequence data contained in the fasta file. Options: "nonspecific", "DNA", "RNA", or "AA"
            #' @importFrom Biostrings readBStringSet readDNAStringSet readRNAStringSet readAAStringSet
            #' @examples
            #' @export
            #' readManyFastas

            readManyFastas <- function( fasta_in_paths, fasta_type = "nonspecific" ) {

        		if ( fasta_type == "nonspecific" ) {
            		temp <- BStringSet()
            	}

            	if ( fasta_type == "DNA" ) {
            		temp <- DNAStringSet()
            	}

            	if ( fasta_type == "RNA" ) {
            		temp <- RNAStringSet()
            	}

            	if ( fasta_type == "AA" ) {
            		temp <- AAStringSet()
            	}

            	for ( fasta in 1:length(fasta_in_paths) ) {

    				if ( fasta_type == "nonspecific" ) {
    	        		temp <- c(temp, Biostrings::readBStringSet( filepath = fasta_in_paths[fasta] ))
    	        	}

    	        	if ( fasta_type == "DNA" ) {
    	        		temp <- c(temp, Biostrings::readDNAStringSet( filepath = fasta_in_paths[fasta] ))
    	        	}

    	        	if ( fasta_type == "RNA" ) {
    	        		temp <- c(temp, Biostrings::readRNAStringSet( filepath = fasta_in_paths[fasta] ))
    	        	}

    	        	if ( fasta_type == "AA" ) {
    	        		temp <- c(temp, Biostrings::readAAStringSet( filepath = fasta_in_paths[fasta] ))
    	        	}

            	}

            	return( temp )

            }

        #### writeManyFastas

        	#' Writes one or more sequences as individual fasta file(s), each with one sequence
            #'
            #' @param XStringSet The sequence(s) to write out
            #' @param fasta_out_directory_path The path to the directory where the file(s) should be written
            #' @importFrom Biostrings writeXStringSet
            #' @examples
            #' @export
            #' writeManyFastas

            writeManyFastas <- function( XStringSet, fasta_out_directory_path ) {
            	for (sequence in 1:length(XStringSet)) {
            		Biostrings::writeXStringSet( x = XStringSet[sequence], filepath = paste0(fasta_out_directory_path, XStringSet@ranges@NAMES[sequence], ".fa") )
            	}
            }

        #### writeFasta

        	#' Writes one or more sequences as a single fasta file
            #'
            #' @param XStringSet The sequences to write out
            #' @param fasta_out_path The path to where the file should be written
            #' @importFrom Biostrings writeXStringSet
            #' @examples
            #' @export
            #' writeFasta

            writeFasta <- function( XStringSet, fasta_out_path, append = FALSE, type = c("DNA", "AA")) {

                if (type[1] == "DNA") {
            	   Biostrings::writeXStringSet( x = DNAStringSet(XStringSet), filepath = fasta_out_path, append = append )
                }

                if (type[1] == "AA") {
                    Biostrings::writeXStringSet( x = XStringSet, filepath = fasta_out_path, append = append )
                }

            }

        #### writeSupplementalTable

            #' Writes a human-readable supplemental table of quantitative data
            #'
            #' @param supplementalTable The data to write out
            #' @param round_level The number of decimal places summary statistics should have
            #' @param file_out_path The path to the location where the file should be written
            #' @examples
            #' @export
            #' writeSupplementalTable

            writeSupplementalTable <- function(supplementalTable, round_level, file_out_path, replicates) {
        
                # Make new object from input
                    supp_table <- supplementalTable
                    supp_table$abundance <- round(supp_table$abundance, round_level)
                
                # Get number of compound levels
                    number_of_compound_levels <- length(grep("compound", colnames(supp_table)))

                # Assign sample_unique_id
                    supp_table$sample_unique_id <- paste(supp_table$sample, supp_table$replicate, sep = "..")

                # Check that sample levels and compound levels make unique IDs
                    if (number_of_compound_levels == 2) {
                        
                        supp_table$compound_unique_id <- paste(supp_table$compound_level_1, supp_table$compound_level_2, sep = "..")
                    }

                    if ( sd(table(supp_table$compound_unique_id)) != 0 ) {
                        stop("Compound levels do not produce unique sample IDs")
                    }
                
                # Spread a subset of the table "supp_table_min"
                    supp_table_min <- supp_table[colnames(supp_table) %in% c("sample_unique_id", "abundance", "compound_unique_id")]
                    supp_table_min <- tidyr::spread(supp_table_min, sample_unique_id, abundance)

                # Assign supp_table_min the compound levels, then remove compound_unique_id
                    if (number_of_compound_levels == 2) {
                        supp_table_min$compound_level_1 <- gsub("\\.\\..*", "", supp_table_min$compound_unique_id)
                        supp_table_min$compound_level_2 <- gsub(".*\\.\\.", "", supp_table_min$compound_unique_id)
                    }
                    supp_table_min <- supp_table_min[,colnames(supp_table_min) != "compound_unique_id"]

                # Order the table according to the compound levels and then the sample names
                    supp_table_min <- supp_table_min[order(supp_table_min$compound_level_1, supp_table_min$compound_level_2),]
                    supp_table_min <- supp_table_min[,match(c("compound_level_1", "compound_level_2", unique(supp_table$sample_unique_id)), colnames(supp_table_min))]

                # Total of each compound class
                    compound_level_1_totals <- list()
                    for (i in 1:length(unique(supp_table_min$compound_level_1))) {
                        data <- supp_table_min[supp_table_min$compound_level_1 == unique(supp_table_min$compound_level_1)[i],]
                        data <- data[,colnames(data) %in% unique(supp_table$sample_unique_id)]
                        compound_level_1_totals[[i]] <- data.frame(
                            compound_level_1 = paste0("Total ", as.character(unique(supp_table_min$compound_level_1)[i])),
                            compound_level_2 = " ",
                            t(colSums(data))
                        )
                    }
                    compound_level_1_totals <- do.call(rbind, compound_level_1_totals)
                
                # Total abundnace for each sample
                    total_abundance <- cbind(
                        data.frame(
                            compound_level_1 = "Total abundance",
                            compound_level_2 = " "
                        ),
                        t(colSums(supp_table_min[,colnames(supp_table_min) %in% unique(supp_table$sample_unique_id)]))
                    )
                    
                # Space between values and summary stats (total unk, total wax)
                    space <- t(data.frame(rep(" ", dim(supp_table_min)[2])))
                    colnames(space) <- colnames(supp_table_min)

                # Bind it all together, Make everything numeric, fix rownames
                    supp_table_min <- rbind(
                        supp_table_min[supp_table_min$compound_level_1 != "UNK",],
                        space,
                        compound_level_1_totals,
                        space,
                        total_abundance
                    )
                    supp_table_min[,3:dim(supp_table_min)[2]] <- t(apply(supp_table_min[,3:dim(supp_table_min)[2]], 1, function(x) as.numeric(x)))
                    rownames(supp_table_min) <- NULL

                # Percents

                    for (i in 1:length(unique(supp_table$sample_unique_id))) {

                        data <- supp_table_min[,colnames(supp_table_min) == unique(supp_table$sample_unique_id)[i]]
                        supp_table_min$percent <- round(100*data/data[length(data)], round_level)

                        colnames(supp_table_min)[colnames(supp_table_min) == "percent"] <- paste0("percent_", unique(supp_table$sample_unique_id)[i])
                    
                    }

                # Make averages and stdevs

                    colnames(supp_table_min)[colnames(supp_table_min) %in% unique(supp_table$sample_unique_id)] <- paste0(
                        "abundance_",
                        colnames(supp_table_min)[colnames(supp_table_min) %in% unique(supp_table$sample_unique_id)]
                    )
                    
                    for (i in 1:length(unique(supp_table$sample)) ) {

                        supp_table_min$average_1 <- round(
                            apply(supp_table_min[,grep(as.character(paste0("abundance_", unique(supp_table$sample)[i])), colnames(supp_table_min))], 1, mean), round_level
                        )
                        supp_table_min$error_1 <- round(
                            apply(supp_table_min[,grep(as.character(paste0("abundance_", unique(supp_table$sample)[i])), colnames(supp_table_min))], 1, sd), round_level
                        )
                        supp_table_min$error_1 <- round(supp_table_min$error_1/sqrt(replicates), round_level)

                        colnames(supp_table_min)[colnames(supp_table_min) == "average_1"] <- paste0("abundance_", unique(supp_table$sample)[i], "_average")
                        colnames(supp_table_min)[colnames(supp_table_min) == "error_1"] <- paste0("abundance_", unique(supp_table$sample)[i], "_stderror")

                        supp_table_min$average_1 <- round(
                            apply(supp_table_min[,grep(as.character(paste0("percent_", unique(supp_table$sample)[i])), colnames(supp_table_min))], 1, mean), round_level
                        )
                        supp_table_min$error_1 <- round(
                            apply(supp_table_min[,grep(as.character(paste0("percent_", unique(supp_table$sample)[i])), colnames(supp_table_min))], 1, sd), round_level
                        )
                        supp_table_min$error_1 <- round(supp_table_min$error_1/sqrt(replicates), round_level)

                        colnames(supp_table_min)[colnames(supp_table_min) == "average_1"] <- paste0("percent_", unique(supp_table$sample)[i], "_average")
                        colnames(supp_table_min)[colnames(supp_table_min) == "error_1"] <- paste0("percent_", unique(supp_table$sample)[i], "_stderror")

                    }

                # Order the columns
                    column_order <- vector()
                    for ( i in 1:length(unique(supp_table$sample)) ) {
                        column_order <- c(column_order, colnames(supp_table_min)[grep(paste0("abundance_", unique(supp_table$sample)[i]), colnames(supp_table_min))])
                        column_order <- c(column_order, colnames(supp_table_min)[grep(paste0("percent_", unique(supp_table$sample)[i]), colnames(supp_table_min))])
                    }

                    column_order <- match(c("compound_level_1", "compound_level_2", column_order), colnames(supp_table_min))

                    supp_table_min <- supp_table_min[,column_order]

                    # tibble::add_column(supp_table_min, .after = c(2,6,9))
                    # grep("stdev", colnames)

                # Finalization of table
                    # supp_table <- cbind(supp_table, t(space), supp_table_percent)
                    supp_table_min[supp_table_min == 0] <- "n.d."
                    # gsub(" ", "_", colnames(supp_table))
                    write.csv(supp_table_min, file = file_out_path, row.names = FALSE)

            }

        #### renameFile

            renameFile <- function(from, to) {
                todir <- dirname(to)
                if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
                file.rename(from = from,  to = to)
            }

        #### openRGD

                #' Open an RScript stored on Google Drive
                #'
                #' @param drive_share_link
                #' @examples
                #' @export
                #' openRGD

                openRGD <- function(drive_share_link) {

                    file_name <- googledrive::drive_get(
                        id = googledrive::as_id(drive_share_link)
                    )$name

                    googledrive::drive_rename(
                        file = googledrive::as_id(drive_share_link),
                        name = paste0("IN_USE___", file_name)
                    )

                    googledrive::drive_download(
                        file = googledrive::as_id(drive_share_link),
                        path = paste0("IN_USE___", file_name),
                        overwrite = TRUE
                    )

                    file.edit(paste0("IN_USE___", file_name))

                }

        #### closeRGD

            #' Save and close the current RScript to Google Drive
            #'
            #' @param drive_share_link
            #' @examples
            #' @export
            #' closeRGD

            closeRGD <- function(drive_share_link){

                rstudioapi::documentSave()

                file_name <- googledrive::drive_get(
                    id = googledrive::as_id(drive_share_link)
                )$name

                googledrive::drive_update(
                    file = googledrive::as_id(drive_share_link),
                    media = paste0(file_name)
                )

                googledrive::drive_rename(
                    file = googledrive::as_id(drive_share_link),
                    name = paste0(gsub("IN_USE___", "", file_name))
                )

                rstudioapi::documentClose()

                file.remove(paste0(file_name))

            }

    ##### Polylist construction and manipulation

        #### buildPolylist

            #' Build a human-readable matrix of samples, analytes, and measurements
            #'
            #' @param samples_monolist_in_path The monolist of samples to be used in building the polylist. Should contain a column "sample_unique_ID".
            #' @param analytes_monolist_in_path The monolist of analytes to be used in building the polylist. Should contain a column "analyte_unique_ID".
            #' @param measurements_monolist_in_path The monolist of measurements to be used in building the polylist. Should contain columns "sample_unique_ID", "analyte_unique_ID", and a column of values corresponding to measurements
            #' @param centroid The characters to use in the centroid
            #' @param polylist_out_path The path to which the polylist should be written
            #' @keywords lists
            #' @examples
            #' @export
            #' buildPolylist

            buildPolylist <-    function(  
                                    samples_monolist_in_path,
                                    analytes_monolist_in_path,
                                    measurements_monolist_in_path = NULL,
                                    centroid = "~~~",
                                    polylist_out_path
                                ) {

                # Read in the monolists
                    samples_monolist <- readMonolist(samples_monolist_in_path)
                    analytes_monolist <- readMonolist(analytes_monolist_in_path)
                    centroid_location <- c((dim(analytes_monolist)[2]+1), (dim(samples_monolist)[2]+1))

                # Build polylist
                    samples_monolist <- as.data.frame(t(samples_monolist))
                    samples_monolist <- cbind(rownames(samples_monolist), samples_monolist)

                    for (i in 1:(dim(analytes_monolist)[2])) {
                        samples_monolist <- cbind("NA", samples_monolist)
                    }
                    colnames(samples_monolist) <- as.character(seq(1,dim(samples_monolist)[2]))
                    samples_monolist <- apply(samples_monolist, 2, FUN = as.character)

                    for (i in 1:(dim(samples_monolist)[2]-dim(analytes_monolist)[2])) {
                        analytes_monolist <- cbind(analytes_monolist, "NA")
                    }
                    analytes_monolist <- apply(analytes_monolist, 2, FUN = as.character)
                    analytes_monolist <- rbind(colnames(analytes_monolist), analytes_monolist)
                    colnames(analytes_monolist) <- as.character(seq(1,dim(analytes_monolist)[2]))
                    
                    polylist <- rbind(samples_monolist,analytes_monolist)

                # Fill in with measurements if specified
                    if ( !is.null(measurements_monolist_in_path) ) {
                        measurements_monolist <- readMonolist(monolist = measurements_monolist_in_path)
                        colnames(measurements_monolist)[!colnames(measurements_monolist) %in% c("sample_unique_ID", "analyte_unique_ID")] <- "value"

                        analyte_unique_ID_col <- grep("analyte_unique_ID", polylist[centroid_location[2],])
                        sample_unique_ID_row <- grep("sample_unique_ID", polylist[,centroid_location[1]])

                        pb <- progress::progress_bar$new(total = dim(measurements_monolist)[1])
                        for (i in 1:dim(measurements_monolist)[1]) {
                            data_point <- measurements_monolist[i,]
                            polylist[
                                match(data_point$analyte_unique_ID, polylist[,analyte_unique_ID_col]),
                                match(data_point$sample_unique_ID, polylist[sample_unique_ID_row,])
                            ] <- as.character(data_point$value)
                            pb$tick()
                        }
                    }
                    
                # Clean up and write out polylist
                    polylist[polylist %in% c("NA", "\"NA\"")] <- ""
                    polylist[centroid_location[2], centroid_location[1]] <- centroid
                    writePolylist(polylist = polylist, polylist_out_path = polylist_out_path)
            }

    ##### Tree and taxa manipulation

        #### gblocks

            #' Subset an alignment
            #'
            #' @param alignment_in_path The path to the information that should be used to build the tree.
            #' @param members The tips of the tree that should be included. Default is to include everything.
            #' @examples
            #' @export
            #' gblocks

            gblocks <- function(
                alignment_in_path = NULL,
                max_gap_percent = 0.5,
                min_conservation_percent = 0.3,
                sequence_type = c("DNA", "AA")
            ) {

                ## Read in the alignment and make it a matrix
                    if (sequence_type == "DNA") {
                        nucl_seqs_aligned <- phangorn::read.phyDat(file = paste(alignment_in_path), format = "fasta", type = "DNA")
                        nucl_seqs_aligned_matrix <- t(as.matrix(as.data.frame(nucl_seqs_aligned)))
                        cat(paste0("Alignment is ", dim(nucl_seqs_aligned_matrix)[2], " positions long.\n"))
                    }
                    if (sequence_type == "AA") {
                        nucl_seqs_aligned <- phangorn::read.phyDat(file = paste(alignment_in_path), format = "fasta", type = "AA")
                        nucl_seqs_aligned_matrix <- t(as.matrix(as.data.frame(nucl_seqs_aligned)))
                        cat(paste0("Alignment is ", dim(nucl_seqs_aligned_matrix)[2], " positions long.\n"))
                    }

                ## Analyze composition of each position
                    position_scores <- list()
                    for (i in 1:dim(nucl_seqs_aligned_matrix)[2]) {
                        position_scores[[i]] <- data.frame(
                            position = i,
                            gap_percent = sum(nucl_seqs_aligned_matrix[,i] == "-", na.rm = TRUE) / dim(nucl_seqs_aligned_matrix)[1],
                            conservation_percent = suppressWarnings(sum(nucl_seqs_aligned_matrix[,i] == mode(nucl_seqs_aligned_matrix[,i]), na.rm = TRUE) / dim(nucl_seqs_aligned_matrix)[1])
                        )
                    }
                    position_scores <- do.call(rbind, position_scores)

                    print(
                        ggplot() +
                            geom_line(data = position_scores, aes(x = position, y = -gap_percent+1.5-max_gap_percent), color = "gold", alpha = 0.8) +
                            geom_line(data = position_scores, aes(x = position, y = conservation_percent+0.5-min_conservation_percent), color = "maroon", alpha = 0.4) +
                            geom_hline(yintercept = 0.5) +
                            coord_cartesian(xlim = c(
                                which(movingAverage(position_scores$conservation_percent) < 0.9)[1],
                                (dim(position_scores)[1] - which(rev(movingAverage(position_scores$conservation_percent) < 0.9))[1])
                            )) +
                            geom_text(aes(x = dim(position_scores)[1]*0.45, y = 1, label = "Positions with gold and red above the line are kept.")) +
                            theme_bw()
                    )

                ## Ask user if these settings are okay

                    if ( !if (interactive()) askYesNo("Do these gblocks parameters look okay?") ) {
                        stop()
                    }

                ## Find all positions with less than threshold conservation and remove them
                    conservation_reject <- position_scores$position[position_scores$conservation_percent < min_conservation_percent]
                    cat(paste0("Removing ", length(conservation_reject), " positions that are too conserved.\n"))
                    gap_reject <- position_scores$position[position_scores$gap_percent > max_gap_percent]
                    cat(paste0("Removing ", length(gap_reject), " positions that contain too many gaps.\n"))
                    
                    reject_positions <- unique(c(conservation_reject, gap_reject))
                    reject_positions[order(reject_positions)]

                    nucl_seqs_aligned_matrix <- nucl_seqs_aligned_matrix[,-c(reject_positions)]
                    cat(paste0("Remaining positions: ", dim(nucl_seqs_aligned_matrix)[2]), "\n")

                ## Write out as fasta
                    if (sequence_type == "DNA") {
                        blocked_alignment <- DNAStringSet()
                        for (i in 1:dim(nucl_seqs_aligned_matrix)[1]) {
                            blocked_alignment <- c(blocked_alignment, DNAStringSet(paste0(nucl_seqs_aligned_matrix[i,], collapse = "")))
                        }
                        blocked_alignment@ranges@NAMES <- rownames(nucl_seqs_aligned_matrix)
                        writeFasta(blocked_alignment, paste0(alignment_in_path, "_blocked"))
                    }
                    if (sequence_type == "AA") {
                        blocked_alignment <- AAStringSet()
                        for (i in 1:dim(nucl_seqs_aligned_matrix)[1]) {
                            blocked_alignment <- c(blocked_alignment, AAStringSet(paste0(nucl_seqs_aligned_matrix[i,], collapse = "")))
                        }
                        blocked_alignment@ranges@NAMES <- rownames(nucl_seqs_aligned_matrix)
                        writeFasta(blocked_alignment, paste0(alignment_in_path, "_blocked"))
                    }

                    blocked_alignment <- t(as.matrix(as.data.frame(blocked_alignment)))
                    position_scores <- list()
                    for (i in 1:dim(blocked_alignment)[2]) {
                        position_scores[[i]] <- data.frame(
                            position = i,
                            gap_percent = sum(blocked_alignment[,i] == "-", na.rm = TRUE) / dim(blocked_alignment)[1],
                            conservation_percent = suppressWarnings(sum(blocked_alignment[,i] == mode(blocked_alignment[,i]), na.rm = TRUE) / dim(nucl_seqs_aligned_matrix)[1])
                        )
                    }
                    position_scores <- do.call(rbind, position_scores)

            }

        #### buildTree

            #' Construct various types of phylogenetic trees from alignments or other trees
            #'
            #' @param scaffold_type The type of information that should be used to build the tree. One of "amin_alignment", "nucl_alignment", or "newick"
            #' @param scaffold_in_path The path to the information that should be used to build the tree.
            #' @param members The tips of the tree that should be included. Default is to include everything.
            #' @param gblocks TRUE/FALSE whether to use gblocks on the alignment
            #' @param ml TRUE/FALSE whether to use maximum likelihood when constructing the tree.
            #' @param model_test TRUE/FALSE whether to test various maximum likelihood models while constructing the tree
            #' @param bootstrap TRUE/FALSE whether to calculate bootstrap values for tree nodes.
            #' @param ancestral_states TRUE/FALSE whether to calculate ancestral states at nodes in the tree. Requires specifying a root via the 'root' parameter
            #' @param root The tree tip to use as the root of the tree
            #' @examples
            #' @export
            #' buildTree

            buildTree <-    function(
                                scaffold_type = c("amin_alignment", "nucl_alignment", "newick"),
                                scaffold_in_path,
                                members = NULL,
                                ml = FALSE, 
                                model_test = FALSE,
                                bootstrap = FALSE,
                                ancestral_states = FALSE,
                                root = NULL
                            ) {

                if ( scaffold_type == "nucl_alignment" ) {

                    ## Create distrance tree
                        cat(paste("Making tree with ", scaffold_in_path," ...\n", sep = ""))
                            nucl_seqs_aligned <- phangorn::read.phyDat(file = paste(scaffold_in_path), format = "fasta", type = "DNA")
                            
                        cat("Creating neighbor-joining tree...\n")
                            dm <- phangorn::dist.ml(nucl_seqs_aligned, "F81")
                            NJ_tree <- phangorn::NJ(dm)
                            output <- NJ_tree

                    ## Make ml tree
                        if ( ml == TRUE ) {
                            ## Test all available nucl models, use the best one to optimize for ml
                                if ( model_test == TRUE ) {
                                    cat(paste("Testing 24 maximum liklihood models... \n"))
                                    mt <- phangorn::modelTest(nucl_seqs_aligned, tree = NJ_tree, multicore = TRUE)
                                    best_nucl_model <- gsub("\\+.*$", "", mt$Model[which.max(mt$logLik)])
                                    cat(paste("Tested 24 models, using best model:", as.character(gsub("\\+.*$","",best_nucl_model)), "\n", sep = " "))
                                } else {
                                    best_nucl_model <- "GTR"
                                }
                                cat("Creating maximum liklihood tree...\n")
                                ML_tree_start <- phangorn::pml(NJ_tree, nucl_seqs_aligned, k = 4)
                                ML_tree_optimized <- phangorn::optim.pml(ML_tree_start, rearrangement = "stochastic", optInv = TRUE, optGamma = TRUE, model = as.character(best_nucl_model))
                                output <- ML_tree_optimized$tree
                        }

                    ## Run bootstrap analysis
                        if ( bootstrap == TRUE ) {
                            if ( ml == FALSE ) {
                                stop("To calculate bootstrap values, please also run maximum likelihood estimation (ml = TRUE).\n")
                            }
                            bootstraps <- phangorn::bootstrap.pml(ML_tree_optimized, bs = 100, optNni = TRUE, multicore = FALSE)
                            ML_tree_optimized$tree$node.label <- phangorn::plotBS(ML_tree_optimized$tree, bootstraps)$node.label
                            output <- ML_tree_optimized$tree
                        }

                    ## Root the tree and run ancestral states
                        if ( ancestral_states == TRUE ) {
                            if ( ml == FALSE ) {
                                cat("To enable ancestral state reconstruction, please also run maximum likelihood estimation (ml = TRUE).\n")
                            }
                            if ( ml == TRUE ) {
                                # Root the tree, then calculate ancestral_states
                                    ML_tree_optimized$tree <- ape::root(ML_tree_optimized$tree, as.character(root))
                                    output <- list()
                                    output$tree <- ML_tree_optimized$tree
                                    output$ancestral_states <- phangorn::ancestral.pml(ML_tree_optimized)
                            }
                        }

                    ## Root the tree if no ancestral_states were requested
                        if ( ancestral_states == FALSE ) {
                            if ( length(root) > 0 ) {
                                output <- ape::root(output, as.character(root))
                            }
                        }

                    ## Return the tree
                        cat("Pro tip: most tree read/write functions reset node numbers.\nFortify your tree and save it as a csv file to preserve node numbering.\nDo not save your tree as a newick or nexus file.\n")
                        return( output )
                }

                if ( scaffold_type == "amin_alignment" ) {

                    ## Make distance tree
                        cat(paste("Making tree with ", scaffold_in_path," ...\n", sep = ""))
                            amin_seqs_aligned <- phangorn::read.phyDat(file = paste(scaffold_in_path), format = "fasta", type = "AA")
                                
                        # cat("Creating neighbor-joining tree...\n")
                        dm = phangorn::dist.ml(amin_seqs_aligned, model = "JTT")
                        NJ_tree = phangorn::NJ(dm)
                        output <- NJ_tree

                    ## Make ml tree
                        if ( ml == TRUE ) {
                            if ( model_test == TRUE ) {
                                ## Test all available amino acid models and extract the best one
                                    cat(paste("Testing 24 maximum liklihood models... \n"))
                                    mt <- phangorn::modelTest(amin_seqs_aligned, tree = NJ_tree, model = "all", multicore = TRUE)
                                    best_amin_model <- gsub("\\+.*$", "", mt$Model[which.max(mt$logLik)])
                                    cat(paste("Tested 24 models, using best model:", as.character(gsub("\\+.*$", "", best_amin_model)), "\n", sep = " "))
                            } else {
                                best_amin_model <- "JTT"
                            }
                            cat("Creating maximum liklihood tree...\n")
                            ML_tree_start <- phangorn::pml(NJ_tree, amin_seqs_aligned, model = as.character(best_amin_model), k = 4, inv = .2)
                            ML_tree_optimized <- phangorn::optim.pml(ML_tree_start, rearrangement = "stochastic", optInv = TRUE, optGamma = TRUE)
                            output <- ML_tree_optimized$tree
                        }

                    ## Run bootstrap analysis
                        if ( bootstrap == TRUE ) {
                            if ( ml == FALSE ) {
                                stop("To calculate bootstrap values, please also run maximum likelihood estimation (ml = TRUE)")
                            }
                            bootstraps <- phangorn::bootstrap.pml(ML_tree_optimized, bs = 100, optNni = TRUE, multicore = FALSE)
                            ML_tree_optimized$tree$node.label <- phangorn::plotBS(ML_tree_optimized$tree, bootstraps)$node.label
                        }

                    ## Root the tree and run ancestral states
                        if ( ancestral_states == TRUE ) {
                            if ( ml == FALSE ) {
                                cat("To enable ancestral state reconstruction, please also run maximum likelihood estimation (ml = TRUE).\n")
                            }
                            if ( ml == TRUE ) {
                                if ( length(root) > 0 ) {
                                    ML_tree_optimized$tree <- ape::root(ML_tree_optimized$tree, as.character(root))
                                }
                                output <- list()
                                output$tree <- ML_tree_optimized$tree
                                ancestral_states <- phangorn::ancestral.pml(ML_tree_optimized, return = "prob")

                                ancestral_states_output <- data.frame()
                                for( i in 1:length(ancestral_states) ) {
                                    temp <- ancestral_states[[i]]
                                    colnames(temp) <- attr(ancestral_states, "levels")
                                    rownames(temp) <- NULL
                                    ancestral_states_output <- rbind(
                                        ancestral_states_output,
                                        cbind(
                                            data.frame(
                                                name = names(ancestral_states)[i],
                                                position = seq(1,dim(temp)[1],1)
                                            ),
                                            temp
                                        )
                                    )
                                }
                                
                                output$ancestral_states <- ancestral_states_output

                            }
                        }

                    ## Return tree
                        cat("Pro tip: most tree read/write functions reset node numbers.\nFortify your tree and save it as a csv file to preserve node numbering.\nDo not save your tree as a newick or nexus file.\n")
                        return ( output )
                }

                if ( scaffold_type == "newick" ) {

                    ## Read in the newick scaffold
                        if (length(grep("google", scaffold_in_path)) > 0) {
                            temp_tree <- tempfile(fileext = ".newick")
                            suppressMessages(googledrive::drive_download(
                                file = googledrive::as_id(scaffold_in_path),
                                path = temp_tree,
                                overwrite = TRUE
                            ))
                            newick <- readTree(temp_tree)
                        } else {
                            newick <- readTree( tree_in_path = scaffold_in_path )    
                        }

                    ## Are the Genus_species in your members in the newick? Are the genera in your members in the newick?
                        compatibility <- data.frame( Genus_species = unique(members), is_species_in_tree = NA, is_genus_in_tree = NA )
                        compatibility$is_species_in_tree <- compatibility$Genus_species %in% newick$tip.label
                        compatibility$is_genus_in_tree <- gsub("_.*$", "", compatibility$Genus_species) %in% gsub("_.*$", "", as.character(newick$tip.label))

                    if ( all(compatibility$is_species_in_tree) == FALSE ) {
                        ## For Genus_species in members whose genus is missing from the tree (orphans), remove them
                            orphans <- as.character(dplyr::filter(compatibility, is_species_in_tree == FALSE & is_genus_in_tree == FALSE)$Genus_species)
                            members <- members[!(members %in% orphans)]
                            if ( length(orphans) > 0 ) {
                                cat("The following species belong to a genus not found in the newick scaffold and were removed: ")
                                for ( orphan in 1:length(orphans) ) {
                                    cat("\n")
                                    cat(orphans[orphan])
                                }
                                cat("\n")
                                cat("\n")
                            }

                        ## Check compatibility again
                            compatibility <- data.frame(Genus_species = unique(members), is_species_in_tree = NA, is_genus_in_tree = NA)
                            compatibility$is_species_in_tree <- compatibility$Genus_species %in% newick$tip.label
                            compatibility$is_genus_in_tree <- gsub("_.*$", "", compatibility$Genus_species) %in% gsub("_.*$", "", as.character(newick$tip.label))

                        ## For unique(members$Genus_species) in members not in the tree but whose genus in the tree (adoptees), substitute
                            adoptees <- as.character(dplyr::filter(compatibility, is_species_in_tree == FALSE & is_genus_in_tree == TRUE)$Genus_species)

                            for ( i in 1:length(adoptees) ) {
                                potential_foster_species <- newick$tip.label[gsub("_.*$", "", newick$tip.label) %in% gsub("_.*$", "", adoptees[i])] # all species in tree of the adoptees genus
                                available_foster_species <- potential_foster_species[!potential_foster_species %in% unique(members)] # potential_foster_species not already in the quantities
                                if ( length(available_foster_species) == 0) {
                                    members <- members[!members %in% adoptees[i]]
                                    cat(paste("There aren't enough fosters to include the following species in the tree so it was removed:", adoptees[i], "\n", sep = " "))
                                } else {
                                    cat(paste("Scaffold newick tip", available_foster_species[1], "substituted with", adoptees[i], "\n", sep = " "))
                                    newick$tip.label[newick$tip.label == as.character(available_foster_species[1])] <- as.character(adoptees[i])
                                }
                            }
                    }

                    ## Drop tree tips not in desired members
                        newick <- ape::drop.tip(newick, newick$tip.label[!newick$tip.label %in% unique(members)])

                    ## Sort members according to the tree
                        ordered_tip_labels <- subset(ggtree::fortify(newick), isTip)$label[order(subset(ggtree::fortify(newick), isTip)$y, decreasing = TRUE)]
                        members <- factor(members, levels = rev(ordered_tip_labels))

                    ## Return tree
                        cat("Pro tip: most tree read/write functions reset node numbers.\nFortify your tree and save it as a csv file to preserve node numbering.\nDo not save your tree as a newick or nexus file.\n")
                        return ( newick )
                }

            }

        #### pruneTree

            #' Prune a tree so it only contains user-specified tips
            #'
            #' @param tree The tree to manipulate
            #' @param tips_to_keep The tips to keep
            #' @importFrom ape drop.tip
            #' @examples
            #' @export
            #' pruneTree

            pruneTree <- function( tree_in_path, tips_to_keep ) {

                ## Read tree, drop all tips but those specified by user and return
                    tree <- readTree(tree_in_path)
                    tree <- ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% tips_to_keep])
                    return( tree )
            }

        #### collapseTree

            #' Create a high-level cladogram from a low-level tree
            #'
            #' @param tree The tree to manipulate, must be a phylo object
            #' @param associations A two-column set of relationships between tree tip names (col 1) and the level on which to summarize (col 2)
            #' @export
            #' @examples
            #' updateCountData()

            collapseTree <- function(tree, associations) {

                ## Convert assocaitions into data.frame
                    associations <- as.data.frame(associations)

                ## Drop all species names not in association list
                    tree <- ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% associations[,1]])

                ## Subset the association list to one memeber per family
                    tip_labels <- associations[associations[,1] %in% tree$tip.label,]

                ## Collapse tree, rename tips, return
                    if ( length(duplicated(tip_labels[,2])) > 0 ) {
                        collapsed_tree <- ape::drop.tip(tree, tip_labels[,1][duplicated(tip_labels[,2])])
                    } else {
                        collapsed_tree <- tree
                    }
                    
                    collapsed_tree$tip.label <- associations[,2][match(collapsed_tree$tip.label, associations[,1])]
                    return(collapsed_tree)
            }

    ##### Other functions

        #### OsDirectoryPathCorrect

            #' OS-aware correction of directory paths
            #'
            #' @param directory_path The path to correct
            #' @examples
            #' @export
            #' OsDirectoryPathCorrect

            OsDirectoryPathCorrect <- function( directory_path ) {

                ## Detect OS

                    OS <- .Platform$OS.type

                ## Add terminal slash if missing

                    if (OS == "unix"){

                        if ( 
                            substr(
                                directory_path, 
                                nchar(directory_path), 
                                nchar(directory_path)
                            ) == "/" ) {
                        } else {
                            directory_path <- paste(directory_path, "/", sep = "")
                        }
                        directory_path_corrected <- directory_path

                    } else if (OS == "windows"){

                        if ( 
                            substr(
                                directory_path, 
                                nchar(directory_path), 
                                nchar(directory_path)
                            ) == "\\\\" ) {
                        } else {
                            directory_path <- paste(directory_path, "\\\\", sep = "")
                        }
                        directory_path_corrected <- directory_path

                    } else {

                        warning("ERROR: OS could not be identified")

                    }

                ## Return
                    
                    return(directory_path_corrected)

            }

        #### OsPathCorrect

            #' OS-aware correction of paths
            #'
            #' @param path The path to correct
            #' @examples
            #' @export
            #' OsPathCorrect

            OsPathCorrect <- function( path ) {

                ## Detect OS

                    OS <- .Platform$OS.type

                ## Replace single windows slashes with doubles

                    if (OS == "unix") {

                        path_corrected <- path

                    } else if (OS == "windows"){

                        path_corrected <- gsub("\\", "\\\\", path)

                    } else {

                        warning("ERROR: OS could not be identified")

                    }

                ## Return
                    
                    return(path_corrected)

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

        #### rotate_coord

            #' rotates a set of x,y coordinates
            #'
            #' @param x The x values
            #' @param y The y values
            #' @param angle The angle to rotate (in degrees)
            #' @param center The center about which to rotate
            #' @examples
            #' @export
            #' rotate_coord

            rotate_coord <- function( x, y, angle, center = c(0,0) ){
  
                x <- x - center[1]
                y <- y - center[2]
                angle <- angle*pi/180

                conversionmatrix <- matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)), ncol=2, nrow=2)
                xy <- cbind(x,y)%*%conversionmatrix

                xy[,1] <- xy[,1]+center[1]
                xy[,2] <- xy[,2]+center[2]

                return(xy)

            }

        #### question

            question <- function(question, distractors, correct, no, fb = "", print_question = TRUE){
                allanswers <- c(distractors, correct)[sample.int(length(distractors)+1)]
                correctanswer <- which(allanswers == correct)
                
                answercode <- paste0(sapply(1:length(allanswers), function(i){
                
                    x <- allanswers[i]
                    paste0('<div class="radio">\n  <label>\n    <input type="radio" name="question', no, '" id="opt', i,'" value="', i, '" onchange="check_answer', no, '()">\n    ', x, '\n  </label>\n</div>')
                
                }), collapse = "\n\n")

                out <- paste0(question, "\n\n", answercode,
                    '<div class="collapse" id="collapseExample', no,'">
                    <div class="card card-body" id="answerFeedback', no, '">
                    </div>
                    </div>',
                    paste0('<script type="text/javascript">
                    function check_answer', no, '()
                    {
                    var radioButtons', no, ' = document.getElementsByName("question', no, '");
                    document.getElementById("answerFeedback', no, '").innerHTML = "Try selecting an answer!!";
                    for(var i = 0; i < radioButtons', no, '.length; i++)
                    {
                    if(radioButtons', no, '[i].checked == true)
                    {
                    var feedback', no, ' = "<p style=\'color:red\'>Wrong', ifelse(fb == "", ".", paste0("; ", fb)),  '</p>";
                    if(radioButtons', no, '[i].value == "', correctanswer, '") {
                    feedback', no, ' = "<p style=\'color:green\'>Correct!</p>"
                    }
                    document.getElementById("answerFeedback', no, '").innerHTML = feedback', no, ';
                    return true;
                    }
                    }
                    }
                    </script>
                    '))
                
                if(print_question){
                    cat(out)
                } else {
                    return(out)
                }
            }

        #### questionnaire

            questionnaire <- function(x, shuffle = TRUE, print_question = TRUE){
                
                if(inherits(x, "character")) x <- read.csv(x, stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM")
                
                if(!all(names(x) == c("question", "distractors", "correct", "fb"))){
                    stop("Incorrect column names")
                }
                
                if(shuffle){
                    x <- x[sample.int(nrow(x)), ]
                }
                
                out <- ""
                for(i in 1:nrow(x)){
                    out <- paste0(out,
                    "**Question ",
                    i,
                    ":**\n",
                    question(question = x$question[i],
                    distractors = eval(parse(text = x$distractors[i])),
                    correct = x$correct[i],
                    no = i,
                    fb = ifelse(is.na(x$fb[i]), "", x$fb[i]),
                    print_question = FALSE),
                    "\n\n"
                    )
                }
                if(print_question){
                    cat(out)
                } else {
                    return(out)
                }
            }

        #### numbersOnly

            numbersOnly <- function(x) { all(suppressWarnings(!is.na(as.numeric(as.character(x))))) }

    ##### Sequence data handling

        #### extractORFs

            #' Extract open reading frames from multifasta files
            #'
            #' @param file_in_path The multifasta file to analyze
            #' @param write_out_ORFs TRUE/FALSE whether to write out a new fasta that contains the ORFs
            #' @param overwrite TRUE/FALSE Optionally overwrite the input file with the ORF file
            ####### Searching for ORFs in the reverse direction leads to issues during codon alignment...
            #' @examples
            #' @export
            #' extractORFs

            extractORFs <- function(
                file_in_path,
                write_out_ORFs = FALSE,
                overwrite = FALSE,
                alternative_ORFs = FALSE
            ) {

                fasta <- seqinr::read.fasta(file_in_path)
                ORFs <- list()
                ORFs_to_write_out <- Biostrings::DNAStringSet()

                if ( write_out_ORFs == TRUE ) {
                    if (file.exists(paste(file_in_path, "_orfs", sep = ""))) { file.remove(paste(file_in_path, "_orfs", sep = "")) }
                }

                for (j in 1:length(fasta)) {
                    orf_coordinates <- data.frame()
                    mRNA <- unlist(fasta[j])

                    ## Search for ORFs in the forward sequence
                        data <- seqinr::translate(mRNA, frame = 0, sens = "F")
                            S <- grep("\\*", data) # find all stop codons
                            M <- grep("M", data) # find all start codons
                            if (length(S)>0 & length(M)>0) { # if there are stop codons, remove all start codons after the last stop codon
                                M <- M[!M > max(S)]
                            }
                            if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                    F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                    orf_coordinates <-  rbind(orf_coordinates,  data.frame(
                                                                                    start_codon = (M[i]*3-2), 
                                                                                    stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                    direction = "forward"
                                                                                )
                                                        )
                                }
                            }
                        data <- seqinr::translate(mRNA, frame = 1, sens = "F")
                            S <- grep("\\*", data) # find all stop codons
                            M <- grep("M", data) # find all start codons
                            if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                M <- M[!M > max(S)] # remove all start codons after the last stop codon
                            }
                            if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                    F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                    orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                    start_codon = (M[i]*3-2+1), 
                                                                                    stop_codon = ((M[i]*3-2+1)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                    direction = "forward"
                                                                                )
                                                        )
                                }
                            }
                        data <- seqinr::translate(mRNA, frame = 2, sens = "F")
                            S <- grep("\\*", data) # find all stop codons
                            M <- grep("M", data) # find all start codons
                            if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                M <- M[!M > max(S)] # remove all start codons after the last stop codon
                            }
                            if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                    F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                    orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                    start_codon = (M[i]*3-2+2),
                                                                                    stop_codon = ((M[i]*3-2+2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                    direction = "forward"
                                                                                )
                                                        )
                                }
                            }

                    if (alternative_ORFs) {

                        # Search for ORFs in the compliment mRNA
                            mRNA_comp <- chartr("ATGC","TACG", mRNA)
                            mRNA_comp <- chartr("atcg","tacg", mRNA_comp)
                            data <- seqinr::translate(mRNA_comp, frame = 0, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "compliment"
                                                                                    )
                                                            )
                                    }
                                }
                            data <- seqinr::translate(mRNA_comp, frame = 1, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "compliment"
                                                                                    )
                                                            )
                                    }
                                }
                            data <- seqinr::translate(mRNA_comp, frame = 2, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "compliment"
                                                                                    )
                                                            )
                                    }
                                }

                        # Search for ORFs in the reverse mRNA
                            mRNA_rev <- rev(mRNA)
                            data <- seqinr::translate(mRNA_rev, frame = 0, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "reverse"
                                                                                    )
                                                            )
                                    }
                                }
                            data <- seqinr::translate(mRNA_rev, frame = 1, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "reverse"
                                                                                    )
                                                            )
                                    }
                                }
                            data <- seqinr::translate(mRNA_rev, frame = 2, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "reverse"
                                                                                    )
                                                            )
                                    }
                                }

                        # Search for ORFs in the reverse compliment mRNA
                            mRNA_rev_comp <- chartr("ATGC","TACG", mRNA)
                            mRNA_rev_comp <- chartr("atcg","tacg", mRNA_rev_comp)
                            mRNA_rev_comp <- rev(mRNA_rev_comp)
                            data <- seqinr::translate(mRNA_rev_comp, frame = 0, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "reverse_compliment"
                                                                                    )
                                                            )
                                    }
                                }
                            data <- seqinr::translate(mRNA_rev_comp, frame = 1, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "reverse_compliment"
                                                                                    )
                                                            )
                                    }
                                }
                            data <- seqinr::translate(mRNA_rev_comp, frame = 2, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "reverse_compliment"
                                                                                    )
                                                            )
                                    }
                                }

                    }

                    ## Process discovered ORFs

                        if ( dim(orf_coordinates)[1] == 0 ) {
                            ORFs[[j]] <-    data.frame(
                                                accession = names(fasta[j]),
                                                start_codon = 0, 
                                                stop_codon = 0,
                                                direction = "forward",
                                                orf_length = 0
                                            )
                        } else {
                            rownames(orf_coordinates) <- NULL # Reset rownames because of previous step
                            orf_coordinates$orf_length <- (orf_coordinates$stop_codon - orf_coordinates$start_codon + 1)
                            orf_coordinates <- orf_coordinates[order(orf_coordinates$orf_length, decreasing = TRUE),]
                            ORFs[[j]] <- cbind(data.frame(accession = names(fasta[j])), orf_coordinates)
                        }

                        if ( write_out_ORFs == TRUE) {
                            if ( dim(orf_coordinates)[1] == 0 ) {
                                ORFs_to_write_out[j] <- c("a","t","g","t","a","a") ## Weird, but necessary so errors are not thrown during codon alignment
                                ORFs_to_write_out@ranges@NAMES[j] <- attr(fasta[j], "name")
                            } else {
                                if ( orf_coordinates$direction[1] == "forward" ) {
                                    ORFs_to_write_out <- c(ORFs_to_write_out, Biostrings::DNAStringSet(paste(mRNA[orf_coordinates$start_codon[1]:orf_coordinates$stop_codon[1]], collapse = ""))) # extract longest orf from the mRNA
                                    ORFs_to_write_out@ranges@NAMES[j] <- attr(fasta[j], "name")
                                }
                                # if ( orf_coordinates$direction[1] == "reverse" ) {
                                #     orf <- rev(mRNA)[orf_coordinates$start_codon[1]:orf_coordinates$stop_codon[1]] # extract longest orf from the mRNA
                                #     seqinr::write.fasta(orf, names = attr(fasta[j], "name"), file = paste(file, "_orfs", sep = ""), open = "a") # append ORF to the list
                                # }
                            }
                        }
                }

                # seqinr::write.fasta(fake_orf, names = attr(fasta[j], "name"), file = paste(file, "_orfs", sep = ""), open = "a") # append fake_ORF to the list
                if ( write_out_ORFs == TRUE) {
                    if (overwrite == TRUE) {
                        Biostrings::writeXStringSet(ORFs_to_write_out, file = file_in_path)
                        # seqinr::write.fasta(orf, names = attr(fasta[j], "name"), file = paste(file, "_orfs", sep = ""), open = "w") # overwrite transcript file
                    } else {
                        Biostrings::writeXStringSet(ORFs_to_write_out, file = paste(file_in_path, "_orfs", sep = ""))
                        # seqinr::write.fasta(orf, names = attr(fasta[j], "name"), file = paste(file, "_orfs", sep = ""), open = "a") # append ORF to the list        
                    }
                }

                ORFs <- do.call(rbind, ORFs)
                rownames(ORFs) <- NULL
                return(ORFs)
            }

        #### polyBlast

            #' BLAST search local sequence collections
            #'
            #' Search locally stored sequence collections for query sequences. Download Blast+ here: https://www.ncbi.nlm.nih.gov/books/NBK279671/
            #' @param named_subjects_list A list of paths to the transcriptomes that should be searched, named by taxonomic identifier (e.g. species names)
            #' @param query_in_path Path to a fasta file containing the query or queries
            #' @param sequences_of_interest_directory_path Path to a directory where blast hits should be written out as fasta files
            #' @param blast_module_directory_path Path to directory containing the BLAST+ module (perhaps something like "/usr/local/ncbi/blast/bin/")
            #' @param blast_mode One of "blastn" or "dc-megablast". "blastn" is a traditional BLASTN requiring an exact match of 11. "dc-megablast" is a discontiguous megablast used to find more distant (e.g., interspecies) sequences.
            #' @param e_value_cutoff e-value cutoff to apply to results. Defaults to 1. Set to 1e-10 for heavy filtering.
            #' @param queries_in_output Should the queries be included in the fasta output?
            #' @param monolist_out_path Path to where the output monolist should be written
            #' @examples
            #' @export
            #' polyBlast

            polyBlast <- function(
                                    named_subjects_list,
                                    query_in_path,
                                    sequences_of_interest_directory_path,
                                    blast_module_directory_path,
                                    blast_mode = c("nnblastn", "ntblastp", "pnblastp"), 
                                    e_value_cutoff = 1,
                                    queries_in_output = TRUE,
                                    monolist_out_path
                                ) {

                ### Make sure transcriptomes object is character and has names

                    if( length(names(transcriptomes)) != length(transcriptomes) ) {
                        stop("Please provide a set of named paths to the argument `transcriptomes`")
                    }
                    
                    names <- names(transcriptomes)
                    transcriptomes <- as.character(transcriptomes)
                    names(transcriptomes) <- names

                ### If *tblast*, then translate all transcriptomes and rename the translations with "_trans" suffix

                    if ( substr(blast_mode, 2, 2) == "t" ) {
                        cat("Translating transcriptomes...\n\n")
                        for ( transcriptome in 1:length(transcriptomes) ) {
                            writeFasta(
                                XStringSet = translate(readDNAStringSet(transcriptomes[transcriptome]), if.fuzzy.codon = "solve"),
                                fasta_out_path = paste(transcriptomes[transcriptome], "trans", sep = "_"),
                                type = "AA"
                            )
                            transcriptomes[transcriptome] <- paste(transcriptomes[transcriptome], "trans", sep = "_")
                        }
                    } else if ( substr(blast_mode, 2, 2) == "n" ) {

                    } else {
                        stop("Please specify t or n to indicate whether the subjects should be translated.")
                    }

                ### Build blast database(s)

                    if ( substr(blast_mode, 1, 1) %in% c("n") ) {
                        dbtype = "nucl"
                    } else if ( substr(blast_mode, 1, 1) %in% c("p") ) {
                        dbtype = "prot"
                    } else {
                        stop("Please specify a subject type of n or p")
                    }

                    if ( .Platform$OS == "unix" ) {

                        for (transcriptome in 1:length(transcriptomes)) {
                            cat(paste0("\nSpecies: ", names(transcriptomes)[transcriptome]))
                            system(
                                paste(
                                    blast_module_directory_path,
                                    "makeblastdb -in ",
                                    transcriptomes[transcriptome],
                                    " -dbtype ",
                                    dbtype,
                                    sep = ""
                                ) 
                            )
                        }

                    } else if ( .Platform$OS == "windows" ) {

                        for (transcriptome in 1:length(transcriptomes)) {
                            cat(paste0("\nSpecies: ", names(transcriptomes)[transcriptome]))
                            shell(
                                paste(
                                    blast_module_directory_path,
                                    "makeblastdb -in ",
                                    transcriptomes[transcriptome],
                                    " -dbtype ",
                                    dbtype,
                                    sep = ""
                                )
                            )
                        }
                    }
                    
                ### Start BLAST process

                    query_seqs <- Biostrings::readBStringSet(filepath = query_in_path, format = "fasta")
                    if ( substr(blast_mode, 8, 8) == "p") { blast_type <- "blastp" }
                    if ( substr(blast_mode, 8, 8) == "n") { blast_type <- "blastn" }

                        ## Loop over each member of the query and use it to blast each transcriptome
                            
                            monolist <- data.frame()
                            cat("\n\n")
                            
                            for ( query_seq in 1:length(query_seqs) ) {

                                ## Write individual files for each member of the query
                                    
                                    Biostrings::writeXStringSet(
                                        query_seqs[query_seq], 
                                        filepath = paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""), 
                                        append = FALSE
                                    )

                                ## Loop over the transcriptomes, run the blast on each, add hits to monolist
                                    
                                    for (transcriptome in 1:length(transcriptomes)) {

                                        ## Run BLASTs on unix system

                                            if ( .Platform$OS == "unix") {

                                                system(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -outfmt '6 sallacc'",
                                                        sep = ""
                                                    )
                                                )

                                                system(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_length", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -outfmt '6 length'",
                                                        sep = ""
                                                    )
                                                )

                                                system(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_pident", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -outfmt '6 pident'",
                                                        sep = ""
                                                    )
                                                )

                                                system(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_evalue", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -outfmt '6 evalue'",
                                                        sep = ""
                                                    )
                                                )

                                                system(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_bitscore", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -outfmt '6 bitscore'",
                                                        sep = ""
                                                    )
                                                )
                                            }

                                        ## Run BLASTs on windows system

                                            if ( .Platform$OS == "windows") {

                                                shell(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -outfmt \"6 sallacc\"",
                                                        sep = ""
                                                    )
                                                )

                                                shell(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_length", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -outfmt \"6 length\"",
                                                        sep = ""
                                                    )
                                                )

                                                shell(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_pident", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -outfmt \"6 pident\"",
                                                        sep = ""
                                                    )
                                                )

                                                shell(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_evalue", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -outfmt \"6 evalue\"",
                                                        sep = ""
                                                    )
                                                )

                                                shell(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_bitscore", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -outfmt \"6 bitscore\"",
                                                        sep = ""
                                                    )
                                                )
                                            }

                                        ## Extract BLAST hits from transcriptome, add them to the monolist, write them to individual files

                                            # Read in whole transcriptome
                                                if ( substr(blast_mode, 8, 8) == "p" ) {
                                                    temp_seqs <- Biostrings::readAAStringSet( 
                                                        filepath = as.character(transcriptomes)[transcriptome], 
                                                        format = "fasta"
                                                    )
                                                } else if ( substr(blast_mode, 8, 8) == "n" ) {
                                                    temp_seqs <- Biostrings::readDNAStringSet( 
                                                        filepath = as.character(gsub("_trans", "", transcriptomes)[transcriptome]), 
                                                        format = "fasta"
                                                    )
                                                }
                                        
                                            # Attempt to read in hits list, subset transcriptome according to that list, add hits to monolist
                                                cat(paste0("Query ", query_seq, ": "))
                                                if ( inherits( try( read.table(file = paste(transcriptomes[transcriptome], ".out", sep = "")), silent = TRUE), "try-error") == TRUE ) { # Skips over empty files
                                                    cat(paste("No BLAST hits found for ", query_seqs[query_seq]@ranges@NAMES, " in ", names(transcriptomes)[transcriptome], "\n", sep = ""))
                                                } else {

                                                    ## Remove duplicate hits from all output files
                                                        temp_hits <- as.character(read.table(file = paste(transcriptomes[transcriptome], ".out", sep = ""))[,1])
                                                        temp_hits_length <- as.character(read.table(file = paste(transcriptomes[transcriptome], ".out_length", sep = ""))[,1])
                                                        temp_hits_pident <- as.character(read.table(file = paste(transcriptomes[transcriptome], ".out_pident", sep = ""))[,1])
                                                        temp_hits_evalue <- as.character(read.table(file = paste(transcriptomes[transcriptome], ".out_evalue", sep = ""))[,1])
                                                        temp_hits_bitscore <- as.character(read.table(file = paste(transcriptomes[transcriptome], ".out_bitscore", sep = ""))[,1])
                                                        
                                                        duplicate_indeces <- duplicated(temp_hits)

                                                        temp_hits <- temp_hits[!duplicate_indeces]
                                                        temp_hits_length <- temp_hits_length[!duplicate_indeces]
                                                        temp_hits_pident <- temp_hits_pident[!duplicate_indeces]
                                                        temp_hits_evalue <- temp_hits_evalue[!duplicate_indeces]
                                                        temp_hits_bitscore <- temp_hits_bitscore[!duplicate_indeces]
                                                        
                                                        writeMonolist(temp_hits, paste(transcriptomes[transcriptome], ".out", sep = ""))
                                                        writeMonolist(temp_hits_length, paste(transcriptomes[transcriptome], ".out_length", sep = ""))
                                                        writeMonolist(temp_hits_pident, paste(transcriptomes[transcriptome], ".out_pident", sep = ""))
                                                        writeMonolist(temp_hits_evalue, paste(transcriptomes[transcriptome], ".out_evalue", sep = ""))
                                                        writeMonolist(temp_hits_bitscore, paste(transcriptomes[transcriptome], ".out_bitscore", sep = ""))

                                                    temp_hits <- readMonolist(paste(transcriptomes[transcriptome], ".out", sep = ""))
                                                    cat(paste(
                                                        "Found ", 
                                                        dim(temp_hits)[1], 
                                                        " BLAST hits found for ", 
                                                        query_seqs[query_seq]@ranges@NAMES, 
                                                        " in ", 
                                                        names(transcriptomes)[transcriptome], 
                                                        "\n", 
                                                        sep = "")
                                                    )

                                                    accession <- gsub(" .*", "", temp_seqs@ranges@NAMES)
                                                    temp_seqs <- temp_seqs[accession %in% temp_hits[,1]]
                                                    accession <- gsub(" .*", "", temp_seqs@ranges@NAMES)
                                                    temp_seqs@ranges@NAMES <- gsub(" .*", "", temp_seqs@ranges@NAMES)

                                                    if ( all(temp_hits[,1] %in% accession) != TRUE ) {
                                                        warning("Couldn't find some BLAST hits within the transcriptome!")
                                                    }

                                                    # If there are hits, add them to the monolist
                                                        if (length(temp_seqs) > 0) {

                                                            # Write out blast hits
                                                                longest_ORFs <- vector()
                                                                for (temp_seq in 1:length(temp_seqs)) {
                                                                    Biostrings::writeXStringSet(
                                                                        temp_seqs[temp_seq], 
                                                                        filepath = paste(sequences_of_interest_directory_path, accession[temp_seq], ".fa", sep = ""), 
                                                                        append = FALSE
                                                                    )
                                                                    longest_ORFs <- c(longest_ORFs, extractORFs(paste(sequences_of_interest_directory_path, accession[temp_seq], ".fa", sep = ""))$orf_length[1])
                                                                }

                                                            # Add hit information to the monolist
                                                                monolist <- rbind(monolist, data.frame(
                                                                    accession = gsub(" .*", "", temp_seqs@ranges@NAMES),
                                                                    Genus = as.character(gsub("_.*$", "", names(transcriptomes)[transcriptome])),
                                                                    species = gsub(".*_", "", names(transcriptomes)[transcriptome]),
                                                                    annotation = temp_seqs@ranges@NAMES,
                                                                    length = temp_seqs@ranges@width,
                                                                    longestORF = longest_ORFs,
                                                                    length_aligned_with_query = readMonolist(paste(transcriptomes[transcriptome], ".out_length", sep = ""))[,1],
                                                                    percent_identity = readMonolist(paste(transcriptomes[transcriptome], ".out_pident", sep = ""))[,1],
                                                                    e_value = readMonolist(paste(transcriptomes[transcriptome], ".out_evalue", sep = ""))[,1],
                                                                    bitscore = readMonolist(paste(transcriptomes[transcriptome], ".out_bitscore", sep = ""))[,1],
                                                                    subset_all = TRUE,
                                                                    query = query_seqs@ranges@NAMES[query_seq],
                                                                    query_length = query_seqs@ranges@width[query_seq],
                                                                    query_longestORF = extractORFs(paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""))$orf_length[1]
                                                                ))
                                                        }
                                                }
                                    }

                                ## Optionally add the queries to the output

                                    if (queries_in_output == TRUE) {

                                        if (blast_type == "blastp") {
                                                longestORF <- 3*query_seqs@ranges@width[query_seq] 
                                            } else {
                                                longestORF <- extractORFs(paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""))$orf_length[1]
                                            }
                                        
                                        monolist <- rbind(monolist, data.frame(
                                            accession = query_seqs@ranges@NAMES[query_seq],
                                            Genus = "Query",
                                            species = "query",
                                            annotation = "query",
                                            length = query_seqs@ranges@width[query_seq],
                                            longestORF = longestORF,
                                            length_aligned_with_query = query_seqs@ranges@width[query_seq],
                                            percent_identity = 100,
                                            e_value = 0,
                                            bitscore = 1000,
                                            subset_all = TRUE,
                                            query = query_seqs@ranges@NAMES[query_seq],
                                            query_length = query_seqs@ranges@width[query_seq],
                                            query_longestORF = extractORFs(paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""))$orf_length[1]
                                        ))

                                        Biostrings::writeXStringSet(
                                            query_seqs[query_seq], 
                                            filepath = paste(sequences_of_interest_directory_path, query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""), 
                                            append = FALSE
                                        )
                                    }

                                ## Remove individual query files
                                    if (file.exists(paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = "")) ) {
                                        file.remove(paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""))
                                    }
                            }

                        ## Write out the monolist

                            monolist <- unique(monolist)
                            rownames(monolist) <- NULL
                            if (file.exists(monolist_out_path)) {file.remove(monolist_out_path)}
                            writeMonolist(monolist = monolist, monolist_out_path = monolist_out_path)
               
                    cat("\nDone!\n\n")
            }

        #### alignSequences

            #' Align sequences in a seqlist
            #'
            #' @param monolist Monolist of the sequences to be aligned. First column should be "accession"
            #' @param subset TRUE/FALSE column in monolist that specifies which sequences should be included in the alignment
            #' @param alignment_directory_path Path to where the alignment should be written
            #' @param sequences_of_interest_directory_path Path to a directory where blast hits should be written out as fasta files
            #' @param input_sequence_type One of "nucl" or "amin"
            #' @param mode One of "nucl_align", "amin_align", or "codon_align"
            #' @param base_fragment
            #' TROUBLESHOOTING: 
            #'      "ERROR: inconsistency between the following pep and nuc seqs" - usually means there are duplicate accessions numbers in the input monolist
            #' @import msa
            #' @examples
            #' @export
            #' alignSequences

            alignSequences <-   function(
                                    monolist, 
                                    subset, 
                                    alignment_directory_path, 
                                    sequences_of_interest_directory_path,
                                    input_sequence_type = c("nucl", "amin"), 
                                    mode = c("nucl_align", "amin_align", "codon_align", "fragment_align"),
                                    base_fragment = NULL
                                ){  

                ## Check directory_paths
                    sequences_of_interest_directory_path <- OsDirectoryPathCorrect(sequences_of_interest_directory_path)
                    alignment_directory_path <- OsDirectoryPathCorrect(alignment_directory_path)

                ## Get appropriate subset of the monolist
                    monolist_subset <- monolist[unlist(monolist[,as.character(colnames(monolist)) == as.character(subset)]),]

                ## If starting with nucleotide sequences
                    
                    if ( input_sequence_type == "nucl" ) {

                        ## Remove existing version of this alignment and it's files
                            if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs.fa", sep = ""))) {
                                file.remove(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs.fa", sep = ""))}
                            if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs_aligned.fa", sep = ""))) {
                                file.remove(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs_aligned.fa", sep = ""))}

                        ## Import the subset's sequences, then write out *_nucl_seqs.fa
                            nucl_seqs_set <- Biostrings::DNAStringSet()
                            for (i in 1:dim(monolist_subset)[1]) {
                                nucl_seqs_set <- c(nucl_seqs_set, Biostrings::readDNAStringSet(paste(sequences_of_interest_directory_path, as.character(monolist_subset$accession[i]), ".fa", sep = "")))
                            }
                            Biostrings::writeXStringSet(nucl_seqs_set, filepath = paste(alignment_directory_path, subset, "_nucl_seqs.fa", sep = ""))

                        ## Run plain nucleotide alignment, if requested
                            if ( mode == "nucl_align") {
                                nucl_seqs_set_aligned <- msa::msa(nucl_seqs_set, order = c("input"))
                                msa <- msa::msaConvert(nucl_seqs_set_aligned, type = "seqinr::alignment")
                                n <- dim(as.data.frame(msa$seq))[1]
                                for (i in 1:n){msa$seq[i] <- substr(msa$seq[i],0,nchar(msa$seq[1]))}
                                seqinr::write.fasta(as.list(msa$seq),as.list(msa$nam), file.out = paste(alignment_directory_path, subset, "_nucl_seqs_aligned.fa", sep = ""))
                            }

                        ## Run codon alignment, if requested
                            if ( mode == "codon_align") {

                                ## Remove existing codon alignment
                                    if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs_codon_aligned.fa", sep = ""))) {
                                    file.remove(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs_codon_aligned.fa", sep = ""))}
                                
                                ## Extarct ORFs
                                    extractORFs(file = paste(alignment_directory_path, subset, "_nucl_seqs.fa", sep = ""), write_out_ORFs = TRUE)
                                
                                ## Translate the nucleotide sequences and write out *_amin_seqs.fa
                                    nucl_seqs_set <- Biostrings::readDNAStringSet(paste(alignment_directory_path, subset, "_nucl_seqs.fa_orfs", sep = ""))
                                    amin_seqs_set <- Biostrings::translate(nucl_seqs_set, if.fuzzy.codon = "solve")
                                    if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))) {
                                        file.remove(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))}
                                    Biostrings::writeXStringSet(amin_seqs_set, filepath = paste(alignment_directory_path, subset, "_amin_seqs.fa", sep = ""))

                                ## Run amino acid alignment and write out *_amin_seqs.fa
                                    amin_seqs_set_aligned <- msa::msa(amin_seqs_set, order = c("input"))
                                    msa <- msa::msaConvert(amin_seqs_set_aligned, type = "seqinr::alignment")
                                    n <- dim(as.data.frame(msa$seq))[1]
                                    for (i in 1:n){msa$seq[i] <- substr(msa$seq[i],0,nchar(msa$seq[1]))}
                                    if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs_aligned.fa", sep = ""))) {
                                    file.remove(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs_aligned.fa", sep = ""))}
                                    seqinr::write.fasta(as.list(msa$seq),as.list(msa$nam), file.out = paste(alignment_directory_path, subset, "_amin_seqs_aligned.fa", sep = ""))

                                ## Run codon alignment and write to nucl_seqs_aligned.fa
                                    nucl_seqs_codon_aligned <-  orthologr::codon_aln(
                                                                    file_aln = paste(alignment_directory_path, subset, "_amin_seqs_aligned.fa", sep = ""), 
                                                                    file_nuc = paste(alignment_directory_path, subset, "_nucl_seqs.fa", sep = ""), 
                                                                    get_aln = TRUE
                                                                )
                                    msa <- nucl_seqs_codon_aligned
                                    n <- dim(as.data.frame(msa$seq))[1]
                                    for (i in 1:n){msa$seq[i] <- substr(msa$seq[i],0,nchar(msa$seq[1]))}
                                    seqinr::write.fasta(as.list(msa$seq),as.list(msa$nam), file.out = paste(alignment_directory_path, subset, "_nucl_seqs_codon_aligned.fa", sep = ""))
                            }

                        ## Run fragment alignment, if requested
                            if ( mode == "fragment_align") {

                                ## Remove existing fragment alignment
                                    if (file.exists(paste(alignment_directory_path, as.character(subset), "_fragment_seqs_aligned.fa", sep = ""))) {
                                    file.remove(paste(alignment_directory_path, as.character(subset), "_fragment_seqs_aligned.fa", sep = ""))}

                                ## Extarct ORFs
                                    extractORFs(file = paste(alignment_directory_path, subset, "_nucl_seqs.fa", sep = ""), write_out_ORFs = TRUE)
                                
                                ## Translate the nucleotide sequences and write out *_amin_seqs.fa
                                    nucl_seqs_set <- Biostrings::readDNAStringSet(paste(alignment_directory_path, subset, "_nucl_seqs.fa_orfs", sep = ""))
                                    amin_seqs_set <- Biostrings::translate(nucl_seqs_set, if.fuzzy.codon = "solve")
                                    if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))) {
                                        file.remove(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))}
                                    Biostrings::writeXStringSet(amin_seqs_set, filepath = paste(alignment_directory_path, subset, "_amin_seqs.fa", sep = ""))

                                ## Define base fragment
                                    base_fragment_seq <- amin_seqs_set[amin_seqs_set@ranges@NAMES == base_fragment]
                                    fragments_seq_set <- amin_seqs_set[amin_seqs_set@ranges@NAMES != base_fragment]
                                    current_base_fragment_seq <- base_fragment_seq
                                
                                ## Align each fragment on its own with the base fragment
                                    aligned_fragments <- AAStringSet()
                                    for ( fragment in 1:length(fragments_seq_set) ) {
                                        fragment_pair_seq_set <- c(base_fragment_seq, fragments_seq_set[fragment])
                                        fragment_pair_seq_set_aligned <- msa::msa(fragment_pair_seq_set, order = c("input"))
                                        fragment_pair_seq_set_aligned <- AAStringSet(fragment_pair_seq_set_aligned)
                                        current_base_fragment_seq <- fragment_pair_seq_set_aligned[fragment_pair_seq_set_aligned@ranges@NAMES == base_fragment]
                                        aligned_fragments <- c(aligned_fragments, fragment_pair_seq_set_aligned[fragment_pair_seq_set_aligned@ranges@NAMES != base_fragment])
                                    }

                                ## Perform final alignment and write it out
                                    fragment_seq_set <- c(current_base_fragment_seq, aligned_fragments)
                                    fragment_seq_set <- AAStringSet(gsub("-", "Z", fragment_seq_set))
                                    fragment_seq_set_aligned <- msa::msa(fragment_seq_set, order = c("input"))
                                    fragment_seq_set_aligned <- AAStringSet(fragment_seq_set_aligned)
                                    writeXStringSet(fragment_seq_set_aligned, filepath = paste(alignment_directory_path, subset, "_fragment_seqs_aligned.fa", sep = ""))
                            }

                        ## Run amin alignment, if requested
                            if ( mode == "amin_align") {

                                ## Remove existing amin alignment
                                    if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs_amin_aligned.fa", sep = ""))) {
                                    file.remove(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs_amin_aligned.fa", sep = ""))}
                                
                                ## Extarct ORFs
                                    extractORFs(file = paste(alignment_directory_path, subset, "_nucl_seqs.fa", sep = ""), write_out_ORFs = TRUE)
                                
                                ## Translate the nucleotide sequences and write out *_amin_seqs.fa
                                    nucl_seqs_set <- Biostrings::readDNAStringSet(paste(alignment_directory_path, subset, "_nucl_seqs.fa_orfs", sep = ""))
                                    amin_seqs_set <- Biostrings::translate(nucl_seqs_set, if.fuzzy.codon = "solve")
                                    if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))) {
                                        file.remove(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))}
                                    Biostrings::writeXStringSet(amin_seqs_set, filepath = paste(alignment_directory_path, subset, "_amin_seqs.fa", sep = ""))

                                ## Run amino acid alignment and write out *_amin_seqs.fa
                                    amin_seqs_set_aligned <- msa::msa(amin_seqs_set, order = c("input"))
                                    msa <- msa::msaConvert(amin_seqs_set_aligned, type = "seqinr::alignment")
                                    n <- dim(as.data.frame(msa$seq))[1]
                                    for (i in 1:n){msa$seq[i] <- substr(msa$seq[i],0,nchar(msa$seq[1]))}
                                    if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs_aligned.fa", sep = ""))) {
                                    file.remove(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs_aligned.fa", sep = ""))}
                                    seqinr::write.fasta(as.list(msa$seq),as.list(msa$nam), file.out = paste(alignment_directory_path, subset, "_amin_seqs_aligned.fa", sep = ""))
                            }
                    }

                ## If starting with amino acid sequences

                    if ( input_sequence_type == "amin" ) {

                        if ( mode == "amin_align" ) {

                            ## Remove existing version of this alignment and it's files
                                if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))) {
                                    file.remove(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))}
                                if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs_aligned.fa", sep = ""))) {
                                    file.remove(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs_aligned.fa", sep = ""))}

                            ## Import the subset's sequences, then write out *_amin_seqs.fa
                                amin_seqs_set <- Biostrings::AAStringSet()
                                for (i in 1:dim(monolist_subset)[1]) {
                                    amin_seqs_set <- c(amin_seqs_set, Biostrings::readAAStringSet(paste(sequences_of_interest_directory_path, as.character(monolist_subset$accession[i]), ".fa", sep = "")))
                                }
                                Biostrings::writeXStringSet(amin_seqs_set, filepath = paste(alignment_directory_path, subset, "_amin_seqs.fa", sep = ""))

                            ## Run the alignment
                                amin_seqs_set_aligned <- msa::msa(amin_seqs_set, order = c("input"))
                                msa <- msa::msaConvert(amin_seqs_set_aligned, type = "seqinr::alignment")
                                n <- dim(as.data.frame(msa$seq))[1]
                                for (i in 1:n){msa$seq[i] <- substr(msa$seq[i],0,nchar(msa$seq[1]))}
                                seqinr::write.fasta(as.list(msa$seq),as.list(msa$nam), file.out = paste(alignment_directory_path, subset, "_amin_seqs_aligned.fa", sep = ""))
                        }
                    }
            }

        #### analyzeAlignment
        
            #' Analyze a multiple sequence alignment and subset it for phylogeny building
            #'
            #' @param alignment_in_path The alignment to analyze, in fasta format
            #' @param type DNA or AA
            #' @examples
            #' @export
            #' extractORFs

            analyzeAlignment <- function(
                alignment_in_path,
                type = c("DNA", "AA"),
                jupyter = FALSE
            ) {

                ## Prelim stuff

                    ## Read in and process alignment
                        
                        alignment <- readAlignment(
                            alignment_in_path = alignment_in_path,
                            type = type
                        )

                        if (type == "DNA") {scaffold_type = "nucl_alignment"}
                        if (type == "AA") {scaffold_type = "amin_alignment"}

                        tree <- fortify(buildTree(
                            scaffold_type = scaffold_type,
                            scaffold_in_path = alignment_in_path,
                        ))

                        alignment$name <- factor(alignment$name, levels = arrange(filter(tree, isTip == TRUE), y)$label)

                    ## Analyze composition of each position
                                
                        alignment %>%
                            group_by(position) %>%
                            summarize(
                                gap_percent = sum(str_count(state, "-"))/n()*100,
                                conservation_percent = suppressWarnings(sum(state == mode(state))/n()*100)
                            ) %>%
                            pivot_longer(cols = c(2,3), names_to = "metric", values_to = "value") -> alignment_stats

                ## User interface

                    ui <- fluidPage(

                        sidebarLayout(

                            sidebarPanel(

                                sliderInput(
                                    inputId = "x_start",
                                    label = "x_start",
                                    min = 0,
                                    max = max(alignment$position),
                                    value = 0
                                ),

                                sliderInput(
                                    inputId = "x_end",
                                    label = "x_end",
                                    min = 0,
                                    max = max(alignment$position),
                                    value = max(alignment$position)
                                ),

                                sliderInput(
                                    inputId = "gap_threshold",
                                    label = "Maximum Gap Percent",
                                    min = 0,
                                    max = 100,
                                    value = 50
                                ),

                                sliderInput(
                                    inputId = "conservation_threshold",
                                    label = "Minimum Conservation Percent",
                                    min = 0,
                                    max = 100,
                                    value = 50
                                ),

                                "Here is a NJ representation of the original alignment:",

                                plotOutput(
                                    outputId = "tree",
                                    height = "300px"
                                ),

                                "Here is a NJ representation of the alignment AFTER trimming:",

                                plotOutput(
                                    outputId = "filtered_tree",
                                    height = "300px"
                                ),

                                tags$head(
                                    HTML(
                                        "
                                        <script>
                                            var socket_timeout_interval
                                            var n = 0
                                            $(document).on('shiny:connected', function(event) {
                                            socket_timeout_interval = setInterval(function(){
                                            Shiny.onInputChange('count', n++)
                                            }, 15000)
                                            });
                                            $(document).on('shiny:disconnected', function(event) {
                                            clearInterval(socket_timeout_interval)
                                            });
                                        </script>
                                        "
                                    )
                                ),

                                textOutput("keepAlive"),

                                width = 3

                            ),

                            mainPanel(

                                plotOutput(
                                    outputId = "plot",
                                    height = "1000px"
                                )
                            )
                        )
                    )

                ## Server

                    server <- function(input, output) {

                        ## Don't let it time out
                            
                            output$keepAlive <- renderText({
                                req(input$count)
                                paste("\nstayin' alive ", input$count)
                            })

                        ## Original tree

                            output$tree <- renderPlot({ ggtree(tree) })

                        ## Plot alignment

                            output$plot <- renderPlot({

                                ## Main plots

                                    stats_plot <- ggplot(alignment_stats, aes(x = position, y = value, color = metric)) +
                                        geom_line(size = 1) +
                                        geom_hline(data = data.frame(
                                            metric = "conservation_percent",
                                            yintercept = input$conservation_threshold
                                            ), aes(yintercept = yintercept)
                                        ) +
                                        geom_hline(data = data.frame(
                                            metric = "gap_percent",
                                            yintercept = input$gap_threshold
                                            ), aes(yintercept = yintercept)
                                        ) +
                                        scale_color_viridis(discrete = TRUE, guide = "none") +
                                        theme_classic() +
                                        scale_y_continuous(name = "Percent") +
                                        scale_x_continuous(expand = c(0,0)) +
                                        coord_cartesian(xlim = c(input$x_start, input$x_end)) +
                                        facet_grid(metric~.) +
                                        theme(
                                            axis.text.y = element_blank()
                                        )
                                
                                    alignment_plot <- ggplot(alignment, aes(x = position, y = name)) +
                                        geom_tile(aes(fill = state)) +
                                        scale_fill_viridis(discrete = TRUE, guide = "none") +
                                        scale_y_discrete(name = "") +
                                        scale_x_continuous(expand = c(0,0)) +
                                        theme_classic() +
                                        coord_cartesian(xlim = c(input$x_start, input$x_end)) +
                                        theme(
                                            axis.text.y = element_blank()
                                        )

                                    alignment_stats %>%
                                        pivot_wider(names_from = "metric", values_from = "value") %>%
                                        filter(
                                            gap_percent < input$gap_threshold &
                                            conservation_percent > input$conservation_threshold
                                        ) %>% 
                                        select(position) %>% as.data.frame() -> positions_to_keep
                                        positions_to_keep <- positions_to_keep[,1]

                                    alignment %>%
                                        filter(position %in% positions_to_keep) %>%
                                        ggplot(aes(x = position, y = name)) +
                                            geom_tile(aes(fill = state)) +
                                            scale_fill_viridis(discrete = TRUE, guide = "none") +
                                            scale_y_discrete(name = "") +
                                            scale_x_continuous(expand = c(0,0)) +
                                            theme_classic() +
                                            coord_cartesian(xlim = c(input$x_start, input$x_end)) +
                                            theme(
                                                axis.text.y = element_blank()
                                            ) -> filtered_alignment_plot

                                    stats_plot + alignment_plot + filtered_alignment_plot + patchwork::plot_layout(ncol = 1)

                            })

                        ## Process and draw tree from filtered_alignment

                            output$filtered_tree <- renderPlot({

                                alignment_stats %>%
                                    pivot_wider(names_from = "metric", values_from = "value") %>%
                                    filter(
                                        gap_percent < input$gap_threshold &
                                        conservation_percent > input$conservation_threshold
                                    ) %>% 
                                    select(position) %>% as.data.frame() -> positions_to_keep
                                    positions_to_keep <- positions_to_keep[,1]

                                alignment %>%
                                    filter(position %in% positions_to_keep) %>%
                                    pivot_wider(names_from = position, values_from = state) -> alignment_wide

                                if (type == "DNA") {
                                    trimmed_alignment <- DNAStringSet()
                                    for (i in 1:dim(alignment_wide)[1]) {
                                        trimmed_alignment <- c(trimmed_alignment, DNAStringSet(paste0(alignment_wide[i,2:(dim(alignment_wide)[2])], collapse = "")))
                                    }
                                    trimmed_alignment@ranges@NAMES <- as.character(as.data.frame(alignment_wide[,1])[,1])
                                    writeFasta(trimmed_alignment, paste0(alignment_in_path, "_trimmed"), type = "DNA")
                                }

                                if (type == "AA") {
                                    trimmed_alignment <- AAStringSet()
                                    for (i in 1:dim(alignment_wide)[1]) {
                                        trimmed_alignment <- c(trimmed_alignment, AAStringSet(paste0(alignment_wide[i,2:(dim(alignment_wide)[2])], collapse = "")))
                                    }
                                    trimmed_alignment@ranges@NAMES <- as.character(as.data.frame(alignment_wide[,1])[,1])
                                    writeFasta(trimmed_alignment, paste0(alignment_in_path, "_trimmed"), type = "AA")
                                }

                                trimmed_tree <- fortify(buildTree(
                                    scaffold_type = scaffold_type,
                                    scaffold_in_path = paste0(alignment_in_path, "_trimmed"),
                                ))

                                ggtree(trimmed_tree)

                            })

                    }

                ## Call the app

                    if ( jupyter == TRUE) {

                        available_ports <- vector()

                        if ( length( unique( c(
                                    randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50)
                        ))) == 2 ) { available_ports <- c(available_ports, 10123) }

                        if ( length( unique( c(
                                    randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50)
                        ))) == 2 ) { available_ports <- c(available_ports, 10456) }

                        if ( length( unique( c(
                                    randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                    randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50)
                        ))) == 2 ) { available_ports <- c(available_ports, 10789) }

                        if (length(available_ports) > 0) {
                            cat(paste0("Connect at: https://shiny", substr(available_ports[1],3,6), ".bustalab.d.umn.edu"))
                            runApp(shinyApp(ui = ui, server = server), host = "127.0.0.1", port = available_ports[1])
                        } else {
                            cat("All available ports are in use. Please try again later.")
                        }
                    
                    }

                    if (jupyter == FALSE) {
                        shinyApp(ui = ui, server = server)
                    }

            }

        #### readGFFs

            #' readGFFs
            #'
            #' Read multiple genome feature files into a tidy dataframe
            #' @param input_frame A dataframe with columns: Genus_species, GFF_in_path, region_reach (20000), concatenation_spacing (10000)
            #' @param query A .fa file containing the query
            #' @param phylochem The phylochem object to which the data should be bound
            #' @examples
            #' @export
            #' readGFFs

                readGFFs <- function(
                                input_frame, # Dataframe with columns
                                subset_mode = c("none", "goi_only", "goi_region", "name_check"), 
                                goi = NULL, # Character vector of gene names. Labelling by species is not necessary
                                concatenate_by_species = TRUE,
                                chromosomes_only = FALSE,
                                omit = NULL
                            ) {

                    framed_GFFs <- list()

                    for ( i in 1:dim(input_frame)[1] ) {

                        # Load GFFs
                            cat(paste("Loading GFF for ", input_frame$Genus_species[i], "...\n", sep = ""))
                            gff_as_GRange <- suppressWarnings(rtracklayer::import.gff(as.character(input_frame$GFF_in_path[i]), version = "3"))
                            chr <- unique(gff_as_GRange@seqnames@values)
                            cat(paste("Number of chromosomes and/or scaffolds: ", length(chr), "\n", sep = ""))

                        # Make the GFF into dataframes "chrom_scaff", and "ranges"
                            library(GenomicRanges)
                            chrom_scaff <- as.data.frame(GenomicRanges::seqinfo(gff_as_GRange))
                            chrom_scaff$chrom_scaff_name <- rownames(chrom_scaff)
                            chrom_scaff$chrom_scaff_number <- seq(1,dim(chrom_scaff)[1],1)
                            rownames(chrom_scaff) <- NULL

                            ranges <- as.data.frame(gff_as_GRange)
                            ranges <- ranges[,colnames(ranges) != "Parent"]

                        # Remove cols to omit
                            if ( length(omit) > 0 ) {
                                ranges <- ranges[, !colnames(ranges) %in% omit]    
                            }
                            
                        # Custom gene name editing on species-by-species basis
                            if ( input_frame$Genus_species[i] == "Zea_mays" ) {
                                # ranges$Name <- gsub("_.*$", "", ranges$Name)
                                ranges$Name <- ranges$locus_tag
                            }

                            if ( input_frame$Genus_species[i] == "Helianthus_annuus" ) {
                                ranges$Name <- ranges$ID
                            }

                            if ( input_frame$Genus_species[i] == "Solanum_lycopersicum" ) {
                                ranges$Name <- gsub("0\\.1", "0", ranges$Name)
                                ranges$Name <- gsub("0\\.2", "0", ranges$Name)
                                ranges$Name <- gsub("0\\.3", "0", ranges$Name)
                            }

                        # No subsetting
                            if ( subset_mode == "none" ) {
                                framed_GFFs <- ranges
                            }

                        # Name check mode
                            if ( subset_mode == "name_check" ) {
                                framed_GFFs[[i]] <- unique(data.frame(
                                                        Genus_species = input_frame$Genus_species[i],
                                                        names = ranges$Name
                                                    ))
                            }

                        # Subset "chrom_scaff" and "ranges"
                            if ( subset_mode == "goi_region" ) {

                                # Set up the genes_of_interest_frame for this species
                                    goi_sp <- goi[goi %in% ranges$Name]
                                    goi_sp_frame <- data.frame(goi_sp = goi_sp, chrom_scaff_of_interest = ranges$seqnames[match(goi_sp, ranges$Name)], to_subset = TRUE)

                                # Remove chrom_scaffs and ranges not on chrom_scaffs that contain goi
                                    chrom_scaff <- chrom_scaff[chrom_scaff$chrom_scaff_name %in% goi_sp_frame$chrom_scaff_of_interest,]
                                    ranges <- ranges[ranges$seqnames %in% goi_sp_frame$chrom_scaff_of_interest,]
                                    
                                # Remove ranges not within region_reach of a goi, if goi are within region_reach of eachother, expand the size of that region iteratively
                                    subsetted_ranges <- data.frame()
                                    goi_region_number <- 1
                                    while ( any(goi_sp_frame$to_subset) == TRUE ) {
                                        j <- which(goi_sp_frame$to_subset == TRUE)[1]
                                        
                                        # Define the start of the goi_start_sites and the goi inside it
                                            ranges_on_this_chrom_scaff <- ranges[ranges$seqnames == goi_sp_frame$chrom_scaff_of_interest[j],]
                                            goi_start_sites <- ranges_on_this_chrom_scaff[ranges_on_this_chrom_scaff$Name %in% goi_sp_frame$goi_sp[j],]$start
                                            known_goi_in_this_region <- vector()
                                            known_goi_in_this_region <- c(known_goi_in_this_region, as.character(goi_sp_frame$goi_sp[j]))
                                            region <- data.frame(start = min(goi_start_sites) - input_frame$region_reach[i], end = max(goi_start_sites) + input_frame$region_reach[i])
                                            region$start[region$start < 0] <- 0

                                        # Expand the range if it includes another goi
                                            finished <- FALSE
                                            while ( finished == FALSE) {
                                                goi_region_bounds <- vector()
                                                for (k in 1:dim(region)[1]) {
                                                    goi_region_bounds <- c(goi_region_bounds,seq(region[k,1], region[k,2], 1))
                                                }
                                                ranges_in_this_goi_region <- ranges_on_this_chrom_scaff[ranges_on_this_chrom_scaff$start %in% goi_region_bounds,]
                                                goi_in_current_goi_region <- goi_sp_frame$goi_sp[goi_sp_frame$goi_sp %in% ranges_in_this_goi_region$Name]
                                                if ( all(goi_in_current_goi_region %in% known_goi_in_this_region) ) {
                                                    finished <- TRUE
                                                } else {
                                                    temp <- ranges_on_this_chrom_scaff$Name %in% goi_in_current_goi_region[!goi_in_current_goi_region %in% known_goi_in_this_region]
                                                    addl_gene <- ranges_on_this_chrom_scaff[temp,]
                                                    addl_gene <- addl_gene[1,]
                                                    goi_start_sites <- c(goi_start_sites, addl_gene$start)
                                                    region <- data.frame(start = min(goi_start_sites) - input_frame$region_reach[i], end = max(goi_start_sites) + input_frame$region_reach[i])
                                                    region$start[region$start < 0] <- 0

                                                    # Mark that the addl_gene is now known in this region
                                                        known_goi_in_this_region <- c(known_goi_in_this_region, as.character(addl_gene$Name))
                                                }
                                            }
                                        goi_sp_frame$to_subset[goi_sp_frame$goi_sp %in% goi_in_current_goi_region] <- FALSE
                                        ranges_in_this_goi_region$goi_region_name <- paste(input_frame$Genus_species[i], goi_region_number, sep="_")
                                        ranges_in_this_goi_region$goi_region_number <- goi_region_number
                                        goi_region_number <- goi_region_number + 1
                                        subsetted_ranges <- rbind(subsetted_ranges, ranges_in_this_goi_region)
                                    }

                                # Reorient each goi_region so it begins at zero
                                    subsetted_ranges <- subsetted_ranges[subsetted_ranges$type != "region",] # remove "ranges" that are whole chrom_scaffs
                                    library(plyr)
                                    subsetted_ranges <- plyr::ddply(subsetted_ranges, .(goi_region_name), mutate, end = end - min(start), start = start - min(start))
                                    subsetted_ranges <- plyr::ddply(subsetted_ranges, .(goi_region_name), mutate, goi_region_length = max(end))

                                # Concatenate range_scaffs of this species
                                    if ( concatenate_by_species == TRUE ) {
                                        if ( length(unique(subsetted_ranges$goi_region_number)) > 1 ) {
                                            end_previous_goi_region <- max(subsetted_ranges[subsetted_ranges$goi_region_number == 1,]$end)
                                            for (this_goi_region_number in 2:length(unique(subsetted_ranges$goi_region_number))) {
                                                amount_to_advance_this_goi_region <- end_previous_goi_region + input_frame$concatenation_spacing[i]
                                                subsetted_ranges[subsetted_ranges$goi_region_number == this_goi_region_number,]$start <- subsetted_ranges[subsetted_ranges$goi_region_number == this_goi_region_number,]$start + amount_to_advance_this_goi_region
                                                subsetted_ranges[subsetted_ranges$goi_region_number == this_goi_region_number,]$end <- subsetted_ranges[subsetted_ranges$goi_region_number == this_goi_region_number,]$end + amount_to_advance_this_goi_region
                                                end_previous_goi_region <- max(subsetted_ranges[subsetted_ranges$goi_region_number == this_goi_region_number,]$end)
                                            }
                                        }
                                    }

                                # Bind to output
                                    subsetted_ranges$Genus_species <- input_frame$Genus_species[i]
                                    subsetted_ranges$Genus <- gsub("_.*$", "", input_frame$Genus_species[i])
                                    subsetted_ranges$species <- gsub(".*_", "", input_frame$Genus_species[i])
                                    subsetted_ranges$y_position <- input_frame$y_position[i]
                                    framed_GFFs[[i]] <- subsetted_ranges
                                    framed_GFFs <- do.call(rbind, framed_GFFs)
                            }
                    }

                    ## Concatenate range_scaffs for single species
                        if ( concatenate_by_species == TRUE & dim(input_frame)[1] == 1) {

                            ## Filter for only main chromosomes if specified, then number seqnames
                                if (chromosomes_only == TRUE) {
                                    chromosome_seqnames <- as.character(framed_GFFs$seqnames[framed_GFFs$chromosome %in% c(dropNA(as.numeric(unique(framed_GFFs$chromosome))), "X", "Y")])
                                    message("NA introduced by coersion here is OKAY!")
                                    framed_GFFs <- framed_GFFs[framed_GFFs$seqnames %in% chromosome_seqnames,]
                                    framed_GFFs$seqname_number <- 0
                                    framed_GFFs$seqname_number <- match(framed_GFFs$seqnames, chromosome_seqnames)
                                } else {
                                    chromosome_seqnames <- unique(framed_GFFs$seqnames)
                                    framed_GFFs$seqname_number <- 0
                                    framed_GFFs$seqname_number <- match(framed_GFFs$seqnames, chromosome_seqnames)
                                }
                            if ( length(unique(framed_GFFs$seqnames)) > 1 ) {
                                end_previous_seqname_number_region <- max(framed_GFFs[framed_GFFs$seqname_number == 1,]$end)
                                for (this_seqname_number in 2:length(unique(framed_GFFs$seqname_number))) {
                                    amount_to_advance_this_seqname_number_region <- end_previous_seqname_number_region + input_frame$concatenation_spacing[i]
                                    framed_GFFs[framed_GFFs$seqname_number == this_seqname_number,]$start <- framed_GFFs[framed_GFFs$seqname_number == this_seqname_number,]$start + amount_to_advance_this_seqname_number_region
                                    framed_GFFs[framed_GFFs$seqname_number == this_seqname_number,]$end <- framed_GFFs[framed_GFFs$seqname_number == this_seqname_number,]$end + amount_to_advance_this_seqname_number_region
                                    end_previous_seqname_number_region <- max(framed_GFFs[framed_GFFs$seqname_number == this_seqname_number,]$end)
                                }
                            }
                        }


                    ## Center the scaffolds, merging the GFFs if there is more than one
                        if (dim(input_frame)[1] > 1) {
                            framed_GFFs <- do.call(plyr::rbind.fill, framed_GFFs)
                            if ( subset_mode != "name_check" ) {
                                center <- (max(framed_GFFs$end) - min(framed_GFFs$start))/2
                                library(plyr)
                                framed_GFFs <- plyr::ddply(framed_GFFs, .(Genus_species), mutate, center_start = (start + (center-(max(end)-min(start))/2)), center_end = (end + (center-(max(end)-min(start))/2)) )
                            }
                        } else {
                            if ( subset_mode != "name_check" ) {
                                center <- (max(framed_GFFs$end) - min(framed_GFFs$start))/2
                                framed_GFFs$center_start = framed_GFFs$start + (center-(max(framed_GFFs$end)-min(framed_GFFs$start))/2)
                                framed_GFFs$center_end = framed_GFFs$end + (center-(max(framed_GFFs$end)-min(framed_GFFs$start))/2)
                            }
                        }

                    return(framed_GFFs)
                }

        #### fastxQC

            #' Run QC on fastx file(s)
            #'
            #' @param paths_to_fastxs A list of the paths to the fastx file(s) you want to analyze
            #' @param max_n_seqs The maximum number of sequences PER FILE to analyze. Default is no limit.
            #' @examples
            #' @export
            #' fastxQC

            fastxQC <-  function(
                            paths_to_fastxs,
                            type = c("fastq", "fasta"),
                            mode = c("fast", "slow"),
                            max_n_seqs = NULL,
                            progress = FALSE
                        ) {

                output <- list()
                for ( path_to_fastx in 1:length(paths_to_fastxs)) {

                    cat(paste("Analyzing file ", paths_to_fastxs[path_to_fastx], "\n", sep = ""))

                    if (length(grep("http", paths_to_fastxs[path_to_fastx])) > 0) {
                        downloadFasta(paths_to_fastxs[path_to_fastx])
                        paths_to_fastxs[path_to_fastx] <- "TEMP___fasta"
                    }

                    ## Determine max number of reads

                        if (length(max_n_seqs) > 0) {
                            n_reads <- max_n_seqs
                        } else {
                            if (type == "fasta") {
                                n_reads <- (dim(data.table::fread(paths_to_fastxs[path_to_fastx], sep = NULL, header = FALSE))[1]/2)
                            }
                            if (type == "fastq") {
                                n_reads <- (dim(data.table::fread(paths_to_fastxs[path_to_fastx], sep = NULL, header = FALSE))[1]/4)
                            }
                        }

                    ## Get length distribution. If fastq, convert to fasta first. If slow, do quality scores

                        if (type == "fasta") {

                            output[[length(output)+1]] <- cbind(
                                data.frame(
                                    file = rep(paths_to_fastxs[path_to_fastx], n_reads)
                                ),
                                read.table(Rsamtools::indexFa(
                                    paths_to_fastxs[path_to_fastx]
                                ))[1:n_reads,1:2]
                            )

                        }

                        if (type == "fastq" & mode == "fast") {
                            
                            cat(paste("Converting to fasta and getting sequence lengths...", "\n", sep = ""))

                            fastqToFasta(
                                file_in_path = paths_to_fastxs[path_to_fastx],
                                file_out_path = paste0(paths_to_fastxs[path_to_fastx], ".fasta")
                            )

                            output[[length(output)+1]] <- cbind(
                                data.frame(
                                    file = rep(paths_to_fastxs[path_to_fastx], n_reads)
                                ),
                                read.table(Rsamtools::indexFa(
                                    paste0(paths_to_fastxs[path_to_fastx], ".fasta")
                                ))[1:n_reads,1:2]
                            )

                        }

                        if (type == "fastq" & mode == "slow") {

                            cat(paste("Using Phred_ASCII_33 score conversions...", "\n", sep = ""))
                            cat(paste("Analyzing quality scores...", "\n", sep = ""))

                            if( progress ) { pb <- progress::progress_bar$new(total = n_reads) }
                            for (i in 1:n_reads) {

                                data <- data.table::fread(
                                    paths_to_fastxs[path_to_fastx], skip = (4*(i-1)),
                                    nrows = 4, sep = NULL, header = FALSE
                                )

                                # if (!length(strsplit(as.character(unlist(data[2,])), split = "")[[1]]) == length(strsplit(as.character(unlist(data[4,])), split = "")[[1]])) {
                                #     stop(paste0("Problem with quality score encoding on line ", i))
                                # }

                                read_score <- mean(as.numeric(phred33_lookup$X2[match(
                                                    strsplit(as.character(unlist(data[4,])), split = "")[[1]],
                                                    phred33_lookup$X1
                                                )]))

                                if (is.na(read_score)) {

                                    bad_characters <- strsplit(as.character(unlist(data[4,])), split = "")[[1]][which(is.na(match(
                                        strsplit(as.character(unlist(data[4,])), split = "")[[1]],
                                        phred33_lookup$X1
                                    )))]

                                    stop(paste0(
                                        "Problem with quality score encoding on line ", i, ".", "\n",
                                        "Unknown characters: ", bad_characters))

                                }
                                
                                output[[length(output)+1]] <- data.frame(
                                    file = paths_to_fastxs[path_to_fastx],
                                    name = data[1],
                                    length = nchar(data[2,]),
                                    score = read_score,
                                    accuracy = (1-10^(-read_score/10))*100
                                )

                                if( progress ) { pb$tick() }
                            
                            }
                        }
                }

                output <- do.call(rbind, output)
                rownames(output) <- NULL
                colnames(output) <- c("file", "name", "length", "score", "accuracy")[1:dim(output)[2]]

                return(as_tibble(output))
            }

        #### reverseTranslate

            #' Reverse translate a protein sequence
            #'
            #' @param fasta_in_path
            #' @param fasta_out_path
            #' @param organism The codon table to use. One of "Arabidopsis_thaliana", "Zea_mays", "Nicotiana_tabacum", or "Saccharomyces_cerevisiae".
            #' @examples
            #' @export
            #' reverseTranslate

            reverseTranslate <- function(
                
                fasta_in_path,
                fasta_out_path,
                organism = c("Arabidopsis_thaliana", "Zea_mays", "Nicotiana_tabacum", "Saccharomyces_cerevisiae")) {

                codon_tables <- readMonolist("https://thebustalab.github.io/R_For_Chemists/sample_data/codon_tables.csv")
                codon_tables$Genus_species <- paste(codon_tables$Genus, codon_tables$species, sep = "_")
                codon_tables <- dplyr::filter(codon_tables, Genus_species == organism)
                    
                substitution_table <- list()
                for (i in 1:length(unique(codon_tables$amino_acid))) {
                    this_codon <- dplyr::filter(codon_tables, amino_acid == unique(codon_tables$amino_acid)[i])
                    substitution_table[[i]] <- this_codon[which.max(this_codon$frequency),]
                }
                substitution_table <- do.call(rbind, substitution_table)
                rownames(substitution_table) <- NULL

                amin_seqs <- Biostrings::readAAStringSet(fasta_in_path)
                output <- Biostrings::DNAStringSet()
                for (i in 1:length(amin_seqs)) {
                    this_amin_seq <- strsplit(as.character(amin_seqs[i]), split = NULL)[[1]]

                    ## Add terminal stop codon if not present
                        if( this_amin_seq[length(this_amin_seq)] != "*" ) {
                            this_amin_seq <- c(this_amin_seq, "*")
                        }
                    
                    ## Determine codons and concatenate them
                        output[[i]] <- paste(substitution_table$codon[match(this_amin_seq, substitution_table$amino_acid)], collapse = "")

                }
                output@ranges@NAMES <- amin_seqs@ranges@NAMES
                Biostrings::writeXStringSet(output, fasta_out_path)

            }

        #### fastqToFasta

            #' Convert fastq to fasta
            #'
            #' @param file_in_path
            #' @param file_out_path
            #' @examples
            #' @export
            #' fastqToFasta

            fastqToFasta <- function(file_in_path, file_out_path) {
                data <- readLines(file_in_path)
                fasta <- DNAStringSet(data[seq(1,length(data),4)+1])
                names(fasta) <- data[seq(1,length(data),4)]
                writeXStringSet(fasta, file_out_path)
            }

        #### canuHistogramToKmerTable

            #' Convert a kmer report file from Canu into a table that can be accepted by genomescope
            #' Use http://qb.cshl.edu/genomescope/genomescope2.0/
            #'
            #' @param file_in_path
            #' @param file_out_path
            #' @examples
            #' @export
            #' canuHistogramToKmerTable

            canuHistogramToKmerTable <- function(file_in_path, file_out_path) {

                dat <- readLines(file_in_path)
                dat <- dat[11:length(dat)]
                dat <- gsub("1\\.0.*$", "", gsub("0\\..*$", "", dat))
                dat <- strsplit(paste(dat, collapse = " "), " ")[[1]]
                dat <- dat[dat != ""]

                dat <- data.frame(
                    x = as.numeric(dat[seq(1,length(dat), 2)]),
                    y = as.numeric(dat[seq(2,length(dat), 2)])
                )

                head(dat)
                tail(dat)

                write.table(
                    dat, file = file_out_path,
                    sep = " ", row.names = FALSE, col.names = FALSE
                )

            }

        #### genomeScope

            #' Run genomescope
            #' See: https://github.com/tbenavi1/genomescope2.0
            #'
            #' @param file_in_path
            #' @param file_out_path
            #' @examples
            #' @export
            #' genomescope

            genomeScope <- function( input_histogram_file, output_dir, ploidy, kmer_length ) {

                ## General set up stuff

                    # library('argparse') ## Load libraries for non-linear least squares and argument parser
                    # library('genomescope') ## Load the genomescope library

                    ## Define default arguments

                        arguments <- vector()
                        arguments$input <- input_histogram_file # parser$add_argument("-i", "--input", help = "input histogram file")
                        arguments$output <- output_dir # parser$add_argument("-o", "--output", help = "output directory name")
                        arguments$ploidy <- ploidy # parser$add_argument("-p", "--ploidy", type = "integer", default = 2, help = "ploidy (1, 2, 3, 4, 5, or 6) for model to use [default 2]")
                        arguments$kmer_length <- kmer_length # parser$add_argument("-k", "--kmer_length", type = "integer", default = 21, help = "kmer length used to calculate kmer spectra [default 21]")
                        arguments$version <- FALSE
                        arguments$name_prefix <- "test" # parser$add_argument("-n", "--name_prefix", default = "", help = "optional name_prefix for output files")
                        arguments$fitted_hist <- FALSE # parser$add_argument("--fitted_hist", action="store_true", default=FALSE, help = "ADVANCED: generates a fitted histogram for kmer multiplicity 0-4 and a lookup table of probabilities")
                        arguments$initial_heterozygosities <- -1 # parser$add_argument("--initial_heterozygosities", type="character", default = -1, help = "ADVANCED: flag to set initial values for nucleotide heterozygosity rates")
                        arguments$initial_repetitiveness <- -1 # parser$add_argument("--initial_repetitiveness", type="character", default = -1, help = "ADVANCED: flag to set initial value for repetitiveness")
                        arguments$lambda <- -1 # parser$add_argument("-l", "--lambda", "--kcov", "--kmercov", type = "integer", default=-1, help = "optional initial kmercov estimate for model to use")
                        arguments$max_kmercov <- -1 # parser$add_argument("-m", "--max_kmercov", type = "integer", default=-1, help = "optional maximum kmer coverage threshold (kmers with coverage greater than max_kmercov are ignored by the model)")
                        arguments$name_prefix <- ""
                        arguments$no_unique_sequence <- FALSE # parser$add_argument("--no_unique_sequence", action="store_true", default=FALSE, help = "optional flag to turn off yellow unique sequence line in plots")
                        arguments$num_rounds <- 4 # parser$add_argument("--num_rounds", type = "integer", default = 4, help = "ADVANCED: parameter for the number of optimization rounds")
                        arguments$start_shift <- 5 # parser$add_argument("--start_shift", type = "integer", default=START_SHIFT, help = "ADVANCED: coverage shifts to exclude between fitting rounds")
                        arguments$testing <- FALSE # parser$add_argument("--testing", action="store_true", default=FALSE, help = "ADVANCED: flag to create testing.tsv file with model parameters")
                        arguments$topology <- 0 # parser$add_argument("-t", "--topology", type = "integer", default = 0, help = "ADVANCED: flag for topology for model to use")
                        arguments$trace_flag <- FALSE # parser$add_argument("--trace_flag", action="store_true", default=FALSE, help = "ADVANCED: flag to turn on printing of iteration progress of nlsLM function")
                        arguments$transform_exp <- 1 # parser$add_argument("--transform_exp", type="integer", default=1, help = "ADVANCED: parameter for the exponent when fitting a transformed (x**transform_exp*y vs. x) kmer histogram [default 1]")
                        arguments$true_params <- -1 # parser$add_argument("--true_params", type="character", default = -1, help = "ADVANCED: flag to state true simulated parameters for testing mode")
                        arguments$typical_error <- 15 # parser$add_argument("--typical_error", type = "integer", default=TYPICAL_ERROR, help = "ADVANCED: typical level of sequencing error")
                        arguments$verbose <- FALSE # parser$add_argument("--verbose", action="store_true", default=FALSE, help = "optional flag to print messages during execution")
                        arguments$version <- FALSE # parser$add_argument("-v", "--version", action="store_true", default=FALSE, help="print the version and exit")
                        version_message <- "GenomeScope 2.0\n"

                        if (is.null(arguments$input) | is.null(arguments$output)) { stop("Please provide input and output") }

                    ## Transfer arguments to variable names
                        
                        histfile <- arguments$input
                        foldername  <- arguments$output
                        p           <- arguments$ploidy
                        k           <- arguments$kmer_length
                        estKmercov  <- arguments$lambda
                        max_kmercov <- arguments$max_kmercov
                        VERBOSE     <- arguments$verbose
                        NO_UNIQUE_SEQUENCE <- arguments$no_unique_sequence
                        topology    <- arguments$topology
                        d_init      <- arguments$initial_repetitiveness
                        r_inits     <- arguments$initial_heterozygosities
                        transform_exp <- arguments$transform_exp
                        TESTING     <- arguments$testing
                        TRUE_PARAMS <- arguments$true_params
                        TRACE_FLAG <- arguments$trace_flag
                        NUM_ROUNDS <- arguments$num_rounds
                        FITTED_HIST <- arguments$fitted_hist
                        START_SHIFT <- arguments$start_shift
                  
                    ## Define some other variables
                        
                        NUM_ROUNDS = 4 ## Number of rounds before giving up
                        START_SHIFT = 5 ## Coverage steps to trim off between rounds
                        TYPICAL_ERROR = 15 ## Typical cutoff for sequencing error
                        MAX_ITERATIONS=200 ## Max rounds on NLS
                        SCORE_CLOSE = 0.20 ## Overrule if two scores are within this percent (0.05 = 5%) but larger difference in het
                        SCORE_HET_FOLD_DIFFERENCE = 10 ## Overrule heterozygosity if there is a large difference in het rate
                    
                    ## Suppress the warnings if the modeling goes crazy, those are in try/catch blocks anyways
                        options(warn = -1)

                    ## Colors for plots
                        COLOR_BGCOLOR  = "light grey"
                        COLOR_HIST     = "#56B4E9"
                        COLOR_2pPEAK   = "black"
                        COLOR_pPEAK    = "#F0E442"
                        COLOR_ERRORS   = "#D55E00"
                        COLOR_KMERPEAK = "black"
                        COLOR_RESIDUAL = "purple"
                        COLOR_COVTHRES = "red"

                    ## Given mean +/- stderr, report min and max value within 2 SE

                        min_max <- function(table) {
                          return (c(max(0,table[1] - 2*table[2]), table[1]+ 2*table[2]))
                        }

                        min_max1 <- function(table) {
                          return (c(max(0,table[1] - 2*table[2]), min(1, table[1]+ 2*table[2])))
                        }

                ## Program

                    cat(paste("GenomeScope analyzing ", histfile, " p=", p, " k=", k, " outdir=", foldername, "\n", sep=""))

                    dir.create(foldername, showWarnings=FALSE)

                    ## Initialize the status
                    progressFilename <- paste(foldername,"/", arguments$name_prefix, "progress.txt",sep="")
                    cat("starting", file=progressFilename, sep="\n")

                    kmer_prof <- read.csv(file=histfile,sep="", header=FALSE,colClasses=c("numeric","numeric"))

                    minkmerx = 1;
                    if (kmer_prof[1,1] == 0) {
                    if (VERBOSE) {cat("Histogram starts with zero, reseting minkmerx\n")}
                    minkmerx = 2;
                    }

                    kmer_prof <- kmer_prof[c(minkmerx:(length(kmer_prof[,2])-1)),] #get rid of the last position
                    kmer_prof_orig <- kmer_prof

                    ## try to find the local minimum between errors and the first (heterozygous) peak
                        kmer_trans = as.numeric(kmer_prof[,1])**transform_exp*as.numeric(kmer_prof[,2])
                        start <- tail(which(kmer_trans[1:TYPICAL_ERROR]==min(kmer_trans[1:TYPICAL_ERROR])),n=1)
                        start_max <- start + which(kmer_trans[start:length(kmer_trans)]==max(kmer_trans[start:length(kmer_trans)])) - 1

                    maxCovIndex = -1

                    ## Figure out which kmers to exclude, if any
                        if(max_kmercov == -1) {
                            maxCovIndex <- length(kmer_prof[,1])
                            max_kmercov <- kmer_prof[maxCovIndex,1]
                        } else {
                    
                    ## Figure out the index we should use for this coverage length
                            x <- kmer_prof[,1]
                            maxCovIndex <- length(x[x<=max_kmercov])
                        }

                    if (VERBOSE) {cat(paste("using max_kmercov:", max_kmercov, " with index:", maxCovIndex, "\n"))}

                    # terminate after NUM_ROUND iterations, store best result so far in container
                    round <- 0
                    best_container <- list(NULL,0)

                ## The loop

                    while (round < NUM_ROUNDS) {
                        
                        cat(paste("round", round, "trimming to", start, "trying 2p peak model... "), file=progressFilename, sep="", append=TRUE)
                        if (VERBOSE) {
                            cat(paste("round", round, "trimming to", start, "trying 2p peak model... \n"))
                        }

                        ## Reset the input trimming off low frequency error kmers
                            kmer_prof=kmer_prof_orig[1:maxCovIndex,]
                            x <- kmer_prof[start:maxCovIndex,1]
                            y <- kmer_prof[start:maxCovIndex,2]

                                ## Below is a non-subfunction implementation of estimate_Genome_peakp, to simplify code

                                    kmer_hist_orig = kmer_prof

                                    if (topology==-1) {
                                        p_to_num_topologies = c(1, 1, 1, 2, 5, 16)
                                        num_topologies = p_to_num_topologies[p]
                                        topologies = 1:num_topologies
                                    } else {
                                        num_topologies = 1
                                        topologies = c(topology)
                                    }
                                    
                                    numofKmers = sum(as.numeric(x)*as.numeric(y))
                                    
                                    if (estKmercov==-1) {
                                        #In situations with low heterozygosity, the peak with highest amplitude typically corresponds to the homozygous peak (i.e. the p-th peak).
                                        #However, with increasing heterozygosity, the highest amplitude peak may be an earlier peak.
                                        #Thus, when setting the estimated kmer coverage, we will need to iterate through these possibilities.
                                        #num_peak_indices indicates how many possibilities we need to iterate through.
                                        num_peak_indices = p
                                        y_transform = as.numeric(x)**transform_exp*as.numeric(y)
                                        estKmercov1 = x[which(y_transform==max(y_transform))][1]
                                    } else {
                                        # When the user sets the estimated kmer coverage, we only need to iterate through one possibility
                                        num_peak_indices = 1
                                        ## We set the estimated kmer coverage to be the user specified value
                                        estKmercov1 = estKmercov
                                    }
                                    
                                    estLength1 = numofKmers/estKmercov1

                                    nls00 = NULL
                                    peak_indices = 1:num_peak_indices
                                    for (i in peak_indices) {
                                        nls0 = NULL
                                        top_count = 0
                                        ## We see what happens when we set the estimated kmer coverage to be 1/i times the x-coordinate where the max peak occurs (1 <= i <= p if the user doesn't set the estimated kmer coverage, and i=1 if they do)
                                        estKmercov2 = estKmercov1 / i
                                        estLength2 = numofKmers/estKmercov2

                                        if (VERBOSE) {cat(paste("trying with kmercov: ", estKmercov2, "\n"))}

                                        for (top in topologies) {
                                            if (VERBOSE) {cat(paste("trying with topology: ", top, "\n"))}
                                            top_count = top_count + 1
                                            

                                            ## Here is a non-subfunction implementation of nls_peak, to simplify code
    
                                                # Initiate variables
                                                    
                                                    model = NULL
                                                    best_deviance = Inf
                                                    d_min = 0
                                                    if (d_init!=-1) {
                                                      d_initial = d_init
                                                    } else {
                                                      d_initial = 0.10
                                                    }
                                                    d_max = 1
                                                    r_min = 0.00001
                                                    if (top==0) {
                                                      p_to_num_r = c(0, 1, 2, 4, 6, 10)
                                                    } else {
                                                      p_to_num_r = c(0, 1, 2, 3, 4, 5)
                                                    }
                                                    num_r = p_to_num_r[p]
                                                    r_max = 1
                                                    kmercov_min = 0
                                                    kmercov_initial = estKmercov2
                                                    kmercov_max = Inf
                                                    bias_min = 0
                                                    bias_initial = 0.5
                                                    bias_max = Inf
                                                    length_min = 0
                                                    length_initial = estLength2/p
                                                    length_max = Inf

                                                # Determine what formula to use, based on p
                                                    
                                                    if (p==1) {
                                                      r_text = ""
                                                    } else {
                                                      r_text = paste(paste(lapply(1:(num_r), function(x) paste("r", as.character(x), sep="")), collapse=", "), ", ")
                                                    }
                                                    
                                                    x = x[1:min(2000,length(x))]
                                                    y = y[1:min(2000,length(y))]
                                                    y_transform = as.numeric(x)**transform_exp*as.numeric(y)
                                                    formula = as.formula(paste("y_transform ~ x**transform_exp*length*predict",p,"_",top,"(",r_text, "k, d, kmercov, bias, x)",sep=""))

                                                    if (VERBOSE) {cat("trying nlsLM algorithm (Levenberg-Marquardt)\n")}

                                                    if (r_inits!=-1) {
                                                        r_initials = unlist(lapply(strsplit(r_inits,","),as.numeric))
                                                        if (length(r_initials)!=num_r) {
                                                            stop("Incorrect number of initial rates supplied.")
                                                        }
                                                        r_initials_list = list(r_initials)
                                                    } else {
                                                        r_initials_list = list(rep(0.001, num_r), 0.001*(1:num_r), 0.001*(num_r:1), rep(0.01, num_r), 0.01*(1:num_r), 0.01*(num_r:1))
                                                    }

                                                    for (r_initials in r_initials_list) {

                                                        model1 = NULL
                                                        r_start = vector("list", num_r)
                                                        if (p > 1) {
                                                            names(r_start) = paste("r", 1:(num_r), sep="")
                                                            for (i in 1:(num_r)) {
                                                                r_start[[paste("r",i,sep="")]] = r_initials[i]
                                                            }
                                                        }

                                                        try(model1 <- nlsLM(formula = formula,
                                                            start   = c(list(d = d_initial), r_start, list(kmercov = kmercov_initial, bias = bias_initial, length = length_initial)),
                                                            lower   = c(c(d_min), rep(r_min, num_r), c(kmercov_min, bias_min, length_min)),
                                                            upper   = c(c(d_max), rep(r_max, num_r), c(kmercov_max, bias_max, length_max)),
                                                            control = list(minFactor=1e-12, maxiter=MAX_ITERATIONS, factor=0.1), trace=TRACE_FLAG), silent = TRUE
                                                        )

                                                        if (!is.null(model1)) {
                                                            current_deviance = model1$m$deviance()
                                                            #cat("Model deviance: ", current_deviance, "\n")
                                                            if (current_deviance < best_deviance) {
                                                              model = model1
                                                              best_deviance = current_deviance
                                                            }
                                                        } else {
                                                            #print("Model did not converge.")
                                                        }

                                                    }

                                                    if (!is.null(model)) {
                                                        model_sum    = summary(model)
                                                        model$p      = p
                                                        model$top = top
                                                        
                                                        if (p==1) {
                                                            model$hets = list(c(0, 0))
                                                        } else {
                                                            model$hets = lapply(1:(num_r), function(x) min_max1(model_sum$coefficients[paste('r', x, sep=""),]))
                                                        }

                                                        #model$het = c(1-Reduce("*", 1-unlist(lapply(model$hets, '[[', 1))), 1-Reduce("*", 1-unlist(lapply(model$hets, '[[', 2))))
                                                        model$het = c(sum(sapply(model$hets, '[[', 1)), sum(sapply(model$hets, '[[', 2)))
                                                        model$homo = 1-model$het
                                                        model$dups   = min_max(model_sum$coefficients['bias',])
                                                        model$kcov   = min_max(model_sum$coefficients['kmercov',])
                                                        model$mlen   = min_max(model_sum$coefficients['length',])
                                                        model$md     = min_max1(model_sum$coefficients['d',])
                                                        
                                                        if (p==1) {
                                                            model$ahets = list(c(0))
                                                        } else {
                                                            model$ahets = lapply(1:(num_r), function(x) model_sum$coefficients[paste('r', x, sep=""),][[1]])
                                                        }
                                                        
                                                        #model$ahet = 1-Reduce("*", 1-unlist(model$ahets))
                                                        model$ahet = Reduce("+", model$ahets)
                                                        model$ahomo = 1-model$ahet
                                                        model$adups = model_sum$coefficients['bias',][[1]]
                                                        model$akcov = model_sum$coefficients['kmercov',][[1]]
                                                        model$amlen = model_sum$coefficients['length',][[1]]
                                                        model$amd   = model_sum$coefficients['d',][[1]]
                                                    }

                                            nls1 = model ## here, model is the output of the non-subfunction version of nls_peak
                                            nls0 = eval_model(kmer_hist_orig, nls0, nls1, p, round, foldername, arguments, score_close = SCORE_CLOSE)[[1]]
                                        }
                                        if (i < num_peak_indices) { #if this is not the last evaluation
                                            nls00 = eval_model(kmer_hist_orig, nls00, nls0, p, round, foldername, arguments, score_close = SCORE_CLOSE)[[1]]
                                        }
                                    }

                                    model_peaks <- eval_model(kmer_hist_orig, nls00, nls0, p, round, foldername, arguments, score_close = SCORE_CLOSE)


                        if (!is.null(model_peaks[[1]])) {
                            cat(paste("converged. score: ", model_peaks[[2]]$all[[1]]), file=progressFilename, sep="\n", append=TRUE)

                            if (VERBOSE) {
                                # mdir = paste(foldername, "/round", round, sep="")
                                # dir.create(mdir, showWarnings = FALSE)
                                # report_results(kmer_prof,kmer_prof_orig, k, p, model_peaks, mdir, arguments, TRUE)
                            }
                        } else { cat(paste("unconverged"), file = progressFilename, sep="\n", append=TRUE) }

                        #check if this result is better than previous
                        if (!is.null(model_peaks[[1]])) {
                            if (is.null(best_container[[1]])) {
                                if (VERBOSE) {cat("no previous best, updating best\n")}
                                best_container = model_peaks
                            }
                            else {
                                best_container_score = best_container[[1]]$m$deviance()
                                model_peaks_score = model_peaks[[1]]$m$deviance()
                                pdiff = abs(model_peaks_score - best_container_score) / max(model_peaks_score, best_container_score)

                                if (pdiff < SCORE_CLOSE) {
                                  hetm = model_peaks[[1]]$ahet
                                  hetb = best_container[[1]]$ahet

                                  #if (hetb * SCORE_HET_FOLD_DIFFERENCE < hetm) {
                                  if (hetb + 0.01 < hetm) {
                                    if (VERBOSE) {cat("model has significantly higher heterozygosity but similar score, overruling\n")}
                                  }
                                  #else if (hetm * SCORE_HET_FOLD_DIFFERENCE < hetb) {
                                  else if (hetm + 0.01 < hetb) {
                                    if (VERBOSE) {cat("previous best has significantly higher heterozygosity and similar score, keeping\n")}
                                    best_container = model_peaks
                                  }
                                  else if (model_peaks_score < best_container_score) {
                                    if (VERBOSE) {cat("score is marginally better but het rate is not extremely different, updating\n")}
                                    best_container = model_peaks
                                  }
                                }
                                else if (model_peaks_score < best_container_score) {
                                      if (VERBOSE) {cat("score is significantly better, updating\n")}
                                      best_container = model_peaks
                                }
                            }
                        }

                        ## Ignore a larger number of kmers as errors
                            start <- start + START_SHIFT
                            round <- round + 1
                    }

                ## Report the results, note using the original full profile
                    ## This is a non-subfunction implementation of report_results
                    
                    ## arguments to report_results
                        kmer_hist = kmer_prof
                        container = best_container
                        IN_VERBOSE = FALSE
  
                    output <- list()

                    ## Some initial stuff

                        kmer_hist = kmer_prof
                        kmer_hist_orig = kmer_prof_orig
                        k = k
                        p = p
                        container = best_container
                        foldername = foldername
                        arguments = arguments
                        IN_VERBOSE = FALSE

                        x=kmer_hist_orig[[1]]
                        y_orig=kmer_hist_orig[[2]]
                        y = as.numeric(x)**transform_exp*as.numeric(y_orig)
                        kmer_hist_transform = kmer_hist_orig
                        kmer_hist_transform$V2 = as.numeric(kmer_hist_transform$V1)**transform_exp * as.numeric(kmer_hist_transform$V2)
                        model = container[[1]]

                    # automatically zoom into the relevant regions of the plot, ignore first 15 positions
                    ## I think this is just for plotting and can be turned off
                        
                        xmax = length(x)
                        start_orig = which(y_orig == min(y_orig[1:TYPICAL_ERROR]))
                        start = which(y == min(y[1:TYPICAL_ERROR]))
                        zoomx = x[start:(xmax-1)]
                        zoomy_orig = y_orig[start_orig:(xmax-1)]
                        zoomy = y[start:(xmax-1)]
            
                    ## allow for a little space above max value past the noise
                    ## I think this is just for plotting and can be turned off
                        
                        y_limit_orig = max(zoomy_orig[start_orig:length(zoomy_orig)])*1.1
                        y_limit = max(zoomy[start:length(zoomy)])*1.1

                        x_limit_orig = which(y_orig == max(y_orig[start_orig:length(zoomx)])) * 3
                        x_limit = which(y == max(y[start:length(zoomx)])) * 3

                        if (min(zoomy_orig) > zoomy_orig[1]){
                          x_limit_orig=max(which(zoomy_orig<zoomy_orig[1])[2],600)
                        }
                        if (min(zoomy) > zoomy[1]){
                          x_limit=max(which(zoomy<zoomy[1])[2],600)
                        }

                        if (!is.null(model))
                        {
                          model_sum=summary(model)
                          kcov = min_max(model_sum$coefficients['kmercov',])[1]
                          x_limit_orig = max(kcov*(2*p+1.1), x_limit_orig)
                          x_limit = max(kcov*(2*p+1.1), x_limit)
                          if (model$top==0) {
                            p_to_num_r = c(0, 1, 2, 4, 6, 10)
                          } else {
                            p_to_num_r = c(0, 1, 2, 3, 4, 5)
                          }
                        } else {
                          if (topology==0) {
                            p_to_num_r = c(0, 1, 2, 4, 6, 10)
                          } else {
                            p_to_num_r = c(0, 1, 2, 3, 4, 5)
                          }
                        }
            
                    ## Uncomment this to enforce a specific number
                    ## I think this is just for plotting and can be turned off
                    # x_limit=150
            
                    ## Define some initial values

                        het=c(-1,-1)
                        homo=c(-1,-1)
                        num_r = p_to_num_r[p]
                        if (p > 1) {
                          hets = lapply(1:(num_r), function(x) c(-1, -1))
                          ahets = lapply(1:(num_r), function(x) -1)
                        }
                        amd = -1
                        akcov = -1
                        adups = -1
                        amlen = -1
                        atotal_len = -1
                        top = -1
                        total_len=c(-1,-1)
                        repeat_len=c(-1,-1)
                        unique_len=c(-1,-1)
                        dups=c(-1,-1)
                        error_rate=c(-1,-1)
                        model_status="fail"

                        model_fit_unique      = c(0,0,0)
                        model_fit_full        = c(0,0,0)
                        model_fit_all         = c(0,0,0)
                        model_fit_allscore    = c(0,0,0)
                        model_fit_fullscore   = c(0,0,0)
                        model_fit_uniquescore = c(0,0,0)

                        plot_size=2000
                        font_size=1.2
                        resolution=300

                    ## Plot the distribution, and hopefully with the model fit (currently turned off)
                        
                        # ylabel_orig = "Frequency"
                        # if (transform_exp == 1) {
                        # ylabel_transform = "Coverage*Frequency"
                        # } else {
                        # ylabel_transform = paste("Coverage^", transform_exp, "*Frequency", sep="")
                        # }
                        # png(paste(foldername, "/", arguments$name_prefix, "linear_plot.png", sep=""),
                        # width=plot_size, height=plot_size, res=resolution)
                        # par(mar = c(5.1,4.1,6.1,2.1))
                        # plot(kmer_hist_orig, type="n", main="GenomeScope Profile\n\n\n",
                        # xlab="Coverage", ylab=ylabel_orig, ylim=c(0,y_limit_orig), xlim=c(0,x_limit_orig),
                        # cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
                        # #rect(0, 0, max(kmer_hist_orig[[1]])*1.1 , max(kmer_hist_orig[[2]])*1.1, col=COLOR_BGCOLOR)
                        # rect(0, 0, x_limit_orig*1.1 , y_limit_orig*1.1, col=COLOR_BGCOLOR)
                        # points(kmer_hist_orig, type="h", col=COLOR_HIST, lwd=2)
                        # #  if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
                        # #    abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
                        # #  }
                        # box(col="black")

                        # png(paste(foldername, "/", arguments$name_prefix, "transformed_linear_plot.png", sep=""),
                        # width=plot_size, height=plot_size, res=resolution)
                        # par(mar = c(5.1,4.1,6.1,2.1))
                        # plot(kmer_hist_transform, type="n", main="GenomeScope Profile\n\n\n",
                        # xlab="Coverage", ylab=ylabel_transform, ylim=c(0,y_limit), xlim=c(0,x_limit),
                        # cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
                        # #rect(0, 0, max(kmer_hist_orig[[1]])*1.1 , max(kmer_hist_orig[[2]])*1.1, col=COLOR_BGCOLOR)
                        # rect(0, 0, x_limit*1.1 , y_limit*1.1, col=COLOR_BGCOLOR)
                        # points(kmer_hist_transform, type="h", col=COLOR_HIST, lwd=2)
                        # #  if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
                        # #    abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
                        # #  }
                        # box(col="black")

                        # ## Make a second plot in log space over entire range
                        # png(paste(foldername, "/", arguments$name_prefix, "log_plot.png", sep=""),
                        # width=plot_size, height=plot_size, res=resolution)
                        # par(mar = c(5.1,4.1,6.1,2.1))
                        # plot(kmer_hist_orig, type="n", main="GenomeScope Profile\n\n\n",
                        # xlab="Coverage", ylab=ylabel_orig, log="xy",
                        # cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
                        # rect(1e-10, 1e-10, max(kmer_hist_orig[,1])*10 , max(kmer_hist_orig[,2])*10, col=COLOR_BGCOLOR)
                        # points(kmer_hist_orig, type="h", col=COLOR_HIST, lwd=2)
                        # if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
                        # abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
                        # }
                        # box(col="black")

                        # png(paste(foldername, "/", arguments$name_prefix, "transformed_log_plot.png", sep=""),
                        # width=plot_size, height=plot_size, res=resolution)
                        # par(mar = c(5.1,4.1,6.1,2.1))
                        # plot(kmer_hist_transform, type="n", main="GenomeScope Profile\n\n\n",
                        # xlab="Coverage", ylab=ylabel_transform, log="xy",
                        # cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
                        # rect(1e-10, 1e-10, max(kmer_hist_transform[,1])*10 , max(kmer_hist_transform[,2])*10, col=COLOR_BGCOLOR)
                        # points(kmer_hist_transform, type="h", col=COLOR_HIST, lwd=2)
                        # if(length(kmer_hist[,1])!=length(kmer_hist_transform[,1])){
                        # abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
                        # }
                        # box(col="black")

                    ## Now the modeling
                            
                            if(!is.null(model)) {
                                
                                x=kmer_hist[[1]]
                                y=kmer_hist[[2]]
                                y_transform = as.numeric(x)**transform_exp*as.numeric(y)

                                ## The model converged!
                                    pred=predict(model, newdata=data.frame(x))

                                ## Compute the genome characteristics
                                    model_sum=summary(model)
                                    #print(model_sum)

                                ## save the model to a file
                                    capture.output(model_sum, file=paste(foldername,"/", arguments$name_prefix, "model.txt", sep=""))

                                ## Identify key values
                                    top   = model$top
                                    hets  = model$hets
                                    het   = model$het
                                    ahets = model$ahets
                                    ahet  = model$ahet
                                    homo  = model$homo
                                    ahomo = model$ahomo

                                    dups = model$dups
                                    kcov = model$kcov
                                    mlen = model$mlen
                                    md   = model$md

                                    adups = model$adups
                                    akcov = model$akcov
                                    amlen = model$amlen
                                    amd   = model$amd

                                ## Compute error rate, by counting kmers unexplained by model through first peak
                                ## truncate errors as soon as it goes to zero, dont allow it to go back up
                                    error_xcutoff = max(1, floor(kcov[1]))
                                    error_xcutoff_ind = tail(which(x<=error_xcutoff),n=1)
                                    if (length(error_xcutoff_ind)==0) {error_xcutoff_ind=1}

                                    error_kmers = x[1:error_xcutoff_ind]**(-transform_exp)*(y_transform[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind])

                                    first_zero = -1

                                    for (i in 1:error_xcutoff_ind) {
                                      if (first_zero == -1)
                                      {
                                        if (error_kmers[i] < 1.0)
                                        {
                                          first_zero = i
                                          if (VERBOSE) {cat(paste("Truncating errors at", i, "\n"))}
                                        }
                                      }
                                      else
                                      {
                                        error_kmers[i] = 0
                                      }
                                    }

                                    if (first_zero == -1) {
                                      first_zero = error_xcutoff_ind
                                      if (VERBOSE) {cat(paste("Truncating errors at", error_xcutoff_ind, "\n"))}
                                    }

                                ## Rather than "0", set to be some very small number so log-log plot looks okay
                                    error_kmers = pmax(error_kmers, 1e-10)

                                    total_error_kmers = sum(as.numeric(error_kmers) * as.numeric(x[1:error_xcutoff_ind]))

                                    total_kmers = sum(as.numeric(x)*as.numeric(y))

                                    error_rate = 1-(1-(total_error_kmers/total_kmers))**(1/k)
                                    error_rate = c(error_rate, error_rate)

                                    total_len = (total_kmers-total_error_kmers)/(p*kcov)
                                    atotal_len = (total_kmers-total_error_kmers)/(p*akcov)

                                ## find kmers that fit the p peak model (no repeats)
                                if (p==1) {
                                  unique_hist = amlen*predict1_1_unique(k, amd, akcov, adups, x)
                                }
                                if (p==2) {
                                  unique_hist = amlen*predict2_1_unique(ahets[[1]], k, amd, akcov, adups, x)
                                }
                                if (p==3) {
                                  unique_hist = amlen*predict3_1_unique(ahets[[1]], ahets[[2]], k, amd, akcov, adups, x)
                                }
                                if (p==4) {
                                  if (top==0) {
                                    unique_hist = amlen*predict4_0_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], k, amd, akcov, adups, x)
                                  } else {
                                    unique_hist = eval(parse(text = paste("amlen*predict4_", top, "_unique(ahets[[1]], ahets[[2]], ahets[[3]], k, amd, akcov, adups, x)", sep="")))
                                  }
                                }
                                if (p==5) {
                                  if (top==0) {
                                    unique_hist = amlen*predict5_0_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], k, amd, akcov, adups, x)
                                  } else {
                                    unique_hist = eval(parse(text = paste("amlen*predict5_", top, "_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], k, amd, akcov, adups, x)", sep="")))
                                  }
                                }
                                if (p==6) {
                                  if (top==0) {
                                    unique_hist = amlen*predict6_0_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], ahets[[7]], ahets[[8]], ahets[[9]], ahets[[10]], k, amd, akcov, adups, x)
                                  } else {
                                    unique_hist = eval(parse(text = paste("amlen*predict6_", top, "_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], k, amd, akcov, adups, x)", sep="")))
                                  }
                                }

                                r0 = 1-ahet #aa
                                t0 = r0**k #AA
                                s0 = t0 #AA
                                s1 = 1-t0 #AB
                                alpha_1 = (1-amd)*(2*s1) + amd*(2*s0*s1 + 2*s1**2)
                                alpha_2 = (1-amd)*(s0) + amd*(s1**2)
                                alpha_3 = amd*(2*s0*s1)
                                alpha_4 = amd*(s0**2)

                                one_hist = alpha_1 * dnbinom(x, size = akcov*1 / adups, mu = akcov*1)
                                two_hist = alpha_2 * dnbinom(x, size = akcov*p / adups, mu = akcov*p)
                                thr_hist = alpha_3 * dnbinom(x, size = akcov*3 / adups, mu = akcov*3)
                                fou_hist = alpha_4 * dnbinom(x, size = akcov*2*p / adups, mu = akcov*2*p)

                                unique_hist_transform = x**transform_exp*unique_hist

                                unique_kmers = sum(as.numeric(x)*as.numeric(unique_hist))
                                repeat_kmers = max(0, total_kmers - unique_kmers - total_error_kmers)

                                repeat_len=repeat_kmers/(p*kcov)
                                if (repeat_kmers == 0) {
                                  unique_len = total_len
                                } else {
                                  unique_len=unique_kmers/(p*kcov)
                                }

                                score = container[[2]]

                                model_fit_allscore    = score$allscore
                                model_fit_fullscore   = score$fullscore
                                model_fit_uniquescore = score$uniquescore

                                model_fit_all    = score$all
                                model_fit_full   = score$full
                                model_fit_unique = score$unique

                                residual_transform = y_transform - pred
                                residual = x**(-transform_exp)*residual_transform

                                hetline_simple = paste0("heterozygosity: ", format(100*ahet, digits=3), "%")

                                if (p==1) {
                                  hetline = paste0("a:", format(100*ahomo, digits=3), "%")
                                }
                                if (p==2) {
                                  hetline = paste0("aa:", format(100*ahomo,      digits=3), "% ",
                                                   "ab:", format(100*ahets[[1]], digits=3), "%")
                                }
                                if (p==3) {
                                  hetline = paste0("aaa:", format(100*ahomo,      digits=3), "% ",
                                                   "aab:", format(100*ahets[[1]], digits=3), "% ",
                                                   "abc:", format(100*ahets[[2]], digits=3), "%")
                                }
                                if (p==4) {
                                  if (top==0) {
                                    hetline = paste0("aaaa:", format(100*ahomo,      digits=3), "% ",
                                                     "aaab:", format(100*ahets[[1]], digits=3), "% ",
                                                     "aabb:", format(100*ahets[[2]], digits=3), "% ",
                                                     "aabc:", format(100*ahets[[3]], digits=3), "% ",
                                                     "abcd:", format(100*ahets[[4]], digits=3), "%")
                                  } else {
                                    hetline = paste0("aaaa:",                       format(100*ahomo,      digits=3), "% ",
                                                     switch(top, "aaab:", "aabb:"), format(100*ahets[[1]], digits=3), "% ",
                                                     "aabc:",                       format(100*ahets[[2]], digits=3), "% ",
                                                     "abcd:",                       format(100*ahets[[3]], digits=3), "%")
                                  }
                                }
                                if (p==5) {
                                  if (top==0) {
                                    hetline = paste0("aaaaa:", format(100*ahomo,      digits=3), "% ",
                                                     "aaaab:", format(100*ahets[[1]], digits=3), "% ",
                                                     "aaabb:", format(100*ahets[[2]], digits=3), "% ",
                                                     "aaabc:", format(100*ahets[[3]], digits=3), "% ",'\n',
                                                     "aabbc:", format(100*ahets[[4]], digits=3), "% ",
                                                     "aabcd:", format(100*ahets[[5]], digits=3), "% ",
                                                     "abcde:", format(100*ahets[[6]], digits=3), "%")
                                  } else {
                                    hetline = paste0("aaaaa:",                                                      format(100*ahomo,      digits=3), "% ",
                                                     switch(top, "aaaab:", "aaaab:", "aaabb:", "aaabb:", "aaabb:"), format(100*ahets[[1]], digits=3), "% ",
                                                     switch(top, "aaabc:", "aabbc:", "aaabc:", "aabcc:", "aabcc:"), format(100*ahets[[2]], digits=3), "% ",'\n',
                                                     switch(top, "aabcd:", "aabcd:", "aabcd:", "aabcd:", "abcdd:"), format(100*ahets[[3]], digits=3), "% ",
                                                     "abcde:",                                                      format(100*ahets[[4]], digits=3), "%")
                                  }
                                }
                                if (p==6) {
                                  if (top==0) {
                                    hetline = paste0("aaaaaa:", format(100*ahomo, digits=3), "% ",
                                                     "aaaaab:", format(100*ahets[[1]], digits=3), "% ",
                                                     "aaaabb:", format(100*ahets[[2]], digits=3), "% ",
                                                     "aaabbb:", format(100*ahets[[3]], digits=3), "% ",'\n',
                                                     "aaaabc:", format(100*ahets[[4]], digits=3), "% ",
                                                     "aaabbc:", format(100*ahets[[5]], digits=3), "% ",
                                                     "aabbcc:", format(100*ahets[[6]], digits=3), "% ",
                                                     "aaabcd:", format(100*ahets[[7]], digits=3), "% ",'\n',
                                                     "aabbcd:", format(100*ahets[[8]], digits=3), "% ",
                                                     "aabcde:", format(100*ahets[[9]], digits=3), "% ",
                                                     "abcdef:", format(100*ahets[[10]], digits=3), "%")
                                  } else {
                                    hetline = paste0("aaaaaa:", format(100*ahomo, digits=3), "% ",
                                                     switch(top, "aaaaab:", "aaaaab:", "aaaaab:", "aaaaab:", "aaaaab:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaabbb:", "aaabbb:", "aaabbb:"), format(100*ahets[[1]], digits=3), "% ",
                                                     switch(top, "aaaabc:", "aaaabc:", "aaabbc:", "aaabbc:", "aaabbc:", "aaaabc:", "aaaabc:", "aaabcc:", "aaabcc:", "aaabcc:", "aabbcc:", "aabbcc:", "aabbcc:", "aaabbc:", "aaabbc:", "aaabbc:"), format(100*ahets[[2]], digits=3), "% ",'\n',
                                                     switch(top, "aaabcd:", "aabbcd:", "aaabcd:", "aabccd:", "aabccd:", "aaabcd:", "aabbcd:", "aaabcd:", "aabcdd:", "aabcdd:", "aabbcd:", "aabcdd:", "aabcdd:", "aaabcd:", "aabccd:", "aabccd:"), format(100*ahets[[3]], digits=3), "% ",
                                                     switch(top, "aabcde:", "aabcde:", "aabcde:", "aabcde:", "abcdde:", "aabcde:", "aabcde:", "aabcde:", "aabcde:", "abcdee:", "aabcde:", "aabcde:", "abcdee:", "aabcde:", "aabcde:", "abcdde:"), format(100*ahets[[4]], digits=3), "% ",
                                                     "abcdef:", format(100*ahets[[5]], digits=3), "%")
                                  }
                                }

                                if (p >= 5) {
                                  hetline = hetline_simple
                                }

                                if (!IN_VERBOSE) {
                                  cat(paste0(hetline,"\n"))
                                }

                                # dev.set(dev.next())

                                ### THIS IS ALL NOW TURNED OFF

                            # ## finish plots (off)

                            #   ## Finish Linear Plot
                            #     title(paste("\n\nlen:",  prettyNum(total_len[1], big.mark=","),
                            #     "bp",
                            #     " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
                            #     "% ", "\n",
                            #     hetline, "\n",
                            #     " kcov:", format(akcov, digits=3),
                            #     " err:",   format(100*error_rate[1], digits=3),
                            #     "% ",
                            #     " dup:",  format(adups, digits=3),
                            #     " ",
                            #     " k:",   format(k, digits=3),
                            #     " p:",   format(p, digits=3),
                            #     sep=""),
                            #     cex.main=.85)

                            #     ## Mark the modes of the peaks
                            #     abline(v=akcov * (1:(2*p)), col=COLOR_KMERPEAK, lty=2)

                            #     ## Draw just the unique portion of the model
                            #     if (!NO_UNIQUE_SEQUENCE) {
                            #       lines(x, unique_hist, col=COLOR_pPEAK, lty=1, lwd=3)
                            #     }
                            #     lines(x, x**(-transform_exp)*pred, col=COLOR_2pPEAK, lwd=3)
                            #     lines(x[1:error_xcutoff_ind], error_kmers, lwd=3, col=COLOR_ERRORS)

                            #     if (VERBOSE) {
                            #       lines(x, residual, col=COLOR_RESIDUAL, lwd=3)
                            #     }

                            #     ## Add legend
                            #     if (NO_UNIQUE_SEQUENCE) {
                            #       legend(.62 * x_limit_orig, 1.0 * y_limit_orig,
                            #       legend=c("observed", "full model", "errors", "kmer-peaks"),
                            #       lty=c("solid", "solid", "solid", "dashed"),
                            #       lwd=c(3,3,3,2),
                            #       col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #       bg="white")
                            #     } else {
                            #       legend(.62 * x_limit_orig, 1.0 * y_limit_orig,
                            #       legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
                            #       lty=c("solid", "solid", "solid", "solid", "dashed"),
                            #       lwd=c(3,3,3,3,2),
                            #       col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #       bg="white")
                            #     }

                            #     dev.set(dev.next())

                            #   ## Finish Linear Plot
                            #     title(paste("\n\nlen:",  prettyNum(total_len[1], big.mark=","),
                            #     "bp",
                            #     " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
                            #     "% ", "\n",
                            #     hetline, "\n",
                            #     " kcov:", format(akcov, digits=3),
                            #     " err:",   format(100*error_rate[1], digits=3),
                            #     "% ",
                            #     " dup:",  format(adups, digits=3),
                            #     " ",
                            #     " k:",   format(k, digits=3),
                            #     " p:",   format(p, digits=3),
                            #     sep=""),
                            #     cex.main=.85)

                            #     ## Mark the modes of the peaks
                            #     abline(v=akcov * (1:(2*p)), col=COLOR_KMERPEAK, lty=2)

                            #     ## Draw just the unique portion of the model
                            #     if (!NO_UNIQUE_SEQUENCE) {
                            #       lines(x, unique_hist_transform, col=COLOR_pPEAK, lty=1, lwd=3)
                            #     }
                            #     lines(x, pred, col=COLOR_2pPEAK, lwd=3)
                            #     lines(x[1:error_xcutoff_ind], (x[1:error_xcutoff_ind]**transform_exp)*error_kmers, lwd=3, col=COLOR_ERRORS)

                            #     if (VERBOSE) {
                            #       lines(x, residual_transform, col=COLOR_RESIDUAL, lwd=3)
                            #     }

                            #     ## Add legend
                            #     if (NO_UNIQUE_SEQUENCE) {
                            #       legend(.62 * x_limit, 1.0 * y_limit,
                            #       legend=c("observed", "full model", "errors", "kmer-peaks"),
                            #       lty=c("solid", "solid", "solid", "dashed"),
                            #       lwd=c(3,3,3,2),
                            #       col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #       bg="white")
                            #     } else {
                            #       legend(.62 * x_limit, 1.0 * y_limit,
                            #       legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
                            #       lty=c("solid", "solid", "solid", "solid", "dashed"),
                            #       lwd=c(3,3,3,3,2),
                            #       col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #       bg="white")
                            #     }

                            #     dev.set(dev.next())

                            #   ## Finish Log plot
                            #     title(paste("\n\nlen:",  prettyNum(total_len[1], big.mark=","),
                            #     "bp",
                            #     " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
                            #     "% ", "\n",
                            #     hetline, "\n",
                            #     " kcov:", format(akcov, digits=3),
                            #     " err:",   format(100*error_rate[1], digits=3),
                            #     "% ",
                            #     " dup:",  format(adups, digits=3),
                            #     " ",
                            #     " k:",   format(k, digits=3),
                            #     " p:",   format(p, digits=3),
                            #     sep=""),
                            #     cex.main=.85)

                            #     ## Mark the modes of the peaks
                            #     abline(v=akcov * (1:(2*p)), col=COLOR_KMERPEAK, lty=2)

                            #     ## Draw just the unique portion of the model
                            #     if (!NO_UNIQUE_SEQUENCE) {
                            #       lines(x, unique_hist, col=COLOR_pPEAK, lty=1, lwd=3)
                            #     }
                            #     lines(x, x**(-transform_exp)*pred, col=COLOR_2pPEAK, lwd=3)
                            #     lines(x[1:error_xcutoff_ind], error_kmers, lwd=3, col=COLOR_ERRORS)

                            #     if (VERBOSE) {
                            #       lines(x, residual, col=COLOR_RESIDUAL, lwd=3)
                            #     }

                            #     ## Add legend
                            #     if(length(kmer_hist[,1])==length(kmer_hist_orig[,1]))
                            #     {
                            #       if (NO_UNIQUE_SEQUENCE) {
                            #         legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "errors", "kmer-peaks"),
                            #         lty=c("solid", "solid", "solid", "dashed"),
                            #         lwd=c(3,3,3,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #         bg="white")
                            #       } else {
                            #         legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
                            #         lty=c("solid", "solid", "solid", "solid", "dashed"),
                            #         lwd=c(3,3,3,3,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #         bg="white")
                            #       }
                            #     }
                            #     else
                            #     {
                            #       if (NO_UNIQUE_SEQUENCE) {
                            #         legend("topright",
                            #         ##legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "errors", "kmer-peaks","cov-threshold"),
                            #         lty=c("solid", "solid", "solid", "dashed", "dashed"),
                            #         lwd=c(3,3,3,2,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK, COLOR_COVTHRES),
                            #         bg="white")
                            #       } else {
                            #         legend("topright",
                            #         ##legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks","cov-threshold"),
                            #         lty=c("solid", "solid", "solid", "solid", "dashed", "dashed"),
                            #         lwd=c(3,3,3,3,2,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK, COLOR_COVTHRES),
                            #         bg="white")
                            #       }
                            #     }

                            #     dev.set(dev.next())

                            #   ## Finish Log plot
                            #     title(paste("\n\nlen:",  prettyNum(total_len[1], big.mark=","),
                            #     "bp",
                            #     " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
                            #     "% ", "\n",
                            #     hetline, "\n",
                            #     " kcov:", format(akcov, digits=3),
                            #     " err:",   format(100*error_rate[1], digits=3),
                            #     "% ",
                            #     " dup:",  format(adups, digits=3),
                            #     " ",
                            #     " k:",   format(k, digits=3),
                            #     " p:",   format(p, digits=3),
                            #     sep=""),
                            #     cex.main=.85)

                            #     ## Mark the modes of the peaks
                            #     abline(v=akcov * (1:(2*p)), col=COLOR_KMERPEAK, lty=2)

                            #     ## Draw just the unique portion of the model
                            #     if (!NO_UNIQUE_SEQUENCE) {
                            #       lines(x, unique_hist_transform, col=COLOR_pPEAK, lty=1, lwd=3)
                            #     }
                            #     lines(x, pred, col=COLOR_2pPEAK, lwd=3)
                            #     lines(x[1:error_xcutoff_ind], (x[1:error_xcutoff_ind]**transform_exp)*error_kmers, lwd=3, col=COLOR_ERRORS)

                            #     if (VERBOSE) {
                            #       lines(x, residual_transform, col=COLOR_RESIDUAL, lwd=3)
                            #     }

                            #     ## Add legend
                            #     if(length(kmer_hist[,1])==length(kmer_hist_orig[,1]))
                            #     {
                            #       if (NO_UNIQUE_SEQUENCE) {
                            #         legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "errors", "kmer-peaks"),
                            #         lty=c("solid", "solid", "solid", "dashed"),
                            #         lwd=c(3,3,3,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #         bg="white")
                            #       } else {
                            #         legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
                            #         lty=c("solid", "solid", "solid", "solid", "dashed"),
                            #         lwd=c(3,3,3,3,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #         bg="white")
                            #       }
                            #     }
                            #     else
                            #     {
                            #       if (NO_UNIQUE_SEQUENCE) {
                            #         legend("topright",
                            #         ##legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "errors", "kmer-peaks","cov-threshold"),
                            #         lty=c("solid", "solid", "solid", "dashed", "dashed"),
                            #         lwd=c(3,3,3,2,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK, COLOR_COVTHRES),
                            #         bg="white")
                            #       } else {
                            #         legend("topright",
                            #         ##legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks","cov-threshold"),
                            #         lty=c("solid", "solid", "solid", "solid", "dashed", "dashed"),
                            #         lwd=c(3,3,3,3,2,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK, COLOR_COVTHRES),
                            #         bg="white")
                            #       }
                            #     }

                                ## other things (also turned off)

                            # model_status="done"

                            # if (!IN_VERBOSE) {
                            #   cat(paste("Model converged het:", format(ahet, digits=3),
                            #   " kcov:", format(akcov, digits=3),
                            #   " err:", format(error_rate[1], digits=3),
                            #   " model fit:", format(adups, digits=3),
                            #   " len:", round(total_len[1]), "\n", sep=""))
                            # }
                            # }
                            # else
                            # {
                            #   title("\nFailed to converge")
                            #   dev.set(dev.next())
                            #   title("\nFailed to converge")
                            #   cat("Failed to converge.", file=paste(foldername,"/", arguments$name_prefix, "model.txt", sep=""))
                            #   cat("Failed to converge.\n")
                            # }

                            # dev.off()
                            # dev.off()
                            # dev.off()
                            # dev.off()

                                ## Write key values to summary file

                                    summaryFile <- paste(foldername,"/", arguments$name_prefix, "summary.txt",sep="")

                                    format_column_1 = "%-30s"
                                    format_column_2 = "%-18s"
                                    format_column_3 = "%-18s"

                                    cat(paste("GenomeScope version 2.0", sep=""), file=summaryFile, sep="\n")
                                    cat(paste("input file = ", arguments$input, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    cat(paste("output directory = ", arguments$output, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    cat(paste("p = ", p,sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    cat(paste("k = ", k,sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    if (arguments$name_prefix!="") {
                                      cat(paste("name prefix = ", substring(arguments$name_prefix,1,nchar(arguments$name_prefix)-1), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (arguments$lambda!=-1) {
                                      cat(paste("initial kmercov estimate = ", arguments$lambda, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (arguments$max_kmercov!=-1) {
                                      cat(paste("max_kmercov = ", arguments$max_kmercov, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (VERBOSE) {
                                      cat(paste("VERBOSE set to TRUE", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (NO_UNIQUE_SEQUENCE) {
                                      cat(paste("NO_UNIQUE_SEQUENCE set to TRUE", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (topology!=0) {
                                      cat(paste("topology = ", topology, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (d_init!=-1) {
                                      cat(paste("initial repetitiveness = ", d_init, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (r_inits!=-1) {
                                      cat(paste("initial heterozygosities = ", r_inits, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (transform_exp != 1) {
                                      cat(paste("TRANSFORM_EXP = ", transform_exp, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (TESTING) {
                                      cat(paste("TESTING set to TRUE", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (TRUE_PARAMS != -1) {
                                      cat(paste("TRUE_PARAMS = ", TRUE_PARAMS, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (TRACE_FLAG) {
                                      cat(paste("TRACE_FLAG set to TRUE", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (NUM_ROUNDS != 4) {
                                      cat(paste("NUM_ROUNDS = ", NUM_ROUNDS, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }

                                    property <- vector()
                                    min <- vector()
                                    max <- vector()

                                    cat(paste("\n",sprintf(format_column_1,"property"),               sprintf(format_column_2,"min"),                              sprintf(format_column_3,"max"), sep=""),                                     file=summaryFile, sep="\n", append=TRUE)
                                    if (p==1) {
                                      cat(paste(sprintf(format_column_1,"Homozygous (a)"),            sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)

                                    }
                                    if (p==2) {
                                      cat(paste(sprintf(format_column_1,"Homozygous (aa)"),           sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1,"Heterozygous (ab)"),         sprintf(format_column_2,percentage_format(hets[[1]][1])),         sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      property <- c(property, "Homozygous (aa)", "Heterozygous (ab)")
                                      min <- c(min, homo[2]*100, hets[[1]][1]*100)
                                      max <- c(max, homo[1]*100, hets[[1]][2]*100)
                                    }
                                    if (p==3) {
                                      cat(paste(sprintf(format_column_1,"Homozygous (aaa)"),          sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1,"Heterozygous (not aaa)"),    sprintf(format_column_2,percentage_format(het[1])),          sprintf(format_column_3,percentage_format(het[2])),  sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1,"aab"),                       sprintf(format_column_2,percentage_format(hets[[1]][1])),         sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1,"abc"),                       sprintf(format_column_2,percentage_format(hets[[2]][1])),         sprintf(format_column_3,percentage_format(hets[[2]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (p==4) {
                                      cat(paste(sprintf(format_column_1,"Homozygous (aaaa)"),                    sprintf(format_column_2,percentage_format(homo[2])),      sprintf(format_column_3,percentage_format(homo[1])),      sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1,"Heterozygous (not aaaa)"),              sprintf(format_column_2,percentage_format(het[1])),       sprintf(format_column_3,percentage_format(het[2])),       sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1, switch(top+1, "aaab", "aaab", "aabb")), sprintf(format_column_2,percentage_format(hets[[1]][1])), sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1, switch(top+1, "aabb", "aabc", "aabc")), sprintf(format_column_2,percentage_format(hets[[2]][1])), sprintf(format_column_3,percentage_format(hets[[2]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1, switch(top+1, "aabc", "abcd", "abcd")), sprintf(format_column_2,percentage_format(hets[[3]][1])), sprintf(format_column_3,percentage_format(hets[[3]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      if (top == 0) {
                                        cat(paste(sprintf(format_column_1,"abcd"),                               sprintf(format_column_2,percentage_format(hets[[4]][1])), sprintf(format_column_3,percentage_format(hets[[4]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      }
                                    }
                                    if (p==5) {
                                      cat(paste(sprintf(format_column_1,"Homozygous (aaaaa)"),                                                sprintf(format_column_2,percentage_format(ahomo)), sep=""),      file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1,"Heterozygous (not aaaaa)"),                                          sprintf(format_column_2,percentage_format(ahet)),  sep=""),      file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1,"aaaab", "aaaab", "aaaab", "aaabb", "aaabb", "aaabb")), sprintf(format_column_2,percentage_format(hets[[1]][1])), sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1,"aaabb", "aaabc", "aabbc", "aaabc", "aabcc", "aabcc")), sprintf(format_column_2,percentage_format(hets[[2]][1])), sprintf(format_column_3,percentage_format(hets[[2]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1,"aaabc", "aabcd", "aabcd", "aabcd", "aabcd", "abcdd")), sprintf(format_column_2,percentage_format(hets[[3]][1])), sprintf(format_column_3,percentage_format(hets[[3]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1,"aabbc", "abcde", "abcde", "abcde", "abcde", "abcde")), sprintf(format_column_2,percentage_format(hets[[4]][1])), sprintf(format_column_3,percentage_format(hets[[4]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      if (top == 0) {
                                        #cat(paste(sprintf(format_column_1,"aabcd"),                     sprintf(format_column_2,percentage_format(hets[[5]][1])),         sprintf(format_column_3,percentage_format(hets[[5]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                        #cat(paste(sprintf(format_column_1,"abcde"),                     sprintf(format_column_2,percentage_format(hets[[6]][1])),         sprintf(format_column_3,percentage_format(hets[[6]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      }
                                    }
                                    if (p==6) {
                                      cat(paste(sprintf(format_column_1,"Homozygous (aaaaaa)"),       sprintf(format_column_2,percentage_format(ahomo)), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1,"Heterozygous (not aaaaaa)"), sprintf(format_column_2,percentage_format(ahet)), sep=""),                 file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1, "aaaaab", "aaaaab:", "aaaaab:", "aaaaab:", "aaaaab:", "aaaaab:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaabbb:", "aaabbb:", "aaabbb:")),                    sprintf(format_column_2,percentage_format(hets[[1]][1])),         sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1, "aaaabb", "aaaabc:", "aaaabc:", "aaabbc:", "aaabbc:", "aaabbc:", "aaaabc:", "aaaabc:", "aaabcc:", "aaabcc:", "aaabcc:", "aabbcc:", "aabbcc:", "aabbcc:", "aaabbc:", "aaabbc:", "aaabbc:")),                    sprintf(format_column_2,percentage_format(hets[[2]][1])),         sprintf(format_column_3,percentage_format(hets[[2]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1, "aaabbb", "aaabcd:", "aabbcd:", "aaabcd:", "aabccd:", "aabccd:", "aaabcd:", "aabbcd:", "aaabcd:", "aabcdd:", "aabcdd:", "aabbcd:", "aabcdd:", "aabcdd:", "aaabcd:", "aabccd:", "aabccd:")),                    sprintf(format_column_2,percentage_format(hets[[3]][1])),         sprintf(format_column_3,percentage_format(hets[[3]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1, "aaaabc", "aabcde:", "aabcde:", "aabcde:", "aabcde:", "abcdde:", "aabcde:", "aabcde:", "aabcde:", "aabcde:", "abcdee:", "aabcde:", "aabcde:", "abcdee:", "aabcde:", "aabcde:", "abcdde:")),                    sprintf(format_column_2,percentage_format(hets[[4]][1])),         sprintf(format_column_3,percentage_format(hets[[4]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1, "aaabbc", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:")),                    sprintf(format_column_2,percentage_format(hets[[5]][1])),         sprintf(format_column_3,percentage_format(hets[[5]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      if (top==0) {
                                        #cat(paste(sprintf(format_column_1,"aabbcc"),                    sprintf(format_column_2,percentage_format(hets[[6]][1])),         sprintf(format_column_3,percentage_format(hets[[6]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                        #cat(paste(sprintf(format_column_1,"aaabcd"),                    sprintf(format_column_2,percentage_format(hets[[7]][1])),         sprintf(format_column_3,percentage_format(hets[[7]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                        #cat(paste(sprintf(format_column_1,"aabbcd"),                    sprintf(format_column_2,percentage_format(hets[[8]][1])),         sprintf(format_column_3,percentage_format(hets[[8]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                        #cat(paste(sprintf(format_column_1,"aabcde"),                    sprintf(format_column_2,percentage_format(hets[[9]][1])),         sprintf(format_column_3,percentage_format(hets[[9]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                        #cat(paste(sprintf(format_column_1,"abcdef"),                    sprintf(format_column_2,percentage_format(hets[[10]][1])),        sprintf(format_column_3,percentage_format(hets[[10]][2])), sep=""),               file=summaryFile, sep="\n", append=TRUE)
                                      }
                                    }
                                    cat(paste(sprintf(format_column_1,"Genome Haploid Length"), sprintf(format_column_2,bp_format(total_len[2])),                  sprintf(format_column_3,bp_format(total_len[1])), sep=""),                   file=summaryFile, sep="\n", append=TRUE)
                                    cat(paste(sprintf(format_column_1,"Genome Repeat Length"),  sprintf(format_column_2,bp_format(repeat_len[2])),                 sprintf(format_column_3,bp_format(repeat_len[1])), sep=""),                  file=summaryFile, sep="\n", append=TRUE)
                                    cat(paste(sprintf(format_column_1,"Genome Unique Length"),  sprintf(format_column_2,bp_format(unique_len[2])),                 sprintf(format_column_3,bp_format(unique_len[1])), sep=""),                  file=summaryFile, sep="\n", append=TRUE)
                                    cat(paste(sprintf(format_column_1,"Model Fit "),            sprintf(format_column_2,percentage_format(model_fit_allscore[1])), sprintf(format_column_3,percentage_format(model_fit_fullscore[1])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    cat(paste(sprintf(format_column_1,"Read Error Rate"),       sprintf(format_column_2,percentage_format(error_rate[1])),         sprintf(format_column_3,percentage_format(error_rate[2])), sep=""),          file=summaryFile, sep="\n", append=TRUE)

                                    property <- c(property, "Genome Haploid Length", "Genome Repeat Length", "Genome Unique Length", "Model Fit", "Read Error Rate")
                                    min <- c(min, total_len[2], repeat_len[2], unique_len[2], model_fit_allscore[1], error_rate[1])
                                    max <- c(max, total_len[1], repeat_len[1], unique_len[1], model_fit_fullscore[1], error_rate[2])

                                    output$summary <- data.frame(
                                      property = property,
                                      min = min,
                                      max = max
                                    )

                                    # if (VERBOSE)
                                    # {
                                    #   cat(paste("\nPercent Kmers Modeled (All Kmers) = ",  percentage_format(model_fit_allscore[1]),    " [", model_fit_allscore[2],    ", ", model_fit_allscore[3],    "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    #   cat(paste("Percent Kmers Modeled (Full Model) = ",   percentage_format(model_fit_fullscore[1]),   " [", model_fit_fullscore[2],   ", ", model_fit_fullscore[3],   "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    #   cat(paste("Percent Kmers Modeled (Unique Kmers) = ", percentage_format(model_fit_uniquescore[1]), " [", model_fit_uniquescore[2], ", ", model_fit_uniquescore[3], "]", sep=""), file=summaryFile, sep="\n", append=TRUE)

                                    #   cat(paste("\nModel RSSE (All Kmers) = ",  model_fit_all[1],    " [", model_fit_all[2],    ", ", model_fit_all[3],    "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    #   cat(paste("Model RSSE (Full Model) = ",   model_fit_full[1],   " [", model_fit_full[2],   ", ", model_fit_full[3],   "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    #   cat(paste("Model RSSE (Unique Model) = ", model_fit_unique[1], " [", model_fit_unique[2], ", ", model_fit_unique[3], "]", sep=""), file=summaryFile, sep="\n", append=TRUE) 
                                    # }

                                    ## Merfin probabilities (off)

                                # colors <- c("black","red","green","purple","blue")

                                # if (p<=2 & FITTED_HIST){

                                #   if (!(p==1)){

                                #   fitted_hist=data.frame(cbind(one_hist,two_hist,thr_hist,fou_hist))

                                #   } else {

                                #   fitted_hist=data.frame(cbind(two_hist,fou_hist))

                                #   }

                                #   ## Fitted histogram  

                                #     png(paste(foldername, "/fitted_hist.png", sep=""), height = plot_size, width = plot_size, res=resolution)
                                #     layout(matrix(c(1,2), nrow=2, byrow = TRUE),heights=lcm(c(11,5.5)))
                                #     par(mar=c(0,5,1,1))
                                  
                                #     plot(kmer_hist_orig, type="n", xlab=NULL, ylab=ylabel_orig, ylim=c(0,y_limit_orig), xlim=c(0,akcov*(p*2.5)),
                                #       cex.main=font_size, cex.sub=font_size, cex.axis=font_size, cex.lab=font_size,
                                #        xaxt='n', yaxt='n')
                                #     myTicks = axTicks(4)
                                #     myTicks[1]<-0
                                #     axis(2, at=myTicks, labels=myTicks, cex.axis=font_size,
                                #        cex.lab=font_size)
                                  
                                #     peaks<-akcov * 1:6

                                #     points(kmer_hist, type="h", col=COLOR_HIST, lwd=2)
                                #     lines(error_kmers, col="black", lwd=1.5)
                                  
                                #     abline(v=peaks, col="black", lty=2, lwd=0.3)
                                  
                                #     for (i in 1:ncol(fitted_hist)) {
                                #     lines(fitted_hist[,i]*amlen, type="l", lwd=1.5, col=colors[i+1])
                                #     }

                                #     lines(rowSums(fitted_hist*amlen), col="darkgray", lwd=2, lty=2)
                                  
                                #     legend_names=c("0-copy", "1-copy", "2-copy", "3-copy", "4-copy")
                                #     legend_lty=c("solid", "solid", "solid", "solid", "solid")
                                #     legend_lwd=c(3,3,3,3,3)
                                  
                                #     legend("topright",
                                #        ##legend(exp(.65 * log(max(x))), 1.0 * max(y),
                                #        legend=c("Observed",legend_names[1:(ncol(fitted_hist)+1)],"Full model"),
                                #        lty=c("solid",legend_lty[1:(ncol(fitted_hist)+1)], "dashed"),
                                #        lwd=c(2, legend_lwd[1:(ncol(fitted_hist)+1)], 3),
                                #        col=c(COLOR_HIST,colors[1:(ncol(fitted_hist)+1)], "darkgray"),
                                #        bg="white")

                                #     ## Generate lookup_table
                                  
                                #       lookup_table <- NULL
                                #       plot_table <- NULL
                                    
                                #       fitted_hist[(akcov*(p*2.5)):length(one_hist),1:ncol(fitted_hist)] <- 0
                                    
                                #       fitted_hist <- na.zero(fitted_hist)
                                    
                                #       for (i in 1:(akcov*(p*2.5)-1)){
                                    
                                #       totalP<-sum(na.zero(error_kmers[i]), rowSums(fitted_hist[i,]*amlen))
                                   
                                #       prob<-c(na.zero(error_kmers[i]/totalP),as.numeric(fitted_hist[i,]*amlen)/totalP)
                                    
                                #       plot_table<-rbind(plot_table,prob)
                                      
                                #       max.p<-max(prob)  
                                #       readK<-which.max(prob)-1
                                    
                                #       lookup_table<-rbind(lookup_table,c(readK,max.p))
                                    
                                #       }
                                    
                                #       lookup_table<-data.frame(lookup_table)
                                #       rownames(plot_table) <- make.names(plot_table[,1], unique = TRUE)
                                    
                                #       write.table(lookup_table, paste(foldername, "/lookup_table.txt", sep=""), sep=",", row.names = FALSE, col.names = FALSE)

                                #     ## Plot lookup values
                                  
                                #       par(mar=c(5,5,0,1))
                                    
                                #       plot_table<-data.frame(plot_table)
                                #       plot(plot_table$X1,type="n", xlab="Coverage", ylab="Probability",xlim=c(0,akcov*(p*2.5)), ylim=c(0,1), cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size, yaxt='n')

                                #       abline(v=peaks, col="black", lty=2, lwd=0.3)
                                    
                                #       for (i in 1:(ncol(fitted_hist)+1)) {
                                #       lines(plot_table[,i], type="l", lwd=1.5, col=colors[i])
                                #       }
                                #       axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), cex.axis=font_size,
                                #          cex.lab=font_size)
                                  
                                #     dev.off()
                                  
                                # }

                                    ## Some sort of advanced testing (turned off)

                            # if (TESTING) {

                            #   if (TRUE_PARAMS!=-1) {
                            #     testingFile <- paste(foldername, "/", arguments$name_prefix, "SIMULATED_testing.tsv",sep="")
                            #   } else {
                            #     testingFile <- paste(foldername,"/SIMULATED_testing.tsv",sep="")
                            #   }
                            #   if (p==1) {
                            #     cat(paste(amd, akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #   }
                            #   if (p==2) {
                            #     cat(paste(amd, ahets[[1]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #   }
                            #   if (p==3) {
                            #     if (TRUE_PARAMS!=-1) {
                            #       true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
                            #       cat(paste(amd, ahets[[1]], ahets[[2]], akcov, adups, atotal_len, top, true_params[1], true_params[2], true_params[3], true_params[4], sep="\t"), file=testingFile, sep="\n", append=FALSE)
                            #     } else {
                            #       cat(paste(amd, ahets[[1]], ahets[[2]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #     }
                            #   }
                            #   if (p==4) {
                            #     if (topology==0) {
                            #       if (TRUE_PARAMS!=-1) {
                            #         true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], akcov, adups, atotal_len, top, true_params[1], true_params[2], true_params[3], true_params[4], true_params[5], true_params[6], sep="\t"), file=testingFile, sep="\n", append=FALSE)
                            #       } else {
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #       }
                            #     } else {
                            #       if (TRUE_PARAMS!=-1) {
                            #         true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
                            #         cat(paste(amd, switch(top, ahets[[1]], 0), switch(top, 0, ahets[[1]]), ahets[[2]], ahets[[3]], akcov, adups, atotal_len, top, true_params[1], switch(true_params[5], true_params[2], 0), switch(true_params[5], 0, true_params[2]), true_params[3], true_params[4], true_params[5], sep="\t"), file=testingFile, sep="\n", append=FALSE)
                            #       } else {
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #       }
                            #     }
                            #   }
                            #   if (p==5) {
                            #     if (topology==0) {
                            #       if (TRUE_PARAMS!=-1) {
                            #         true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], akcov, adups, atotal_len, top, true_params[1], true_params[2], true_params[3], true_params[4], true_params[5], true_params[6], true_params[7], true_params[8], sep="\t"), file=testingFile, sep="\n", append=FALSE)
                            #       } else {
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #       }
                            #     } else {
                            #       if (TRUE_PARAMS!=-1) {
                            #         true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
                            #         cat(paste(amd, switch(top, ahets[[1]], ahets[[1]], 0, 0, 0), switch(top, 0, 0, ahets[[1]], ahets[[1]], ahets[[1]]), switch(top, ahets[[2]], 0, ahets[[2]], 0, 0), switch(top, 0, ahets[[2]], 0, 0, 0), switch(top, 0, 0, 0, ahets[[2]], ahets[[2]]), switch(top, ahets[[3]], ahets[[3]], ahets[[3]], ahets[[3]], 0), switch(top, 0, 0, 0, 0, ahets[[3]]), ahets[[4]], akcov, adups, atotal_len, top, true_params[1], switch(true_params[6], true_params[2], true_params[2], 0, 0, 0), switch(true_params[6], 0, 0, true_params[2], true_params[2], true_params[2]), switch(true_params[6], true_params[3], 0, true_params[3], 0, 0), switch(true_params[6], 0, true_params[3], 0, 0, 0), switch(true_params[6], 0, 0, 0, true_params[3], true_params[3]), switch(true_params[6], true_params[4], true_params[4], true_params[4], true_params[4], 0), switch(true_params[6], 0, 0, 0, 0, true_params[4]), true_params[5], true_params[6], sep="\t"), file=testingFile, sep="\n", append=FALSE)
                            #       } else {
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #       }
                            #     }
                            #   }
                            #   if (p==6) {
                            #     if (topology==0) {
                            #       if (TRUE_PARAMS!=-1) {
                            #         true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], ahets[[7]], ahets[[8]], ahets[[9]], ahets[[10]], akcov, adups, atotal_len, top, true_params[1], true_params[2], true_params[3], true_params[4], true_params[5], true_params[6], true_params[7], true_params[8], true_params[9], true_params[10], true_params[11], true_params[12], sep="\t"), file=testingFile, sep="\n", append=FALSE)
                            #       } else {
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], ahets[[7]], ahets[[8]], ahets[[9]], ahets[[10]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #       }
                            #     } else {
                            #       if (TRUE_PARAMS!=-1) {
                            #         true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
                            #         cat(paste(amd, ifelse(top %in% c(1,2,3,4,5), ahets[[1]], 0), ifelse(top %in% c(6,7,8,9,10,11,12,13), ahets[[1]], 0), ifelse(top %in% c(14,15,16), ahets[[1]], 0), ifelse(top %in% c(1,2,6,7), ahets[[2]], 0), ifelse(top %in% c(3,4,5,14,15,16), ahets[[2]], 0), ifelse(top %in% c(8,9,10), ahets[[2]], 0), ifelse(top %in% c(11,12,13), ahets[[2]], 0), ifelse(top %in% c(1,3,6,8,14), ahets[[3]], 0), ifelse(top %in% c(2,7,11), ahets[[3]], 0), ifelse(top %in% c(4,5,15,16), ahets[[3]], 0), ifelse(top %in% c(9,10,12,13), ahets[[3]], 0), ifelse(top %in% c(1,2,3,4,6,7,8,9,11,12,14,15), ahets[[4]], 0), ifelse(top %in% c(5,16), ahets[[4]], 0), ifelse(top %in% c(10,13), ahets[[4]], 0), ahets[[5]], akcov, adups, atotal_len, top, true_params[1], ifelse(true_params[7] %in% c(1,2,3,4,5), true_params[2], 0), ifelse(true_params[7] %in% c(6,7,8,9,10,11,12,13), true_params[2], 0), ifelse(true_params[7] %in% c(14,15,16), true_params[2], 0), ifelse(true_params[7] %in% c(1,2,6,7), true_params[3], 0), ifelse(true_params[7] %in% c(3,4,5,14,15,16), true_params[3], 0), ifelse(true_params[7] %in% c(8,9,10), true_params[3], 0), ifelse(true_params[7] %in% c(11,12,13), true_params[3], 0), ifelse(true_params[7] %in% c(1,3,6,8,14), true_params[4], 0), ifelse(true_params[7] %in% c(2,7,11), true_params[4], 0), ifelse(true_params[7] %in% c(4,5,15,16), true_params[4], 0), ifelse(true_params[7] %in% c(9,10,12,13), true_params[4], 0), ifelse(true_params[7] %in% c(1,2,3,4,6,7,8,9,11,12,14,15), true_params[5], 0), ifelse(true_params[7] %in% c(5,16), true_params[5], 0), ifelse(true_params[7] %in% c(10,13), true_params[5], 0), true_params[6], true_params[7], sep="\t"), file=testingFile, sep="\n", append=FALSE)
                            #       } else {
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #       }
                            #     }
                            #   }
                            # }

                            }

                    ## Make a ggplot-friendly version of the plots

                        output$plots <- data.frame(
                            kmer_hist_x = kmer_hist_orig[,1],
                            kmer_hist_y = kmer_hist_orig[,2],
                            kmer_hist_unique = unique_hist,
                            kmer_hist_predicted = x**(-transform_exp)*pred,
                            kmer_hits_error = c(error_kmers, rep(0,length(x)-length(error_kmers))),
                            kmer_hist_transformed_x = kmer_hist_transform[,1],
                            kmer_hist_transformed_y = kmer_hist_transform[,2],
                            kmer_hist_unique_transformed = unique_hist_transform,
                            kmer_hist_predicted_transformed = x**(-transform_exp)*pred,
                            kmer_hits_error = c(((x[1:error_xcutoff_ind]**transform_exp)*error_kmers), rep(0,length(x)-length(error_kmers)))
                        )

                    ## the ggplots

                        # thing1 <- ggplot() +
                        #   geom_col(data = kmer_hist_orig, aes(x = V1, y = V2), color = "black", fill = "white") +
                        #   geom_line(data = data.frame(
                        #     x = x, y = unique_hist
                        #   ), aes(x = x, y = y), color = "darkgreen", size = 2, alpha = 0.6) +
                        #   geom_line(data = data.frame(
                        #     x = x, y = x**(-transform_exp)*pred
                        #   ), aes(x = x, y = y), color = "darkblue", size = 2, alpha = 0.6) +
                        #   geom_line(data = data.frame(
                        #     x = x[1:error_xcutoff_ind], y = error_kmers
                        #   ), aes(x = x, y = y), color = "red", size = 2, alpha = 0.6) +
                        #   scale_x_continuous(limits = c(0,x_limit_orig)) +
                        #   scale_y_continuous(limits = c(0,y_limit_orig)) +
                        #   theme_bw()

                        # thing2 <- ggplot() +
                        #   geom_col(data = kmer_hist_transform, aes(x = V1, y = V2), color = "black", fill = "white") +
                        #   geom_line(data = data.frame(
                        #     x = x, y = unique_hist_transform
                        #   ), aes(x = x, y = y), color = "darkgreen", size = 2, alpha = 0.6) +
                        #   geom_line(data = data.frame(
                        #     x = x, y = pred
                        #   ), aes(x = x, y = y), color = "darkblue", size = 2, alpha = 0.6) +
                        #   geom_line(data = data.frame(
                        #     x = x[1:error_xcutoff_ind], y = (x[1:error_xcutoff_ind]**transform_exp)*error_kmers
                        #   ), aes(x = x, y = y), color = "red", size = 2, alpha = 0.6) +
                        #   scale_x_continuous(limits = c(0,x_limit)) +
                        #   scale_y_continuous(limits = c(0,y_limit)) +
                        #   theme_bw()

                        # thing3 <- ggplot() +
                        #   geom_col(data = kmer_hist_orig, aes(x = log(V1), y = log(V2)), color = "black", fill = "white") +
                        #   # geom_line(data = data.frame(
                        #   #   x = x, y = unique_hist
                        #   # ), aes(x = x, y = y), color = "darkgreen", size = 2, alpha = 0.6) +
                        #   # scale_x_continuous(limits = c(0,x_limit_orig)) +
                        #   # scale_y_continuous(limits = c(0,y_limit_orig)) +
                        #   theme_bw()

                        # thing4 <- ggplot() +
                        #   geom_col(data = kmer_hist_transform, aes(x = log(V1), y = log(V2)), color = "black", fill = "white") +
                        #   # geom_line(data = data.frame(
                        #   #   x = x, y = unique_hist_transform
                        #   # ), aes(x = x, y = y), color = "darkgreen", size = 2, alpha = 0.6) +
                        #   # scale_x_continuous(limits = c(0,x_limit)) +
                        #   # scale_y_continuous(limits = c(0,y_limit)) +
                        #   theme_bw()

                        output
            }

    ##### Chemical data handling

        #### convertCDFstoCSVs

            #' Convert mass spectral datafiles (CDF) into a csv file
            #'
            #' @param paths_to_cdfs Paths to CDF files
            #' @param min_mz Smallest m/z value to process (acts like a filter)
            #' @param max_mz Largest m/z value to process (acts like a filter)
            #' @param min_rt (optional) Smallest retention time to process (acts like a filter)
            #' @param max_rt (optional) Largest retention time to process (acts like a filter)
            #' @param force Convert all CDF files, even if they've already been converted previously
            #' @examples
            #' @export
            #' convertCDFstoCSVs

                convertCDFstoCSVs <- 	function(
    						            	paths_to_cdfs, 
    						            	min_mz = 50, 
    						            	max_mz = 800, 
    						            	min_rt = NULL, 
    						            	max_rt = NULL, 
    						            	force = FALSE
    						            ) {

                    if ( force == FALSE ) {
                        paths_to_cdfs <- paths_to_cdfs[!file.exists(paste0(paths_to_cdfs, ".csv"))]
                    }

                    if (length(paths_to_cdfs) > 0) {

    	                for (file in 1:length(paths_to_cdfs)) {

    	                    cat(paste("Reading data file ", paths_to_cdfs[file], sep = ""))
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

    	                    cat("   Filling in blank m/z values ...\n")
    	                        spreadDataFile <- tidyr::spread(framedDataFile, rt, intensity, fill = 0)
    	                        
    	                        numbers <- round(seq(min_mz, max_mz, 0.1), digits = 1)
    	                        add <- numbers[!numbers %in% spreadDataFile$mz]
    	                        addition <- spreadDataFile[1:length(add),]
    	                        
    	                        if ( length(add) > 0 ) {
    	                            addition$mz <- add     
    	                        }
    	                        
                                addition <- tidyr::pivot_longer(addition, cols = 2:dim(addition)[2], names_to = "rt", values_to = "intensity")
    	                        addition$intensity <- 0
    	                        addition$rt <- as.numeric(addition$rt)
    	                        
                                framedDataFile <- tidyr::pivot_longer(spreadDataFile, cols = 2:dim(spreadDataFile)[2], names_to = "rt", values_to = "intensity")

                            cat("   still working ...\n")
    	                        framedDataFile <- rbind(framedDataFile, addition)

                            cat("   still working ... patience please\n")
    	                        framedDataFile$rt <- as.numeric(framedDataFile$rt)
    	                        framedDataFile <- framedDataFile[sort.list(framedDataFile$mz),]

                            cat("   still working ... almost done\n")
    	                        framedDataFile <- framedDataFile[sort.list(framedDataFile$rt),]
    	                        rownames(framedDataFile) <- NULL

    	                    if ( length(min_rt) > 0 ) {
    	                        cat("   Filtering by minimum retention time thresholds ...\n")
    	                            framedDataFile <- dplyr::filter(framedDataFile, rt > min_rt)
    	                    }

    	                    if ( length(max_rt) > 0 ) {
    	                        cat("   Filtering by maximum retention time thresholds ...\n")
    	                            framedDataFile <- dplyr::filter(framedDataFile, rt < max_rt)    
    	                    }

    	                    cat("   Writing out data file as CSV... \n\n")
    	                        data.table::fwrite(framedDataFile, file = paste(paths_to_cdfs[file], ".csv", sep = ""), col.names = TRUE, row.names = FALSE)
    	                }
    	            }
                }

        #### extractChromatogramsFromCSVs

            #' Extract total ion chromatograms from csv files containing mass spectral data
            #'
            #' @param paths_to_cdf_csvs Paths to the cdf.csvs
            #' @param path_to_existing_chromatograms_csv Path to an existing chromatograms.csv file. The function will append new chromatograms to this file that are not already in it.
            #' @examples
            #' @export
            #' extractChromatogramsFromCSVs

                extractChromatogramsFromCSVs <- function( paths_to_cdf_csvs, path_to_existing_chromatograms_csv = NULL ) {

                    ## If the chromatograms.csv already exists, figure out which chromatograms are already in there and don't analyze those
                        if ( is.null(path_to_existing_chromatograms_csv) == FALSE ) {
                            paths_to_cdf_csvs <- paths_to_cdf_csvs[!paths_to_cdf_csvs %in% unique(readMonolist(path_to_existing_chromatograms_csv)$path_to_cdf_csv)]
                        }

                    ## Extract the chromatograms from the full dataset
                        chromatograms <- list()
                        for ( file in 1:length(paths_to_cdf_csvs) ) {

                            cat(paste("Reading data file ", paths_to_cdf_csvs[file], "\n", sep = ""))
                                framedDataFile <- as.data.frame(data.table::fread(paths_to_cdf_csvs[file]))

                            cat("   Extracting the total ion chromatogram...\n")
                                # library(plyr)
                                # chromatogram <- plyr::ddply(framedDataFile, .(rt), summarize, tic = sum(intensity))
                                framedDataFile %>% group_by(rt) %>% summarize(tic = sum(intensity)) -> chromatogram
                                chromatogram$rt <- as.numeric(chromatogram$rt)
                                chromatogram$path_to_cdf_csv <- paths_to_cdf_csvs[file]

                            cat("   Appending chromatogram to list...\n")
                                chromatograms[[file]] <- chromatogram
                        }

                        chromatograms <- do.call(rbind, chromatograms)

                    ## If the chromatograms file already exists, append to it, otherwise, return the chromatograms 
                        if ( is.null(path_to_existing_chromatograms_csv) == FALSE ) {
                            writeMonolist(
                                monolist = rbind(
                                    readMonolist(path_to_existing_chromatograms_csv),
                                    chromatograms
                                ),
                                monolist_out_path = path_to_existing_chromatograms_csv
                            )
                        } else {
                            return( chromatograms )        
                        }                
                }

        #### mergePeakLists

    	    #' Merge multiple peaklists
    	    #'
    	    #' @param analysis_directory_path Path to the directory above all the peak lists
    	    #' @examples
    	    #' @export
    	    #' mergePeakLists

    	    mergePeakLists <- function(analysis_directory_path) {

    	        paths_to_peak_monolists <- paste(
    	            analysis_directory_path,
    	            dir(analysis_directory_path),
    	            "peaks_monolist.csv",
    	            sep = "/"
    	        )

    	        peak_list <- list()
    	        for (i in 1:length(paths_to_peak_monolists)) {
                    
                    if ( file.exists(paths_to_peak_monolists[i]) ) {

                        data <- readMonolist(paths_to_peak_monolists[i])
                        data$area <- as.numeric(data$area)

                        if (dim(data)[1] > 0) {
                            peak_list[[i]] <- data[,1:9]
                        }
                    
                    }
    	        
                }

    	        do.call(rbind, peak_list)

    	    }

        #### analyzeGCMSdata

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
            #' analyzeGCMSdata

            analyzeGCMSdata <- function(
                    CDF_directory_path = getwd(),
                    zoom_and_scroll_rate = 100,
                    baseline_window = 400,
                    x_axis_start_default = NULL,
                    x_axis_end_default = NULL,
                    path_to_reference_library = busta_spectral_library,
                    samples_monolist_subset = NULL,
                    ions = 0,
                    jupyter = FALSE
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
                                        missing_ions <- as.numeric(as.character(ions[!ions %in% ions_for_this_cdf_csv])) ## Here as numeric mess with TIC
                                        missing_ions <- dropNA(missing_ions)

                                    } else {

                                        missing_ions <- ions

                                    }

                                    if (length(missing_ions) > 0) {
                                        
                                        cat(paste("Chromatogram extraction. Reading data file ", paths_to_cdf_csvs[file], "\n", sep = ""))    
                                            framedDataFile <- as.data.frame(data.table::fread(paths_to_cdf_csvs[file]))
                                        
                                        cat("   Extracting chromatograms...\n")
                                            
                                            if (0 %in% ions) {
                                                framedDataFile$row_number <- seq(1,dim(framedDataFile)[1],1)
                                                framedDataFile %>% 
                                                    group_by(rt) %>% summarize(
                                                    abundance = sum(intensity),
                                                    ion = 0,
                                                    rt_first_row_in_raw = min(row_number),
                                                    rt_last_row_in_raw = max(row_number)
                                                ) -> chromatogram
                                                chromatogram <- as.data.frame(chromatogram)
                                                chromatogram$rt <- as.numeric(chromatogram$rt)
                                                chromatogram$path_to_cdf_csv <- paste(paths_to_cdfs[file], ".csv", sep = "")
                                                chromatograms_to_add <- rbind(chromatograms_to_add, chromatogram)
                                            }

                                            if ( length(ions[ions != 0]) > 0 ) {

                                                numeric_ions <- as.numeric(as.character(ions[ions != 0]))
                                                for ( ion in 1:length(numeric_ions) ){
                                                    framedDataFile$row_number <- seq(1,dim(framedDataFile)[1],1)
                                                    framedDataFile %>% 
                                                        group_by(rt) %>% 
                                                        filter(mz > (numeric_ions[ion] - 0.6)) %>%
                                                        filter(mz < (numeric_ions[ion] + 0.6)) -> signal
                                                        summarize(signal,
                                                            abundance = sum(intensity),
                                                            ion = numeric_ions[ion],
                                                            rt_first_row_in_raw = if (dim(signal)[1] > 0) { min(row_number) } else { 0 },
                                                            rt_last_row_in_raw = if (dim(signal)[1] > 0) { max(row_number) } else { 0 }
                                                        ) -> chromatogram
                                                    chromatogram <- as.data.frame(chromatogram)
                                                    chromatogram$rt <- as.numeric(chromatogram$rt)

                                                    if (dim(signal)[1] > 0) { 
                                                        chromatogram$path_to_cdf_csv <- paste(paths_to_cdfs[file], ".csv", sep = "")
                                                        chromatograms_to_add <- rbind(chromatograms_to_add, chromatogram)   
                                                    }
                                                }
                                            }
                                    }

                                ## If the chromatograms file already exists, append to it and re-write out, else create it

                                    if ( file.exists("chromatograms.csv") ) {
                                        
                                        # print("writing it")

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
                                    ),

                                    tags$head(
                                        HTML(
                                            "
                                            <script>
                                                var socket_timeout_interval
                                                var n = 0
                                                $(document).on('shiny:connected', function(event) {
                                                socket_timeout_interval = setInterval(function(){
                                                Shiny.onInputChange('count', n++)
                                                }, 15000)
                                                });
                                                $(document).on('shiny:disconnected', function(event) {
                                                clearInterval(socket_timeout_interval)
                                                });
                                            </script>
                                            "
                                        )
                                    ),

                                    textOutput("keepAlive")
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

                            ## Don't let it time out
                                output$keepAlive <- renderText({
                                    req(input$count)
                                    paste("\nstayin' alive ", input$count)
                                })

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
                                                tic <- filter(chromatogram, ion == 0)

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
                                                    y_axis_end <<- max(chromatograms$abundance)
                                                }

                                            ## If brush is not null, assign brush values to start and end

                                                if ( !is.null(input$chromatogram_brush) ) {
                                                    peak_points <<- isolate(brushedPoints(chromatograms_updated, input$chromatogram_brush))
                                                    x_axis_start <<- min(peak_points$rt)
                                                    x_axis_end <<- max(peak_points$rt)
                                                    y_axis_end <<- max(peak_points$abundance)
                                                    # x_axis_end <<- max(peak_points$rt)
                                                }
                                            
                                            ## Filter chromatogram
                                                
                                                chromatograms_updated_filtered <- dplyr::filter(
                                                    chromatograms_updated, rt_rt_offset > x_axis_start & rt_rt_offset < x_axis_end
                                                )

                                        ## Plot
                                                
                                            facet_labels <- gsub(".CDF.csv", "", gsub(".*/", "", chromatograms_updated_filtered$path_to_cdf_csv))
                                            names(facet_labels) <- chromatograms_updated_filtered$path_to_cdf_csv

                                            chromatogram_plot <- ggplot() +
                                                geom_line(
                                                    data = filter(chromatograms_updated_filtered, ion == "baseline"),
                                                    mapping = aes(x = rt_rt_offset, y = abundance), color = "grey"
                                                ) +
                                                # geom_line(
                                                #     data = filter(chromatograms_updated_filtered, ion != "baseline"),
                                                #     mapping = aes(x = rt_rt_offset, y = abundance, color = ion),
                                                #     alpha = 0.8
                                                # ) +
                                                scale_x_continuous(limits = c(x_axis_start, x_axis_end), name = "Retention (Scan number)") +
                                                scale_y_continuous(limits = c(0, y_axis_end), name = "Abundance (counts)", oob = scales::squish) +
                                                facet_grid(path_to_cdf_csv~., scales = "free_y", labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                                theme_classic() +
                                                guides(fill = "none") +
                                                scale_fill_continuous(type = "viridis") +
                                                scale_color_manual(values = discrete_palette)

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
                                                            ion == 0)$abundance
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
                                                    # print("filter passed")

                                                    for (peak in 1:dim(peak_table)[1]) {
                                                        # print(peak)
                                                        signal_for_this_peak <- dplyr::filter(
                                                            chromatograms_updated[chromatograms_updated$path_to_cdf_csv == peak_table[peak,]$path_to_cdf_csv,], 
                                                            rt_rt_offset > peak_table[peak,]$peak_start_rt_offset, 
                                                            rt_rt_offset < peak_table[peak,]$peak_end_rt_offset
                                                        )
                                                        # print("\n")
                                                        # print("thing1")
                                                        if (dim(signal_for_this_peak)[1] > 0) {

                                                            signal_for_this_peak$peak_number_within_sample <- peak_table$peak_number_within_sample[peak]
                                                        # print("\n")
                                                        # print("thing2")
                                                            ribbon <- filter(signal_for_this_peak, ion == 0)
                                                            ribbon$baseline <- filter(signal_for_this_peak, ion == "baseline")$abundance
                                                        # print("\n")
                                                        # print("thing3")
                                                            chromatogram_plot <- chromatogram_plot +
                                                                geom_vline(data = signal_for_this_peak[1,], mapping = aes(xintercept = rt_rt_offset), alpha = 0.3) +
                                                                geom_ribbon(
                                                                    data = ribbon,
                                                                    mapping = aes(x = rt_rt_offset, ymax = abundance, ymin = baseline, fill = peak_number_within_sample),
                                                                    alpha = 0.8
                                                                ) +
                                                                geom_text(
                                                                    data = filter(signal_for_this_peak, ion == 0),
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
                                            area = sum(peak_points$abundance[peak_points$ion == 0])
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
                        
                        if ( jupyter == TRUE) {

                            available_ports <- vector()

                            if ( length( unique( c(
                                        randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10122, max = 10123, host = "127.0.0.1", n = 50)
                            ))) == 2 ) { available_ports <- c(available_ports, 10123) }

                            if ( length( unique( c(
                                        randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10455, max = 10456, host = "127.0.0.1", n = 50)
                            ))) == 2 ) { available_ports <- c(available_ports, 10456) }

                            if ( length( unique( c(
                                        randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50),
                                        randomPort(min = 10788, max = 10789, host = "127.0.0.1", n = 50)
                            ))) == 2 ) { available_ports <- c(available_ports, 10789) }

                            if (length(available_ports) > 0) {
                                cat(paste0("Connect at: https://shiny", substr(available_ports[1],3,6), ".bustalab.d.umn.edu"))
                                runApp(shinyApp(ui = ui, server = server), host = "127.0.0.1", port = available_ports[1])
                            } else {
                                cat("All available ports are in use. Please try again later.")
                            }
                        
                        }

                        if (jupyter == FALSE) {
                            shinyApp(ui = ui, server = server)
                        }
                    

            }

        #### readChromatograms

            #' Import chromatograms stored as ChemStation exported .csv files
            #'
            #' Allows the user to import chromatograms stored as .csv files. Files should have originated from an Agilent GC system running ChemStation
            #' Files should be in the format of character~var1-var2.csv
            #' "character" is the name of the chromatogram, var1 and var2 are NUMERIC variables you want associated with that chromatogram.
            #' @param dir A directory containing (ONLY) the chromatogram(s) that are to be plotted
            #' @param normalize_level Level at which to normalize chromatogram
            #' @param normalize_range y range to consider during normalization. Useful for normalizing but excluding solvent peak(s), for example
            #' @export
            #' @examples
            #' readChromatograms()

            readChromatograms <- function( dir, normalize_level = c("none", "var1", "var2"), normalize_range = NULL ) {

                ## Get paths to the .csv files in the directory "dir"
                    list <- paste(dir, dir(dir), sep = "/")

                ## Read in time scale = retention times
                    temp <- read.csv(list[1])
                    s <- dim(temp)[1]
                    ret <- as.numeric(as.character(temp[3:s,1]))

                ## Read in and process data into a gathered dataframe
                    data <- data.frame(V1 = as.numeric(1), value = as.numeric(1), character = as.character("1"), var1 = as.character("1"), var2 = as.character("1"))

                        for (i in 1:length(list)){
                            temp <- read.csv(list[i])
                            s <- dim(temp)[1]
                            temp <- temp[3:s,1:2]
                            value <- as.numeric(as.character(temp[,2]))
                            value <- c(rep(0,as.numeric(gsub("-.*$", "", gsub(".*~", "", list[i])))),value)
                            
                            # character is whatever comes before the "~" in the filename
                                character <- as.character(rep(gsub("~.*$", "", gsub(".*/", "", list[i])),length(value)))
                            
                            # var1 is whatever is in between the "~" and the "-" in the filename
                                var1 <- as.character(rep(gsub("-.*$","", gsub(".*~", "", list[i])),length(value)))

                            # var2 is whatever comes after the "-" in the filename
                                var2 <- as.character(rep(substr(list[i],regexpr("-",list[i])[1]+1,nchar(list[i])-4),length(value)))
                            
                            data <- rbind(data,cbind(ret[1:length(value)],value, character, var1, var2))
                        }

                        data <- data[-1,]
                        colnames(data) <- c("ret", "value", "character", "var1", "var2")
                        rownames(data) <- NULL
                        data$ret <- as.numeric(as.character(data$ret))
                        data$value <- as.numeric(as.character(data$value))

                ## Remove entries with var2 = NA
                    data <- data[!is.na(data$ret),]

                ## Normalize by var1
                    if (normalize_level[1] == "var1") {
                        
                        if (length(normalize_range) > 0) {
                            xmin <- normalize_range[1]
                            xmax <- normalize_range[2]
                            data <- data[data$ret >= xmin & data$ret <= xmax,]
                        }

                        data$value_norm <- data$value
                        for (i in 1:length(unique(data$var1))) {
                            data$value_norm[data$var1==unique(data$var1)[i]] <- phylochemistry::normalize(
                                data$value[data$var1==unique(data$var1)[i]],
                                old_min = min(as.numeric(data$value[data$var1==unique(data$var1)[i]])),
                                old_max = max(as.numeric(data$value[data$var1==unique(data$var1)[i]])),
                                new_min = 0, new_max = 100)
                        }
                    }

                ## Normalize by var2
                    if (normalize_level[1] == "var2") {
                        
                        if (length(normalize_range) > 0) {
                            xmin <- normalize_range[1]
                            xmax <- normalize_range[2]
                            data <- data[data$ret >= xmin & data$ret <= xmax,]
                        }

                        data$value_norm <- data$value
                        for (i in 1:length(unique(data$var2))) {
                            data$value_norm[data$var2==unique(data$var2)[i]] <- phylochemistry::normalize(
                                data$value[data$var2==unique(data$var2)[i]],
                                old_min = min(as.numeric(data$value[data$var2==unique(data$var2)[i]])),
                                old_max = max(as.numeric(data$value[data$var2==unique(data$var2)[i]])),
                                new_min = 0, new_max = 100)
                        }
                    }

                return(data)
            }

        #### readSpectra

            #' Import mass spectral data
            #'
            #' Used to import one or more mass spectra and collect the data in a dataframe.
            #' The names of the files containing the mass spectra must have the following format: "var2~var3-var4.csv", with var1 being the file's position in the directory
            #' @param dir The directory containing the mass spectral (.csv) files of interest
            #' @export
            #' @examples
            #' readSpectra()

            readSpectra <- function( dir, source = c("phylochemistry", "agilent") ) {

                ## Set working directory and prepare the data frame
                    setwd(dir)
                    data <- data.frame(
                        mz = as.numeric(1), abu = as.numeric(1), var1 = as.character("1"),
                        var2 = as.character("1"), var3 = as.character("1"), var4 = as.character("1")
                    )

                ## Read in all data
                    for (i in 1:length(dir())) {
                        temp <- readMonolist(dir()[i])
                        if (source == "phylochemistry" ) {
                            mz <- as.numeric(as.character(temp$mz))
                            abu <- 100*as.numeric(as.character(temp$abu))/(max(as.numeric(as.character(temp$abu))))
                        }
                        if (source == "agilent") {
                            temp <- temp[3:dim(temp)[1],1:2]
                            mz <- as.numeric(as.character(temp[,1]))
                            abu <- 100*as.numeric(as.character(temp[,2]))/(max(as.numeric(as.character(temp[,2]))))    
                        }
                        var1 <- as.character(rep(as.character(i))) # Var1 is the position in the directory
                        var2 <- as.character(rep(gsub("~.*$", "", dir()[i])))
                        var3 <- as.character(rep(gsub("-.*$", "", gsub(".*~", "", dir()[i]))))
                        var4 <- as.character(rep(gsub("\\..*$", "", gsub(".*-", "", dir()[i]))))
                        data <- rbind(data,cbind(mz,abu,var1,var2,var3,var4))
                    }

                ## Clean up the data frame and correct variable types
                    data <- data[-1,]
                    data$mz <- as.numeric(data$mz)
                    data$abu <- as.numeric(data$abu)
                    data$var1 <- factor(data$var1, levels=unique(data$var1))
                    data$var2 <- factor(data$var2, levels=unique(data$var2))
                    data$var3 <- factor(data$var3, levels=unique(data$var3))
                    data$var4 <- factor(data$var4, levels=unique(data$var4))

                ## Return data
                    return( data )
            }

        #### analyzeMassSpectralImages

            #' Use image analysis to obtain csv of mass spectrum from an image.
            #'
            #' Read a png with left-most pixel low mass, right-most pixel high mass
            #' Filename: reference~compoundName_lowmz-rightmz
            #' @param image_directory_in_path The directory containing the images to be read.
            #' @export
            #' @examples
            #' analyzeMassSpectralImages()

            analyzeMassSpectralImages <- function(image_directory_in_path) {

                images_to_analyze <- dir(image_directory_in_path)[grep(".png", dir(image_directory_in_path))]
                
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

                        image <- imager::load.image(paste0(image_directory_in_path, "/", images_to_analyze[image_to_analyze]))
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
                    write_csv(processed_image, paste0(image_directory_in_path, "/", gsub(".png", "", images_to_analyze[image_to_analyze]), ".csv"))

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
                    ggsave(filename = paste0(image_directory_in_path, "/", gsub(".png", "", images_to_analyze[image_to_analyze]), "_comparison.png"), plot = plot_output)

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
                        write_csv(merged_spectral_output, paste0(image_directory_in_path, "/_merged_spectral_output.csv"))
            }

        #### crystallinity

            #' Extract 720 and 720 wavenumber intensities from FTIR data
            #'
            #' Filenames: 123~123-1.csv
            #' To use plots: library(gridExtra); nCol <- floor(sqrt(length(plots))); do.call("grid.arrange", c(plots, ncol=nCol))
            #' @param data_directory_in_path The directory containing the data to be read and processed.
            #' @export
            #' @examples
            #' crystallinity()

            crystallinity <- function(data_directory_in_path, make_plots = TRUE) {

                data_directory_in_path <- OsDirectoryPathCorrect(data_directory_in_path)

                files <- dir(data_directory_in_path)
                out <- list()

                if (make_plots == TRUE) {
                    plots <- list()
                }
                
                for (i in 1:length(files)) {

                    cat(paste0("File ", files[i], "\n"))
                    spectrum <- read_csv(paste0(data_directory_in_path, files[[i]]), col_names = FALSE, col_types = cols())
                    colnames(spectrum) <- c("x", "y")

                    spectrum <- dplyr::filter(spectrum, x < 800 & x > 600)

                    ## Baseline the signal using windows. Doesn't work as well as window-defined line below
                        # baseline <- findBaseline(
                        #   x = spectrum$x,
                        #   y = spectrum$y,
                        #   window_width = 80
                        # ) 
                        # spectrum <- spectrum[spectrum$x %in% baseline$x,]
                        # spectrum$b <- baseline$y

                    ## Baseline the signal using a window-defined line
                        baseline_guide <- rbind(
                            dplyr::filter(spectrum, x < 760 & x > 740),
                            dplyr::filter(spectrum, x < 710 & x > 680)
                        )
                        spectrum$b <- lm(y~x, data = baseline_guide)$coefficients[2]*spectrum$x + lm(y~x, data = baseline_guide)$coefficients[1]
                    
                    ## Subtract the baseline, smooth, and plot the signal
                    ## If you dont filter, you will get extra little crappy peaks when you try to model!!
                        spectrum$y <- signal::sgolayfilt(spectrum$y, p = 1, n = 11)
                        spectrum$yb <- spectrum$y - spectrum$b

                    ## filter and just find maxima
                        peak1 <- max(dplyr::filter(spectrum, x > 728 & x < 740)$yb)
                        peak2 <- max(dplyr::filter(spectrum, x > 715 & x < 722)$yb)
                        wave1 <- dplyr::filter(spectrum, x > 728 & x < 740)[dplyr::filter(spectrum, x > 728 & x < 740)$yb == peak1,]$x
                        wave2 <- dplyr::filter(spectrum, x > 715 & x < 722)[dplyr::filter(spectrum, x > 715 & x < 722)$yb == peak2,]$x
                        ratio <- peak1/peak2
                        print(ratio)

                    ## Plot to inspect

                        if (make_plots == TRUE) {
                            plots[[i]] <- ggplot(data = spectrum) +
                                # geom_line(aes(x = x, y = y)) +
                                # geom_point(aes(x = x, y = y)) +
                                # geom_line(aes(x = x, y = b)) +
                                geom_vline(xintercept = wave1) +
                                geom_vline(xintercept = wave2) +
                                scale_x_continuous(breaks = seq(0,1000,10)) +
                                geom_line(aes(x = x, y = yb)) +
                                coord_cartesian(xlim = c(800, 600))
                        }

                    ## Filter down and fit gaussians. Just use heights. It's more reliable
                        # spectrum <- dplyr::filter(spectrum, x > 705 & x < 745)
                        # gaussians <- fitGaussians(
                        #   x = spectrum$x, 
                        #   y = spectrum$yb,
                        #   peak_detection_threshold = as.numeric(gsub(".CSV", "", gsub(".*~", "", files[i])))
                        # )

                    ## Get areas
                        # gaussians$file <- files[i]
                        # gaussians %>%
                        #   as_tibble() %>%
                        #   group_by(peak_number) %>%
                        #   summarize(
                        #       wavenumber = x[which.max(y)],
                        #       area = unique(peak_area),
                        #       file = unique(file),
                        #       height = max(y)
                        #   ) -> out[[i]]

                    cat("\n")
                    # readline(prompt="Press [enter] to continue")

                    out[[i]] <- 
                        data.frame(
                            sample = files[i],
                            ratio = ratio
                        )
                }

                out <- do.call(rbind, out)

                out$rep <- gsub(".CSV", "", gsub(".*-", "", out$sample))
                out$name <- gsub("-.*$", "", gsub("", "", out$sample))

                out %>%
                    group_by(name) %>%
                    dplyr::summarize(mean_ratio = mean(ratio)) -> ratio_order

                out$name <- factor(
                    out$name, 
                    levels = ratio_order$name[order(ratio_order$mean_ratio, decreasing = TRUE)]
                )

                return(out)

            }

    ##### Text data handling

        #### searchField

            searchField <- function(subject_list, query_list) {

                ## Make query list unique, remove all commas from the subjects so output can be written as csv

                    query_list <- unique(query_list)

                    subject_list <- gsub(",", "", subject_list)

                ## Start the output dataframe
                
                    output <- data.frame(
                        subject = subject_list,
                        hit_indeces = NA,
                        hits = NA
                    )

                    pb <- progress::progress_bar$new(
                        total = length(query_list),
                        format = "  analyzing [:bar] :percent eta: :eta"
                    )

                ## Loop over each query and record what subjects it's found in
                    
                    for( i in 1:length(query_list) ) {
                        
                        ## Identify where hits are 
                            
                            positive_subject_index <- grep(query_list[i], subject_list)
                        
                        ## If this query results in hits, start or add them to the list

                            if( length(positive_subject_index) > 0 ) {
                                
                                ## If it's the first entry in the list, start the list.

                                    if( any(is.na(output$hit_indeces[positive_subject_index])) ) {
                                        output$hit_indeces[positive_subject_index] <- i
                                        output$hits[positive_subject_index] <- as.character(query_list[i])
                                
                                ## If it's not the first entry in the list, append to the list.

                                    } else {

                                        # print(i)

                                        ## If adding to the list, run unique on it to remove duplicate names
                                            
                                            output$hit_indeces[positive_subject_index] <- 
                                                paste(
                                                    unique(
                                                        stringr::str_split(
                                                            paste(c(output$hit_indeces[positive_subject_index], i), collapse = ","),
                                                            pattern = ","
                                                        )[[1]]
                                                    ), collapse = ","
                                                )
                                            
                                            output$hits[positive_subject_index] <- 
                                                paste(
                                                    unique(
                                                        stringr::str_split(
                                                            paste(c(output$hits[positive_subject_index], as.character(query_list[i])), collapse = ","),
                                                            pattern = ","
                                                        )[[1]]
                                                    ), collapse = ","
                                                )
                                    }
                            }

                        ## Move progress bar along
                            
                            pb$tick()
                    }

                return(as_tibble(output))
            }

    ##### Networks

        #### buildNetwork

            #' Build a network from an edgelist
            #'
            #' @param infile The multiple alignment to (analyze and) plot.
            #' @param roi Plot a subset of the alignment. 
            #' @param consensus Include in the output graphic a line plot corresponding to the degree of conservation at each site in the alignment.
            #' @param funct_assoc Include in the output graphic a line plot of the association between AA residues and 
            #' @importFrom phangorn read.aa
            #' @import tidyr
            #' @import ggplot2
            #' @import ggtree  
            #' @export
            #' @examples
            #' buildNetwork()

            buildNetwork <- function( edgelist, node_attributes = NULL, facet_variable = NULL ) {

                ## Make sure the edgelist is a dataframe and that its first two columns are characters
                    edgelist <- as.data.frame(edgelist)
                    edgelist[,1] <- as.character(edgelist[,1])
                    edgelist[,2] <- as.character(edgelist[,2])

                ## If edgelist has more than two columns, assume that 3rd column is edge weights and build matrix from adjacency matrix
                    # if (dim(edgelist)[2] > 2) {
                    #     adjacency_matrix <- as.matrix(as.data.frame(tidyr::pivot_wider(edgelist, names_from = 1, values_from = 3)))
                    #     rownames(adjacency_matrix) <- adjacency_matrix[,1]
                    #     adjacency_matrix <- adjacency_matrix[,-1]
                    #     network_object <- network::network(adjacency_matrix, matrix.type = "adjacency")
                    # }

                ## Just build a network with the edgelist
                    if (dim(edgelist)[2] == 3) {
                        network_object <- network::network(edgelist, matrix.type = "edgelist", ignore.eval = FALSE)
                    }

                ## Fortify the network object into a dataframe, modify it with new colnames and additional information, crack into a list of edges and nodes
                    combined_network_frame <- ggnetwork::ggnetwork(network_object, arrow.gap = 0, layout = "kamadakawai")
                    # combined_network_frame <- ggnetwork::ggnetwork(network_object, arrow.gap = 0)
                    combined_network_frame$type <- "edge"
                    colnames(combined_network_frame)[colnames(combined_network_frame) == "vertex.names"] <- "start_node"
                    combined_network_frame$type[apply(cbind(combined_network_frame$x == combined_network_frame$xend, combined_network_frame$y == combined_network_frame$yend), 1, all)] <- "node"
                    combined_network_frame$end_node <- combined_network_frame$start_node[match(combined_network_frame$xend, combined_network_frame$x)]

                    network_frame <- list()
                    network_frame$edges <- combined_network_frame[combined_network_frame$type == "edge",]
                    network_frame$edges <- network_frame$edges[,colnames(network_frame$edges) %in% c("x", "y", "start_node", "xend", "yend", "end_node")]
                    network_frame$nodes <- combined_network_frame[combined_network_frame$type == "node",]
                    network_frame$nodes <- network_frame$nodes[,colnames(network_frame$nodes) %in% c("x", "y", "start_node")]
                    colnames(network_frame$nodes)[3] <- "node_name"

                    combined_network_frame$end_node <- as.character(combined_network_frame$end_node)

                ## Bind edge attributes to the network_frame
                    if (dim(edgelist)[2] > 2) {
                        edge_attributes <- list()
                        for ( edge in 1:dim(edgelist)[1] ) {
                            network_frame_entry <- which(apply(cbind(combined_network_frame$start_node == edgelist[edge,1], combined_network_frame$end_node == edgelist[edge,2]), 1, all))
                            edge_attributes[[network_frame_entry]] <- as.character(edgelist[edge, 3:dim(edgelist)[2]])
                        }
                        edge_attributes <- as.data.frame(do.call(rbind, edge_attributes))
                        colnames(edge_attributes) <- colnames(edgelist)[3:dim(edgelist)[2]]
                        network_frame$edges <- cbind(network_frame$edges, edge_attributes)
                    }

                ## Merge node_attributes frame with network_frame$nodes
                    if (!is.null(node_attributes)) {
                        network_frame$nodes <- cbind(network_frame$nodes, node_attributes[match(network_frame$nodes$node_name, node_attributes[,1]),])
                    }

                ## Create network_frame_expanded, which contains entries for each facet_variable
                    # network_frame_expanded <- data.frame()
                    # facet_variables <- as.character(unique(annotation_frame[,which(colnames(annotation_frame) == facet_variable)]))

                    # if ( length(facet_variables) > 0 ) {
                    #     for (i in 1:length(facet_variables)) {
                    #         frame_subset <- annotation_frame[annotation_frame[,which(colnames(annotation_frame) == facet_variable)] == facet_variables[i],]
                    #         frame_subset <- frame_subset[frame_subset$abundance > 0,]
                            
                    #         # Show only nodes present in that facet_variable
                    #             # nodes <- network_frame[apply(cbind(network_frame$type == "node", network_frame$start_node %in% frame_subset$compound_name), 1, all),]

                    #         # Show all nodes
                    #             nodes <- network_frame[network_frame$type == "node",]
                            
                    #         # Show only edges present in that facet variable
                    #             edges <- network_frame[apply(   cbind(
                    #                                                 network_frame$type == "edge", 
                    #                                                 network_frame$start_node %in% frame_subset$compound_name,
                    #                                                 network_frame$end_node %in% frame_subset$compound_name
                    #                                             ), 1, all),]

                    #         temp <- rbind(nodes, edges)
                    #         temp$facet_variable <- facet_variables[i]
                    #         colnames(temp)[colnames(temp) == "facet_variable"] <- facet_variable

                    #         network_frame_expanded <- rbind(network_frame_expanded, temp)
                    #     }
                    # } else {
                    #     network_frame_expanded <- network_frame
                    # }

                return(network_frame)
            }

    ##### Data Visualization
        
        #### drawAlignment

            #' Generate a multiple alignment graphic using ggplot
            #'
            #' Allows the user to plot a multiple alignent using ggplot's grammar of graphics
            #' @param infile The alignment to use (fasta file)
            #' @param monolist A dataframe with metadata for the alignment
            #' @param alignment_labels The name of the column to label the entries in the alignment with
            #' @param wrap TRUE/FALSE whether to use a column of plots to show the alignment
            #' @param wrap_length Length at which to cut the wrap
            #' @param roi TRUE/FALSE whether to only show regions of interest in the alignment
            #' @param roi_data Dataframe defining the regions to show
            #' @param consensus TRUE/FALSE whether to include a consensus line in the alignment
            #' @param consensus_height Height of the consensus line plot
            #' @param funct_assoc TRUE/FALSE whether to show functional associations in the alignment
            #' @param funct_assoc_data Columns in the metadata dataframe to use in searching for functional association (e.g. "c(sterol~other, sterol~cyclo)")
            #' @param funct_assoc_height Height of the functional association lines
            #' @param highlights TRUE/FALSE whether to highlight certain sites in the alignment
            #' @param highlights_data Dataframe specifying which site to highlight
            #' @param hlines TRUE/FALSE whether to divide the alignment with horizontal lines
            #' @param hlines_data Dataframe specifying where to divide alignment
            #' @param tick_spacing Spacing between x-axis ticks
            #' @param ticks_text_size Size of x-axis ticks text
            #' @param order Order in which to display the members of the alignment
            #' @param color_pal Color palette to use
            #' @importFrom phangorn read.aa
            #' @import tidyr
            #' @import ggplot2
            #' @import ggtree  
            #' @export
            #' @examples
            #' drawAlignment()

            drawAlignment <- function(   
                                infile,
                                monolist,
                                alignment_labels,

                                wrap = FALSE,
                                wrap_length = NULL,

                                roi = FALSE,
                                roi_data = NULL,

                                consensus = FALSE,
                                consensus_height = 5,

                                funct_assoc = FALSE,
                                funct_assoc_data = NULL,
                                funct_assoc_height = 5,

                                highlights = FALSE,
                                highlights_data = NULL,

                                hlines = FALSE,
                                hlines_data = NULL,

                                tick_spacing = 5,
                                ticks_text_size = 30,

                                order = NULL,
                                color_pal = NULL
                            ) {

                if (wrap == FALSE) {
                    if (roi == FALSE) {
                        print("Please specify either wrap OR roi")
                        stop()
                    }
                }

                if (wrap == TRUE) {
                    if (roi == TRUE) {
                        print("Please specify either wrap OR roi")
                        stop()
                    }
                }

                # Read alignment, create basic plottable
                    AA_phydat <- as.data.frame(as.character(phangorn::read.aa(file = infile, format = "fasta")))
                    AA_phydat <- as.data.frame(t(AA_phydat))
                    colnames(AA_phydat) <- as.character(seq(1,dim(AA_phydat)[2],1))
                    AA_alignment_df <-  cbind(  
                                            data.frame(protein = paste(
                                                monolist[,which(colnames(monolist) == alignment_labels)][match(rownames(as.matrix(AA_phydat)), monolist$accession)]
                                            )),
                                            # data.frame(funct = monolist[,colnames(monolist)=="Function"][match(rownames(as.matrix(AA_phydat)), monolist[,1])]),
                                            data.frame(y = rep(0,dim(AA_phydat)[1])),
                                            AA_phydat
                                        )
                    # AA_alignment_df_plottable <- tidyr::gather(AA_alignment_df
                    first_number_col <- min(which(is.na(suppressWarnings(as.numeric(colnames(AA_alignment_df)))) == FALSE))
                    AA_alignment_df_plottable <- as.data.frame(tidyr::pivot_longer(AA_alignment_df, cols = first_number_col:dim(AA_alignment_df)[2], names_to = "position", values_to = "residue"))
                    AA_alignment_df_plottable$position <- as.numeric(as.character(AA_alignment_df_plottable$position))
                    AA_alignment_df_plottable <- AA_alignment_df_plottable[order(AA_alignment_df_plottable$protein, AA_alignment_df_plottable$position),]
                    head(AA_alignment_df_plottable)

                # Make wrapping list
                    if (wrap == TRUE) {
                        wrap_break_points <- list()
                        for (i in 1:(dim(AA_phydat)[2]%/%wrap_length)) {
                            wrap_break_points <- c(wrap_break_points, 1+(wrap_length*(i)))
                        }
                        head(wrap_break_points)
                    }

                    if (roi == TRUE) {
                        wrap_break_points <- dim(AA_phydat)[2]
                        head(wrap_break_points)
                        wrap_length = dim(AA_phydat)[2]
                    }

                # ROI
                    if (roi == TRUE) {
                        #Define the ROI and Modify the alignment plottable
                            roi_ranges <- rep("NA",dim(AA_phydat)[2])
                                range_names <- roi_data[,1]
                                for (i in 1:dim(roi_data)[1]) {
                                    roi_ranges[as.numeric(as.character(roi_data[i,2])):as.numeric(as.character(roi_data[i,3]))] <- as.character(roi_data[i,1])
                                }
                            AA_alignment_df_plottable$roi_ranges <- rep(roi_ranges, length(unique(AA_alignment_df_plottable$protein)))
                            AA_alignment_df_plottable <- subset(AA_alignment_df_plottable, roi_ranges %in% as.character(range_names))
                    }

                # ORDER the sequences to optionally match the alignment with a tree or something
                    if (!is.null(order)) {
                        AA_alignment_df_plottable$protein <- factor(AA_alignment_df_plottable$protein, levels=order)
                    }
                    
                # Initiate the plot(s)
                    plot_list <- list()
                    for (i in 1:length(wrap_break_points)) {
                        plot_list[[i]] <-   ggplot2::ggplot(data = AA_alignment_df_plottable[AA_alignment_df_plottable$position %in% 
                                                            as.numeric(as.character(unlist(wrap_break_points)[i]-wrap_length)):
                                                            as.numeric(as.character(unlist(wrap_break_points)[i])),],
                                                        aes(x = position, y = y)
                                            )   
                    }

                # HIGHLIGHTS
                    if (highlights == TRUE) {
                        if (roi == TRUE) {
                            highlights_data <- cbind(highlights_data, data.frame(roi_ranges=roi_ranges[highlights_data$xint]))
                        }
                        for (i in 1:length(wrap_break_points)) {
                            plot_list[[i]] <- plot_list[[i]] +  geom_vline(
                                                                    data=highlights_data[highlights_data$xint %in% 
                                                                        as.numeric(as.character(unlist(wrap_break_points)[i]-wrap_length)):
                                                                        as.numeric(as.character(unlist(wrap_break_points)[i])),],
                                                                    aes(xintercept=xint, color=as.character(xint)),
                                                                    size=3
                                                                )
                        }
                    }
                    
                # CONSENSUS
                    if (consensus == TRUE) {
                        # Calculate consensus scores
                            #Calculate max possible score:
                                AA_phydat_char_distrib <- t(as.data.frame(lapply(apply(AA_phydat,2,table), sd)))
                                AA_phydat_char_distrib[is.na(AA_phydat_char_distrib)] <- 0
                                max_consensus_score <- ceiling(max(AA_phydat_char_distrib))

                            # Calculate consensus scores:
                                AA_phydat_char_distrib <- t(as.data.frame(lapply(apply(AA_phydat,2,table), sd)))
                                AA_phydat_char_distrib[is.na(AA_phydat_char_distrib)] <- max_consensus_score
                                colnames(AA_phydat_char_distrib) <- "consensus_score"
                                AA_phydat_char_distrib <- normalize(AA_phydat_char_distrib)*consensus_height
                                
                        # Bind consensus score
                            consensus_protein <- data.frame(
                                protein = "consensus",
                                funct="none",
                                y=AA_phydat_char_distrib[,1],
                                position=seq(1,dim(AA_phydat)[2]),
                                residue="+"
                            )
                            head(consensus_protein)

                        # Truncate the consensus_protein to fit roi
                            if (roi == TRUE) {
                                consensus_protein <- cbind(
                                    consensus_protein, 
                                    data.frame(roi_ranges=  rep(    roi_ranges,
                                                                    length(unique(consensus_protein$protein))
                                                            )
                                    )
                                )
                                consensus_protein <- subset(consensus_protein, roi_ranges %in% range_names)
                            }

                        # Add consensus_protein to the plot(s) in plot_list
                            for (i in 1:length(wrap_break_points)) {
                                plot_list[[i]] <-  plot_list[[i]] + layer(  data=consensus_protein[consensus_protein$position %in% 
                                                                                    as.numeric(as.character(unlist(wrap_break_points)[i]-wrap_length)):
                                                                                    as.numeric(as.character(unlist(wrap_break_points)[i])),],
                                                                                aes(x=position, y=y), geom="line", stat="identity", position="identity"
                                                                            )
                            }
                    }

                # FUNCTIONAL ASSOCIATION
                    if (funct_assoc == TRUE) {
                        for (j in 1:length(funct_assoc_data)) {
                            # Subset alignment matrices for the two variables
                                alignment_matrix <- toupper(t(as.data.frame(as.character(phangorn::read.aa(file = infile, format = "fasta")))))
                                
                                spec1 <- list()
                                spec1[[1]] <- monolist[,colnames(monolist) == funct_assoc_data[j]] == as.character(gsub("~.*$", "", funct_assoc_data[j]))
                                alignment_matrix_dim1 <- alignment_matrix[Reduce("&", spec1),]
                                dim(alignment_matrix_dim1)

                                spec2 <- list()
                                spec2[[1]] <- monolist[,colnames(monolist)==funct_assoc_data[j]] == as.character(gsub(".*~", "", funct_assoc_data[j]))
                                alignment_matrix_dim2 <- alignment_matrix[Reduce("&", spec2),]
                                dim(alignment_matrix_dim2)

                            #Set up fixed substitution table
                                fixed_sub_table <- data.frame(funct = c(as.character(gsub("~.*$", "", funct_assoc_data[j])), as.character(gsub(".*~", "", funct_assoc_data[j]))), code=c(1,2))

                            #Get frequency stats from the first alignment subset
                                alignment_matrix_dim1_char_distrib <- lapply(apply(alignment_matrix_dim1,2,table), sort, decreasing=TRUE)
                                alignment_matrix_dim2_char_distrib <- lapply(apply(alignment_matrix_dim2,2,table), sort, decreasing=TRUE)

                            #Set up empty assoc matrix
                                association_list <- list()

                                for (i in 1:dim(alignment_matrix)[2]){

                                    # Set up dim1 residue substitution table
                                        base_sub_table_dim1 <- data.frame(letter=LETTERS, number=rep(2,26))
                                        base_sub_table_dim1 <- rbind(base_sub_table_dim1, data.frame(letter="-", number=2))
                                        sub_position_dim1 <- data.frame(consensus=names(alignment_matrix_dim1_char_distrib[[i]][1]), number=1)
                                        positional_sub_table_dim1 <- base_sub_table_dim1
                                        positional_sub_table_dim1$number <- sub_position_dim1$number[match(positional_sub_table_dim1$letter, sub_position_dim1$consensus)]
                                        positional_sub_table_dim1[is.na(positional_sub_table_dim1$number),]$number <- 2
                                        positional_sub_table_dim1

                                    # Substitute in the first alignment subset
                                        positional_data_frame_dim1 <- data.frame(
                                            enzyme=names(alignment_matrix_dim1[,i]),
                                            funct=as.character(gsub("~.*$", "", funct_assoc_data[j])),
                                            x_num=as.character(gsub("~.*$", "", funct_assoc_data[j])),
                                            position=i,
                                            residue=alignment_matrix_dim1[,i],
                                            y_num=alignment_matrix_dim1[,i]
                                        )

                                        positional_data_frame_dim1$x_num <- fixed_sub_table$code[match(positional_data_frame_dim1$x_num, fixed_sub_table$funct)]
                                        positional_data_frame_dim1$y_num <- positional_sub_table_dim1$number[match(positional_data_frame_dim1$y_num, positional_sub_table_dim1$letter)]
                                        positional_data_frame_dim1

                                    # Set up dim2 residue substitution table
                                        base_sub_table_dim2 <- data.frame(letter=LETTERS, number=rep(1,26))
                                        base_sub_table_dim2 <- rbind(base_sub_table_dim2, data.frame(letter="-", number=1))
                                        sub_position_dim2 <- data.frame(consensus=names(alignment_matrix_dim2_char_distrib[[i]][1]), number=2)
                                        if (as.character(sub_position_dim2$consensus) == as.character(sub_position_dim1$consensus)) {sub_position_dim2$number <- 1}
                                        positional_sub_table_dim2 <- base_sub_table_dim2
                                        positional_sub_table_dim2$number <- sub_position_dim2$number[match(positional_sub_table_dim2$letter, sub_position_dim2$consensus)]
                                        positional_sub_table_dim2[is.na(positional_sub_table_dim2$number),]$number <- 1
                                        positional_sub_table_dim2

                                    #Substitute in the second alignment subset
                                        positional_data_frame_dim2 <- data.frame(
                                            enzyme=names(alignment_matrix_dim2[,i]),
                                            funct=as.character(gsub(".*~", "", funct_assoc_data[j])),
                                            x_num=as.character(gsub(".*~", "", funct_assoc_data[j])),
                                            position=i,
                                            residue=alignment_matrix_dim2[,i],
                                            y_num=alignment_matrix_dim2[,i]
                                        )

                                        positional_data_frame_dim2$x_num <- fixed_sub_table$code[match(positional_data_frame_dim2$x_num, fixed_sub_table$funct)]
                                        positional_data_frame_dim2$y_num <- positional_sub_table_dim2$number[match(positional_data_frame_dim2$y_num, positional_sub_table_dim2$letter)]

                                    #Add the positional_data_frame to the list
                                        association_list[[i]] <- rbind(positional_data_frame_dim1, positional_data_frame_dim2)
                                }

                                # USE THIS LINE TO MANUALLY INSPECT SITE ASSOCIATAIONS
                                    # association_list[[273]]

                                cor_list <- list()
                                for (i in 1:length(association_list)) {
                                    cor_list[i] <- cor(cbind(association_list[[i]]$x_num,association_list[[i]]$y_num))[1,2]
                                }

                                cor_list[is.na(cor_list)] <- 0
                                cor_list[cor_list < 0] <- 0
                                correlation_protein <-  data.frame(
                                                            protein = as.character(funct_assoc_data[j]), 
                                                            funct = "NA", 
                                                            y = t(as.data.frame(cor_list))*funct_assoc_height, 
                                                            position = seq(1,length(cor_list),1), 
                                                            residue = "NA"
                                                        )
                                head(correlation_protein)

                                # Truncate the correlation_protein plottable so it only contains things in the roi
                                    if (roi == TRUE) {
                                        correlation_protein <- cbind(
                                            correlation_protein, 
                                            data.frame(roi_ranges=rep(
                                                roi_ranges,
                                                length(unique(correlation_protein$protein)))
                                            )
                                        )
                                        correlation_protein <- subset(correlation_protein, roi_ranges %in% range_names)
                                    }

                                # Assign correlation protein to global environment
                                    assign("correlation_protein", correlation_protein, envir = .GlobalEnv)

                            # Add the correlation_protein to the plot(s) in plot_list
                                for (i in 1:length(wrap_break_points)) {
                                    plot_list[[i]] <-  plot_list[[i]] + layer(data=correlation_protein[correlation_protein$position %in% 
                                                                                as.numeric(as.character(unlist(wrap_break_points)[i]-wrap_length)):
                                                                                as.numeric(as.character(unlist(wrap_break_points)[i])),],
                                                                            aes(x=position, y=y), geom="line", stat="identity", position="identity"
                                                                        )
                                }
                        }
                    }

                # HLINES
                    if (hlines == TRUE) {
                        for (i in 1:length(wrap_break_points)) {
                            plot_list[[i]] <- plot_list[[i]] +  geom_hline(
                                                                    data = hlines_data,
                                                                    aes(yintercept = yint),
                                                                    size = 1,
                                                                    linetype = 2
                                                                )
                        }
                    }

                # Build the rest of the alignment
                    #Specify facets based on presence/absence of ROIs
                        if (roi == TRUE) { 
                            for (i in 1:length(wrap_break_points)) {
                                plot_list[[i]] <- plot_list[[i]] + facet_grid(protein~roi_ranges, scales="free", space="free", switch="y") 
                            }
                        }
                        if (wrap == TRUE) { 
                            for (i in 1:length(wrap_break_points)) {
                                plot_list[[i]] <- plot_list[[i]] + facet_grid(protein~., scales="free", space="free", switch="y") 
                            }
                        }
                        
                    # Make rest of plot
                        for (i in 1:length(wrap_break_points)) {
                            plot_list[[i]] <- plot_list[[i]] + theme_classic()
                            plot_list[[i]] <- plot_list[[i]] +  scale_x_continuous(
                                                                    expand = c(0.005,0.005),
                                                                    name = ""
                                                                    # breaks = as.numeric(generateTicks(seq(0,5000,tick_spacing))$all_breaks),
                                                                    # labels = as.character(generateTicks(seq(0,5000,tick_spacing))$all_labels)
                                                                )
                            plot_list[[i]] <- plot_list[[i]] + scale_y_continuous(name = "", breaks = 0)
                            plot_list[[i]] <- plot_list[[i]] + scale_color_manual(values = color_pal)
                            plot_list[[i]] <- plot_list[[i]] + theme(
                                        panel.spacing.y = unit(0.01, "lines"),
                                        panel.spacing.x = unit(0.2, "cm"),
                                        panel.border = element_blank(),
                                        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                                        axis.text.y = element_blank(),
                                        text = element_text(size = ticks_text_size),
                                        legend.position = "none",
                                        axis.title = element_text(color = "#737373", face = "bold", size = 40),
                                        axis.ticks.length = unit(0.2, "cm"),
                                        axis.ticks = element_line(color = "#737373", size = 1, lineend = 6),
                                        axis.text = element_text(color = "#737373", face = "bold", size = ticks_text_size),
                                        axis.line = element_line(color = "#737373", size = 1),
                                        strip.text = element_text(hjust = 1, size = ticks_text_size),
                                        strip.text.y.left = element_text(angle = 0),
                                        strip.background = element_blank(),
                                        strip.placement = "outside"
                                    )
                        }

                # Return the Alignment
                    if (wrap == TRUE) { ## IF USING WRAP, RETURNS A LIST THAT MUST BE PLOTTED USING do.call(gridExtra::grid.arrange,c(alignment, ncol=1))
                        plots <- list()
                        for (i in 1:length(wrap_break_points)) {
                            plots[[i]] <- ggplot2::ggplot_gtable(
                                ggplot2::ggplot_build(plot_list[[i]] + layer(
                                    data = AA_alignment_df_plottable[AA_alignment_df_plottable$position %in% 
                                        as.numeric(as.character(unlist(wrap_break_points)[i]-wrap_length)):
                                        as.numeric(as.character(unlist(wrap_break_points)[i])),],
                                    geom = "text",
                                    mapping = aes(label = residue),
                                    stat = "identity",
                                    position = "identity")
                                )
                            )
                        }
                        return(plots)
                    }

                    if (roi == TRUE) { ## IF USING ROI, RETURNS A GGPLOT THAT CAN BE PLOTTED USING NORMAL GGPLOT PLOTTING TECHNIQUES
                        plot_list[[i]] <- plot_list[[i]] + layer(
                            data = AA_alignment_df_plottable[AA_alignment_df_plottable$position %in% 
                                as.numeric(as.character(unlist(wrap_break_points)[i]-wrap_length)):
                                as.numeric(as.character(unlist(wrap_break_points)[i])),],
                            geom = "text", 
                            mapping = aes(label = residue),
                            stat = "identity",
                            position = "identity")
                        return(plot_list[[1]])
                    }           
            }

        #### generateTicks

            #' Create major and minor axes ticks and labels
            #'
            #' @param major_ticks Values at which major ticks should be generated
            #' @param minor_freq Number of minor ticks between each major tick 
            #' @examples
            #' @export
            #' generateTicks

            generateTicks <-    function(
                                    major_ticks,
                                    minor_freq = 4,
                                    major_tick_size = 2,
                                    minor_tick_size = 1
                                ) {
                                    
                major_labels <- vector()
                
                for (tick in 1:(length(major_ticks)-1)) {
                    major_labels <- c(major_labels, major_ticks[tick])
                    major_labels <- c(major_labels, rep("", (minor_freq)))
                }
                
                all_ticks <- vector()
                
                for (tick in 1:(length(major_ticks)-1)) {
                    all_ticks <- c(all_ticks, major_ticks[tick])
                    minor_tick_values <- vector()
                    for (minor_tick in 1:minor_freq) {
                        minor_tick_values <- c(minor_tick_values, major_ticks[tick] + minor_tick*(((major_ticks[tick+1] - major_ticks[tick])/(minor_freq+1))))
                    }
                    all_ticks <- c(all_ticks, minor_tick_values)
                }
                
                major_labels <- c(major_labels, major_ticks[length(major_ticks)])
                all_ticks <- c(all_ticks, major_ticks[length(major_ticks)])

                return <- data.frame(all_ticks = as.numeric(all_ticks), major_labels = as.character(major_labels))
                return$tick_size <- minor_tick_size
                return$tick_size[return$major_labels != ""] <- major_tick_size

                return(return)
            }

        #### drawMolecules

            #' Draw chemical structures from CSV
            #'
            #' @param csv_in_path Path to the csv with the molecular coordinates.
            #' @param data Data that the functions should plot.
            #' @export
            #' @examples
            #' drawMolecules()

            drawMolecules <-    function(
                                    csv_in_path = NULL, 
                                    data = NULL, 
                                    atom_color_column = "default",
                                    bond_color_column = "default"
                                ) {

                ## Check input data / path

                    if (is.null(csv_in_path) & is.null(data)) {
                        stop("Please provide a data source using either the `data` argument OR the `csv_in_path` argument.")
                    }

                    if (!is.null(csv_in_path) & !is.null(data)) {
                        stop("Please provide a SINGLE data source using either the `data` argument OR the `csv_in_path` argument.")
                    }

                    if (is.null(csv_in_path) & !is.null(data)) {}

                    if (!is.null(csv_in_path) & is.null(data)) {
                        data <- read_csv(csv_in_path, col_types = cols())
                    }

                ## Clean up data
                    output <- list()
                    for (i in 1:length(unique(data$molecule_name))) {

                        temp <- filter(data, molecule_name == unique(data$molecule_name)[i])
                        temp$bond_start_x <- temp$x[match(temp$bond_start_atom, temp$atom_number)]
                        temp$bond_start_y <- temp$y[match(temp$bond_start_atom, temp$atom_number)]
                        temp$bond_end_x <- temp$x[match(temp$bond_end_atom, temp$atom_number)]
                        temp$bond_end_y <- temp$y[match(temp$bond_end_atom, temp$atom_number)]
                        temp$bond_start_element_or_group <- temp$element_or_group[match(temp$bond_start_atom, temp$atom_number)]
                        temp$bond_end_element_or_group <- temp$element_or_group[match(temp$bond_end_atom, temp$atom_number)]
                        output[[i]] <- temp
                    
                    }

                    plot_data <- do.call(rbind, output)

                ## Define color schemes

                    atom_colors_pre <- as.data.frame(rbind(
                        c("C", "black"),
                        c("CH2", "grey"),
                        c("CH2OH", "#4daf4a"), # green
                        c("CH3", "#ff7f00"), # orange
                        c("COOH", "#ffff33"), # yellow
                        c("H", "white"),
                        c("O", "#e41a1c"), # red
                        c("OH", "#377eb8"), # blue
                        c("CHO", "grey20"),
                        c("OCOCH3", "grey30"),
                        c("COOCH3", "grey40"),
                        c("CH2CH3", "grey50")
                    ))
                    atom_colors <- atom_colors_pre[,2]
                    names(atom_colors) <- atom_colors_pre[,1]

                    #984ea3 purple
                    #a65628 brown

                ## Plotting

                    plot <- ggplot()

                        ## Add bonds under default color scheme

                            if (bond_color_column == "default") {

                                ## Add achiral single bonds
                                    plot <- plot + geom_link(
                                        data = filter(plot_data, 
                                            molecule_component == "bond" &  bond_type == "single" & bond_direction == "flat"
                                        ),
                                        aes(
                                            x = bond_start_x, y = bond_start_y, xend = bond_end_x, yend = bond_end_y
                                        ), size = 2, color = "grey30"
                                    )

                                ## Add down chiral single bonds
                                    plot <- plot + geom_link(
                                        data = filter(plot_data, 
                                            molecule_component == "bond" & bond_type == "single" & bond_direction == "down"
                                        ),
                                        aes(
                                            x = bond_start_x, y = bond_start_y, xend = bond_end_x, yend = bond_end_y,
                                            size = stat(index)
                                        ), color = "grey60"
                                    )
                                
                                ## Add up chiral single bonds
                                    plot <- plot + geom_link(
                                        data = filter(plot_data, 
                                            molecule_component == "bond" & bond_type == "single" & bond_direction == "up"
                                        ),
                                        aes(
                                            x = bond_start_x, y = bond_start_y, xend = bond_end_x, yend = bond_end_y,
                                            size = stat(index)
                                        ), color = "grey0"
                                    )

                                ## Add steep non-horizontal double bonds
                                    plot <- plot + geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x != bond_end_x & bond_start_y != bond_end_y & abs(bond_start_y- bond_end_y) != 1),
                                        aes(x = bond_start_x-0.25, y = bond_start_y-0.15, xend = bond_end_x-0.25, yend = bond_end_y-0.15),
                                        size = 1.5, alpha = 0.7
                                    ) +
                                    geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x != bond_end_x & bond_start_y != bond_end_y & abs(bond_start_y- bond_end_y) != 1),
                                        aes(x = bond_start_x+0.25, y = bond_start_y+0.15, xend = bond_end_x+0.25, yend = bond_end_y+0.15),
                                        size = 1.5, alpha = 0.7
                                    )

                                ## Add nearly horizontal double bonds
                                    plot <- plot + geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x != bond_end_x & abs(bond_start_y- bond_end_y) == 1),
                                        aes(x = bond_start_x, y = bond_start_y-0.3, xend = bond_end_x, yend = bond_end_y-0.3),
                                        size = 1.5, alpha = 0.7
                                    ) +
                                    geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x != bond_end_x & abs(bond_start_y- bond_end_y) == 1),
                                        aes(x = bond_start_x, y = bond_start_y+0.3, xend = bond_end_x, yend = bond_end_y+0.3),
                                        size = 1.5, alpha = 0.7
                                    )

                                ## Add horizontal double bonds
                                    plot <- plot + geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_y == bond_end_y),
                                        aes(x = bond_start_x, y = bond_start_y-0.2, xend = bond_end_x, yend = bond_end_y-0.2),
                                        size = 1.5, alpha = 0.7
                                    ) +
                                    geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_y == bond_end_y),
                                        aes(x = bond_start_x, y = bond_start_y+0.2, xend = bond_end_x, yend = bond_end_y+0.2),
                                        size = 1.5, alpha = 0.7
                                    )

                                ## Add vertical double bonds
                                    plot <- plot + geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x == bond_end_x),
                                        aes(x = bond_start_x-0.25, y = bond_start_y, xend = bond_end_x-0.25, yend = bond_end_y),
                                        size = 1.5, alpha = 0.7
                                    ) +
                                    geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x == bond_end_x),
                                        aes(x = bond_start_x+0.25, y = bond_start_y, xend = bond_end_x+0.25, yend = bond_end_y),
                                        size = 1.5, alpha = 0.7
                                    )

                                ## Color scheme
                                    # plot <- plot + scale_color_gradient(low = "black", high = "black")
                           
                            }

                        ## Add bonds under custom color scheme

                            if (bond_color_column != "default") {

                                ## Add achiral single bonds
                                    plot <- plot + geom_link(
                                        data = filter(plot_data, 
                                            molecule_component == "bond" &  bond_type == "single" & bond_direction == "flat"
                                        ),
                                        aes_string(
                                            x = "bond_start_x", y = "bond_start_y", xend = "bond_end_x", yend = "bond_end_y",
                                            color = bond_color_column
                                        ), size = 2
                                    )

                                ## Add down chiral single bonds
                                    plot <- plot + geom_link(
                                        data = filter(plot_data, 
                                            molecule_component == "bond" & bond_type == "single" & bond_direction == "down"
                                        ),
                                        aes_string(
                                            x = "bond_start_x", y = "bond_start_y", xend = "bond_end_x", yend = "bond_end_y",
                                            size = "stat(index)", color = bond_color_column
                                        )
                                    )
                                
                                ## Add up chiral single bonds
                                    plot <- plot + geom_link(
                                        data = filter(plot_data, 
                                            molecule_component == "bond" & bond_type == "single" & bond_direction == "up"
                                        ),
                                        aes_string(
                                            x = "bond_start_x", y = "bond_start_y", xend = "bond_end_x", yend = "bond_end_y",
                                            size = "stat(index)", color = bond_color_column
                                        )
                                    )

                                ## Add nearly horizontal double bonds
                                    plot <- plot + geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x != bond_end_x),
                                        aes_string(x = "bond_start_x", y = "bond_start_y-0.2", xend = "bond_end_x", yend = "bond_end_y-0.2"),
                                        size = 1.5, alpha = 0.7, color = "blue"
                                    ) +
                                    geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x != bond_end_x),
                                        aes_string(x = "bond_start_x", y = "bond_start_y+0.2", xend = "bond_end_x", yend = "bond_end_y+0.2"),
                                        size = 1.5, alpha = 0.7, color = "blue"
                                    )

                                ## Add vertical double bonds
                                    plot <- plot + geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x == bond_end_x),
                                        aes_string(x = "bond_start_x-0.18", y = "bond_start_y", xend = "bond_end_x-0.18", yend = "bond_end_y"),
                                        size = 1.5, alpha = 0.7, color = "blue"
                                    ) +
                                    geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x == bond_end_x),
                                        aes_string(x = "bond_start_x+0.18", y = "bond_start_y", xend = "bond_end_x+0.18", yend = "bond_end_y"),
                                        size = 1.5, alpha = 0.7, color = "blue"
                                    )
                            }

                        
                        ## Add atom circles 
                          
                            if (atom_color_column == "default") {
                                plot <- plot + geom_point(
                                    data = filter(plot_data, molecule_component == "atom"),
                                    aes(x = x, y = y, fill = element_or_group),
                                    shape = 21, size = 4
                                ) +
                                scale_fill_manual(values = atom_colors, name = "")
                            }

                            if (atom_color_column != "default") {
                                plot <- plot + geom_point(
                                    data = filter(plot_data, molecule_component == "atom"),
                                    aes_string(x = "x", y = "y", fill = atom_color_column),
                                    shape = 21, size = 4
                                ) +
                                scale_fill_viridis(na.value = "#440154FF")
                            }

                        ## Add atom number labels
                          
                            plot <- plot + geom_text(
                                data = filter(plot_data, molecule_component == "atom", atom_number %in% c(10, 50, 71)),
                                aes(x = x, y = y, label = atom_number),
                                size = 2, color = "black"
                            ) +

                            geom_text(
                                data = filter(plot_data, molecule_component == "atom", !atom_number %in% c(10, 50, 71)),
                                aes(x = x, y = y, label = atom_number),
                                size = 2, color = "white"
                            ) +

                        ## Add molecule name as label
                            geom_text(
                                data = drop_na(unique(select(plot_data, molecule_name))),
                                aes(x = 1, y = 13, label = molecule_name),
                                size = 4, color = "black", hjust = 0
                            ) +

                        ## Scales and theme
                            scale_size_continuous(range = c(0,4)) +
                            scale_x_continuous(breaks = seq(0,20,1)) +
                            scale_color_viridis(na.value = "#440154FF") +
                            # scale_linetype_manual(values = bond_line_types) +
                            scale_y_continuous(breaks = seq(0,20,1)) +
                            facet_wrap(.~molecule_name, ncol = 3) +
                            # facet_wrap(.~molecule_name) +
                            theme_void() +
                            guides(size = "none", alpha = "none") +
                            coord_fixed() +
                            theme(
                                strip.text = element_blank()
                            )

                    return(plot)
            }           

        #### annotatePDB

            #### annotatePDB

                #' Import chemcial data
                #'
                #' Allows the user to plot a multiple alignent using ggplot's grammar of graphics
                #' @param infile The multiple alignment to (analyze and) plot.
                #' @param roi Plot a subset of the alignment. 
                #' @param consensus Include in the output graphic a line plot corresponding to the degree of conservation at each site in the alignment.
                #' @param funct_assoc Include in the output graphic a line plot of the association between AA residues and 
                #' @importFrom phangorn read.aa
                #' @import phangorn
                #' @import ggtree  
                #' @export
                #' @examples
                #' annotatePDB()

                    annotatePDB <- function(pdb_path, alignment, pdb_name_in_alignment, highlights_data, out_file) {

                        ## Read in PDB file
                            pdb <- bio3d::read.pdb(pdb_path)
                            pdb$atom
                            pdb_sequence <- bio3d::aa321(pdb$atom$resid[bio3d::atom.select(pdb, "calpha")$atom])
                            pdb_frame <- data.frame(pdb_sequence=pdb_sequence, pdb_position=seq(1,length(pdb_sequence),1))
                            pdb_frame$pdb_position <- as.numeric(pdb_frame$pdb_position)
                            head(pdb_frame)

                        ## Read in alignment
                            AA_phydat <- as.data.frame(as.character(phangorn::read.aa(file=alignment, format="fasta")))
                            alignment_sequence <- as.data.frame(AA_phydat)[,grep(pdb_name_in_alignment, colnames(as.data.frame(AA_phydat)))]
                            alignment_frame <- data.frame(alignment_sequence = alignment_sequence, alignment_position = seq(1,length(alignment_sequence),1))
                            head(alignment_frame)

                        ## Add spaces at the front of the PDB frame until its first AA lines up with alignment

                            head(pdb_frame)
                            head(alignment_frame)

                            space <- data.frame(pdb_sequence = "-", pdb_position = "NA")
                            space$pdb_position <- as.numeric(as.character(space$pdb_position))

                            for (i in 1:length(pdb_sequence)) {
                                if ( as.character(alignment_frame[i,1]) == as.character(pdb_frame[i,1]) ) {
                                } else {
                                    pdb_frame <- rbind(pdb_frame[0:(i-1),],space,pdb_frame[i:dim(pdb_frame)[1],])
                                }
                            }

                            head(pdb_frame, 30)
                            head(alignment_frame, 30)

                        ## Add spaces to the end of the PDB frame until it lines up with alignment

                            tail(pdb_frame)
                            tail(alignment_frame)

                            final_space <- data.frame(pdb_sequence="-", pdb_position="NA")
                            final_space$pdb_position <- as.numeric(as.character(final_space$pdb_position))        
                            
                            if (((dim(alignment_frame)[1]-dim(pdb_frame)[1])-1) > 0) {
                                for (i in 1:((dim(alignment_frame)[1]-dim(pdb_frame)[1])-1)) {
                                    final_space <- rbind(final_space, space)
                                }
                                pdb_frame <- rbind(pdb_frame, final_space)
                            } else {
                                pdb_frame <- rbind(pdb_frame, final_space)
                            }            

                            tail(pdb_frame, 20)
                            tail(alignment_frame, 20)
                            
                        ## Bind them together - they should have the same dimensions at this point
                            
                            dim(pdb_frame)
                            dim(alignment_frame)
                            translational_frame <- cbind(alignment_frame, pdb_frame)

                        ## Write it out

                            if (file.exists(out_file)) { file.remove(out_file) }

                            write("## set background", file=out_file, append=TRUE)
                            write("set bg_rgb, white", file=out_file, append=TRUE)
                            write("", file=out_file, append=TRUE)
                            write("## Set base color", file=out_file, append=TRUE)
                            write("set_color myblue = [55,126,184]", file=out_file, append=TRUE)
                            write("select base, resi 1-800", file=out_file, append=TRUE)
                            write("color myblue, base", file=out_file, append=TRUE)
                            write("", file=out_file, append=TRUE)
                            write("## General image improvement", file=out_file, append=TRUE)
                            write("set spec_reflect, 0.2", file=out_file, append=TRUE)

                            for (i in 1:dim(highlights_data)[1]) {
                                color <- paste("[", col2rgb(highlights_data$color[i])[1], ",", col2rgb(highlights_data$color[i])[2], ",", col2rgb(highlights_data$color[i])[3], "]", sep = "") 
                                write(paste("set_color ", paste("COLOR", i, sep = ""), " = ", color, sep = ""), file=out_file, append=TRUE)
                                write(paste("select ", paste("RES", i, sep = ""), paste(", resi", translational_frame$pdb_position[match(highlights_data$xint[i], translational_frame$alignment_position)]), sep = ""), file=out_file, append=TRUE)
                                write(paste("color ", paste("COLOR", i, sep = ""), paste(", RES", i, sep = ""), sep = ""), file=out_file, append=TRUE)
                                write(paste("show sticks", paste(", RES", i, sep = ""), sep = ""), file=out_file, append=TRUE)
                                write(paste(""), file=out_file, append=TRUE)
                            }

                            return(translational_frame)
                    }

        #### vennAnalysis

            #' Get coordinates for quantitative Venn Diagram
            #'
            #' @param data Data frame of TRUE FALSE for whether each observation belongs to each category
            #' @examples
            #' vennAnalysis()

            vennAnalysis <- function(data) {
                plot(euler(data, shape = "circle"), quantities = TRUE)
                venn_circle_data <- plot(eulerr::euler(data, shape = "circle"), quantities = TRUE)$data$ellipses[,1:3]
                colnames(venn_circle_data) <- c("x", "y", "r")
                venn_circle_data$category <- rownames(venn_circle_data)
                return(venn_circle_data)
            }

        #### points3D

            #' Get coordinates for a 3D scatter plot
            #'
            #' @param data Data frame with columns x, y, z, and sample_unique_ID
            #' @param angle Angle of the z axis
            #' @param tick_round Integer to which tick labels should be rounded
            #' @param x_tick_interval Interval of ticks on x axis
            #' @param y_tick_interval Interval of ticks on y axis
            #' @param z_tick_interval Interval of ticks on z axis
            #' @examples
            #' points3D()

            points3D <- function( data, angle, tick_round, x_tick_interval, y_tick_interval, z_tick_interval ) {

                ## Initial calculations

                    output <- list()

                    ymin = plyr::round_any(min(data$y), tick_round, f = floor)
                    ymax = plyr::round_any(max(data$y), tick_round, f = ceiling)
                    ylength = ymax - ymin

                    zmin = plyr::round_any(min(data$z), tick_round, f = floor)
                    zmax = plyr::round_any(max(data$z), tick_round, f = ceiling)
                    zlength = zmax - zmin

                    xmin = plyr::round_any(min(data$x), tick_round, f = floor)
                    xmax = plyr::round_any(max(data$x), tick_round, f = ceiling)
                    xend = xmax+(zlength*cos(angle))
                    xlength = xmax - xmin

                    xintervals = seq(xmin, xmax, x_tick_interval)
                    yintervals = seq(ymin, ymax, y_tick_interval)
                    zintervals = seq(zmin, zmax, z_tick_interval)

                ## Draw grid

                    grid_x <- data.frame(
                        y = ((zintervals-zmin) * sin(angle)),
                        yend = ((zintervals-zmin) * sin(angle)),
                        x = ((zintervals-zmin) * cos(angle)) + xmin,
                        xend = (((zintervals-zmin) * cos(angle)) + xmax)
                    )

                    grid_z <- data.frame(
                        y = 0,
                        yend = zlength*sin(angle),
                        x = xintervals,
                        xend = zlength*cos(angle) + xintervals
                    )

                    output$grid <- rbind(grid_z, grid_x)

                ## Draw ticks

                    ticks_x <- data.frame(
                        y = 0,
                        yend = zlength*0.02*sin(angle),
                        x = xintervals,
                        xend = zlength*0.02*cos(angle) + xintervals
                    )

                    ticks_x2 <- data.frame(
                        y = 0,
                        yend = -zlength*0.02*sin(angle),
                        x = xintervals,
                        xend = -zlength*0.02*cos(angle) + xintervals
                    )

                    ticks_y <- data.frame( ### HERE
                        y = yintervals-ymin,
                        yend = yintervals-ymin,
                        x = -0.02*xlength + xmin,
                        xend = 0.02*xlength + xmin
                    )
                    
                    ticks_z <- data.frame(
                        y = ((zintervals-zmin) * sin(angle)),
                        yend = ((zintervals-zmin) * sin(angle)),
                        x = (((zintervals-zmin) * cos(angle)) + xmax) - 0.02*xlength,
                        xend = (((zintervals-zmin) * cos(angle)) + xmax) + 0.02*xlength
                    )

                    output$ticks <- rbind(ticks_x, ticks_x2, ticks_y, ticks_z)

                ## Draw tick labels

                    labels_x <- data.frame(
                        y = -ylength*0.04,
                        x = xintervals-0.04*xlength,
                        label = xintervals
                    )

                    labels_y <- data.frame(
                        y = yintervals-ymin,
                        x = -xlength*0.08 + xmin,
                        label = yintervals
                    )

                    labels_z <- data.frame(
                        y = ((zintervals-zmin) * sin(angle)),
                        x = (((zintervals-zmin) * cos(angle)) + xmax) + xlength*0.08,
                        label = zintervals
                    )

                    output$labels <- rbind(labels_x, labels_y, labels_z)

                ## Draw axes

                    output$axes <- data.frame(
                        x = c(xmin, xmin, xmax),
                        xend = c(xmax, xmin, xend),
                        y = c(0, 0, 0),
                        yend = c(0, ylength, zlength*sin(angle))
                    )
                
                ## Draw points and their dashed lines

                    output$point_segments <- data.frame(
                        x = data$x + (data$z-zmin)*cos(angle),
                        xend = data$x + (data$z-zmin)*cos(angle),
                        y = data$y - ymin + (data$z-zmin)*sin(angle),
                        yend = (data$z-zmin)*sin(angle)
                    )

                    output$points <- data.frame(
                        x = data$x + (data$z-zmin)*cos(angle),
                        y = data$y - ymin + (data$z-zmin)*sin(angle),
                        sample_unique_ID = data$sample_unique_ID
                    )

                return(output)
            }

    ##### Image Analysis

        #### analyzeImages

            #' Create major and minor axes ticks and labels
            #'
            #' @param major_ticks Values at which major ticks should be generated
            #' @param minor_freq Number of minor ticks between each major tick 
            #' @examples
            #' @export
            #' analyzeImages

            analyzeImages <- function(
                                folder_URL,
                                monolist_out_path
                            ) {

                files <- googledrive::drive_ls(folder_URL)
                file <- 1
                prop <- 0.4
                corner_points_x <<- vector()
                corner_points_y <<- vector()

                ui <- fluidPage(
                    verticalLayout(
                        splitLayout(
                            verticalLayout(
                                splitLayout(
                                    sliderInput("file", "File Number:",
                                        min = 1, max = dim(files)[1],
                                        value = 1, step = 1
                                    ),
                                    actionButton("download_image", "Download Image")
                                ),
                                verbatimTextOutput("photo_unique_id", placeholder = TRUE),
                                verbatimTextOutput("photo_metadata1", placeholder = TRUE),
                                verbatimTextOutput("photo_metadata2", placeholder = TRUE),
                                splitLayout(
                                    verticalLayout(
                                        sliderInput("rotation", "Rotation:",
                                            min = 0, max = 360,
                                            value = 0, step = 90
                                        ),
                                        verbatimTextOutput('selected')
                                    ),
                                    verticalLayout(
                                        actionButton("draw_image", "Draw Image"),
                                        actionButton("colorchecker_select", "ColorChecker Selected. Calibrate."),
                                        actionButton("rgb_select", "RGB Selected. Calibrate."),
                                        actionButton("reset", "Reset.")
                                    )
                                )
                            ),
                            splitLayout(
                                verticalLayout(
                                    actionButton("write_out", "Append to Spreadsheet")
                                ),
                                plotOutput(
                                    outputId = "sample_window",
                                    height = 400
                                )
                            )
                        ),
                        splitLayout(
                            plotOutput(
                                outputId = "photograph",
                                brush = brushOpts(
                                    id = "photograph_brush"
                                ),
                                click = "photograph_click", 
                                dblclick = "photograph_double_click",
                                height = 500
                            ),
                            pickerOutput(
                                outputId = "calibrated_photograph",
                                width = '100%',
                                height = "500px"
                            )
                        )
                    )
                )

                server <- function(input, output, session) {

                    ## Download image

                        observeEvent(input$download_image, {

                            googledrive::drive_download(
                                file = files$id[input$file],
                                path = "temp.JPEG",
                                overwrite = TRUE
                            )

                            path_to_image <- "temp.JPEG"

                            imRed <<- imager::load.image(path_to_image)

                            ## Extract and print QR code data and timestamp

                                photo_metadata <- list()
                                qr_location_data <- quadrangle::qr_scan_cpp(magick::image_read(path_to_image), lighten = TRUE, darken = TRUE)
                                photo_metadata$photo_unique_id <- as.character(files$id[input$file])
                                # photo_metadata$photo_unique_id <- as.character(files$id[file])
                                photo_metadata$name <- quadrangle::qr_scan_js_from_corners(magick::image_read(path_to_image), qr_location_data$points)$data
                                if( length(exifr::read_exif(path_to_image)$SubSecCreateDate) > 0 ){ photo_metadata$time_stamp <- exifr::read_exif(path_to_image)$SubSecCreateDate } else { photo_metadata$time_stamp <- exifr::read_exif(path_to_image)$FileModifyDate }
                                output$photo_unique_id <- renderText({ photo_metadata$photo_unique_id })
                                output$photo_metadata1 <- renderText({ photo_metadata$name })
                                output$photo_metadata2 <- renderPrint({ photo_metadata$time_stamp })
                                photo_metadata$date <- gsub(":", "/", gsub(" .*$", "", photo_metadata$time_stamp))
                                photo_metadata$year <- lubridate::year(gsub(":", "/", gsub(" .*$", "", photo_metadata$time_stamp)))
                                photo_metadata$month <- lubridate::month(gsub(":", "/", gsub(" .*$", "", photo_metadata$time_stamp)))
                                photo_metadata$day <- lubridate::day(gsub(":", "/", gsub(" .*$", "", photo_metadata$time_stamp)))
                                photo_metadata <<- photo_metadata
                        })

                    ## Draw or rotate photograph when asked

                        observeEvent(input$draw_image, {

                            imRed <<- imager::imrotate(imRed, input$rotation)
                            
                            image_as_df <- as_tibble(as.data.frame(as.cimg(imRed)))
                            image_as_df$y <- -as.numeric(image_as_df$y)
                            image_as_df$y <- image_as_df$y + -min(image_as_df$y)
                            image_as_df <- pivot_wider(image_as_df, names_from = cc, values_from = value)
                            names(image_as_df)[3:5] <- c("R", "G", "B")
                            image_as_df$hex <- rgb(image_as_df$R, image_as_df$G, image_as_df$B)
                            image_as_df <<- image_as_df

                            photograph <- ggplot() +
                                geom_tile(data = filter(image_as_df), aes(x = x, y = y), fill = image_as_df$hex) + theme_void()
                                photograph
                            output$photograph <- renderPlot({photograph})

                            corner_points_x <<- vector()
                            corner_points_y <<- vector()

                            # output$selected <- isolate(renderPrint({
                            #     corner_points_x
                            # }))

                        })

                    ## Extract sample spot

                        output$sample_window <- renderPlot({ 

                            if ( !is.null(input$calibrated_photograph_selected_points) ) {

                                sample_pixel_data <<- image_calibrated[input$calibrated_photograph_selected_points,]

                                ggplot() +
                                    geom_tile(
                                        data = sample_pixel_data,
                                        aes(x = x, y = y),
                                        fill = sample_pixel_data$hex
                                    ) +
                                    theme_void()

                                # sample_pixel_output <<- data.frame(
                                #     photo_unique_id = photo_metadata$photo_unique_id,
                                #     name = photo_metadata$name,
                                #     time_stamp = photo_metadata$time_stamp,
                                #     x = sample_pixel_data$x,
                                #     y = sample_pixel_data$y,
                                #     color = sample_pixel_data$hex
                                # )

                            } else { NULL}

                        })
                                    
                    ## Run ColorChecker correction

                        output$selected <- renderPrint({
                            corner_points_x <<- c(corner_points_x, input$photograph_click$x)
                            corner_points_y <<- c(corner_points_y, input$photograph_click$y)
                            corner_points_x
                        })

                        # observeEvent(input$reset, {
                            

                        #     output$selected <- renderPrint({
                        #         corner_points_x
                        #     })
                        # })

                        observeEvent(input$colorchecker_select, {

                            ## Pre-calculations
                                
                                mR <- raster::as.matrix(imager::R(imRed)) * 255
                                mG <- raster::as.matrix(imager::G(imRed)) * 255
                                mB <- raster::as.matrix(imager::B(imRed)) * 255
                                rR <- raster::raster(mR)
                                rG <- raster::raster(mG)
                                rB <- raster::raster(mB)
                                raster::extent(rR) <- c(0, dim(imRed)[2], 0, dim(imRed)[1])
                                raster::extent(rG) <- c(0, dim(imRed)[2], 0, dim(imRed)[1])
                                raster::extent(rB) <- c(0, dim(imRed)[2], 0, dim(imRed)[1])
                                rR <- raster::flip(t(rR), "y")
                                rG <- raster::flip(t(rG), "y")
                                rB <- raster::flip(t(rB), "y")

                            ## Get the border of the colorchecker and swatch positions

                                xy <- list()
                                xy$x <- corner_points_x
                                xy$y <- corner_points_y

                                xyDF <- as.data.frame(xy)
                                line1 <- xyDF[1:2, ]
                                line2 <- xyDF[2:3, ]
                                line3 <- xyDF[3:4, ]
                                line4 <- xyDF[c(4, 1), ]
                                xdiff <- (line1$x[2] - line1$x[1])/6
                                ydiff <- (line1$y[2] - line1$y[1])/6
                                xmin <- line1$x[1]
                                ymin <- line1$y[1]
                                xySubA <- list(x = c((xmin + xmin + xdiff)/2, (xmin + 
                                    xdiff + xmin + xdiff * 2)/2, (xmin + xdiff * 
                                    2 + xmin + xdiff * 3)/2, (xmin + xdiff * 3 + 
                                    xmin + xdiff * 4)/2, (xmin + xdiff * 4 + xmin + 
                                    xdiff * 5)/2, (xmin + xdiff * 5 + xmin + xdiff * 
                                    6)/2), y = c((ymin + ymin + ydiff)/2, (ymin + 
                                    ydiff + ymin + ydiff * 2)/2, (ymin + ydiff * 
                                    2 + ymin + ydiff * 3)/2, (ymin + ydiff * 3 + 
                                    ymin + ydiff * 4)/2, (ymin + ydiff * 4 + ymin + 
                                    ydiff * 5)/2, (ymin + ydiff * 5 + ymin + ydiff * 
                                    6)/2))
                                xySubDF_1A <- as.data.frame(xySubA)
                                xySubBa <- list(x = c(xmin + xdiff - (xdiff/2) * 
                                    prop, xmin + xdiff * 2 - (xdiff/2) * prop, xmin + 
                                    xdiff * 3 - (xdiff/2) * prop, xmin + xdiff * 
                                    4 - (xdiff/2) * prop, xmin + xdiff * 5 - (xdiff/2) * 
                                    prop, xmin + xdiff * 6 - (xdiff/2) * prop), y = c(ymin + 
                                    ydiff - (ydiff/2) * prop, ymin + ydiff * 2 - 
                                    (ydiff/2) * prop, ymin + ydiff * 3 - (ydiff/2) * 
                                    prop, ymin + ydiff * 4 - (ydiff/2) * prop, ymin + 
                                    ydiff * 5 - (ydiff/2) * prop, ymin + ydiff * 
                                    6 - (ydiff/2) * prop))
                                xySubBb <- list(x = c(xmin + (xdiff/2) * prop, xmin + 
                                    xdiff + (xdiff/2) * prop, xmin + xdiff * 2 + 
                                    (xdiff/2) * prop, xmin + xdiff * 3 + (xdiff/2) * 
                                    prop, xmin + xdiff * 4 + (xdiff/2) * prop, xmin + 
                                    xdiff * 5 + (xdiff/2) * prop), y = c(ymin + (ydiff/2) * 
                                    prop, ymin + ydiff + (ydiff/2) * prop, ymin + 
                                    ydiff * 2 + (ydiff/2) * prop, ymin + ydiff * 
                                    3 + (ydiff/2) * prop, ymin + ydiff * 4 + (ydiff/2) * 
                                    prop, ymin + ydiff * 5 + (ydiff/2) * prop))
                                xySubDF_1Ba <- as.data.frame(xySubBa)
                                xySubDF_1Bb <- as.data.frame(xySubBb)
                                xdiff <- (line3$x[2] - line3$x[1])/6
                                ydiff <- (line3$y[2] - line3$y[1])/6
                                xmin <- line3$x[1]
                                ymin <- line3$y[1]
                                xySubA <- list(x = c((xmin + xmin + xdiff)/2, (xmin + 
                                    xdiff + xmin + xdiff * 2)/2, (xmin + xdiff * 
                                    2 + xmin + xdiff * 3)/2, (xmin + xdiff * 3 + 
                                    xmin + xdiff * 4)/2, (xmin + xdiff * 4 + xmin + 
                                    xdiff * 5)/2, (xmin + xdiff * 5 + xmin + xdiff * 
                                    6)/2), y = c((ymin + ymin + ydiff)/2, (ymin + 
                                    ydiff + ymin + ydiff * 2)/2, (ymin + ydiff * 
                                    2 + ymin + ydiff * 3)/2, (ymin + ydiff * 3 + 
                                    ymin + ydiff * 4)/2, (ymin + ydiff * 4 + ymin + 
                                    ydiff * 5)/2, (ymin + ydiff * 5 + ymin + ydiff * 
                                    6)/2))
                                xySubDF_3A <- as.data.frame(xySubA)
                                xySubBa <- list(x = c(xmin + xdiff - (xdiff/2) * 
                                    prop, xmin + xdiff * 2 - (xdiff/2) * prop, xmin + 
                                    xdiff * 3 - (xdiff/2) * prop, xmin + xdiff * 
                                    4 - (xdiff/2) * prop, xmin + xdiff * 5 - (xdiff/2) * 
                                    prop, xmin + xdiff * 6 - (xdiff/2) * prop), y = c(ymin + 
                                    ydiff - (ydiff/2) * prop, ymin + ydiff * 2 - 
                                    (ydiff/2) * prop, ymin + ydiff * 3 - (ydiff/2) * 
                                    prop, ymin + ydiff * 4 - (ydiff/2) * prop, ymin + 
                                    ydiff * 5 - (ydiff/2) * prop, ymin + ydiff * 
                                    6 - (ydiff/2) * prop))
                                xySubBb <- list(x = c(xmin + (xdiff/2) * prop, xmin + 
                                    xdiff + (xdiff/2) * prop, xmin + xdiff * 2 + 
                                    (xdiff/2) * prop, xmin + xdiff * 3 + (xdiff/2) * 
                                    prop, xmin + xdiff * 4 + (xdiff/2) * prop, xmin + 
                                    xdiff * 5 + (xdiff/2) * prop), y = c(ymin + (ydiff/2) * 
                                    prop, ymin + ydiff + (ydiff/2) * prop, ymin + 
                                    ydiff * 2 + (ydiff/2) * prop, ymin + ydiff * 
                                    3 + (ydiff/2) * prop, ymin + ydiff * 4 + (ydiff/2) * 
                                    prop, ymin + ydiff * 5 + (ydiff/2) * prop))
                                xySubDF_3Ba <- as.data.frame(xySubBa)
                                xySubDF_3Bb <- as.data.frame(xySubBb)
                                xdiff <- (line2$x[2] - line2$x[1])/4
                                ydiff <- (line2$y[2] - line2$y[1])/4
                                xmin <- line2$x[1]
                                ymin <- line2$y[1]
                                xySubA <- list(x = c((xmin + xmin + xdiff)/2, (xmin + 
                                    xdiff + xmin + xdiff * 2)/2, (xmin + xdiff * 
                                    2 + xmin + xdiff * 3)/2, (xmin + xdiff * 3 + 
                                    xmin + xdiff * 4)/2), y = c((ymin + ymin + ydiff)/2, 
                                    (ymin + ydiff + ymin + ydiff * 2)/2, (ymin + 
                                      ydiff * 2 + ymin + ydiff * 3)/2, (ymin + ydiff * 
                                      3 + ymin + ydiff * 4)/2))
                                xySubDF_2A <- as.data.frame(xySubA)
                                xySubBa <- list(x = c(xmin + xdiff - (xdiff/2) * 
                                    prop, xmin + xdiff * 2 - (xdiff/2) * prop, xmin + 
                                    xdiff * 3 - (xdiff/2) * prop, xmin + xdiff * 
                                    4 - (xdiff/2) * prop), y = c(ymin + ydiff - (ydiff/2) * 
                                    prop, ymin + ydiff * 2 - (ydiff/2) * prop, ymin + 
                                    ydiff * 3 - (ydiff/2) * prop, ymin + ydiff * 
                                    4 - (ydiff/2) * prop))
                                xySubBb <- list(x = c(xmin + (xdiff/2) * prop, xmin + 
                                    xdiff + (xdiff/2) * prop, xmin + xdiff * 2 + 
                                    (xdiff/2) * prop, xmin + xdiff * 3 + (xdiff/2) * 
                                    prop), y = c(ymin + (ydiff/2) * prop, ymin + 
                                    ydiff + (ydiff/2) * prop, ymin + ydiff * 2 + 
                                    (ydiff/2) * prop, ymin + ydiff * 3 + (ydiff/2) * 
                                    prop))
                                xySubDF_2Ba <- as.data.frame(xySubBa)
                                xySubDF_2Bb <- as.data.frame(xySubBb)
                                xdiff <- (line4$x[2] - line4$x[1])/4
                                ydiff <- (line4$y[2] - line4$y[1])/4
                                xmin <- line4$x[1]
                                ymin <- line4$y[1]
                                xySubA <- list(x = c((xmin + xmin + xdiff)/2, (xmin + 
                                    xdiff + xmin + xdiff * 2)/2, (xmin + xdiff * 
                                    2 + xmin + xdiff * 3)/2, (xmin + xdiff * 3 + 
                                    xmin + xdiff * 4)/2), y = c((ymin + ymin + ydiff)/2, 
                                    (ymin + ydiff + ymin + ydiff * 2)/2, (ymin + 
                                      ydiff * 2 + ymin + ydiff * 3)/2, (ymin + ydiff * 
                                      3 + ymin + ydiff * 4)/2))
                                xySubDF_4A <- as.data.frame(xySubA)
                                xySubBa <- list(x = c(xmin + xdiff - (xdiff/2) * 
                                    prop, xmin + xdiff * 2 - (xdiff/2) * prop, xmin + 
                                    xdiff * 3 - (xdiff/2) * prop, xmin + xdiff * 
                                    4 - (xdiff/2) * prop), y = c(ymin + ydiff - (ydiff/2) * 
                                    prop, ymin + ydiff * 2 - (ydiff/2) * prop, ymin + 
                                    ydiff * 3 - (ydiff/2) * prop, ymin + ydiff * 
                                    4 - (ydiff/2) * prop))
                                xySubBb <- list(x = c(xmin + (xdiff/2) * prop, xmin + 
                                    xdiff + (xdiff/2) * prop, xmin + xdiff * 2 + 
                                    (xdiff/2) * prop, xmin + xdiff * 3 + (xdiff/2) * 
                                    prop), y = c(ymin + (ydiff/2) * prop, ymin + 
                                    ydiff + (ydiff/2) * prop, ymin + ydiff * 2 + 
                                    (ydiff/2) * prop, ymin + ydiff * 3 + (ydiff/2) * 
                                    prop))
                                xySubDF_4Ba <- as.data.frame(xySubBa)
                                xySubDF_4Bb <- as.data.frame(xySubBb)

                                sub_data <<- rbind(xySubDF_1A, xySubDF_1Ba, xySubDF_1Bb, xySubDF_2A, xySubDF_2Ba, xySubDF_2Bb, xySubDF_3A, xySubDF_3Ba, xySubDF_3Bb, xySubDF_4A, xySubDF_4Ba, xySubDF_4Bb)

                                labels <- list(
                                    c(1, 7, 13, 19), c(2, 8, 14, 20), 
                                    c(3, 9, 15, 21), c(4, 10, 16, 22),
                                    c(5, 11, 17, 23), c(6, 12, 18, 24)
                                )

                            ## Get swatch positions and labels
                                
                                swatch_labels <- list()
                                xyTot <- c()
                                for (e in 1:nrow(xySubDF_1A)) {
                                    xyLine1 <- xySubDF_1A[e, ]
                                    xyLine2 <- xySubDF_3A[6:1, ][e, ]
                                    xdiff <- (xyLine2$x - xyLine1$x)/4
                                    ydiff <- (xyLine2$y - xyLine1$y)/4
                                    xmin <- xyLine1$x
                                    ymin <- xyLine1$y
                                    xySub <- list(x = c((xmin + xmin + xdiff)/2, 
                                      (xmin + xdiff + xmin + xdiff * 2)/2, (xmin + 
                                        xdiff * 2 + xmin + xdiff * 3)/2, (xmin + 
                                        xdiff * 3 + xmin + xdiff * 4)/2), y = c((ymin + 
                                      ymin + ydiff)/2, (ymin + ydiff + ymin + ydiff * 
                                      2)/2, (ymin + ydiff * 2 + ymin + ydiff * 3)/2, 
                                      (ymin + ydiff * 3 + ymin + ydiff * 4)/2))
                                    xySubDF <- as.data.frame(xySub)

                                    swatch_labels[[e]] <- xySubDF
                                    
                                    xySubDFLabel <- cbind(xySubDF, label = labels[[e]])
                                    xyLine1a <- xySubDF_1Ba[e, ]
                                    xyLine1b <- xySubDF_1Bb[e, ]
                                    xyLine3a <- xySubDF_3Bb[6:1, ][e, ]
                                    xyLine3b <- xySubDF_3Ba[6:1, ][e, ]
                                    xdiffa <- (xyLine3a$x - xyLine1a$x)/4
                                    ydiffa <- (xyLine3a$y - xyLine1a$y)/4
                                    xdiffb <- (xyLine3b$x - xyLine1b$x)/4
                                    ydiffb <- (xyLine3b$y - xyLine1b$y)/4
                                    xmina <- xyLine1a$x
                                    ymina <- xyLine1a$y
                                    xminb <- xyLine1b$x
                                    yminb <- xyLine1b$y
                                    xySubAa <- list(x = c(xmina + xdiffa - (xdiffa/2) * 
                                      prop, xmina + xdiffa * 2 - (xdiffa/2) * prop, 
                                      xmina + xdiffa * 3 - (xdiffa/2) * prop, xmina + 
                                        xdiffa * 4 - (xdiffa/2) * prop), y = c(ymina + 
                                      ydiffa - (ydiffa/2) * prop, ymina + ydiffa * 
                                      2 - (ydiffa/2) * prop, ymina + ydiffa * 3 - 
                                      (ydiffa/2) * prop, ymina + ydiffa * 4 - (ydiffa/2) * 
                                      prop))
                                    xySubAb <- list(x = c(xmina + (xdiffa/2) * prop, 
                                      xmina + xdiffa + (xdiffa/2) * prop, xmina + 
                                        xdiffa * 2 + (xdiffa/2) * prop, xmina + xdiffa * 
                                        3 + (xdiffa/2) * prop), y = c(ymina + (ydiffa/2) * 
                                      prop, ymina + ydiffa + (ydiffa/2) * prop, ymina + 
                                      ydiffa * 2 + (ydiffa/2) * prop, ymina + ydiffa * 
                                      3 + (ydiffa/2) * prop))
                                    xySubBa <- list(x = c(xminb + xdiffb - (xdiffb/2) * 
                                      prop, xminb + xdiffb * 2 - (xdiffb/2) * prop, 
                                      xminb + xdiffb * 3 - (xdiffb/2) * prop, xminb + 
                                        xdiffb * 4 - (xdiffb/2) * prop), y = c(yminb + 
                                      ydiffb - (ydiffb/2) * prop, yminb + ydiffb * 
                                      2 - (ydiffb/2) * prop, yminb + ydiffb * 3 - 
                                      (ydiffb/2) * prop, yminb + ydiffb * 4 - (ydiffb/2) * 
                                      prop))
                                    xySubBb <- list(x = c(xminb + (xdiffb/2) * prop, 
                                      xminb + xdiffb + (xdiffb/2) * prop, xminb + 
                                        xdiffb * 2 + (xdiffb/2) * prop, xminb + xdiffb * 
                                        3 + (xdiffb/2) * prop), y = c(yminb + (ydiffb/2) * 
                                      prop, yminb + ydiffb + (ydiffb/2) * prop, yminb + 
                                      ydiffb * 2 + (ydiffb/2) * prop, yminb + ydiffb * 
                                      3 + (ydiffb/2) * prop))
                                    xySubAaDF <- as.data.frame(xySubAa)
                                    xySubAbDF <- as.data.frame(xySubAb)
                                    xySubBaDF <- as.data.frame(xySubBa)
                                    xySubBbDF <- as.data.frame(xySubBb)
                                    xySubRow <- cbind(xySubAaDF, xySubAbDF, xySubBaDF, 
                                      xySubBbDF, label = labels[[e]])
                                    colnames(xySubRow) <- c("x1", "y1", "x2", "y2", 
                                      "x3", "y3", "x4", "y4", "label")
                                    xyTot <- rbind(xyTot, xySubRow)
                                }

                                swatch_labels <- do.call(rbind, swatch_labels)
                                swatch_labels$label <- unlist(labels)

                            ## Get swatch colors

                                xyTot$imR <- NA
                                xyTot$imG <- NA
                                xyTot$imB <- NA
                                xyTot <- xyTot[order(xyTot$label), ]
                                polygons <- list()
                                for (e in 1:nrow(xyTot)) {
                                    
                                    polygons[[e]] <- data.frame(
                                        xmax = xyTot$x1[e],
                                        xmin = xyTot$x2[e],
                                        ymax = xyTot$y1[e],
                                        ymin = xyTot$y3[e],
                                        label = e
                                    )

                                    image_as_df %>%
                                        filter(
                                            x < xyTot$x1[e] &
                                            x > xyTot$x2[e] &
                                            y < xyTot$y1[e] &
                                            y > xyTot$y3[e]
                                        ) -> polygon_data

                                    xyTot$imR[e] <- mode(polygon_data$R) * 255
                                    xyTot$imG[e] <- mode(polygon_data$G) * 255
                                    xyTot$imB[e] <- mode(polygon_data$B) * 255

                                }
                                polygons <- do.call(rbind, polygons)


                                l1 <- c(1, 115, 82, 68)
                                l2 <- c(2, 194, 150, 130)
                                l3 <- c(3, 98, 122, 157)
                                l4 <- c(4, 87, 108, 67)
                                l5 <- c(5, 133, 128, 177)
                                l6 <- c(6, 103, 189, 170)
                                l7 <- c(7, 214, 126, 44)
                                l8 <- c(8, 80, 91, 166)
                                l9 <- c(9, 193, 90, 99)
                                l10 <- c(10, 94, 60, 108)
                                l11 <- c(11, 157, 188, 64)
                                l12 <- c(12, 224, 163, 46)
                                l13 <- c(13, 56, 61, 150)
                                l14 <- c(14, 70, 148, 73)
                                l15 <- c(15, 175, 54, 60)
                                l16 <- c(16, 231, 199, 31)
                                l17 <- c(17, 187, 86, 149)
                                l18 <- c(18, 8, 133, 161)
                                l19 <- c(19, 243, 243, 242)
                                l20 <- c(20, 200, 200, 200)
                                l21 <- c(21, 160, 160, 160)
                                l22 <- c(22, 122, 122, 121)
                                l23 <- c(23, 85, 85, 85)
                                l24 <- c(24, 52, 52, 52)
                                ColorCheckerRGB <- as.data.frame(rbind(l1, l2, l3, 
                                    l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, 
                                    l15, l16, l17, l18, l19, l20, l21, l22, l23, 
                                    l24)
                                )
                                colnames(ColorCheckerRGB) <- c("label", "sR", "sG", "sB")

                                dat <<- merge(xyTot, ColorCheckerRGB, by = "label")

                                ## Calculate regression for the colors

                                    print("Calculating polynomial regression...")
                                        sR <- dat$sR
                                        sG <- dat$sG
                                        sB <- dat$sB
                                        imR <- dat$imR
                                        imG <- dat$imG
                                        imB <- dat$imB
                                        modelR <- lm(sR ~ imR + imG + imB + imR^2 + imG^2 + 
                                            imB^2)
                                        modelG <- lm(sG ~ imR + imG + imB + imR^2 + imG^2 + 
                                            imB^2)
                                        modelB <- lm(sB ~ imR + imG + imB + imR^2 + imG^2 + 
                                            imB^2)
                                        dfIm = data.frame(imR = matrix(mR, ncol = 1), imG = matrix(mG, 
                                            ncol = 1), imB = matrix(mB, ncol = 1))

                                ## Calibrate colors

                                    print("Calibrating colors...")
                                        prR <- predict(modelR, dfIm)
                                        prG <- predict(modelG, dfIm)
                                        prB <- predict(modelB, dfIm)
                                        dfCal <- as.data.frame(cbind(prR, prG, prB))
                                    
                                    print("Rebuilding image...")
                                        Ri = matrix(dfCal$prR, nrow = dim(imRed)[1])
                                        Gi = matrix(dfCal$prG, nrow = dim(imRed)[1])
                                        Bi = matrix(dfCal$prB, nrow = dim(imRed)[1])
                                        imCal = array(dim = dim(imRed))
                                        imCal[, , , 1] = Ri
                                        imCal[, , , 2] = Gi
                                        imCal[, , , 3] = Bi
                                        imCal <- imager::as.cimg(imCal)
                                        
                                        image_calibrated <<- as_tibble(as.data.frame(as.cimg(imCal)))
                                        image_calibrated$y <- -as.numeric(image_calibrated$y)
                                        image_calibrated$y <- image_calibrated$y + -min(image_calibrated$y)
                                        image_calibrated <- pivot_wider(image_calibrated, names_from = cc, values_from = value)
                                        names(image_calibrated)[3:5] <- c("R", "G", "B")

                                        image_calibrated$R[image_calibrated$R/255 < 0] <- 0
                                        image_calibrated$G[image_calibrated$G/255 < 0] <- 0
                                        image_calibrated$B[image_calibrated$B/255 < 0] <- 0

                                        image_calibrated$R[image_calibrated$R/255 > 1] <- 1
                                        image_calibrated$G[image_calibrated$G/255 > 1] <- 1
                                        image_calibrated$B[image_calibrated$B/255 > 1] <- 1

                                        image_calibrated$hex <- rgb(image_calibrated$R/255, image_calibrated$G/255, image_calibrated$B/255)
                                        image_calibrated <<- image_calibrated

                                output$photograph <- renderPlot({
                                    
                                    ggplot() +
                                        geom_tile(data = filter(image_as_df), aes(x = x, y = y), fill = image_as_df$hex) + theme_void() +
                                        geom_point(
                                            data = sub_data,
                                            aes(x = x, y = y), shape = 21, size = 3, color = "black", fill = "maroon"
                                        ) +
                                        geom_text(
                                            data = swatch_labels,
                                            aes(x = x, y = y, label = label), color = "gold"
                                        ) +
                                        geom_rect(
                                            data = polygons,
                                            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = "gold", alpha = 0
                                        )

                                })

                                output$calibrated_photograph <- renderPicker({

                                    picker(
                                        coords = as.data.frame(image_calibrated[,1:2]),
                                        colors = as.data.frame(image_calibrated[,6])[[1]],
                                        # cluster_colors = image_as_df[,6],
                                        labels = as.character(as.data.frame(image_calibrated[,1])[[1]]) 
                                        # labels = NULL
                                        # label_coords = label_coords,
                                        # polygons = polygons, 
                                        # text_props = text_props,
                                        # point_color_polygons = 'white',
                                        # grid_legend_items = grid_legend_items
                                    )
                                    
                                    # ggplot() +
                                    #     geom_tile(data = filter(image_calibrated), aes(x = x, y = y), fill = image_calibrated$hex) + theme_void()
                                        # geom_point(
                                        #     data = sub_data,
                                        #     aes(x = x, y = y), shape = 21, size = 3, color = "black", fill = "maroon"
                                        # ) +
                                        # geom_text(
                                        #     data = swatch_labels,
                                        #     aes(x = x, y = y, label = label), color = "gold"
                                        # ) +
                                        # geom_rect(
                                        #     data = polygons,
                                        #     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = "gold"
                                        # )

                                })

                        })

                    ## Extract rgb_window

                        observeEvent(input$rgb_select, {

                            if ( isolate(!is.null(input$photograph_brush)) ) {
                                
                                image_points <<- isolate(brushedPoints(image_as_df, input$photograph_brush))
                                
                                Rmax <- cbind(
                                    data.frame(color_swatch = "red"),
                                    image_points[order(image_points$R, decreasing = TRUE),][1:50,]
                                )
                                Gmax <- cbind(
                                    data.frame(color_swatch = "green"),
                                    image_points[order(image_points$G, decreasing = TRUE),][1:50,]
                                )
                                Bmax <- cbind(
                                    data.frame(color_swatch = "blue"),
                                    image_points[order(image_points$B, decreasing = TRUE),][1:50,]
                                )
                                
                                rgb_window <- ggplot() +
                                    geom_tile(data = image_points, aes(x = x, y = y), fill = image_points$hex) +
                                    geom_rect(aes(xmin = min(Rmax$x), xmax = max(Rmax$x), ymin = min(Rmax$y), ymax = max(Rmax$y)), alpha = 0.5) +
                                    geom_rect(aes(xmin = min(Gmax$x), xmax = max(Gmax$x), ymin = min(Gmax$y), ymax = max(Gmax$y)), alpha = 0.5) +
                                    geom_rect(aes(xmin = min(Bmax$x), xmax = max(Bmax$x), ymin = min(Bmax$y), ymax = max(Bmax$y)), alpha = 0.5)

                                output$rgb_window <- renderPlot({ rgb_window })

                                swatch_data <<- rbind(
                                    select(Rmax, color_swatch, x, y, R) %>%
                                    set_colnames(c("color_swatch", "x", "y", "value")),
                                    select(Gmax, color_swatch, x, y, G) %>%
                                    set_colnames(c("color_swatch", "x", "y", "value")),
                                    select(Bmax, color_swatch, x, y, B) %>%
                                    set_colnames(c("color_swatch", "x", "y", "value"))
                                )
                                swatch_data$color_swatch <- factor(swatch_data$color_swatch, levels = c("red", "green", "blue"))

                                rgb_histogram <- ggplot(data = swatch_data) + 
                                    geom_histogram(aes(x = value, fill = color_swatch), binwidth = 0.01, color = "black") +
                                    scale_fill_manual(values = c("red", "green", "blue")) +
                                    facet_grid(color_swatch~.) +
                                    scale_x_continuous(limits = c(-0.1,1.1)) +
                                    theme_bw()

                                output$rgb_histogram <- renderPlot({ rgb_histogram })

                                # mode_reference_value_R <<- mode(swatch_data$value[swatch_data$color_swatch == "red"])
                                # mode_reference_value_G <<- mode(swatch_data$value[swatch_data$color_swatch == "green"])
                                # mode_reference_value_B <<- mode(swatch_data$value[swatch_data$color_swatch == "blue"])

                                mode_reference_value_R <<- mean(swatch_data$value[swatch_data$color_swatch == "red"])
                                mode_reference_value_G <<- mean(swatch_data$value[swatch_data$color_swatch == "green"])
                                mode_reference_value_B <<- mean(swatch_data$value[swatch_data$color_swatch == "blue"])

                                output$mode_reference_value_R <- renderText({ mode_reference_value_R })
                                output$mode_reference_value_G <- renderText({ mode_reference_value_G })
                                output$mode_reference_value_B <- renderText({ mode_reference_value_B })

                                reference_pixel_output <<- data.frame(
                                    QR = if(is.null(photo_metadata$name)) {"NA"} else {photo_metadata$name},
                                    Year = photo_metadata$year,
                                    Month = photo_metadata$month,
                                    Day = photo_metadata$day,
                                    Date = photo_metadata$date,
                                    type = "reference",
                                    x = image_points$x,
                                    y = image_points$y,
                                    R = image_points$R,
                                    G = image_points$G,
                                    B = image_points$B
                                )

                            } else { NULL }

                        })

                    ## Write out

                        ## Here be sure to use a unique_photo_ID and have the app check to see if its already in output before writing out!

                        ## Maybe a share link could be the unique photo ID

                        observeEvent(input$write_out, {

                            if (file.exists(monolist_out_path)) {

                                writeMonolist(
                                    cbind(as.data.frame(photo_metadata), sample_pixel_data),
                                    append = TRUE,
                                    monolist_out_path = monolist_out_path,
                                    col.names = FALSE
                                )
                                cat("Appended to Spreadsheet")

                            } else {

                                writeMonolist(
                                    cbind(as.data.frame(photo_metadata), sample_pixel_data),
                                    append = TRUE,
                                    monolist_out_path = monolist_out_path,
                                    col.names = TRUE
                                )
                                cat("Appended to Spreadsheet")

                            }
                                
                        })

                    ## Meta analysis tab

                            # observeEvent(input$metadata_file, {

                            #     ## Read in metadata file
                            #         metadata_file <<- input$metadata_file
                            #         print(metadata_file)
                            #         prelim_data <<- readMonolist(metadata_file$datapath)
                            #         print(head(prelim_data))
                            
                            #     ## Extract swatches
                            #         prelim_data %>%
                            #             pivot_longer(cols = 9:11, names_to = "color", values_to = "value") %>%
                            #             group_by(QR, color) %>%
                            #             arrange(desc(value), .by_group = TRUE) %>%
                            #             dplyr::slice(1:25) -> swatch_data

                            #     ## Determine global means
                            #         swatch_data$global_mean <- 0  
                            #         swatch_data$global_mean[swatch_data$color == "R"] <- mean(swatch_data$value[swatch_data$color == "R"])
                            #         swatch_data$global_mean[swatch_data$color == "G"] <- mean(swatch_data$value[swatch_data$color == "B"])
                            #         swatch_data$global_mean[swatch_data$color == "B"] <- mean(swatch_data$value[swatch_data$color == "G"])

                            #     ## Create dataframe according to which normalization should be performed
                            #         swatch_data %>%
                            #         group_by(QR, color) %>%
                            #         summarize(difference = unique(global_mean) - unique(mean(value))) %>%
                            #         # mutate(unique_id = paste(QR, "_", color)) %>%
                            #         pivot_wider(names_from = color, values_from = difference, names_prefix = "adjust_") -> normalization_data

                            #     ## Modify swatch data according to global means
                            #         swatch_data %>%
                            #         group_by(QR, color) %>%
                            #         mutate(modified_value = value + (global_mean - mean(value))) %>%
                            #         pivot_longer(cols = 10:12, names_to = "value_type", values_to = "value") %>%
                            #         ggplot() +
                            #             geom_histogram(aes(x = value, fill = value_type), bins = 30, alpha = 0.7, color = "black") +
                            #             facet_grid(QR~color) +
                            #             theme_bw() -> global_means_plot

                            #         output$metaplot <- renderPlot({ global_means_plot })

                            # })

                        ## Modify reference data according to global means
                            
                            # prelim_data %>%
                            #     mutate(
                            #         R_adjusted = R + normalization_data$adjust_R[match(QR, normalization_data$QR)],
                            #         G_adjusted = G + normalization_data$adjust_G[match(QR, normalization_data$QR)],
                            #         B_adjusted = B + normalization_data$adjust_B[match(QR, normalization_data$QR)]
                            #     ) -> prelim_data_processed

                            # prelim_data_processed$R_adjusted[prelim_data_processed$R_adjusted < 0] <- 0
                            # prelim_data_processed$G_adjusted[prelim_data_processed$G_adjusted < 0] <- 0
                            # prelim_data_processed$B_adjusted[prelim_data_processed$B_adjusted < 0] <- 0

                            # prelim_data_processed$R_adjusted[prelim_data_processed$R_adjusted > 1] <- 1
                            # prelim_data_processed$G_adjusted[prelim_data_processed$G_adjusted > 1] <- 1
                            # prelim_data_processed$B_adjusted[prelim_data_processed$B_adjusted > 1] <- 1

                            # prelim_data_processed %>%
                            # filter(type == "reference") %>%
                            # mutate(
                            # hex_before = rgb(R, G, B),
                            # hex_after = rgb(R_adjusted, G_adjusted, B_adjusted)
                            # ) %>%
                            # pivot_longer(cols = c("hex_before", "hex_after"), names_to = "state", values_to = "hex") %>%
                            # filter(QR %in% unique(prelim_data_processed$QR)) -> plot_data

                            # ggplot() +
                            # geom_tile(data = plot_data, aes(x = x, y = y), fill = plot_data$hex) +
                            # facet_grid(QR~state, space = "free", scales = "free") -> plot

                            # png(filename = "test.png", width = 10, height = 10, units = "in", res = 300)
                            # plot
                            # dev.off()

                            # ## Plot swatches before and after
                            #     ggplot() +
                            #         geom_tile(
                            #         data = filter(swatch_data, value_type == "value"),
                            #         aes(x = x, y = y, fill = hex)
                            #     )

                            # swatch_data %>%
                            #     group_by(QR, color) %>%
                            #     mutate(modified_value = value + (global_mean - mean(value))) %>%
                            #     pivot_longer(cols = 11:13, names_to = "value_type", values_to = "value") %>%
                            #     filter(value_type == "value")
                            #     ggplot() +
                            #     geom_tile(aes(x = x, y = y, fill = hex)) +
                            #     facet_grid(QR~color) +
                            #     theme_bw()

                            #     out <- list()
                            #     seq <- seq(0.6,0.8,0.001)
                            #     for (i in 1:length(seq)) {
                            #     test$black_v_white <- 0
                            #     test$black_v_white[test$greyscale > seq[i]] <- 1
                            #     out[[i]] <- data.frame(
                            #     ratio = 
                            #     as.data.frame(table(test$black_v_white))$Freq[1]/
                            #     as.data.frame(table(test$black_v_white))$Freq[2],
                            #     i = i
                            #     )
                            #     }
                            #     out <- do.call(rbind, out)

                            #     ggplot(out, aes(x = i, y = ratio)) + geom_col()



                            #     plot_a <- test %>%
                            #     ggplot(aes(x = x, y = y)) +
                            #     geom_tile(
                            #     aes(fill = greyscale),
                            #     show.legend = FALSE
                            #     ) +
                            #     scale_fill_gradient(low = "black", high = "white")

                            #     plot_b <- test %>%
                            #     ggplot(aes(x = x, y = y)) +
                            #     geom_tile(
                            #     aes(fill = black_v_white),
                            #     show.legend = FALSE
                            #     ) +
                            #     scale_fill_gradient(low = "black", high = "white")

                            #     plot_c <- ggplot(test) +
                            #     geom_histogram(aes(x = greyscale), bins = 80) +
                            #     geom_vline(xintercept = threshold)

                            #     d <- as.data.frame(table(test$black_v_white))

                            #     plot_d <- ggplot(data = d) + geom_col(aes(x = 1, y = Freq, fill = Var1))

                            #     cowplot::plot_grid(plot_a, plot_b, plot_c, plot_d)

                }

                shinyApp(ui = ui, server = server)

            }

        #### generateQRcodes

            #' Create major and minor axes ticks and labels
            #'
            #' @param share_link Share link to Google Sheet with QR code metadata
            #' @examples
            #' @export
            #' analyzeImage

            generateQRcodes <- function( share_link ) {

                ## Read in the Google docs spreadsheet

                    data <- read_sheet(share_link)

                ## Create QR codes for the plants

                    plot_list <- list()
                    for (i in 1:dim(data)[1]) {
                        
                        metadata <- unite(data, col = "metadata")$metadata

                        png(paste0(metadata, ".png"))
                            image(
                                1 - qrencoder::qrencode(metadata),
                                asp = 1, xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1),
                                col = c("black", "white")
                            )
                        dev.off()

                        plot_list[[i]] <- cowplot::plot_grid(
                            ggpubr::text_grob(metadata, size = 5),
                            grid::rasterGrob(png::readPNG(paste0(metadata, ".png"))),
                            ncol = 1, rel_heights = c(1,30)
                        )

                        file.remove(paste0(metadata, ".png"))

                    }

                ## return
                    return(plot_list)
                    message("Consider using 'do.call(\"grid.arrange\", c(plot_list))' to plot all the plots")
        }

    ##### Mathematics, Statistical Testing, Modeling, Signal Processing

        #### drawBaseline

            #' Draws a baseline on a signal
            #'
            #' @param data x, y dataframe
            #' @param prelim_baseline_window interval at which to draw baseline
            #' @examples
            #' @export
            #' drawBaseline

            drawBaseline <- function(data, prelim_baseline_window) {
      
                n_prelim_baseline_windows <- floor(length(data$x)/prelim_baseline_window)
                prelim_baseline <- list()
                for ( i in 1:n_prelim_baseline_windows ) {
                  abundances_in_window <- data$y[((prelim_baseline_window*(i-1))+1):(prelim_baseline_window*i)]
                  prelim_baseline[[i]] <- data.frame(
                      x = data$x[(which.min(abundances_in_window)+((i-1)*prelim_baseline_window))],
                      min = min(abundances_in_window)
                  )
                }
                prelim_baseline <- do.call(rbind, prelim_baseline)
                data$in_prelim_baseline <- FALSE
                data$in_prelim_baseline[data$x %in% prelim_baseline$x] <- TRUE

                y = prelim_baseline$min
                x = prelim_baseline$x

                baseline2 <- data.frame(
                  x = data$x,
                  y = approx(x, y, xout = data$x)$y
                )
                baseline2 <- baseline2[!is.na(baseline2$y),]
                data <- data[data$x %in% baseline2$x,]
                data$baseline <- baseline2$y

                return(data)

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

            normalize <- function( x, old_min = NULL, old_max = NULL, new_min = 0, new_max = 1, na_zero = FALSE ) {

                if ( length(old_min) == 0 & length(old_max) == 0 ) {
                    
                    x <- (x - min(x)) * (new_max - new_min) / (max(x) - min(x)) + new_min
                
                } else {

                    x <- (x - old_min) * (new_max - new_min) / (old_max - old_min) + new_min

                }

                if (na_zero == TRUE) {
                    x[is.na(x)] <- 0
                }

                return(x)

            }

        #### movingAverage

            #' Compute a moving average
            #'
            #' movingAverage

            movingAverage <- function(x, n = 5){
            
                stats::filter(x, rep(1 / n, n), sides = 2)
            
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
                                                "pca", "pca_ord", "pca_dim",
                                                "mca", "mca_ord", "mca_dim",
                                                "mds", "mds_ord", "mds_dim",
                                                "tsne", "dbscan", "kmeans",
                                                "hclust", "hclust_phylo"
                                            ),
                                            parameters = NULL,
                                            column_w_names_of_multiple_analytes = NULL,
                                            column_w_values_for_multiple_analytes = NULL,
                                            columns_w_values_for_single_analyte = NULL,
                                            columns_w_additional_analyte_info = NULL,
                                            columns_w_sample_ID_info = NULL,
                                            transpose = FALSE,
                                            distance_method = c("euclidean", "coeff_unlike"),
                                            unknown_sample_ID_info = NULL,
                                            scale_variance = TRUE,
                                            na_replacement = c("mean", "none", "zero", "drop"),
                                            output_format = c("wide", "long"),
                                            ...
                                        ) {

                    # Pre-process data

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
                                    "unknown_sample_ID_info",
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
                                    
                                    cat(paste0("Replacing NAs in your data with ", na_replacement[1]), "\n")

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

                        ## Scale data, unless not requested

                            if( scale_variance == TRUE & !analysis %in% c("mca", "mca_ord", "mca_dim")) {
                                
                                scaled_matrix <- scale(matrix)

                                # for (column in 1:length(matrix)) {
                                #     scaled_matrix <- matrix
                                #     scaled_matrix[,column] <- normalize(scaled_matrix[,column], new_min = -1, na_zero = TRUE)
                                # }

                                if( any(is.na(scaled_matrix)) ) {
                                    cat("Some analytes have zero variance and will be assigned a value of zero in the scaled matrix.")
                                    scaled_matrix[is.na(scaled_matrix)] <- 0
                                }

                            }

                            if( scale_variance == FALSE ) {
                                scaled_matrix <- matrix
                            }

                        ## Get distance matrix

                            if(  !analysis %in% c("mca", "mca_ord", "mca_dim") ) {

                                if( distance_method[1] == "euclidean") {
                                    dist_matrix <- stats::dist(scaled_matrix, method = "euclidean")
                                } else {
                                    stop("Please specify distance method")
                                }

                                if( analysis == "dist") {
                                    return(dist_matrix)
                                    stop()
                                }

                            }

                            if( analysis %in% c("mca", "mca_ord", "mca_dim") ) {
                            
                                scaled_matrix <- matrix
                            
                            }

                    # Run the matrix analysis selected

                        ## Dimensionality reduction

                            ## MDS

                                if( analysis == "mds" ) {
                                    coords <- stats::cmdscale(dist_matrix)
                                    colnames(coords) <- c("Dim.1", "Dim.2")
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
                                    cat("Running Multiple Correspondence Analysis, extracting sample coordinates...\n")
                                    coords <- FactoMineR::MCA(matrix, graph = FALSE)$ind$coord[,c(1:2)]
                                    clustering <- as_tibble(coords)
                                    clustering$sample_unique_ID <- rownames(coords)
                                    colnames(clustering) <- c("Dim_1", "Dim_2", "sample_unique_ID")
                                    cat("Done!\n")
                                }

                                if( analysis == "mca_ord" ) {
                                    cat("Running Multiple Correspondence Analysis, extracting ordination plot...\n")
                                    coords <- FactoMineR::MCA(matrix, graph = FALSE)$var$eta2[,c(1,2)]
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
                                    coords <- FactoMineR::PCA(scaled_matrix, graph = FALSE, scale.unit = FALSE)$ind$coord[,c(1:2)]
                                    clustering <- as_tibble(coords)
                                    clustering$sample_unique_ID <- rownames(coords)
                                }

                                if( analysis == "pca_ord" ) {
                                    coords <- FactoMineR::PCA(scaled_matrix, graph = FALSE, scale.unit = FALSE)$var$coord[,c(1,2)]
                                    clustering <- as_tibble(coords)
                                    clustering$analyte <- rownames(coords)
                                    clustering <- select(clustering, analyte, Dim.1, Dim.2)
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

                            ## DBSCAN

                                if( analysis == "dbscan" ) {

                                    if ( length(parameters) > 0 ) {
                                        cluster_k <- parameters[1]
                                        cluster_threshold <- parameters[2]
                                    }

                                    if ( length(parameters) == 0 ) {
                                        findClusterParameters(dist_matrix = dist_matrix, matrix = matrix, analysis = "dbscan")
                                    }

                                    cat("Using", cluster_k, "as a value for k.\n")
                                    cat("Using", cluster_threshold, "as a value for threshold.\n")
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

                                    cat("Using", n_clusters, "as a value for cluster_number.\n")
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
                                }

                                return( clustering )
                        }
                }

        #### findClusterParameters

            #' Interactively find clustering parameters from a distance matrix and source matrix
            #'
            #' @param dist_matrix
            #' @param matrix
            #' @examples
            #' @export
            #' findClusterParameters

            findClusterParameters <- function(dist_matrix, matrix, analysis = c("dbscan", "kmeans")) {

                ui <-   navbarPage(title = "Cluster Parameter Selection", id = "navbar",
                            
                            tabPanel(title = "Main",

                                sidebarLayout(

                                    sidebarPanel(

                                        conditionalPanel( "input.selection == 'kmeans'",
                                        
                                            sliderInput(
                                                inputId = "n_clusters",
                                                label = "Number of clusters",
                                                min = 0, max = 25, value = 1, step = 1
                                            )
                                        
                                        ),
                                
                                        conditionalPanel( "input.selection == 'dbscan'",

                                            sliderInput(
                                                inputId = "cluster_k", label = "k",
                                                min = 1, max = 40, value = 2
                                            ),

                                            sliderInput(
                                                inputId = "cluster_threshold", label = "threshold",
                                                min = 0, max = 50, value = 2, step = 0.05
                                            )

                                        ),

                                        plotOutput(
                                            outputId = "plot1",
                                            width = 200, height = 200
                                        )
                                    ),

                                    mainPanel(
                                            plotOutput(
                                                outputId = "plot2", width = "100%"
                                                # height = 800
                                            )
                                    )
                                )

                            ),

                            tabPanel(title = "Quit", value = "stop", icon = icon("circle-notch"),
                            
                                selectInput(
                                    inputId = "selection",
                                    label = "Do not touch this",
                                    choices = c("dbscan", "kmeans"),
                                    selected = analysis
                                )

                            )

                )
                
                server <- function(input, output, session) {

                    output$plot1 <- renderPlot({

                        n_clusters <<- input$n_clusters
                        cluster_k <<- input$cluster_k
                        cluster_threshold <<- input$cluster_threshold

                        if (analysis == "dbscan") {
                            distances <- dbscan::kNNdist(dist_matrix, k = input$cluster_k)
                        }

                        if (analysis == "kmeans") {
                            kmeans_results <- list()
                            for( i in 1:(dim(matrix)[1]-1) ) {
                                kmeans_results[[i]] <- sum(stats::kmeans(x = matrix, centers = i, nstart = 25, iter.max = 1000)$withinss)
                            }
                            distances <- do.call(rbind, kmeans_results)
                        }

                        distances <- as.numeric(distances)[rev(order(as.numeric(distances)))]

                        ggplot(
                            data.frame(
                                index = seq(1,length(distances),1),
                                distance = distances
                            ), aes(x = index, y = distance)
                        ) + geom_point() + theme_bw() +
                        geom_hline(yintercept = input$cluster_threshold)

                    })

                    output$plot2 <- renderPlot ({

                        if (analysis == "dbscan") {
                            clustering <- as_tibble(data.frame(
                                sample_unique_ID = colnames(as.matrix(dist_matrix)),
                                cluster = paste0("cluster_", fpc::dbscan(dist_matrix, eps = input$cluster_threshold , MinPts = input$cluster_k, scale = FALSE, method = "dist")[[1]])
                            ))
                            clustering$cluster[clustering$cluster == "cluster_0"] <- NA
                        }

                        if (analysis == "kmeans") {
                            kmeans_clusters <- stats::kmeans(x = matrix, centers = input$n_clusters, nstart = 25, iter.max = 1000)$cluster
                            clustering <- as_tibble(data.frame(sample_unique_ID = names(kmeans_clusters), cluster = paste0("cluster_", kmeans_clusters)))
                        }

                        clustering <- cbind(
                            clustering,
                            matrix[match(rownames(matrix), clustering$sample_unique_ID),]
                        )

                        ggplot(clustering, aes_string(x = colnames(matrix)[1], y = colnames(matrix)[2], fill = "cluster")) +
                            geom_point(shape = 21, size = 3, alpha = 0.6) +
                            theme_bw() +
                            scale_fill_manual(values = discrete_palette) +
                            coord_fixed()

                    })

                    session$onSessionEnded(function() { stopApp() })

                    observe({ if (input$navbar == "stop") stopApp() })

                }

                runApp(shinyApp(ui, server))

            }

        #### pGroups

            #' Determine which samples are identical based on p values
            #'
            #' @param data 
            #' @examples
            #' @export
            #' pGroups

            pGroups <- function(data) {
                
                p <- data$p.adj
                names(p) <- paste(data$group1, data$group2, sep = "-")

                if (is.null(multcompLetters(p)$monospacedLetters)) {
                    spaced_group <- multcompLetters(p)$Letters 
                } else {
                    spaced_group <- multcompLetters(p)$monospacedLetters   
                }

                output <- data.frame(
                    treatment = names(multcompLetters(p)$Letters),
                    group = multcompLetters(p)$Letters,
                    spaced_group = spaced_group
                )
                return(output)
            }

        #### p_groups

            #' Determine which samples are identical based on p values
            #'
            #' @param data 
            #' @examples
            #' @export
            #' p_groups

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

        #### fitGaussians

            #' Fit a Gaussian curve to a set of x y data
            #'
            #' This function may work better if you feed it a smoothed signal, or add smoothing to the derivatives
            #' @param x
            #' @param y
            #' @param smoothing_order
            #' @param smoothing_window
            #' @export
            #' @examples
            #' fitGaussians()

            fitGaussians <- function(
                                x, 
                                y, 
                                peak_detection_threshold = 10,
                                peak_means_stdevs_heights = NULL
                            ) {

                ## Find peaks, if peak number is not provided.

                    if (is.null(peak_means_stdevs_heights)) {

                        ## Find "change points" maxima and minima in first derivative = peak half heights
                            deriv1_y <- as.vector(stats::predict(pspline::sm.spline(x, y), x, 1))
                            deriv1_extreme_indices <- c(cumsum(rle(abs(diff(deriv1_y)) > median(diff(deriv1_y))*2)$lengths)+1) # This from https://stackoverflow.com/questions/46029090/how-to-find-changing-points-in-a-dataset
                            deriv1_extremes <- x[deriv1_extreme_indices]

                        ## Find peaks using second derivative - better at detecting shoulder peaks v. 1st derivative
                            deriv2_y <- stats::predict(pspline::sm.spline(x, deriv1_y), x, 1)

                            peak_means <- list()
                            for (i in peak_detection_threshold:(length(x)-peak_detection_threshold)) {
                                if ( # Find all points in 2nd derivative where values on either side of it are increasing for n points in a row
                                    all(
                                        all(deriv2_y[(i-(peak_detection_threshold-1)):(i-1)] > deriv2_y[i]), 
                                        all(deriv2_y[(i+1):(i+peak_detection_threshold)] > deriv2_y[i])
                                    )
                                ) {
                                    peak_means[[i]] <- x[i]
                                }
                            }
                            peak_means <- unlist(peak_means)

                        ## Find sigma for each peak based on mean and "change points"
                            peak_sigma <- vector()
                            for (peak in 1:length(peak_means)) {
                                peak_sigma[[peak]] <- 1.5*min(abs(deriv1_extremes - peak_means[peak]))
                            }
                            
                        ## Peak heights are just the height at the mean
                            peak_heights <- y[x %in% peak_means]

                        ## Report on peaks found and plot them
                            cat(paste0("Found ", length(peak_heights), " peaks at peak_detection_threshold = ", peak_detection_threshold, "\n"))

                            plot <- ggplot() +
                                geom_line(data = data.frame(x = x, y = y), aes(x = x, y = y), color = "red") +
                                geom_line(data = data.frame(x = x, y = y), aes(x = x, y = deriv1_y), color = "black") +
                                geom_line(data = data.frame(x = x, y = y), aes(x = x, y = deriv2_y), color = "blue") +
                                geom_vline(aes(xintercept = peak_means)) +
                                geom_vline(aes(xintercept = peak_means + peak_sigma), color = "red") +
                                geom_vline(aes(xintercept = peak_means - peak_sigma), color = "red")
                            print(plot)

                    }

                ## If peak_means_stdevs_heights provided, then use those

                    if (!is.null(peak_means_stdevs_heights)) {
                        n_peaks <- length(peak_means_stdevs_heights)/3
                        peak_means <- peak_means_stdevs_heights[1:n_peaks]
                        peak_sigma <- peak_means_stdevs_heights[(n_peaks+1):(n_peaks*2)]
                        peak_heights <- peak_means_stdevs_heights[((n_peaks*2)+1):(n_peaks*3)]
                    }

                ## Fit peaks with gaussians using nls or manually and extract fit data

                    algorithm = "port"
                    fit <- NULL

                    if (length(peak_means) == 1) {
                        try(
                            fit <-  nls(    y ~ (
                                            C1*exp(-(x-mean1)**2/(2 * sigma1**2))
                                        ),
                                    data = data.frame(x = x, y = y),
                                    start = list(
                                                mean1 = peak_means[1],
                                                sigma1 = peak_sigma[1],
                                                C1 = peak_heights[1]
                                            ),
                                    algorithm = algorithm, control = list(maxiter = 5000000))
                        , silent = TRUE)
                    }

                    if (length(peak_means) == 2) {
                        try(
                            fit <- nls( y ~ (
                                            C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                                            C2*exp(-(x-mean2)**2/(2 * sigma2**2))
                                        ),
                                data = data.frame(x = x, y = y),
                                start = list(
                                            mean1 = peak_means[1],
                                            mean2 = peak_means[2],
                                            sigma1 = peak_sigma[1],
                                            sigma2 = peak_sigma[2],
                                            C1 = peak_heights[1],
                                            C2 = peak_heights[2]
                                        ),
                                algorithm = algorithm, control = list(maxiter = 500000000))
                    , silent = TRUE)
                    }
                    
                    if (length(peak_means) == 3) {
                        try(
                            fit <- nls( y ~ (
                                            C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                                            C2*exp(-(x-mean2)**2/(2 * sigma2**2)) +
                                            C3*exp(-(x-mean3)**2/(2 * sigma3**2))
                                        ),
                                    data = data.frame(x = x, y = y),
                                    start = list(
                                                mean1 = peak_means[1],
                                                mean2 = peak_means[2],
                                                mean3 = peak_means[3],
                                                sigma1 = peak_sigma[1],
                                                sigma2 = peak_sigma[2],
                                                sigma3 = peak_sigma[3],
                                                C1 = peak_heights[1],
                                                C2 = peak_heights[2],
                                                C3 = peak_heights[3]
                                            ),
                                    algorithm = algorithm, control = list(maxiter = 5000000))
                        , silent = TRUE)
                    }

                    if (length(peak_means) == 4) {
                        try(
                            fit <-  nls(    y ~ (
                                            C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                                            C2*exp(-(x-mean2)**2/(2 * sigma2**2)) +
                                            C3*exp(-(x-mean3)**2/(2 * sigma3**2)) +
                                            C4*exp(-(x-mean4)**2/(2 * sigma4**2))
                                        ),
                                    data = data.frame(x = x, y = y),
                                    start = list(
                                                mean1 = peak_means[1],
                                                mean2 = peak_means[2],
                                                mean3 = peak_means[3],
                                                mean4 = peak_means[4],
                                                sigma1 = peak_sigma[1],
                                                sigma2 = peak_sigma[2],
                                                sigma3 = peak_sigma[3],
                                                sigma4 = peak_sigma[4],
                                                C1 = peak_heights[1],
                                                C2 = peak_heights[2],
                                                C3 = peak_heights[3],
                                                C4 = peak_heights[4]
                                            ),
                                    algorithm = algorithm, control = list(maxiter = 5000000))
                        , silent = TRUE)
                    }

                    if (length(peak_means) == 5) {
                        try(
                            fit <- nls( y ~ (
                                            C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                                            C2*exp(-(x-mean2)**2/(2 * sigma2**2)) +
                                            C3*exp(-(x-mean3)**2/(2 * sigma3**2)) +
                                            C4*exp(-(x-mean4)**2/(2 * sigma4**2)) +
                                            C5*exp(-(x-mean5)**2/(2 * sigma5**2))
                                        ),
                                    data = data.frame(x = x, y = y),
                                    start = list(
                                                mean1 = peak_means[1],
                                                mean2 = peak_means[2],
                                                mean3 = peak_means[3],
                                                mean4 = peak_means[4],
                                                mean5 = peak_means[5],
                                                sigma1 = peak_sigma[1],
                                                sigma2 = peak_sigma[2],
                                                sigma3 = peak_sigma[3],
                                                sigma4 = peak_sigma[4],
                                                sigma5 = peak_sigma[5],
                                                C1 = peak_heights[1],
                                                C2 = peak_heights[2],
                                                C3 = peak_heights[3],
                                                C4 = peak_heights[4],
                                                C5 = peak_heights[5]
                                            ),
                                    algorithm = algorithm, control = list(maxiter = 5000000))
                        , silent = TRUE)
                    }

                    if (length(peak_means) == 6) {
                        try(
                            fit <- nls( y ~ (
                                            C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                                            C2*exp(-(x-mean2)**2/(2 * sigma2**2)) +
                                            C3*exp(-(x-mean3)**2/(2 * sigma3**2)) +
                                            C4*exp(-(x-mean4)**2/(2 * sigma4**2)) +
                                            C5*exp(-(x-mean5)**2/(2 * sigma5**2)) +
                                            C6*exp(-(x-mean6)**2/(2 * sigma6**2))
                                        ),
                                    data = data.frame(x = x, y = y),
                                    start = list(
                                                mean1 = peak_means[1],
                                                mean2 = peak_means[2],
                                                mean3 = peak_means[3],
                                                mean4 = peak_means[4],
                                                mean5 = peak_means[5],
                                                mean6 = peak_means[6],
                                                sigma1 = peak_sigma[1],
                                                sigma2 = peak_sigma[2],
                                                sigma3 = peak_sigma[3],
                                                sigma4 = peak_sigma[4],
                                                sigma5 = peak_sigma[5],
                                                sigma6 = peak_sigma[6],
                                                C1 = peak_heights[1],
                                                C2 = peak_heights[2],
                                                C3 = peak_heights[3],
                                                C4 = peak_heights[4],
                                                C5 = peak_heights[5],
                                                C6 = peak_heights[6]
                                            ),
                                    algorithm = algorithm, control = list(maxiter = 10000000))
                        , silent = TRUE)
                    }

                    if (length(peak_means) == 7) {
                        try(
                            fit <- nls( y ~ (
                                            C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                                            C2*exp(-(x-mean2)**2/(2 * sigma2**2)) +
                                            C3*exp(-(x-mean3)**2/(2 * sigma3**2)) +
                                            C4*exp(-(x-mean4)**2/(2 * sigma4**2)) +
                                            C5*exp(-(x-mean5)**2/(2 * sigma5**2)) +
                                            C6*exp(-(x-mean6)**2/(2 * sigma6**2)) +
                                            C7*exp(-(x-mean7)**2/(2 * sigma7**2))
                                        ),
                                    data = data.frame(x = x, y = y),
                                    start = list(
                                                mean1 = peak_means[1],
                                                mean2 = peak_means[2],
                                                mean3 = peak_means[3],
                                                mean4 = peak_means[4],
                                                mean5 = peak_means[5],
                                                mean6 = peak_means[6],
                                                mean7 = peak_means[7],
                                                sigma1 = peak_sigma[1],
                                                sigma2 = peak_sigma[2],
                                                sigma3 = peak_sigma[3],
                                                sigma4 = peak_sigma[4],
                                                sigma5 = peak_sigma[5],
                                                sigma6 = peak_sigma[6],
                                                sigma7 = peak_sigma[7],
                                                C1 = peak_heights[1],
                                                C2 = peak_heights[2],
                                                C3 = peak_heights[3],
                                                C4 = peak_heights[4],
                                                C5 = peak_heights[5],
                                                C6 = peak_heights[6],
                                                C7 = peak_heights[7]
                                            ),
                                    algorithm = algorithm, control = list(maxiter = 10000000))
                        , silent = TRUE)
                    }

                    if ( !is.null(fit) ) {
                        
                        cat(paste0("Fit ", length(peak_means), " peak(s) with gaussians at this peak detection threshold.\n\n"))
                        dffit <- data.frame(x = x)
                        dffit$y <- predict(fit, newdata = dffit)
                        fit.sum <- summary(fit)
                        coef.fit <- fit.sum$coefficients[,1]
                        mean.fit <- coef.fit[1:length(peak_means)]
                        sigma.fit <- coef.fit[(length(peak_means)+1):(2*length(peak_means))]
                        C.fit <- coef.fit[((2*length(peak_means))+1):(3*length(peak_means))]
                        
                    } else {
                        
                        cat(paste0("Couldn't fit gaussian(s) at this peak detection threshold\n"))
                        gaussians <- NULL
                        cat(paste0("Trying to fit manually using least squares starting from initial values\n"))

                        # ## Generate gaussian curves, using least squares to optimize

                            mean.fit <- peak_means
                            sigma.fit <- peak_sigma
                            C.fit <- peak_heights

                            init_step <- 0.9
                            step_size <- 0.05
                            final_step <- 1.1

                            grid <- expand.grid(
                                x = seq(sigma.fit[1]*init_step, sigma.fit[1]*final_step, sigma.fit[1]*step_size),
                                y = seq(sigma.fit[2]*init_step, sigma.fit[2]*final_step, sigma.fit[2]*step_size),
                                z = seq(sigma.fit[3]*init_step, sigma.fit[3]*final_step, sigma.fit[3]*step_size),
                                l = seq(mean.fit[1]*init_step, mean.fit[1]*final_step, mean.fit[1]*step_size),
                                m = seq(mean.fit[2]*init_step, mean.fit[2]*final_step, mean.fit[2]*step_size),
                                n = seq(mean.fit[3]*init_step, mean.fit[3]*final_step, mean.fit[3]*step_size),
                                a = seq(C.fit[1]*init_step, C.fit[1]*final_step, C.fit[1]*step_size),
                                b = seq(C.fit[2]*init_step, C.fit[2]*final_step, C.fit[2]*step_size),
                                c = seq(C.fit[3]*init_step, C.fit[3]*final_step, C.fit[3]*step_size)
                            )

                        #     # grid <- expand.grid(
                        #     #     x = seq(sigma.fit[1]*0.999, sigma.fit[1]*1.001, sigma.fit[1]*0.001),
                        #     #     y = seq(sigma.fit[2]*0.999, sigma.fit[2]*1.001, sigma.fit[2]*0.001),
                        #     #     z = seq(sigma.fit[3]*0.999, sigma.fit[3]*1.001, sigma.fit[3]*0.001),
                        #     #     l = seq(mean.fit[1]*0.999, mean.fit[1]*1.001, mean.fit[1]*0.001),
                        #     #     m = seq(mean.fit[2]*0.999, mean.fit[2]*1.001, mean.fit[2]*0.001),
                        #     #     n = seq(mean.fit[3]*0.999, mean.fit[3]*1.001, mean.fit[3]*0.001)
                        #     # )

                        #     # # grid <- expand.grid(
                        #     # #     x = seq(sigma.fit[1]*0.9999, sigma.fit[1]*1.0001, sigma.fit[1]*0.0001),
                        #     # #     y = seq(sigma.fit[2]*0.9999, sigma.fit[2]*1.0001, sigma.fit[2]*0.0001),
                        #     # #     z = seq(sigma.fit[3]*0.9999, sigma.fit[3]*1.0001, sigma.fit[3]*0.0001),
                        #     # #     l = seq(mean.fit[1]*0.9999, mean.fit[1]*1.0001, mean.fit[1]*0.0001),
                        #     # #     m = seq(mean.fit[2]*0.9999, mean.fit[2]*1.0001, mean.fit[2]*0.0001),
                        #     # #     n = seq(mean.fit[3]*0.9999, mean.fit[3]*1.0001, mean.fit[3]*0.0001)
                        #     # # )

                            out <- list()
                            pb <- progress::progress_bar$new(total = dim(grid)[1])
                            for (i in 1:dim(grid)[1]) {
                                pb$tick()
                                gaussians <- list()
                                for (peak in 1:length(peak_means)) {
                                    peak_data <-    data.frame(
                                                        x = x,
                                                        # y = C.fit[peak] * exp(-((x)-mean.fit[peak])**2/(2 * sigma.fit[peak]**2)),
                                                        # y = C.fit[peak] * exp(-((x)-grid[i,peak])**2/(2 * sigma.fit[peak]**2)),
                                                        # y = C.fit[peak] * exp(-((x)-mean.fit[peak])**2/(2 * grid[i,peak]**2)),
                                                        y = grid[i,peak+6] * exp(-((x)-grid[i,peak+3])**2/(2 * grid[i,peak]**2)),
                                                        peak_number = peak
                                                    )
                                    peak_data$peak_area <- sum(peak_data$y)
                                    gaussians[[peak]] <- peak_data
                                }
                                gaussians <- do.call(rbind, gaussians)

                                test_set <- gaussians
                                test_set$ref <- y
                                test_set %>%
                                    group_by(x) %>%
                                    summarize(sum = abs(sum(y) - ref)) %>%
                                    ungroup() %>%
                                    select(sum) %>%
                                    summarize(total = sum(sum)) -> qual

                                out[[i]] <- qual

                            }

                            totals <- as.numeric(as.data.frame(do.call(rbind, out))[,1])

                            out <- data.frame(
                                x = seq(1,length(totals), 1),
                                y = totals
                            )

                            fit_data <- as.numeric(grid[which.min(out$y),])
                            sigma.fit <- fit_data[1:n_peaks]
                            mean.fit <- fit_data[(n_peaks+1):(n_peaks*2)]
                            C.fit <- fit_data[((n_peaks*2)+1):(n_peaks*3)]
                    }

                ## Generate gaussian curves

                    gaussians <- list()
                    for (peak in 1:length(peak_means)) {
                        peak_data <-    data.frame(
                                            x = x,
                                            y = C.fit[peak] * exp(-((x)-mean.fit[peak])**2/(2 * sigma.fit[peak]**2)),
                                            peak_number = peak
                                        )
                        peak_data$peak_area <- sum(peak_data$y)
                        gaussians[[peak]] <- peak_data
                    }
                    gaussians <- do.call(rbind, gaussians)

                ## Plot the results

                    plot <- ggplot() +
                        geom_line(data = data.frame(x = x, y = y), aes(x = x, y = y), color = "red", size = 2) +
                        geom_vline(aes(xintercept = peak_means)) +
                        stat_summary(data = gaussians, geom = "line", fun = "sum", color = "blue", aes(x = x, y = y), size = 2) +
                        geom_line(data = gaussians, aes(x = x, y = y, group = peak_number))

                    print(plot)

                ## Return gaussians and fit parameters

                    output <- list()

                    output$fit_data <- data.frame(
                        peak_height = C.fit,
                        peak_mean = mean.fit,
                        peak_sd = sigma.fit,
                        peak_number = seq(1,length(C.fit),1)
                    )
                    output$gaussians <- gaussians
                    
                    return(output)

            }

        #### buildLinearModel

            #' buildLinearModel
            #'
            #' @param data
            #' @param formula
            #' @export
            #' @examples
            #' buildLinearModel()

                buildLinearModel <- function(data, formula) {

                    ## Correct the formula and the data

                        data <- as.data.frame(data)
                        formula <- gsub("=", "~", formula)
                        formula <- gsub("~", " ~ ", formula)
                        formula <- gsub("  ", " ", formula)

                    ## Run the fit and start the output

                        fit <- lm(
                            data = data,
                            formula = formula,
                            x = TRUE, y = TRUE, model = TRUE, qr = TRUE
                        )
                        output <- fit$model
                        output$residuals <- fit$residuals
                        output$model_y <- fit$model[,1] - output$residuals

                    ## Determine input x values

                        eqn_side <- gsub(".*~ ", "", formula)
                        for (j in 1:length(colnames(data))) {
                            result <- grep(paste0(" ", colnames(data)[j]), paste0(" ", eqn_side), fixed = TRUE)
                            if (length(result) > 0) {
                                # print(j)
                                index <- gregexpr(pattern = paste0(" ", colnames(data)[j]),paste0(" ", eqn_side))[[1]][1]
                                eqn_side <- paste0(
                                    substr(eqn_side, 0, index-1),
                                    "data$",
                                    substr(eqn_side, index, nchar(eqn_side))
                                )
                            }
                        }
                        input_x <- eval(parse(text=eqn_side))

                    ## Determine input y values

                        eqn_side <- gsub(" ~ .*$", "", formula)
                        for (j in 1:length(colnames(data))) {
                            result <- grep(paste0(" ", colnames(data)[j]), paste0(" ", eqn_side), fixed = TRUE)
                            if (length(result) > 0) {
                                # print(j)
                                index <- gregexpr(pattern = paste0(" ", colnames(data)[j]),paste0(" ", eqn_side))[[1]][1]
                                eqn_side <- paste0(
                                    substr(eqn_side, 0, index-1),
                                    "data$",
                                    substr(eqn_side, index, nchar(eqn_side))
                                )
                            }
                        }
                        input_y <- eval(parse(text=eqn_side))

                        output <- cbind(
                            cbind(input_x, input_y)[!apply(apply(cbind(input_x, input_y), 2, is.na), 1, any),],
                            output
                        )
                        # colnames(output)[1]

                    ## Start the master output and put the coefficients and stats in it

                        out <- list()
                        
                        coefficients <- data.frame(
                            variable = names(coefficients(fit)),
                            value = as.numeric(coefficients(fit)),
                            std_err = round(as.numeric(summary(fit)$coefficients[,2]), 4),
                            type = "coefficient",
                            p_value = round(summary(fit)$coefficients[,4], 4)
                        )

                        total_sum_squares <- sum((input_y - mean(input_y, na.rm = TRUE))^2, na.rm = TRUE)
                        residual_sum_squares <- sum((summary(fit)$residuals)^2, na.rm = TRUE)

                        statistics <- rbind(
                            c("median_residual", median(summary(fit)$residuals), NA, "statistic", NA),
                            c("total_sum_squares", total_sum_squares, NA, "statistic", NA),
                            c("residual_sum_squares", residual_sum_squares, NA, "statistic", NA),
                            c("r_squared", summary(fit)$r.squared, NA, "statistic", NA)
                        )
                        colnames(statistics) <- c("variable", "value", "std_err", "type", "p_value")

                        # round(1-(residual_sum_squares/total_sum_squares), 4) == round(summary(fit)$r.squared, 4)

                        metrics <- rbind(coefficients, statistics)
                        metrics$value <- round(as.numeric(metrics$value), 4)
                        rownames(metrics) <- NULL
                        metrics

                        out$metrics <- metrics

                        output$model_x <- output$input_x

                        out$data <- output

                    return(out)
                }

        #### mode

            #' mode
            #'
            #' @param x
            #' @param ignore_zero
            #' @export
            #' @examples
            #' mode()

            mode <- function(x, ignore_zero = FALSE) {
                
                x <- as.data.frame(table(x))
                x <- x[order(x$Freq, decreasing = TRUE),]
                if (ignore_zero) {
                  x <- filter(x, x != 0)
                }
                
                ## If it can be numeric, return numeric, otherwise, return character
                    if(any(is.na(as.numeric(as.character(x$x))))) {
                        return(as.character(x$x[1]))
                    } else {
                        return(as.numeric(as.character(x$x[1])))
                    }
            }

        #### geomSignif

            #' geomSignif
            #'
            #' @param data
            #' @param orientation
            #' @export
            #' @examples
            #' geomSignif()

            geomSignif <- function(data, orientation) {

                if (orientation == "horizontal") {
                    return(
                        list(
                            geom_segment(
                                data = data,
                                aes(x = xmin, xend = xmax, y = y_position, yend = y_position),
                                inherit.aes = FALSE
                            ),
                            geom_segment(
                                data = data,
                                aes(x = xmin, xend = xmin, y = y_position, yend = y_position - tip_length_xmin),
                                inherit.aes = FALSE
                            ),
                            geom_segment(
                                data = data,
                                aes(x = xmax, xend = xmax, y = y_position, yend = y_position - tip_length_xmax),
                                inherit.aes = FALSE
                            ),
                            geom_text(
                                data = data,
                                aes(
                                    x = text_horiz_offset, y = y_position+text_vert_offset,
                                    label = text, size = text_size, hjust = hjust, vjust = vjust
                                ),
                                show.legend = FALSE, inherit.aes = FALSE
                            )
                        )
                    )
                }

                if (orientation == "vertical") {
                    return(
                        list(
                            geom_segment(
                                data = data,
                                aes(x = x_position, xend = x_position, y = ymin, yend = ymax),
                                inherit.aes = FALSE
                            ),
                            geom_segment(
                                data = data,
                                aes(x = x_position, xend = x_position - tip_length_ymin, y = ymin, yend = ymin),
                                inherit.aes = FALSE
                            ),
                            geom_segment(
                                data = data,
                                aes(x = x_position, xend = x_position - tip_length_ymax, y = ymax, yend = ymax),
                                inherit.aes = FALSE
                            ),
                            geom_text(
                                data = data,
                                aes(
                                    x = x_position+text_horiz_offset, y = text_vert_offset,
                                    label = text, size = text_size, hjust = hjust, vjust = vjust
                                ),
                                show.legend = FALSE, inherit.aes = FALSE
                            )
                        )
                    )
                }
                
            }

        #### grubbsFilter

            #' grubbsFilter on grouped data
            #'
            #' @param data
            #' @param threshold
            #' @export
            #' @examples
            #' grubbsFilter()

            grubbsFilter <- function(data, col, threshold = 0.05, ...) {

                if (dplyr::is_grouped_df(data)) {
                  return(dplyr::do(data, grubbsFilter(., col, ...))) 
                }

                test_vector <- as.data.frame(data[,which(colnames(data) == col)])[[1]]

                if (length(test_vector) > 3) {
                    result <- outliers::grubbs.test(test_vector)

                    if (result$p.value < threshold) {
                        # print(result)
                        tail <- gsub(" .*$", "", result$alternative)
                        if (tail == "lowest") { data <- data[-which.min(test_vector),] }
                        if (tail == "highest") { data <- data[-which.max(test_vector),] }
                    }
                }

                return(data)

                # data.frame(
                #     id = c(rep("A", 5), rep("B", 5)),
                #     val = c(-10,2,3,4,5,1,2,3,4,15)
                # ) %>%
                # as_tibble() %>%
                # group_by(id) -> data

                # grubbsFilter(data = data, col = "val")

            }

    ##### Phylogenetic statistical testing

        #### phylogeneticSignal

            #' Compute phylogenetic signal for continuous and discrete traits
            #'
            #' @param trait A data frame of traits where the first column is tree tip names that exactly match the tree tip labels
            #' @param tree A phylogenetic tree (class phylo) with tips that exactly match the names of trait
            #' @param replicates Number of random replications to run
            #' @param cost Optional. A specialized transition matrix
            #' @examples
            #' @export
            #' phylogeneticSignal

                phylogeneticSignal <- function( traits, tree, replicates = 999, cost = NULL ) {

                    results <- list()

                    ## Loop over columns
                    
                        for (i in 2:dim(traits)[2]) {

                            trait <- traits[,i]
                            names(trait) <- traits[,1]

                            if (class(trait) %in% c("factor", "character")) {

                                ### For discrete traits

                                    ## Get the states in which the trait may exist (levels)

                                        levels <- attributes(factor(trait))$levels

                                    ## Chech that the variable is indeed categorical
                                                    
                                        if (length(levels) == length(trait)) {
                                        
                                            warn("Are you sure this variable is categorical?")

                                        }

                                    ## Make the transition matrix
                                        
                                        if (is.null(cost)) {
                                        
                                            cost1 <- 1-diag(length(levels))
                                        
                                        } else {
                                        
                                        if (length(levels) != dim(cost)[1])
                                            
                                            stop("Dimensions of the character state transition matrix do not agree with the number of levels")
                                            
                                            cost1 <- t(cost)
                                        }
                                        dimnames(cost1) <- list(levels, levels)

                                    ## Make the trait numeric
                                    
                                        trait_as_numeric <- as.numeric(as.factor(trait))
                                        names(trait_as_numeric) <- names(trait)
                                    
                                    ## Make the phyDat object and get the parsimony score for the tree with the associated observations

                                        # obs <- t(data.frame(trait))
                                        obs <- phyDat( trait, type = "USER", levels = attributes(factor(trait))$levels )
                                        OBS <- parsimony( tree, obs, method = "sankoff", cost = cost1 )

                                    ## Make "replicates" number of random tree-trait associations and check their parsimony score

                                        null_model <- matrix(NA, replicates, 1)
                                        for (i in 1:replicates){

                                            ## Randomize the traits and get the parsimony score for the random traits on that tree
                                                null <- sample(as.numeric(trait_as_numeric))
                                                attributes(null)$names <- attributes(trait_as_numeric)$names
                                                # null <- t(data.frame(null))
                                                null <- phyDat( null,type = "USER",levels = attributes(factor(null))$levels )
                                                null_model[i,] <- parsimony( tree, null, method = "sankoff", cost = cost1 )

                                        }

                                    ## Assess observed parsimony score in the context of the random ones
                                        
                                        p_value <- sum(OBS >= null_model)/(replicates + 1)
                                        p_value

                                    ## Summarize output and report it
                                        
                                        results[[i-1]] <- data.frame(
                                            trait = colnames(traits)[i],
                                            trait_type = "categorical",
                                            n_species = length(tree$tip.label),
                                            number_of_levels = length(attributes(factor(trait))$levels), 
                                            evolutionary_transitions_observed = OBS,
                                            median_evolutionary_transitions_in_randomization = median(null_model),
                                            minimum_evolutionary_transitions_in_randomization = min(null_model),
                                            evolutionary_transitions_in_randomization = max(null_model),
                                            p_value = p_value
                                        )

                            }

                            if (class(trait) %in% c("numeric", "integer")) {

                                ### For continuous traits

                                    # for (trait in 1:length(cont_traits)) {

                                    #     # Pull out one variable into a new dataframe
                                    #         var <- as.matrix(pdata_continuous[,colnames(pdata_continuous) == cont_traits[trait]])
                                    #         colnames(var) <- cont_traits[trait]
                                    #         rownames(var) <- pdata_continuous$Species
                                    #         var <- na.omit(var)

                                    #     # Drop tips not in data, data not in tree
                                    #         matches <- name.check(nepenthes_tree,var)
                                    #         matches
                                    #         nepenthes_tree_2 <- ape::drop.tip(phy = nepenthes_tree, tip = matches$tree_not_data)
                                    #         name.check(nepenthes_tree_2,var)
                                    #         var <- var[nepenthes_tree_2$tip.label,]

                                    #     # Root the tree and resolve multichotomies
                                    #         # nepenthes_tree_3 <- root(nepenthes_tree_2, "N._pervillei_matK")
                                    #         nepenthes_tree_3 <- nepenthes_tree_2
                                    #         nepenthes_tree_4 <- multi2di(nepenthes_tree_3)

                                        # Determine phylogenetic signal
                                            results[[(i-1)]] <- data.frame(
                                                trait = colnames(traits)[i],
                                                trait_type = "continuous",
                                                n_species = length(tree$tip.label),
                                                number_of_levels = NA,
                                                evolutionary_transitions_observed = NA,
                                                median_evolutionary_transitions_in_randomization = NA,
                                                minimum_evolutionary_transitions_in_randomization = NA,
                                                evolutionary_transitions_in_randomization = NA,
                                                p_value = as.numeric(picante::phylosignal(trait, tree)[4])
                                            )
                            }

                        }

                    ## Return results

                        results <- do.call(rbind, results)
                        results
                        return(results)

                }

###### Datasets

        if ( exists("datasets") ) {

            message("Object 'datasets' exists, not loading phylochemistry datasets....")

        } else {
            
            message("Loading datasets...")

            ## Sample datasets for CHEM5725

                sample_datasets <- as.data.frame(rbind(
                    c("algae_data", "https://thebustalab.github.io/R_For_Chemists_2/sample_data/algae_data.csv"),
                    c("alaska_lake_data", "https://thebustalab.github.io/R_For_Chemists_2/sample_data/alaska_lake_data.csv"),
                    c("solvents", "https://thebustalab.github.io/R_For_Chemists_2/sample_data/solvents.csv"),
                    c("periodic_table", "https://thebustalab.github.io/R_For_Chemists_2/sample_data/per_table.csv"),
                    c("fadb_sample", "https://thebustalab.github.io/R_For_Chemists_2/sample_data/fadb_sample.csv"),
                    c("periodic_table_subset", "https://thebustalab.github.io/R_For_Chemists_2/sample_data/per_table_small.csv"),
                    c("ny_trees", "https://thebustalab.github.io/R_For_Chemists_2/sample_data/ny_trees.csv"),
                    c("metabolomics_data", "https://thebustalab.github.io/R_For_Chemists_2/sample_data/metabolomics_data.csv"),
                    c("wine_grape_data", "https://thebustalab.github.io/R_For_Chemists_2/sample_data/wine_grape_data.csv"),
                    c("hawaii_aquifers", "https://thebustalab.github.io/R_For_Chemists_2/sample_data/hawaii_aquifers.csv"),
                    c("beer_components", "https://thebustalab.github.io/R_For_Chemists_2/sample_data/beer_components.csv"),
                    c("wood_smoke", "https://thebustalab.github.io/R_For_Chemists_2/sample_data/wood_smoke_data.csv"),
                    c("hops_components", "https://thebustalab.github.io/R_For_Chemists_2/sample_data/hops_components.csv")
                ))

                pb <- progress::progress_bar$new(total = dim(sample_datasets)[1])
                for (i in 1:dim(sample_datasets)[1]) {
                    temp_obj <- readr::read_csv(sample_datasets[i,2], show_col_types = FALSE)
                    gdata::mv("temp_obj", as.character(sample_datasets[i,1]))
                    pb$tick()
                }

            ## Busta lab specific datasets

                if (exists("bustalab")) {
                    
                    if (bustalab == TRUE) {

                        busta_spectral_library <- read_csv("https://thebustalab.github.io/R_For_Chemists_2/sample_data/busta_spectral_library_v1.csv", col_types = c(Compound_common_name = "c"))
                        plant_phylogeny <- read.tree("https://thebustalab.github.io/data/plant_phylogeny.newick")
                        plant_species <- readMonolist("https://thebustalab.github.io/data/plant_species.csv")
                    
                    }
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

message("phylochemistry loaded!!")
