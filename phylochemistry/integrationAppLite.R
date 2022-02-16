###### Libraries

    if ( !exists("packages") ) {

        ## Define necessary libraries
            
            CRAN_packages <- c(
                "gridExtra",
                # "ape",
                "multcompView",
                # "imager",
                "shiny",
                "DT",
                "RColorBrewer",
                "data.table",
                "rhandsontable",
                # "ips",
                # "phangorn",
                # "seqinr",
                # "Rfast",
                # "picante",
                # "BiocManager",
                "googlesheets4",
                # "Hmisc",
                # "ggforce",
                # "network",
                # "pracma",
                # "ggnetwork",
                # "FactoMineR",
                "dplyr",
                "stringr", 
                "progress",
                "tidyverse",
                # "ggrepel",
                # "cowplot",
                # "rstatix",
                # "agricolae",
                # "ggpmisc",
                # "exifr",
                # "lubridate",
                # "bio3d",
                # "shipunov",
                # "remotes",
                # "gdata"
            )

            Bioconductor_packages <- c(
                "ggtree",
                # "xcms",
                # "msa",
                # "rtracklayer",
                # "Biostrings",
                # "GenomicRanges",
                # "GenomicFeatures",
                # "Rsamtools"
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

    ## Set up better warning messages for some functions

        shapiroTest <- function( data, ... ) {

            if ( any(summarize(data, size = n())$size < 3) ) {
                stop("One of the groups defined has fewer than 3 members. A Shapiro test cannot be run on such a group. Please filter your data or choose new groups.")
            } else {
                rstatix::shapiro_test(data = data, ...)
            }

        }

        leveneTest <- function( data, ... ) {rstatix::levene_test( data = data, ... )}
        tTest <- function( data, ... ) {rstatix::t_test( data = data, ... )}
        wilcoxTest <- function( data, ... ) {rstatix::wilcox_test( data = data, ... )}
        anovaTest <- function( data, ... ) {rstatix::anova_test( data = data, ... )}
        tukeyTest <- function( data, ... ) {rstatix::tukey_hsd( data = data, ... )}
        kruskalTest <- function( data, ... ) {rstatix::kruskal_test( data = data, ... )}
        dunnTest <- function( data, ... ) {rstatix::dunn_test( data = data, ... )}

    ## Set up lookups

        ## Phred_ascii_33 lookup
            phred33_lookup <- data.frame(rbind(
                c("!","0"),c("\"","1"),c("#","2"),c("$","3"),c("%","4"),c("&","5"),c("'","6"),c("(","7"),c(")","8"),c("*","9"),c("+","10"),c(",","11"),c("-","12"),c(".","13"),c("/","14"),c("0","15"),c("1","16"),c("2","17"),c("3","18"),c("4","19"),c("5","20"),c("6","21"),c("7","22"),c("8","23"),c("9","24"),c(":","25"),c(";","26"),c("<","27"),c("=","28"),c(">","29"),c("?","30"),c("@","31"),c("A","32"),c("B","33"),c("C","34"),c("D","35"),c("E","36"),c("F","37"),c("G","38"),c("H","39"),c("I","40"),c("J","41"),c("K","42"),c("L","43"),c("M","44"),c("N","45"),c("O","46"),c("P","47"),c("Q","48"),c("R","49"),c("S","50"),c("T","51"),c("U","52"),c("V","53"),c("W","54"),c("X","55"),c("Y","56"),c("Z","57"),c("[","58"),c("\\","59"),c("]","60"),c("^","61"),c("_","62"),c("`","63"),c("a","64"),c("b","65"),c("c","66"),c("d","67"),c("e","68"),c("f","69"),c("g","70"),c("h","71"),c("i","72"),c("j","73"),c("k","74"),c("l","75"),c("m","76"),c("n","77"),c("o","78"),c("p","79"),c("q","80"),c("r","81"),c("s","82"),c("t","83"),c("u","84"),c("v","85"),c("w","86"),c("x","87"),c("y","88"),c("z","89")
            ))

###### Datasets

    message("Loading MS Library...")

        busta_spectral_library <- read_csv("https://thebustalab.github.io/R_For_Chemists_2/sample_data/busta_spectral_library_v1.csv", col_types = c(Compound_common_name = "c"))

###### Functions

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

            convertCDFstoCSVs <-    function(
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

                    if (dim(data)[1] > 0) {
                        peak_list[[i]] <- data[,1:9]
                    }
                
                }
            
            }

            do.call(rbind, peak_list)

        }

    #### integrationAppLite

        #' A Shiny app to integrate GC-FID and GC-MS data
        #'
        #' could not find function "." is from a plyr/dplyr issue
        #' @param chromatograms A data frame containing columns: "rt", "tic", and "path_to_cdf_csv", which contain retention time, total ion chromatogram intensities, and paths to CDF.csv files generated by the convertCDFstoCSVs function.
        #' @param x_axis_start A numeric value for the lower x-axis bounds on the plot generated by the app. Defaults to full length.
        #' @param x_axis_end A numeric value for the upper x-axis bounds on the plot generated by the app. Defaults to full length.
        #' @param samples_monolist_path A path to a .csv file containing metadata for the samples you wish to analyze. Requied columns are: "rt_offset", "baseline_window", and "path_to_cdf_csv", which are for aligning chromatograms, adjusting baseline determination, and defining the path to the CDF.csv files for each sample, respectively.
        #' @param samples_monolist_subset Optional, a numeric vector (for example, "c(1:10)"), defining a subset of samples to be loaded.
        #' @param peaks_monolist_path A path to a .csv file containing metadata for all peaks in the sample set. Required columns are: peak_start", "peak_end", "path_to_cdf_csv", "area", "peak_number_within_sample", "rt_offset", "peak_start_rt_offset", "peak_end_rt_offset". This file is automatically generated by the app.
        #' @param zoom_and_scroll_rate Defines intervals of zooming and scrolling movement while running the app
        #' @examples
        #' @export
        #' integrationAppLite

        integrationAppLite <- function(
                CDF_directory_path,
                zoom_and_scroll_rate = 100,
                baseline_window = 400,
                path_to_reference_library = busta_spectral_library
            ) {

                ## Set up monolists and status
                    samples_monolist_subset = NULL

                    if ( .Platform$OS == "unix") {

                        samples_monolist_path <- paste0(CDF_directory_path, "/samples_monolist.csv")
                        peaks_monolist_path <- paste0(CDF_directory_path, "/peaks_monolist.csv")

                    }

                    if ( .Platform$OS == "windows") {

                        samples_monolist_path <- paste0(CDF_directory_path, "\\samples_monolist.csv")
                        peaks_monolist_path <- paste0(CDF_directory_path, "\\peaks_monolist.csv")

                    }
                    
                ## Covert CDF to csv, if necessary

                    if ( .Platform$OS == "unix") {

                         ## Extract paths_to_cdfs

                            paths_to_cdfs <- paste(
                                CDF_directory_path,
                                dir(CDF_directory_path)[grep("*.CDF$", dir(CDF_directory_path), ignore.case = TRUE)],
                                sep = "/"
                            )

                        ## Remove any trailing slashes

                            for ( i in 1:length(paths_to_cdfs)) {

                                if (substr(paths_to_cdfs, nchar(paths_to_cdfs[i]), nchar(paths_to_cdfs[i])) == "/") {
                                    paths_to_cdfs <- substr(paths_to_cdfs, 0, nchar(paths_to_cdfs[i])-1)
                                }
                                if (substr(paths_to_cdfs, nchar(paths_to_cdfs[i]), nchar(paths_to_cdfs[i])) == "/") {
                                    paths_to_cdfs <- substr(paths_to_cdfs, 0, nchar(paths_to_cdfs[i])-1)
                                }
                            
                            }
                    }

                    if ( .Platform$OS == "windows") {

                        ## Extract paths_to_cdfs

                            paths_to_cdfs <- paste(
                                CDF_directory_path,
                                dir(CDF_directory_path)[grep("*.CDF$", dir(CDF_directory_path), ignore.case = TRUE)],
                                sep = "\\"
                            )

                        ## Remove any trailing slashes

                            for ( i in 1:length(paths_to_cdfs)) {

                                if (substr(paths_to_cdfs, nchar(paths_to_cdfs[i]), nchar(paths_to_cdfs[i])) == "\\") {
                                    paths_to_cdfs <- substr(paths_to_cdfs, 0, nchar(paths_to_cdfs[i])-1)
                                }
                                if (substr(paths_to_cdfs, nchar(paths_to_cdfs[i]), nchar(paths_to_cdfs[i])) == "\\") {
                                    paths_to_cdfs <- substr(paths_to_cdfs, 0, nchar(paths_to_cdfs[i])-1)
                                }
                            
                            }
                    }

                    convertCDFstoCSVs(paths_to_cdfs)

                ## Recreate chromatograms.csv, if necessary

                    paths_to_cdf_csvs <- paste0(paths_to_cdfs, ".csv")

                    ## Handle chromatograms file

                        if ( .Platform$OS == "unix") {

                            # If chromatograms doesn't exist, make it 
                                if (!file.exists(paste0(CDF_directory_path, "/chromatograms.csv"))) {
                                    chromatograms <- extractChromatogramsFromCSVs(paths_to_cdf_csvs)
                                    writeMonolist(chromatograms, paste0(CDF_directory_path, "/chromatograms.csv"))
                                } else {
                                    # If it exists, check to see if all cdfs in this folder are in it, if not, recreate
                                    chromatograms <- readMonolist(paste0(CDF_directory_path, "/chromatograms.csv"))
                                    if (!all(paths_to_cdf_csvs %in% unique(chromatograms$path_to_cdf_csv))) {
                                        chromatograms <- extractChromatogramsFromCSVs(paths_to_cdf_csvs)
                                        writeMonolist(chromatograms, paste0(CDF_directory_path, "/chromatograms.csv"))
                                    }
                                }

                        }

                        if ( .Platform$OS == "windows") {

                            # If chromatograms doesn't exist, make it 
                                if (!file.exists(paste0(CDF_directory_path, "\\chromatograms.csv"))) {
                                    chromatograms <- extractChromatogramsFromCSVs(paths_to_cdf_csvs)
                                    writeMonolist(chromatograms, paste0(CDF_directory_path, "\\chromatograms.csv"))
                                } else {
                                    # If it exists, check to see if all cdfs in this folder are in it, if not, recreate
                                    chromatograms <- readMonolist(paste0(CDF_directory_path, "\\chromatograms.csv"))
                                    if (!all(paths_to_cdf_csvs %in% unique(chromatograms$path_to_cdf_csv))) {
                                        chromatograms <- extractChromatogramsFromCSVs(paths_to_cdf_csvs)
                                        writeMonolist(chromatograms, paste0(CDF_directory_path, "\\chromatograms.csv"))
                                    }
                                }

                        }

                ## Set up new samples monolist
                        
                    samples_monolist <- data.frame(
                        Sample_ID = gsub("\\..*$", "", gsub(".*/", "", unique(chromatograms$path_to_cdf_csv))),
                        rt_offset = 0,
                        baseline_window = baseline_window,
                        path_to_cdf_csv = unique(chromatograms$path_to_cdf_csv)
                    )

                    write.table(
                        x = samples_monolist,
                        file = samples_monolist_path,
                        row.names = FALSE,
                        sep = ","
                    )

                ## Set up several variables, plot_height, and x_axis limits if not specified in function call
                    
                    peak_data <- NULL
                    peak_points <- NULL
                    if ( length(samples_monolist_subset) > 0 ) {
                        plot_height <- 200 + 100*length(samples_monolist_subset)    
                    } else {
                        plot_height <- 200 + 100*dim(samples_monolist)[1]
                    }

                    x_axis_start <<- min(chromatograms$rt)
                    x_axis_end <<- max(chromatograms$rt)
                    
                ## Set up new peak monolist if it doesn't exist
                    
                    if ( !file.exists(peaks_monolist_path) ) {
                        
                        peak_data <- data.frame(
                          peak_start = 0,
                          peak_end = 0,
                          peak_ID = "unknown",
                          path_to_cdf_csv = "a",
                          area = 0
                        )

                        write.table(
                          x = peak_data[-1,],
                          file = peaks_monolist_path,
                          append = FALSE,
                          row.names = FALSE,
                          col.names = TRUE,
                          sep = ","
                        )

                    }

                ## Set up user interface

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

                ## Set up the server

                    server <- function(input, output, session) {

                        ## Check keystoke value
                            output$key <- renderPrint({
                                input$keypress
                            })

                        ## Keys to move chromatogram view - zoom in and out, move L and R
                            observeEvent(input$keypress, {
                                if( input$keypress == 70 ) { x_axis_start <<- x_axis_start + zoom_and_scroll_rate } # Forward on "F"
                                if( input$keypress == 70 ) { x_axis_end <<- x_axis_end + zoom_and_scroll_rate } # Forward on "F"
                                if( input$keypress == 68 ) { x_axis_start <<- x_axis_start - zoom_and_scroll_rate } # Backward on "D"
                                if( input$keypress == 68 ) { x_axis_end <<- x_axis_end - zoom_and_scroll_rate } # Backward on "D"
                                if( input$keypress == 86 ) { x_axis_start <<- x_axis_start - zoom_and_scroll_rate } # Wider on "V"
                                if( input$keypress == 86 ) { x_axis_end <<- x_axis_end + zoom_and_scroll_rate } # Wider on "V"
                                if( input$keypress == 67 ) { x_axis_start <<- x_axis_start + zoom_and_scroll_rate } # Closer on "C"
                                if( input$keypress == 67 ) { x_axis_end <<- x_axis_end - zoom_and_scroll_rate } # Closer on "C"
                            })

                        ## Save manual changes to table on "Z" (90) keystroke

                            observeEvent(input$keypress, {

                                if (input$keypress == 90 ) {

                                    ## Write out any modifications to peak table (i.e. sample IDs)
                                        
                                        hot = isolate(input$peak_table)
                                        if (!is.null(hot)) {
                                            writeMonolist(rhandsontable::hot_to_r(input$peak_table), peaks_monolist_path)
                                            print(peaks_monolist_path)
                                        }

                                }

                            })

                        ## Update chromatogram on "Q" (81) keystroke
                            
                            observeEvent(input$keypress, {      
                                
                                if( input$keypress == 81 ) { # Update on "Q"
                                    
                                    output$chromatograms <- renderPlot({

                                        ## Read in samples monolist and put chromatograms into chromatograms_updated
                                            
                                            samples_monolist <- read.csv(samples_monolist_path)
                                            if ( length(samples_monolist_subset) > 0 ) {
                                                samples_monolist <- samples_monolist[samples_monolist_subset,]    
                                            }
                                            chromatograms_updated <- dplyr::filter(chromatograms, path_to_cdf_csv %in% samples_monolist$path_to_cdf_csv)

                                        ## Calculate baseline for each sample

                                            baselined_chromatograms <- list()

                                            for ( chrom in 1:length(unique(chromatograms_updated$path_to_cdf_csv)) ) {
                                      
                                                chromatogram <- dplyr::filter(chromatograms_updated, path_to_cdf_csv == unique(chromatograms_updated$path_to_cdf_csv)[chrom])

                                                prelim_baseline_window <- samples_monolist$baseline_window[match(chromatogram$path_to_cdf_csv[1], samples_monolist$path_to_cdf_csv)]

                                                n_prelim_baseline_windows <- floor(length(chromatogram$rt)/prelim_baseline_window)
                                                prelim_baseline <- list()
                                                for ( i in 1:n_prelim_baseline_windows ) {
                                                  min <- min(chromatogram$tic[((prelim_baseline_window*(i-1))+1):(prelim_baseline_window*i)])
                                                  prelim_baseline[[i]] <-     data.frame(
                                                                          rt = chromatogram$rt[chromatogram$tic == min],
                                                                          min = min
                                                                      )
                                                }
                                                prelim_baseline <- do.call(rbind, prelim_baseline)
                                                chromatogram$in_prelim_baseline <- FALSE
                                                chromatogram$in_prelim_baseline[chromatogram$rt %in% prelim_baseline$rt] <- TRUE

                                                y = prelim_baseline$min
                                                x = prelim_baseline$rt

                                                baseline2 <- data.frame(
                                                              rt = chromatogram$rt,
                                                              y = approx(x, y, xout = chromatogram$rt)$y
                                                          )
                                                baseline2 <- baseline2[!is.na(baseline2$y),]
                                                chromatogram <- chromatogram[chromatogram$rt %in% baseline2$rt,]
                                                chromatogram$baseline <- baseline2$y

                                                baselined_chromatograms[[chrom]] <- data.frame(
                                                    rt = chromatogram$rt,
                                                    tic = chromatogram$tic,
                                                    path_to_cdf_csv = chromatogram$path_to_cdf_csv,
                                                    in_prelim_baseline = chromatogram$in_prelim_baseline,
                                                    baseline = chromatogram$baseline
                                                )

                                            }

                                            chromatograms_updated <- do.call(rbind, baselined_chromatograms)

                                        ## Add rt offset information for all chromatograms

                                            chromatograms_updated$rt_offset <- samples_monolist$rt_offset[match(chromatograms_updated$path_to_cdf_csv, samples_monolist$path_to_cdf_csv)]
                                            chromatograms_updated$rt_rt_offset <- chromatograms_updated$rt + chromatograms_updated$rt_offset
                                            chromatograms_updated <<- chromatograms_updated

                                        ## Plot with peaks, if any
                                            
                                            peak_table <- read.csv(peaks_monolist_path)
                                    
                                            if (dim(peak_table)[1] > 0) {

                                                ## Filter out duplicate peaks and NA peaks
                                                    
                                                    library(plyr)
                                                    peak_table <- plyr::ddply(peak_table, .(path_to_cdf_csv), mutate, duplicated = duplicated(peak_start))
                                                    peak_table <- dplyr::filter(peak_table, duplicated == FALSE)
                                                    peak_table <- plyr::ddply(peak_table, .(path_to_cdf_csv), mutate, duplicated = duplicated(peak_end))
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
                                                            rt >= peaks_in_this_sample$peak_start[peak] & rt <= peaks_in_this_sample$peak_end[peak])$tic
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
                                        
                                                    write.table(peak_table, file = peaks_monolist_path, col.names = TRUE, sep = ",", row.names = FALSE)

                                                    output$peak_table <- rhandsontable::renderRHandsontable(rhandsontable::rhandsontable({
                                                        peak_table2 <- read.csv(peaks_monolist_path)
                                                        # peak_table2[,1:6]
                                                        peak_table2
                                                    }))

                                                ## Create the plot with subsetted data to make it faster

                                                    ## Subset the chromatograms and peaks

                                                        chromatograms_updated <- dplyr::filter(chromatograms_updated, rt_rt_offset > x_axis_start & rt_rt_offset < x_axis_end)
                                                        peak_table <- dplyr::filter(peak_table, peak_start_rt_offset > x_axis_start & peak_end_rt_offset < x_axis_end)

                                                    ## Make chromatogram plot object

                                                        # Make their labels easy to read
                                                            facet_labels <- gsub(".CDF.csv", "", gsub(".*/", "", chromatograms_updated$path_to_cdf_csv))
                                                            names(facet_labels) <- chromatograms_updated$path_to_cdf_csv

                                                            p <-    ggplot() + 
                                                                geom_line(data = chromatograms_updated, mapping = aes(x = rt_rt_offset, y = baseline), color = "grey") +
                                                                geom_line(data = chromatograms_updated, mapping = aes(x = rt_rt_offset, y = tic)) +
                                                                scale_x_continuous(limits = c(x_axis_start, x_axis_end)) +
                                                                facet_grid(path_to_cdf_csv~., scales = "free_y", labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                                                # facet_grid(path_to_cdf_csv~., labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                                                # scale_y_continuous(limits = c(0, max(dplyr::filter(chromatograms_updated, rt > x_axis_start & rt < x_axis_end)$tic))) +
                                                                theme_classic() +
                                                                scale_fill_continuous(type = "viridis") +
                                                                theme(
                                                                  legend.position = 'none',
                                                                  legend.title = element_blank()
                                                                )

                                                    ## Add peaks

                                                        for (peak in 1:dim(peak_table)[1]) {
                                                            
                                                            signal_for_this_peak <- dplyr::filter(
                                                                chromatograms_updated[chromatograms_updated$path_to_cdf_csv == peak_table[peak,]$path_to_cdf_csv,], 
                                                                rt_rt_offset > peak_table[peak,]$peak_start_rt_offset, 
                                                                rt_rt_offset < peak_table[peak,]$peak_end_rt_offset
                                                            )

                                                            if (dim(signal_for_this_peak)[1] > 0) {

                                                                signal_for_this_peak$peak_number_within_sample <- peak_table$peak_number_within_sample[peak]
                                                                
                                                                p <- p +  geom_vline(data = signal_for_this_peak[1,], mapping = aes(xintercept = rt_rt_offset), alpha = 0.3) +
                                                                          geom_ribbon(data = signal_for_this_peak, mapping = aes(x = rt_rt_offset, ymax = tic, ymin = baseline, fill = peak_number_within_sample)) +
                                                                          geom_text(data = signal_for_this_peak, mapping = aes(label = peak_number_within_sample, x = median(rt_rt_offset), y = max(tic)))
                                                            }
                                                        }

                                            } else {

                                                ## Make chromatogram plot object

                                                    chromatograms_updated <- dplyr::filter(chromatograms_updated, rt_rt_offset > x_axis_start & rt_rt_offset < x_axis_end)
                                                
                                                        # Make their labels easy to read
                                                            facet_labels <- gsub(".CDF.csv", "", gsub(".*/", "", chromatograms_updated$path_to_cdf_csv))
                                                            names(facet_labels) <- chromatograms_updated$path_to_cdf_csv

                                                        p <-  ggplot() + 
                                                                geom_line(data = chromatograms_updated, mapping = aes(x = rt_rt_offset, y = baseline), color = "grey") +
                                                                geom_line(data = chromatograms_updated, mapping = aes(x = rt_rt_offset, y = tic)) +
                                                                scale_x_continuous(limits = c(x_axis_start, x_axis_end)) +
                                                                facet_grid(path_to_cdf_csv~., scales = "free_y", labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                                                # facet_grid(path_to_cdf_csv~., labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                                                # scale_y_continuous(limits = c(0, max(dplyr::filter(chromatograms_updated, rt > x_axis_start & rt < x_axis_end)$tic))) +
                                                                theme_classic() +
                                                                scale_fill_continuous(type = "viridis") +
                                                                theme(
                                                                  legend.position = 'none',
                                                                  legend.title = element_blank()
                                                                )
                                            }

                                        ## Draw the plot
                                            
                                            p
                                    })
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
                                        area = sum(peak_points$tic)
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
                                        peak_table <- read.csv(peaks_monolist_path)

                                        peak_table <- peak_table[
                                            !apply(cbind(
                                                peak_table$peak_start > selection_start,
                                                peak_table$peak_end < selection_end,
                                                peak_table$path_to_cdf_csv == as.character(peak_points$path_to_cdf_csv[1])
                                            ), 1, all)
                                        ,]

                                        write.table(
                                            x = peak_table,
                                            file = peaks_monolist_path,
                                            append = FALSE,
                                            row.names = FALSE,
                                            col.names = TRUE,
                                            sep = ","
                                        )
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
                                            file = peaks_monolist_path,
                                            append = TRUE,
                                            row.names = FALSE,
                                            col.names = FALSE,
                                            sep = ","
                                        )

                                        output$peak_table <- rhandsontable::renderRHandsontable(rhandsontable::rhandsontable({
                                            peak_table2 <- read.csv(peaks_monolist_path)
                                            peak_table2
                                        }))
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
                                            file = peaks_monolist_path,
                                            append = TRUE,
                                            row.names = FALSE,
                                            col.names = FALSE,
                                            sep = ","
                                        )

                                        output$peak_table <- DT::renderDataTable(DT::datatable({
                                            peak_table <- read.csv(peaks_monolist_path)
                                            peak_table
                                        }))
                                    }
                            })

                        ## [MS1 extract ("shift+1"), update ("shift+2"), subtract ("shift+3"), library search ("shift+4"), save ("shift+5")]
                            
                            observeEvent(input$keypress, {

                                ## If "shift+1", MS from chromatogram brush -> MS_out_1
                                    if( input$keypress == 33 ) {
                                        
                                        framedDataFile <- isolate(as.data.frame(
                                                            data.table::fread(as.character(
                                                                brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1]
                                                            ))
                                        ))
                                        framedDataFile <- isolate(dplyr::filter(
                                                                framedDataFile, 
                                                                rt > min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
                                                                rt < max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
                                                            ))
                                        library(plyr)
                                        framedDataFile$mz <- round(framedDataFile$mz, 1)
                                        framedDataFile <- plyr::ddply(framedDataFile, .(mz), summarize, intensity = sum(intensity))
                                        MS_out_1 <<- framedDataFile

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
                                            library(plyr)
                                            framedDataFile_to_subtract$mz <- round(framedDataFile_to_subtract$mz, 1)
                                            framedDataFile_to_subtract <- plyr::ddply(framedDataFile_to_subtract, .(mz), summarize, intensity = sum(intensity))
                                            MS_out_1$intensity <- MS_out_1$intensity - framedDataFile_to_subtract$intensity
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
                                                mz_missing <- seq(min(MS_out_1$mz), max(MS_out_1$mz), 1)[!seq(min(MS_out_1$mz), max(MS_out_1$mz), 1) %in% MS_out_1$mz]
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
                                                    t(
                                                        c(
                                                            rep(0,min(MS_out_1$mz)-1),
                                                            MS_out_1$intensity
                                                        )
                                                    )
                                                )
                                                colnames(unknown)[6:805] <- paste("mz_", c(seq(1,min(MS_out_1$mz)-1,1), MS_out_1$mz), sep = "")

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
                                                    scale_variance = TRUE,
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
                                                    scale_variance = TRUE,
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

                        ## [MS1 Subtract] Subtract selected MS as background from MS panel one on "shift+3" (35) keypress
                            # observeEvent(input$keypress, {

                            #     # Do nothing if no selection
                            #         if(is.null(input$chromatogram_brush)) {
                            #             return()
                            #         }

                            #     # If selection and "35" is pressed, extract and print mass spectra
                            #         if( input$keypress == 35 ) {

                            #             output$massSpectra_1 <- renderPlot({
                            #                 framedDataFile_to_subtract <- isolate(as.data.frame(
                            #                                     data.table::fread(as.character(
                            #                                         brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1]
                            #                                     ))
                            #                 ))
                            #                 framedDataFile_to_subtract <- isolate(dplyr::filter(
                            #                                         framedDataFile_to_subtract, 
                            #                                         rt > min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
                            #                                         rt < max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
                            #                                     ))
                            #                 library(plyr)
                            #                 framedDataFile_to_subtract$mz <- round(framedDataFile_to_subtract$mz, 0)
                            #                 framedDataFile_to_subtract <- plyr::ddply(framedDataFile_to_subtract, .(mz), summarize, intensity = sum(intensity))
                            #                 framedDataFile_1$intensity <- framedDataFile_1$intensity - framedDataFile_to_subtract$intensity
                            #                 framedDataFile_1 <<- framedDataFile_1
                            #                 framedDataFile_1$intensity <- framedDataFile_1$intensity*100/max(framedDataFile_1$intensity)

                            #                 # Set ranges on mass spec (allow zooming in with selection)
                            #                     if(isolate(is.null(input$massSpectra_1_brush))) {
                            #                         MS1_low_x_limit <- 0
                            #                         MS1_high_x_limit <- 800
                            #                         MS1_high_y_limit <- 110
                            #                     } else {
                            #                         MS1_low_x_limit <- isolate(min(brushedPoints(framedDataFile_1, input$massSpectra_1_brush)$mz))
                            #                         MS1_high_x_limit <- isolate(max(brushedPoints(framedDataFile_1, input$massSpectra_1_brush)$mz))
                            #                         MS1_high_y_limit <- max(dplyr::filter(framedDataFile_1, intensity > MS1_low_x_limit & intensity < MS1_high_x_limit)$intensity) + 10
                            #                     }

                            #                 ggplot() + 
                            #                     geom_bar(
                            #                         data = framedDataFile_1,
                            #                         mapping = aes(x = mz, y = intensity),
                            #                         stat = "identity", width = 1,
                            #                         color = "black", fill = "grey"
                            #                     ) +
                            #                     theme_classic() +
                            #                     scale_x_continuous(expand = c(0,0)) +
                            #                     scale_y_continuous(expand = c(0,0)) +
                            #                     coord_cartesian(xlim = c(MS1_low_x_limit, MS1_high_x_limit), ylim = c(0, MS1_high_y_limit)) +
                            #                     geom_text(
                            #                         data = 
                            #                             dplyr::filter(framedDataFile_1, mz > MS1_low_x_limit, mz < MS1_high_x_limit)[
                            #                                 order(
                            #                                     dplyr::filter(framedDataFile_1, mz > MS1_low_x_limit & mz < MS1_high_x_limit)$intensity,
                            #                                     decreasing = TRUE
                            #                                 )[1:10]
                            #                             ,],
                            #                         mapping = aes(x = mz, y = intensity + 5, label = mz)
                            #                     )
                            #             })
                            #         }
                            # })
                            
                           

                        ## Extract selected MS to MS panel two on "shift+2" (64) keypress
                            
                        #     observeEvent(input$keypress, {

                        #         # Do nothing if no selection
                        #             if(is.null(input$chromatogram_brush)) {
                        #                 return()
                        #             }

                        #         # If selection and "shift+22" is pressed, extract and print mass spectra
                        #             if( input$keypress == 64 ) {

                        #                 output$massSpectra_2 <- renderPlot({
                        #                     framedDataFile <- isolate(as.data.frame(
                        #                                         data.table::fread(as.character(
                        #                                             brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1]
                        #                                         ))
                        #                 ))
                        #                 framedDataFile <- isolate(dplyr::filter(
                        #                                             framedDataFile, 
                        #                                             rt > min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
                        #                                             rt < max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
                        #                                         ))
                        #                 library(plyr)
                        #                 framedDataFile <- plyr::ddply(framedDataFile, .(mz), summarize, intensity = sum(intensity))
                        #                 framedDataFile_2 <<- framedDataFile
                        #                 framedDataFile$intensity <- framedDataFile$intensity*100/max(framedDataFile$intensity)
                        #                 ggplot() + 
                        #                     geom_bar(data = framedDataFile, mapping = aes(x = mz, y = intensity), stat = "identity", width = 1) +
                        #                     theme_classic() +
                        #                     scale_x_continuous(expand = c(0,0)) +
                        #                     scale_y_continuous(expand = c(0,0), limits = c(0,110)) +
                        #                     geom_text(data = dplyr::filter(framedDataFile, intensity > 10), mapping = aes(x = mz, y = intensity+5, label = mz))
                        #                 })
                        #             }
                        #     })

                        # ## Subtract selected MS as background from MS panel one on "shift+4" (36) keypress
                            
                        #     observeEvent(input$keypress, {

                        #         # Do nothing if no selection
                        #             if(is.null(input$chromatogram_brush)) {
                        #                 return()
                        #             }

                        #             # If selection and "shift+4" is pressed, extract and print mass spectra
                        #                 if( input$keypress == 36 ) {

                        #                     output$massSpectra_2 <- renderPlot({
                        #                         framedDataFile_to_subtract <- isolate(as.data.frame(
                        #                                             data.table::fread(as.character(
                        #                                                 brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1]
                        #                                             ))
                        #                         ))
                        #                         framedDataFile_to_subtract <- isolate(dplyr::filter(
                        #                                                 framedDataFile_to_subtract, 
                        #                                                 rt > min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
                        #                                                 rt < max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
                        #                                             ))
                        #                         library(plyr)
                        #                         framedDataFile_to_subtract <- plyr::ddply(framedDataFile_to_subtract, .(mz), summarize, intensity = sum(intensity))
                        #                         framedDataFile_2$intensity <- framedDataFile_2$intensity - framedDataFile_to_subtract$intensity
                        #                         framedDataFile_2 <<- framedDataFile_2
                        #                         framedDataFile_2$intensity <- framedDataFile_2$intensity*100/max(framedDataFile_2$intensity)
                        #                         ggplot() + 
                        #                             geom_bar(data = framedDataFile_2, mapping = aes(x = mz, y = intensity), stat = "identity", width = 1) +
                        #                             theme_classic() +
                        #                             scale_x_continuous(expand = c(0,0)) +
                        #                             scale_y_continuous(expand = c(0,0), limits = c(0,110)) +
                        #                             geom_text(data = dplyr::filter(framedDataFile_2, intensity > 10), mapping = aes(x = mz, y = intensity+5, label = mz))
                        #                     })
                        #                 }
                        #         })

                        # ## Save MS in MS panel one on "shift+6" (94) keypress

                        #     observeEvent(input$keypress, {

                        #         # Do nothing if no selection
                        #             if(is.null(input$chromatogram_brush)) {
                        #                 return()
                        #             }

                        #         # If selection and "shift+6" is pressed, extract and print mass spectra
                        #             if( input$keypress == 96 ) {
                        #                 framedDataFile_2$intensity[framedDataFile_2$intensity < 0] <- 0
                        #                 write.csv(framedDataFile_2, "integration_app_spectrum_2.csv")
                        #             }
                        #     })

                    }

                ## Call the app
                    
                    shinyApp(ui = ui, server = server)

            }

