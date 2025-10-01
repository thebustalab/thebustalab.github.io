###### Libraries

    ## Define necessary libraries
        
        CRAN_packages <- c(
            "shiny",
            "shinythemes",
            "data.table",
            "rhandsontable",
            "dplyr",
            "ggplot2",
            "readr",
            "tidyr",
            "plyr"
        )

        Bioconductor_packages <- character(0)

        Github_packages <- c(
            # "HajkD/orthologr"
        )

        packages_needed <- c(CRAN_packages, Bioconductor_packages, Github_packages)[!c(CRAN_packages, Bioconductor_packages, gsub(".*/", "", Github_packages)) %in% rownames(installed.packages())]

        skip_package_check <- tolower(Sys.getenv("SKIP_PACKAGE_CHECK", "false")) %in% c("true", "1", "yes")

    ## Determine if anything needs to be installed
            
            if (!skip_package_check && length(packages_needed) > 0) {

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
            } else if (skip_package_check && length(packages_needed) > 0) {
                stop(
                    paste0(
                        "Required packages missing but SKIP_PACKAGE_CHECK is set. Missing: ",
                        paste(packages_needed, collapse = ", ")
                    )
                )
            }

            message("Loading packages...")
            invisible(suppressMessages(suppressWarnings(lapply(c(CRAN_packages, Bioconductor_packages), library, character.only = TRUE))))
    
    ## Set up prioriy functions

        mutate <- dplyr::mutate
        summarize <- dplyr::summarize
        filter <- dplyr::filter

    ## Set up options

        options(readr.show_progress = FALSE)
        options(dplyr.summarise.inform = FALSE)

###### Functions

    message("Loading functions...")

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
                                framedDataFile$row_number <- seq_len(nrow(framedDataFile))
                                chromatogram <- framedDataFile %>%
                                    group_by(rt) %>%
                                    summarize(
                                        tic = sum(intensity),
                                        rt_first_row_in_raw = min(row_number),
                                        rt_last_row_in_raw = max(row_number),
                                        .groups = "drop"
                                    )
                                chromatogram <- as.data.frame(chromatogram)
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

    #### writeCSV

            #' Write a monolist to a CSV file
            #'
            #' This function allows the user to interactively select a location to save a monolist as a CSV file.
            #'
            #' @param monolist A data frame representing the monolist to be written to a CSV file.
            #' @return None. The function writes the monolist to a file and does not return a value.
            #' @examples
            #' \dontrun{
            #' # To write a monolist to a CSV file
            #' writeCSV(monolist)
            #' }
            #' @export

            writeCSV <- function(monolist) {
                writeMonolist(monolist = monolist, monolist_out_path = file.choose(new = TRUE))
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
                # write_csv(monolist, file = monolist_out_path)
            }

    #### readMonolist

        #' Read a monolist from a CSV file or URL
        #'
        #' This function reads a monolist from a specified path or URL. The monolist should be in CSV format.
        #'
        #' @param monolist_in_path A string specifying the path or URL to the monolist CSV file.
        #' @return A data frame containing the contents of the monolist.
        #' @examples
        #' \dontrun{
        #' # To read a monolist from a local file
        #' monolist <- readMonolist("path/to/monolist.csv")
        #'
        #' # To read a monolist from a URL
        #' monolist <- readMonolist("http://example.com/monolist.csv")
        #' }
        #' @export

        readMonolist <- function( monolist_in_path ) {
            if (length(grep("http", monolist_in_path)) > 0) {
                monolist <- read_csv(monolist_in_path, show_col_types = FALSE)        
            } else {
                monolist <- as.data.frame(data.table::fread(file = monolist_in_path))
            }
            return( monolist )
        }

    #### analyzeGCMSdata

        #' A Shiny app to integrate GC-FID and GC-MS data
        #'
        #' could not find function "." is from a plyr/dplyr issue
        #' @param chromatograms A data frame containing columns: "rt", "tic", and "path_to_cdf_csv", which contain retention time, total ion chromatogram intensities, and paths to pre-generated sample CSV files.
        #' @param x_axis_start A numeric value for the lower x-axis bounds on the plot generated by the app. Defaults to full length.
        #' @param x_axis_end A numeric value for the upper x-axis bounds on the plot generated by the app. Defaults to full length.
        #' @param samples_monolist_path A path to a .csv file containing metadata for the samples you wish to analyze. Requied columns are: "rt_offset", "baseline_window", and "path_to_cdf_csv", which are for aligning chromatograms, adjusting baseline determination, and defining the path to the CDF.csv files for each sample, respectively.
        #' @param samples_monolist_subset Optional, a numeric vector (for example, "c(1:10)"), defining a subset of samples to be loaded.
        #' @param peaks_monolist_path A path to a .csv file containing metadata for all peaks in the sample set. Required columns are: peak_start", "peak_end", "path_to_cdf_csv", "area", "peak_number_within_sample", "rt_offset", "peak_start_rt_offset", "peak_end_rt_offset". This file is automatically generated by the app.
        #' @param zoom_and_scroll_rate Defines intervals of zooming and scrolling movement while running the app
        #' @param zooming Logical flag that, when TRUE, zooms the chromatogram view to the active brush range after updates (default FALSE).
        #' @examples
        #' @export
        #' analyzeGCMSdata

        analyzeGCMSdata <- function(
                CDF_csv_directory_path,
                zoom_and_scroll_rate = 100,
                baseline_window = 100,
                zooming = TRUE
            ) {

                ## Set up monolists and status
                    samples_monolist_subset = NULL

                    samples_monolist_path <- file.path(CDF_csv_directory_path, "samples_monolist.csv")
                    peaks_monolist_path <- file.path(CDF_csv_directory_path, "peaks_monolist.csv")
                    chromatograms_path <- file.path(CDF_csv_directory_path, "chromatograms.csv")

                ## Identify pre-generated sample CSV files

                    csv_candidates <- list.files(
                        CDF_csv_directory_path,
                        pattern = "\\.csv$",
                        ignore.case = TRUE,
                        full.names = TRUE
                    )

                    reserved_filenames <- c(
                        basename(samples_monolist_path),
                        basename(peaks_monolist_path),
                        basename(chromatograms_path)
                    )

                    paths_to_cdf_csvs <- csv_candidates[!tolower(basename(csv_candidates)) %in% tolower(reserved_filenames)]

                    if (length(paths_to_cdf_csvs) == 0) {
                        stop(
                            "No sample CSV files found. Provide pre-converted chromatographic CSV files in the target directory.",
                            call. = FALSE
                        )
                    }

                ## Recreate chromatograms.csv, if necessary

                    if (!file.exists(chromatograms_path)) {
                        chromatograms <- extractChromatogramsFromCSVs(paths_to_cdf_csvs)
                        writeMonolist(chromatograms, chromatograms_path)
                    } else {
                        chromatograms <- readMonolist(chromatograms_path)
                        if (!all(paths_to_cdf_csvs %in% unique(chromatograms$path_to_cdf_csv))) {
                            chromatograms <- extractChromatogramsFromCSVs(paths_to_cdf_csvs)
                            writeMonolist(chromatograms, chromatograms_path)
                        }
                    }

                ## Set up new samples monolist
                        
                    samples_monolist <- data.frame(
                        Sample_ID = sub("\\.[^.]+$", "", basename(unique(chromatograms$path_to_cdf_csv))),
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
                    x_axis_start_default <- NULL
                    x_axis_end_default <- NULL
                    
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

                        theme = shinythemes::shinytheme("yeti"),

                        sidebarLayout(

                            sidebarPanel(
                                style = "overflow-y: auto; max-height: 90vh;",

                                h4("GCMS Analysis App"),

                                # img(
                                #     src = "https://raw.githubusercontent.com/thebustalab/thebustalab.github.io/refs/heads/master/logo.png",
                                #     height = "260px", 
                                #     style = "border-radius:8px; margin-bottom:5px;"
                                # ),

                                tags$hr(),
                                strong("Keyboard Shortcuts:"),
                                tags$ul(
                                    tags$li("Shift + Q => Update chromatogram"),
                                    tags$li("------------------"),
                                    tags$li("Shift + 6 => Detect peaks"),
                                    tags$li("Shift + A => Add single peak"),
                                    tags$li("Shift + G => Add global peak"),
                                    tags$li("Shift + E => Excise single peak"),
                                    tags$li("Shift + R => Remove global peak"),
                                    tags$li("Shift + 0 => Remove ALL peaks"),
                                    tags$li("Shift + Z => Save peak table"),
                                    tags$li("------------------"),
                                    tags$li("Shift + 1 => Extract MS from selection"),
                                    tags$li("Shift + 2 => Zoom in on selected MS portion"),
                                    tags$li("Shift + 3 => Subtract selection MS from current MS"),
                                    tags$li("Shift + 4 => Library search"),
                                    tags$li("Shift + 5 => Save current MS")
                                ),

                                verbatimTextOutput("key", placeholder = TRUE),

                                h5("Messages"),
                                verbatimTextOutput("message_window", placeholder = TRUE)
                            ),

                            mainPanel(

                                tags$script('
                                    $(document).on("keypress", function (e) {
                                       Shiny.onInputChange("keypress", e.which);
                                    });
                                '),

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

                                    rhandsontable::rHandsontableOutput("peak_table")
                                ),

                                textOutput("keepAlive")
                            )
                        )
                    )

                ## Set up the server

                    server <- function(input, output, session) {

                        messages <- reactiveVal(character())
                        append_message <- function(text) {
                            timestamp <- format(Sys.time(), "%H:%M:%S")
                            current <- isolate(messages())
                            messages(c(current, paste0("[", timestamp, "] ", text)))
                        }
                        append_message("App ready.")

                        output$keepAlive <- renderText({
                            req(input$count)
                            paste("\nstayin' alive", input$count)
                        })

                        ## Check keystoke value
                            output$key <- renderPrint({
                                input$keypress
                            })

                        output$message_window <- renderText({
                            logs <- messages()
                            if (length(logs) == 0) {
                                "No messages yet."
                            } else {
                                paste(logs, collapse = "\n")
                            }
                        })

                        read_ms_slice <- function(sample_path, sample_chrom_df, rt_start, rt_end) {
                            if (is.null(sample_chrom_df) || nrow(sample_chrom_df) == 0) {
                                return(data.frame())
                            }

                            has_row_lookup <- all(c("rt_first_row_in_raw", "rt_last_row_in_raw") %in% colnames(sample_chrom_df))

                            if (has_row_lookup) {
                                start_idx <- sample_chrom_df$rt_first_row_in_raw[which.min(abs(sample_chrom_df$rt - rt_start))]
                                end_idx <- sample_chrom_df$rt_last_row_in_raw[which.min(abs(sample_chrom_df$rt - rt_end))]

                                if (length(start_idx) == 0 || length(end_idx) == 0 || is.na(start_idx) || is.na(end_idx)) {
                                    has_row_lookup <- FALSE
                                } else {
                                    start_line <- as.integer(floor(start_idx) + 1L)
                                    end_line <- as.integer(floor(end_idx) + 1L)

                                    if (start_line < 2L) {
                                        start_line <- 2L
                                    }
                                    if (end_line < start_line) {
                                        end_line <- start_line
                                    }

                                    header <- names(data.table::fread(sample_path, nrows = 0, showProgress = FALSE))
                                    lines_to_read <- end_line - start_line + 1L

                                    slice <- data.table::fread(
                                        sample_path,
                                        skip = start_line - 1L,
                                        nrows = lines_to_read,
                                        header = FALSE,
                                        col.names = header,
                                        showProgress = FALSE
                                    )

                                    slice <- as.data.frame(slice)
                                }
                            }

                            if (!has_row_lookup) {
                                slice <- as.data.frame(data.table::fread(sample_path, showProgress = FALSE))
                            }

                            dplyr::filter(slice, rt > rt_start, rt < rt_end)
                        }

                        ## Keys to move chromatogram view - zoom in and out, move L and R
                            observeEvent(input$keypress, {

                                if (zooming) {
                                    if (is.null(x_axis_start_default)) {
                                        x_axis_start_default <<- x_axis_start
                                    }
                                    if (is.null(x_axis_end_default)) {
                                        x_axis_end_default <<- x_axis_end
                                    }
                                }

                                adjust_axes <- function(delta_start, delta_end) {
                                    rate_start <- delta_start
                                    rate_end <- delta_end

                                    if (zooming) {
                                        proposed_start <- x_axis_start_default + rate_start
                                        proposed_end <- x_axis_end_default + rate_end
                                        if (proposed_start >= proposed_end) {
                                            return(NULL)
                                        }
                                        x_axis_start_default <<- proposed_start
                                        x_axis_end_default <<- proposed_end
                                        x_axis_start <<- x_axis_start_default
                                        x_axis_end <<- x_axis_end_default
                                    } else {
                                        proposed_start <- x_axis_start + rate_start
                                        proposed_end <- x_axis_end + rate_end
                                        if (proposed_start >= proposed_end) {
                                            return(NULL)
                                        }
                                        x_axis_start <<- proposed_start
                                        x_axis_end <<- proposed_end
                                    }
                                }

                                step <- zoom_and_scroll_rate

                                if (input$keypress == 70) { adjust_axes(step, step) }      # Forward on "F"
                                if (input$keypress == 68) { adjust_axes(-step, -step) }   # Backward on "D"
                                if (input$keypress == 86) { adjust_axes(-step, step) }    # Wider on "V"
                                if (input$keypress == 67) { adjust_axes(step, -step) }    # Closer on "C"
                            })

                        ## Save manual changes to table on "Z" (90) keystroke

                            observeEvent(input$keypress, {

                                if (input$keypress == 90 ) {

                                    ## Write out any modifications to peak table (i.e. sample IDs)
                                        
                                        hot = isolate(input$peak_table)
                                        if (!is.null(hot)) {
                                            writeMonolist(rhandsontable::hot_to_r(input$peak_table), peaks_monolist_path)
                                            print(peaks_monolist_path)
                                            append_message(paste0("Peak table saved to ", peaks_monolist_path))
                                        } else {
                                            append_message("No peak table edits available to save.")
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

                                                y <- prelim_baseline$min
                                                x <- prelim_baseline$rt
                                                valid_points <- which(!is.na(x) & !is.na(y))

                                                fallback_baseline <- if (all(is.na(chromatogram$tic))) 0 else min(chromatogram$tic, na.rm = TRUE)
                                                baseline_values <- rep(fallback_baseline, length(chromatogram$rt))

                                                if (length(valid_points) >= 2) {
                                                    baseline_values <- tryCatch(
                                                        approx(
                                                            x[valid_points],
                                                            y[valid_points],
                                                            xout = chromatogram$rt,
                                                            rule = 2
                                                        )$y,
                                                        error = function(e) {
                                                            baseline_values
                                                        }
                                                    )
                                                } else if (length(valid_points) == 1) {
                                                    baseline_values <- rep(y[valid_points], length(chromatogram$rt))
                                                }

                                                if (all(is.na(baseline_values))) {
                                                    baseline_values <- rep(fallback_baseline, length(chromatogram$rt))
                                                }

                                                baseline2 <- data.frame(
                                                    rt = chromatogram$rt,
                                                    y = baseline_values
                                                )
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

                                        if (zooming) {
                                            full_start <- min(chromatograms_updated$rt_rt_offset, na.rm = TRUE)
                                            full_end <- max(chromatograms_updated$rt_rt_offset, na.rm = TRUE)

                                            if (is.null(x_axis_start_default)) {
                                                x_axis_start_default <<- full_start
                                            }
                                            if (is.null(x_axis_end_default)) {
                                                x_axis_end_default <<- full_end
                                            }

                                            brush_info <- isolate(input$chromatogram_brush)

                                            if (!is.null(brush_info)) {
                                                selected_start <- suppressWarnings(as.numeric(brush_info$xmin))
                                                selected_end <- suppressWarnings(as.numeric(brush_info$xmax))

                                                if (!is.finite(selected_start) || !is.finite(selected_end)) {
                                                    selected_start <- full_start
                                                    selected_end <- full_end
                                                }

                                                if (selected_start == selected_end) {
                                                    selected_start <- selected_start - zoom_and_scroll_rate/2
                                                    selected_end <- selected_end + zoom_and_scroll_rate/2
                                                }

                                                if (selected_start < selected_end) {
                                                    x_axis_start <<- selected_start
                                                    x_axis_end <<- selected_end
                                                    x_axis_start_default <<- selected_start
                                                    x_axis_end_default <<- selected_end
                                                }
                                            } else {
                                                x_axis_start_default <<- full_start
                                                x_axis_end_default <<- full_end
                                                x_axis_start <<- full_start
                                                x_axis_end <<- full_end
                                            }
                                        } else {
                                            x_axis_start <<- min(chromatograms_updated$rt_rt_offset, na.rm = TRUE)
                                            x_axis_end <<- max(chromatograms_updated$rt_rt_offset, na.rm = TRUE)
                                        }

                                        ## Prepare chromatograms for plotting

                                            chromatograms_for_plot <- dplyr::filter(
                                                chromatograms_updated,
                                                rt_rt_offset > x_axis_start & rt_rt_offset < x_axis_end
                                            )

                                            if (nrow(chromatograms_for_plot) == 0) {
                                                append_message("No chromatogram data within the current zoom window. Resetting view.")
                                                x_axis_start <<- min(chromatograms_updated$rt_rt_offset, na.rm = TRUE)
                                                x_axis_end <<- max(chromatograms_updated$rt_rt_offset, na.rm = TRUE)
                                                if (zooming) {
                                                    x_axis_start_default <<- x_axis_start
                                                    x_axis_end_default <<- x_axis_end
                                                }
                                                chromatograms_for_plot <- chromatograms_updated
                                            }

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

                                                        chromatograms_for_plot <- dplyr::filter(chromatograms_for_plot, rt_rt_offset > x_axis_start & rt_rt_offset < x_axis_end)
                                                        peak_table <- dplyr::filter(peak_table, peak_start_rt_offset > x_axis_start & peak_end_rt_offset < x_axis_end)

                                                    ## Make chromatogram plot object

                                                        # Make their labels easy to read
                                                            facet_labels <- sub("\\.csv$", "", basename(chromatograms_for_plot$path_to_cdf_csv), ignore.case = TRUE)
                                                            names(facet_labels) <- chromatograms_for_plot$path_to_cdf_csv

                                                            p <-    ggplot() + 
                                                                geom_line(data = chromatograms_for_plot, mapping = aes(x = rt_rt_offset, y = baseline), color = "grey") +
                                                                geom_line(data = chromatograms_for_plot, mapping = aes(x = rt_rt_offset, y = tic)) +
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
                                                                chromatograms_for_plot[chromatograms_for_plot$path_to_cdf_csv == peak_table[peak,]$path_to_cdf_csv,], 
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

                                                        # Make their labels easy to read
                                                            facet_labels <- sub("\\.csv$", "", basename(chromatograms_for_plot$path_to_cdf_csv), ignore.case = TRUE)
                                                            names(facet_labels) <- chromatograms_for_plot$path_to_cdf_csv

                                                        p <-  ggplot() + 
                                                                geom_line(data = chromatograms_for_plot, mapping = aes(x = rt_rt_offset, y = baseline), color = "grey") +
                                                                geom_line(data = chromatograms_for_plot, mapping = aes(x = rt_rt_offset, y = tic)) +
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
                                    append_message("Chromatograms updated.")
                                }
                            

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

                                        output$peak_table <- rhandsontable::renderRHandsontable(rhandsontable::rhandsontable({
                                            peak_table <- read.csv(peaks_monolist_path)
                                            peak_table
                                        }))
                                    }
                            })

                        ## [MS1 extract ("shift+1"), update ("shift+2"), subtract ("shift+3"), library search ("shift+4" disabled), save ("shift+5")]
                            
                            observeEvent(input$keypress, {

                                ## If "shift+1", MS from chromatogram brush -> MS_out_1
                                    if( input$keypress == 33 ) {

                                        if (is.null(input$chromatogram_brush)) {
                                            cat("No chromatogram selection to extract MS from.\n")
                                            return()
                                        }

                                        ret_start_MS <- min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
                                        ret_end_MS <- max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
                                        sample_name_MS <- as.character(brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1])

                                        chromatogram_updated_MS <- dplyr::filter(chromatograms_updated, path_to_cdf_csv == sample_name_MS)

                                        if (nrow(chromatogram_updated_MS) == 0) {
                                            cat("Could not locate chromatogram rows for the selected sample.\n")
                                            return()
                                        }

                                        framedDataFile <- isolate(read_ms_slice(sample_name_MS, chromatogram_updated_MS, ret_start_MS, ret_end_MS))

                                        if (nrow(framedDataFile) == 0) {
                                            cat("No mass spectrum extracted for the selected region.\n")
                                            return()
                                        }

                                        framedDataFile$mz <- round(framedDataFile$mz, 1)
                                        MS_out_1 <- framedDataFile %>%
                                            dplyr::group_by(mz) %>%
                                            dplyr::summarize(intensity = sum(intensity), .groups = "drop") %>%
                                            as.data.frame()

                                        MS_out_1 <<- MS_out_1

                                    }

                                ## If "shift+3" subtract chromatogram brush from MS_out_1 -> MS_out_1

                                    if ( input$keypress == 35 ) {

                                        if (!exists("MS_out_1")) {
                                            cat("No mass spectrum extracted yet.\n")
                                            return()
                                        }

                                        if (is.null(input$chromatogram_brush)) {
                                            cat("No chromatogram selection to subtract.\n")
                                            return()
                                        }

                                        ret_start_sub <- min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
                                        ret_end_sub <- max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
                                        sample_name_sub <- as.character(brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1])

                                        chromatogram_updated_sub <- dplyr::filter(chromatograms_updated, path_to_cdf_csv == sample_name_sub)

                                        framedDataFile_to_subtract <- isolate(read_ms_slice(sample_name_sub, chromatogram_updated_sub, ret_start_sub, ret_end_sub))

                                        if (nrow(framedDataFile_to_subtract) == 0) {
                                            cat("Nothing to subtract for the selected region.\n")
                                            return()
                                        }

                                        framedDataFile_to_subtract$mz <- round(framedDataFile_to_subtract$mz, 1)
                                        framedDataFile_to_subtract <- framedDataFile_to_subtract %>%
                                            dplyr::group_by(mz) %>%
                                            dplyr::summarize(intensity = sum(intensity), .groups = "drop") %>%
                                            as.data.frame()

                                        subtraction <- dplyr::left_join(
                                            MS_out_1,
                                            framedDataFile_to_subtract,
                                            by = "mz",
                                            suffix = c("_orig", "_sub")
                                        )
                                        subtraction$intensity_sub[is.na(subtraction$intensity_sub)] <- 0

                                        updated_MS <- subtraction %>%
                                            dplyr::transmute(
                                                mz = mz,
                                                intensity = pmax(intensity_orig - intensity_sub, 0)
                                            )

                                        MS_out_1 <<- as.data.frame(updated_MS)

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
                                        message("Mass spectrum library lookup is disabled in this build.\n")
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

    message("Done!")
