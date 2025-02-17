#### analyzeGCMSdata4
            
    analyzeGCMSdata4 <- function(
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
                                framedDataFile <- drop_na(framedDataFile)

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

                theme = shinythemes::shinytheme("yeti"),

                sidebarLayout(

                    sidebarPanel(
                        style = "overflow-y: auto; max-height: 90vh;",

                        h4("Busta Lab GCMS Analysis App"),

                        img(
                            src = "https://raw.githubusercontent.com/thebustalab/thebustalab.github.io/refs/heads/master/logo.png",
                            height = "260px", 
                            style = "border-radius:8px; margin-bottom:5px;"
                        ),

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

                        # Message window
                            h5("Messages"),
                            verbatimTextOutput("message_window", placeholder = TRUE),
                    ),

                    mainPanel(

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

                            tabPanel("Peak Detection",
                              fluidRow(
                                column(3,
                                  sliderInput("peakStartSlope", "Peak Start Slope",
                                              min = 1, max = 100, value = 10, step = 1),
                                  sliderInput("apexSlope", "Apex Slope",
                                              min = -20, max = 20, value = 0, step = 1),
                                  sliderInput("peakEndSlope", "Peak Tail Slope",
                                              min = 1, max = 100, value = 5, step = 1),
                                  sliderInput("minSignalAboveBaseline", "Min Signal Above Baseline",
                                              min = 0, max = 50, value = 2, step = 1),
                                  sliderInput("minPeakArea", "Minimum Peak Area",
                                              min = 0, max = 1e7, value = 1e5, step = 1e4)
                                ),
                                column(9,
                                  h4("Press Shift + 6 to auto-detect peaks using these settings."),
                                  p("peakStartSlope: how steep the slope must be to call a peak start."),
                                  p("Once derivative goes < 0 => apex is reached."),
                                  p("peakEndSlope: how negative the slope can stay before we say the peak is done."),
                                  p("minSignalAboveBaseline: how close to baseline to consider the tail done."),
                                  p("minPeakArea: if the integrated area is below this, ignore the peak."),
                                  plotOutput("peak_detection_plot", height = "600px")
                                )
                              )
                            ),


                            tabPanel("MS Library",

                                verticalLayout(
                                    fluidRow(
                                        column(3,
                                            actionButton("extract_ms", "Extract Mass Spectra"),
                                            br(), br(),
                                            actionButton("classify_ms", "Classify Spectra")
                                        ),
                                        column(9,
                                            DT::dataTableOutput("ms_results")
                                        )
                                    ),
                                    
                                    plotOutput(
                                        outputId = "massSpectrumLookup",
                                        height = 1200
                                    )
                                )
                            )
                        )
                    )
                )
            )

        ## SET UP SERVER

            server <- function(input, output, session) {

                ## Store console info:
                    message_data <- reactiveVal("")

                    # Helper function to append a new line
                    logMessage <- function(msg) {
                        old_val <- message_data()
                        new_val <- paste(old_val, msg, sep = "\n")
                        message_data(new_val)
                    }

                    # Render them in the message window
                    output$message_window <- renderText({
                        message_data()
                    })

                ## Set up MS values
                    ms_data_reactive <- reactiveVal(NULL)
                    ms_results_reactive <- reactiveVal(NULL)

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

                ## Set up null values for peak detection
                    peakDetectionPlotData <- reactiveVal(NULL)

                        output$peak_detection_plot <- renderPlot({
                        # If there's no data yet, do nothing
                        df <- peakDetectionPlotData()
                        validate(
                          need(!is.null(df), "No random peak subset to display yet. Run Shift+6 to detect peaks.")
                        )
                        
                        # Plot them in facet_wrap:
                        ggplot(df, aes(x = rt, y = abundance, color = ion)) +
                          geom_line() +
                          geom_line(aes(y = baseline), color = "gray50") +
                          geom_ribbon(aes(ymin = baseline, ymax = abundance, fill = peak_label, group = peak_label),
                                      alpha = 0.2, color = NA) +
                          facet_wrap(~peak_label, scales = "free") +
                          theme_bw() +
                          labs(title = "Randomly Sampled Peak Regions",
                               x = "Scan (RT)", y = "Abundance (counts)") +
                          guides(fill = "none")
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

                ## Delete peaks (Shift+0)
                    observeEvent(input$keypress, {
                      # SHIFT+0 => `)` => ASCII code 41
                      if (input$keypress == 41) {
                        cat("Removing ALL peaks...\n")
                        if (file.exists("peaks_monolist.csv")) {
                          # Overwrite with empty table
                          empty_peaks <- data.frame(
                            peak_start = numeric(),
                            peak_end   = numeric(),
                            peak_ID    = character(),
                            path_to_cdf_csv = character(),
                            area       = numeric(),
                            stringsAsFactors = FALSE
                          )
                          write.table(
                            empty_peaks,
                            "peaks_monolist.csv",
                            row.names = FALSE,
                            col.names = TRUE,
                            sep = ","
                          )
                          cat("All peaks removed.\n")
                          
                          # ALSO CLEAR THE facet-plot reactiveVal
                          peakDetectionPlotData(NULL)
                          
                        } else {
                          cat("No peaks_monolist.csv file exists, so nothing to remove.\n")
                        }
                      }
                    })

                ## Peak detection
                    observeEvent(input$keypress, {
                      if (input$keypress == 94) {  # '^' = Shift+6
                        peakDetectionPlotData(NULL)
                        cat("Starting custom peak detection...\n")
                        logMessage("Starting custom peak detection...\n")

                        # 1) Read thresholds from your sliders
                        peakStartSlope         <- input$peakStartSlope
                        apexSlope              <- input$apexSlope    # might be zero, used for crossing detection
                        peakTailSlope          <- input$peakEndSlope # how negative it can be before flattening
                        minSignalAboveBaseline <- input$minSignalAboveBaseline
                        minPeakArea            <- input$minPeakArea  # default = 100,000

                        # 2) Load any existing peaks (so we can append to them)
                        if (file.exists("peaks_monolist.csv")) {
                          peak_data_existing <- read.csv("peaks_monolist.csv", stringsAsFactors = FALSE)
                        } else {
                          peak_data_existing <- data.frame(
                            peak_start       = numeric(),
                            peak_end         = numeric(),
                            peak_ID          = character(),
                            path_to_cdf_csv  = character(),
                            area             = numeric(),
                            stringsAsFactors = FALSE
                          )
                        }

                        peak_data_new <- data.frame(
                          peak_start       = numeric(),
                          peak_end         = numeric(),
                          peak_ID          = character(),
                          path_to_cdf_csv  = character(),
                          area             = numeric(),
                          stringsAsFactors = FALSE
                        )

                        # 3) Ensure chromatograms_updated exists
                        if (!exists("chromatograms_updated")) {
                          cat("No chromatogram data loaded yet. Press Q first to update.\n")
                          return()
                        }

                        # 4) If you have sample filters
                        samples_monolist <- read.csv("samples_monolist.csv", stringsAsFactors = FALSE)
                        if (length(samples_monolist_subset) > 0) {
                          samples_monolist <- samples_monolist[samples_monolist_subset, ]
                        }

                        # 5) Loop through each sample’s TIC
                        for (this_file in unique(samples_monolist$path_to_cdf_csv)) {

                          # Pull data for this sample
                          this_chrom <- dplyr::filter(chromatograms_updated, path_to_cdf_csv == this_file)
                          # TIC only
                          tic_data <- dplyr::filter(this_chrom, ion == 0)
                          if (nrow(tic_data) < 3) next

                          # Sort by RT
                          tic_data <- tic_data[order(tic_data$rt), ]

                          # Subtract baseline
                          baseline_data <- dplyr::filter(this_chrom, ion == "baseline")
                          baseline_data <- baseline_data[order(baseline_data$rt), ]

                          merged_data <- merge(
                            tic_data, 
                            baseline_data[, c("rt", "abundance")],
                            by = "rt", 
                            suffixes = c("", ".bl")
                          )
                          merged_data$net_abundance <- merged_data$abundance - merged_data$abundance.bl

                          # Derivative of net_abundance
                          d1 <- c(0, diff(merged_data$net_abundance))

                          # 6) State machine
                          outsidePeak <- 0
                          ascending   <- 1
                          descending  <- 2

                          currentState <- outsidePeak
                          peak_start_idx <- NA

                          for (i in seq_along(d1)) {
                            if (currentState == outsidePeak) {
                              # Start if derivative > peakStartSlope
                              if (d1[i] > peakStartSlope) {
                                currentState <- ascending
                                peak_start_idx <- i
                              }

                            } else if (currentState == ascending) {
                              # If derivative crosses below apexSlope => apex reached => begin descending
                              if (d1[i] < apexSlope) {
                                currentState <- descending
                              }

                            } else if (currentState == descending) {
                              # If slope is now flattening out or net_abundance is near baseline => end peak
                              # i.e. derivative > -peakTailSlope OR net_abundance < minSignalAboveBaseline
                              if (d1[i] > -peakTailSlope || merged_data$net_abundance[i] < minSignalAboveBaseline) {
                                
                                # Record the peak
                                rt_start <- merged_data$rt[peak_start_idx]
                                rt_end   <- merged_data$rt[i]
                                area_val <- sum(merged_data$net_abundance[peak_start_idx:i])
                                if (area_val < 0) area_val <- 0

                                # Apply the minPeakArea filter
                                if (area_val >= minPeakArea) {
                                  peak_data_new <- rbind(
                                    peak_data_new,
                                    data.frame(
                                      peak_start      = rt_start,
                                      peak_end        = rt_end,
                                      peak_ID         = "auto_detected",
                                      path_to_cdf_csv = this_file,
                                      area            = area_val,
                                      stringsAsFactors = FALSE
                                    )
                                  )
                                }
                                
                                # Reset
                                currentState <- outsidePeak
                                peak_start_idx <- NA
                              }
                            }
                          }

                          # If we exit while still in ascending or descending, end at the last data point
                          if (currentState != outsidePeak && !is.na(peak_start_idx)) {
                            last_idx <- nrow(merged_data)
                            rt_start <- merged_data$rt[peak_start_idx]
                            rt_end   <- merged_data$rt[last_idx]
                            area_val <- sum(merged_data$net_abundance[peak_start_idx:last_idx])
                            if (area_val < 0) area_val <- 0

                            if (area_val >= minPeakArea) {
                              peak_data_new <- rbind(
                                peak_data_new,
                                data.frame(
                                  peak_start      = rt_start,
                                  peak_end        = rt_end,
                                  peak_ID         = "auto_detected_unfinished",
                                  path_to_cdf_csv = this_file,
                                  area            = area_val,
                                  stringsAsFactors = FALSE
                                )
                              )
                            }
                          }

                        } # end for each file

                        # 7) Append the new peaks to existing, write out
                            if (nrow(peak_data_new) > 0) {
                              peak_data_combined <- rbind(peak_data_existing, peak_data_new)
                              write.table(
                                peak_data_combined,
                                file = "peaks_monolist.csv",
                                row.names = FALSE,
                                col.names = TRUE,
                                sep = ","
                              )
                              cat("Auto-detected", nrow(peak_data_new), "peaks (passed minArea filter).\n")
                              logMessage(paste("Auto-detected", nrow(peak_data_new), "peaks (passed minArea filter).\n"))
                            } else {
                              cat("No new peaks found above threshold or min area.\n")
                              # No new peaks => might not want to do the random plot:
                              return()
                            }

                            ##
                            ## NEW: Build random subset of ~30 peaks for facet plotting
                            ##
                            # 1) Reload final table from disk (or just reuse `peak_data_combined`)
                            peak_table_final <- read.csv("peaks_monolist.csv")

                            # 2) If there are more than 30 peaks, sample 30
                            peak_table_sampled <- peak_table_final %>%
                              dplyr::arrange(desc(area)) %>%
                              dplyr::slice_head(n=30)

                            # 3) For each sampled peak, gather the chromatogram data
                            all_peaks_data <- list()

                            # Make sure we have 'chromatograms_updated' loaded
                            if (!exists("chromatograms_updated")) {
                              cat("No chromatograms_updated found. Press Q to update.\n")
                              return()
                            }

                            for (i in seq_len(nrow(peak_table_sampled))) {
                              
                              rowi <- peak_table_sampled[i, ]
                              
                              # Subset data from chromatograms_updated
                              # Only the relevant file, plus RT slice
                              df_signal <- dplyr::filter(
                                chromatograms_updated,
                                path_to_cdf_csv == rowi$path_to_cdf_csv,
                                rt >= rowi$peak_start,
                                rt <= rowi$peak_end
                              )
                              
                              # If no data, skip
                              if (nrow(df_signal) == 0) next
                              
                              # Optionally keep only ion==0 and ion=="baseline" if you want
                              # so the plot isn't cluttered. Or keep all if you prefer.
                              # Let's keep both TIC (ion=0) and baseline (ion="baseline"):
                              df_signal <- df_signal[df_signal$ion %in% c(0, "baseline"), ]
                              
                              # Add a facet label, e.g. "Peak #1"
                              df_signal$peak_label <- paste0("Peak #", i)
                              
                              all_peaks_data[[i]] <- df_signal
                            }

                            # 4) Combine into one big data frame
                            if (length(all_peaks_data) > 0) {
                              all_peaks_data <- dplyr::bind_rows(all_peaks_data)
                              
                              # To help the plotting, let's define a separate 'baseline' column
                              # For baseline rows, abundance is already the baseline,
                              # For TIC rows, let's store the baseline in a new column. 
                              # We'll match them by RT if we want. Or we can do a simpler approach:
                              
                              # For convenience, let's do a small join to assign a 'baseline' column
                              # to the entire data, so we can do a ribbon from baseline->abundance.
                              
                              df_baselines <- dplyr::filter(all_peaks_data, ion=="baseline") %>%
                                dplyr::select(rt, path_to_cdf_csv, peak_label, baseline=abundance)
                              
                              # For TIC rows
                              df_signal_only <- dplyr::filter(all_peaks_data, ion==0)
                              
                              # left_join by path_to_cdf_csv, rt, peak_label
                              df_signal_joined <- dplyr::left_join(
                                df_signal_only,
                                df_baselines,
                                by = c("rt", "path_to_cdf_csv", "peak_label")
                              )
                              
                              # For baseline rows themselves, we can keep them separate or combine them,
                              # but the simplest might be: "baseline" is the same as "abundance" for baseline.
                              # Let’s store them in the same df, though you can do separate geoms if you prefer.
                              df_baselines$baseline <- df_baselines$baseline # unchanged
                              df_baselines$ion <- "baseline"
                              # We'll rename abundance to something else to not confuse the ribbon
                              df_baselines$abundance <- df_baselines$baseline
                              
                              # Combine them
                              df_plot <- dplyr::bind_rows(df_signal_joined, df_baselines)
                              
                              # 5) Store in reactiveVal => triggers the plot
                              peakDetectionPlotData(df_plot)
                            } else {
                              # If no peaks or no data slices, set to NULL
                              peakDetectionPlotData(NULL)
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

                                        # 1. Make empty containers
                                        all_ribbons <- list()
                                        all_vlines  <- list()
                                        all_labels  <- list()


                                        # 2. Loop over peaks, but only assemble data frames
                                        if (nrow(peak_table) > 0) {
                                            for (peak in 1:nrow(peak_table)) {
                                                # Filter chromatogram data for this peak
                                                signal_for_this_peak <- dplyr::filter(
                                                    chromatograms_updated,
                                                    path_to_cdf_csv == peak_table[peak, ]$path_to_cdf_csv,
                                                    rt_rt_offset > peak_table[peak, ]$peak_start_rt_offset,
                                                    rt_rt_offset < peak_table[peak, ]$peak_end_rt_offset
                                                )

                                                # Only proceed if there is valid data
                                                if (nrow(signal_for_this_peak) > 0) {
                                                    # Assign peak number
                                                    signal_for_this_peak$peak_number_within_sample <- 
                                                        peak_table$peak_number_within_sample[peak]

                                                    # Extract ribbon data
                                                    ribbon <- dplyr::filter(signal_for_this_peak, ion == 0)
                                                    ribbon$baseline <- dplyr::filter(signal_for_this_peak, ion == "baseline")$abundance

                                                    # Store data in lists
                                                    all_ribbons[[peak]] <- ribbon
                                                    all_vlines[[peak]]  <- signal_for_this_peak[1, ]

                                                    # Create label data
                                                    label_df <- dplyr::filter(signal_for_this_peak, ion == 0) %>%
                                                        dplyr::summarize(
                                                            peak_number_within_sample = peak_number_within_sample[1],
                                                            x = median(rt_rt_offset),
                                                            y = max(abundance),
                                                            path_to_cdf_csv = path_to_cdf_csv[1]
                                                        )
                                                    all_labels[[peak]] <- label_df
                                                }
                                            }
                                        
                                            # 3. Combine all stored data frames
                                            all_ribbons <- dplyr::bind_rows(all_ribbons)
                                            all_vlines  <- dplyr::bind_rows(all_vlines)
                                            all_labels  <- dplyr::bind_rows(all_labels)

                                            # 4. Add a single set of ggplot layers
                                            chromatogram_plot <- chromatogram_plot +
                                                geom_vline(
                                                    data = all_vlines, 
                                                    mapping = aes(xintercept = rt_rt_offset), 
                                                    alpha = 0.3
                                                ) +
                                                geom_ribbon(
                                                    data = all_ribbons,
                                                    mapping = aes(
                                                        x = rt_rt_offset, 
                                                        ymax = abundance, 
                                                        ymin = baseline, 
                                                        fill = peak_number_within_sample,
                                                        group = peak_number_within_sample
                                                    ),
                                                    alpha = 0.8
                                                ) +
                                                geom_text(
                                                    data = all_labels,
                                                    mapping = aes(
                                                        label = peak_number_within_sample, 
                                                        x = x, 
                                                        y = y
                                                    ),
                                                    color = "black"
                                                )
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

                ## Single-peak remove => SHIFT+E => ASCII code 69
                      observeEvent(input$keypress, {
                        if( input$keypress == 69 ) { # Update on "E" for "excise"

                          cat("Excising single peak...\n")
                          if ( !is.null(input$chromatogram_brush )) {

                            peak_points <- brushedPoints(chromatograms_updated, input$chromatogram_brush)
                            selection_start = min(peak_points$rt)
                            selection_end   = max(peak_points$rt)
                            path_to_cdf_csv = peak_points$path_to_cdf_csv[1]
                            
                            peak_table <- read.csv("peaks_monolist.csv")

                            # Remove any peak that fully resides within [selection_start, selection_end]
                            # for that single cdf
                            peak_table <- peak_table[!
                              apply(cbind(
                                peak_table$peak_start > selection_start,
                                peak_table$peak_end < selection_end,
                                peak_table$path_to_cdf_csv == as.character(path_to_cdf_csv)
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
                            cat("Excised single peak.\n")
                          } else {
                            cat("No brush selection to excise.\n")
                          }
                        }
                      })

                ## Global peak remove => SHIFT+R => ASCII code 82
                      observeEvent(input$keypress, {
                        if (input$keypress == 82) { # "R" for "Remove" globally

                          cat("Removing selected peaks GLOBALLY...\n")
                          if ( !is.null(input$chromatogram_brush )) {

                            peak_points <- brushedPoints(chromatograms_updated, input$chromatogram_brush)
                            selection_start = min(peak_points$rt)
                            selection_end   = max(peak_points$rt)

                            # read peak table
                            peak_table <- read.csv("peaks_monolist.csv")

                            # remove ANY peak in ANY file that fully resides in [start, end]
                            peak_table <- peak_table[!
                              apply(cbind(
                                peak_table$peak_start > selection_start,
                                peak_table$peak_end < selection_end
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
                            cat("Removed global peaks within brush.\n")
                          } else {
                            cat("No brush selection to remove globally.\n")
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
                                            MS1_low_x_limit <- 0; MS1_high_x_limit <- 1200; MS1_high_y_limit <- 110
                                        } else {
                                            MS1_low_x_limit <- isolate(min(brushedPoints(MS_out_1, input$massSpectra_1_brush)$mz))
                                            MS1_high_x_limit <- isolate(max(brushedPoints(MS_out_1, input$massSpectra_1_brush)$mz))
                                            MS1_high_y_limit <- max(dplyr::filter(MS_out_1, mz > MS1_low_x_limit & mz < MS1_high_x_limit)$intensity) + 8
                                        }
                                        if (MS1_low_x_limit %in% c(Inf, -Inf) | MS1_high_x_limit %in% c(Inf, -Inf) | MS1_high_y_limit %in% c(Inf, -Inf)) {
                                            MS1_low_x_limit <- 0; MS1_high_x_limit <- 1200; MS1_high_y_limit <- 110
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
                                    
                                # message("Searching reference library for unknown spectrum...\n")

                                #     ## Round to nominal mass spectrum
                                #         MS_out_1$mz <- floor(MS_out_1$mz)
                                #         MS_out_1 %>%
                                #             group_by(mz) %>%
                                #             dplyr::summarize(intensity = sum(intensity)) -> MS_out_1
                                #         MS_out_1$intensity <- MS_out_1$intensity/max(MS_out_1$intensity)*100

                                #     ## Add zeros for mz values missing from unknown spectrum if necessary
                                #         mz_missing <- seq(1, 800, 1)[!seq(1, 800, 1) %in% MS_out_1$mz]
                                #         if (length(mz_missing) > 0) {
                                #             MS_out_1 <- rbind(MS_out_1, data.frame(mz = mz_missing, intensity = 0))
                                #             MS_out_1 <- MS_out_1[order(MS_out_1$mz),]
                                #         }

                                #     ## Bind it with metadata
                                #         unknown <- cbind(
                                #             data.frame(
                                #                 Accession_number = "unknown",
                                #                 Compound_systematic_name = "unknown",
                                #                 Compound_common_name = "unknown",
                                #                 SMILES = NA,
                                #                 Source = "unknown"
                                #             ),
                                #             t(MS_out_1$intensity)
                                #         )
                                #         colnames(unknown)[6:805] <- paste("mz_", seq(1,800, 1), sep = "")

                                #     ## Bind it to the library and run the lookup
                                #         lookup_data <- rbind(busta_spectral_library, unknown)
                                        
                                #         hits <- runMatrixAnalysis(
                                #             data = lookup_data,
                                #             analysis = c("hclust"),
                                #             column_w_names_of_multiple_analytes = NULL,
                                #             column_w_values_for_multiple_analytes = NULL,
                                #             columns_w_values_for_single_analyte = colnames(lookup_data)[6:805],
                                #             columns_w_additional_analyte_info = NULL,
                                #             columns_w_sample_ID_info = c("Accession_number", "Compound_systematic_name"),
                                #             transpose = FALSE,
                                #             unknown_sample_ID_info = c("unknown_unknown"),
                                #             scale_variance = FALSE,
                                #             kmeans = "none",
                                #             na_replacement = "drop",
                                #             output_format = "long"
                                #         )

                                #     ## Find distance between tips

                                #         phylo <- runMatrixAnalysis(
                                #             data = hits[!is.na(hits$sample_unique_ID),],
                                #             analysis = c("hclust_phylo"),
                                #             column_w_names_of_multiple_analytes = "analyte_name",
                                #             column_w_values_for_multiple_analytes = "value",
                                #             columns_w_values_for_single_analyte = NULL,
                                #             columns_w_additional_analyte_info = NULL,
                                #             columns_w_sample_ID_info = c("Accession_number", "Compound_systematic_name"),
                                #             transpose = FALSE,
                                #             unknown_sample_ID_info = NULL,
                                #             scale_variance = FALSE,
                                #             kmeans = "none",
                                #             na_replacement = "drop",
                                #             output_format = "long"
                                #         )

                                #         distances_1 <- as.data.frame(cophenetic.phylo(phylo)[,colnames(cophenetic.phylo(phylo)) == "unknown_unknown"])
                                #         distances_2 <- data.frame(
                                #             sample_unique_ID = rownames(distances_1),
                                #             distance = distances_1[,1]
                                #         )
                                #         distances_2$distance <- normalize(distances_2$distance, old_min = min(distances_2$distance), old_max = max(distances_2$distance), new_min = 100, new_max = 0)

                                #     ## Make the bar data

                                #         hits[!is.na(hits$sample_unique_ID),] %>%
                                #             select(analyte_name, value, sample_unique_ID) %>%
                                #             unique() -> bars

                                #         bars$analyte_name <- gsub(".*_", "", bars$analyte_name)

                                #         bars %>%
                                #             group_by(sample_unique_ID) %>%
                                #             arrange(desc(value)) -> bar_labels

                                #         bar_labels <- bar_labels[1:110,]

                                #     ## Order bar data

                                #         bars$sample_unique_ID <- factor(
                                #             bars$sample_unique_ID, levels = distances_2$sample_unique_ID[order(distances_2$distance, decreasing = TRUE)]
                                #         )

                                #         distances_2$sample_unique_ID <- factor(
                                #             distances_2$sample_unique_ID, levels = distances_2$sample_unique_ID[order(distances_2$distance, decreasing = TRUE)]
                                #         )

                                #         bar_labels$sample_unique_ID <- factor(
                                #             bar_labels$sample_unique_ID, levels = distances_2$sample_unique_ID[order(distances_2$distance, decreasing = TRUE)]
                                #         )

                                #     ## Make the bar plot
                                        
                                #         plot <- ggplot() +
                                #             geom_col(data = bars, aes(x = as.numeric(as.character(analyte_name)), y = value)) +
                                #             geom_text(data = unique(select(bars, sample_unique_ID)), aes(label = sample_unique_ID, x = 400, y = 90), hjust = 0.5, size = 4) +
                                #             facet_grid(sample_unique_ID~.) +
                                #             theme_bw() +
                                #             scale_x_continuous(name = "m/z") +
                                #             scale_y_continuous(name = "Relative intensity (%)") +
                                #             geom_text(
                                #                     data = bar_labels,
                                #                     mapping = aes(
                                #                         x = as.numeric(as.character(analyte_name)),
                                #                         y = value + 5, label = analyte_name
                                #                     )
                                #                 ) +
                                #             geom_text(
                                #                     data = distances_2,
                                #                     mapping = aes(
                                #                         x = 800,
                                #                         y = 75,
                                #                         label = paste0(
                                #                             "Relative similarity to unknown: ",
                                #                             round(distance, 2), "%"
                                #                         ), hjust = 1
                                #                     )
                                #                 )

                                #         output$massSpectrumLookup <- renderPlot({plot})

                                # message("Done.\n")

                            }
                    
                        ## If "shift+5" save mass spectrum

                            # if( input$keypress == 37 ) {

                            #     # Do nothing if no MS extracted
                            #         if (!exists("MS_out_1")) {
                            #             cat("No mass spectrum extracted yet.\n")
                            #             return()
                            #         } else {
                            #             MS_out_1_to_write <- MS_out_1
                            #             MS_out_1_to_write$mz <- round(MS_out_1_to_write$mz)
                            #             MS_out_1_to_write %>% 
                            #                 group_by(mz) %>%
                            #                 dplyr::summarize(intensity = sum(intensity)) -> MS_out_1_to_write
                            #             MS_out_1_to_write$intensity <- MS_out_1_to_write$intensity*100/max(MS_out_1_to_write$intensity)
                            #             MS_out_1_to_write <- data.frame(
                            #                 Compound_common_name = NA,
                            #                 Compound_systematic_name = NA,
                            #                 SMILES = NA,
                            #                 Source = "Busta",
                            #                 mz = MS_out_1_to_write$mz,
                            #                 abu = MS_out_1_to_write$intensity
                            #             )
                            #             write_csv(MS_out_1_to_write, paste0(CDF_directory_path, "/selected_MS.csv"))
                            #         }
                                
                            # }

                        ## Extract and classify mass spectra on button press (shift + 4)
                                
                            if (input$keypress == 36) {
                                message("UPDATED!")
                                # Read in the peak table (ensure your file is in the correct format)
                                    peak_table <<- read.csv("peaks_monolist.csv", stringsAsFactors = FALSE)
                                    if(nrow(peak_table) == 0){
                                        cat("No peaks available.", type = "error")
                                        return()
                                    }

                                ### somethign to do with the duplicate names
                                ### also all predicted peaks ID should be prefixed with "prediction"
                                ### perhaps it should only output predictions on things labelled "unknown"?

                                # Loop over each detected peak
                                    ms_list <- list()
                                    for(i in 1:nrow(peak_table)){
                                        # message(paste0("\n*getting spectrum for peak ", i))
                                        peak <- peak_table[i, ]
                                        # Get the corresponding .CDF.csv file
                                        cdf_file <- as.character(peak$path_to_cdf_csv)
                                        if(file.exists(cdf_file)){
                                            # Read the CSV file (using data.table::fread for speed)
                                            df <- data.table::fread(cdf_file)
                                            # Filter rows that fall within the peak’s retention time window
                                            df_subset <- df[df$rt >= peak$peak_start & df$rt <= peak$peak_end, ]
                                            if(nrow(df_subset) > 0){
                                                # Round m/z values and average intensities by m/z
                                                df_agg <- df_subset %>%
                                                    dplyr::mutate(mz = round(mz, 0)) %>%
                                                    dplyr::group_by(mz) %>%
                                                    dplyr::summarize(intensity = mean(intensity, na.rm = TRUE))
                                                # Tag the result with a peak identifier
                                                df_agg$peak_id <- peak$peak_ID
                                                df_agg$peak_unique_id <- i
                                                ms_list[[length(ms_list) + 1]] <- df_agg
                                            }
                                        }
                                    }
                                    
                                    if(length(ms_list) == 0){
                                        message("No mass spectra extracted from peaks.")
                                        return()
                                    }
                                
                                # Combine all peaks into one data frame
                                    ms_data <<- dplyr::bind_rows(ms_list)
                                    ms_data_reactive(ms_data)
                                    message("Mass spectra extracted successfully.")

                                # When the user clicks "Classify Spectra"
                                    ms_data <- ms_data_reactive()
                                    if(is.null(ms_data)){
                                        message("Please extract mass spectra first.")
                                        return()
                                    }
                                
                                # Convert the long-format mass spectral data into a wide format
                                # so that each row corresponds to one peak and columns are m/z values.
                                ms_data %>%
                                    mutate(mz = round(mz)) %>%
                                    group_by(peak_unique_id, mz) %>%
                                    summarize(intensity = sum(intensity)) %>%
                                    bind_rows(data.frame(peak_unique_id = 0, mz = 0:1000, intensity = 0)) %>%
                                    pivot_wider(names_from = mz, values_from = intensity, values_fill = 0) ->> ms_wide
                                message("ms_wide made.")

                                predictions <<- as.character(predictWithModel(
                                    data = ms_wide[ms_wide$peak_unique_id != 0,][,-1],
                                    model_type = "random_forest_classification",
                                    model = readRDS("/project_data/shared/mass_spectral_library/busta_lab_rfc_model_v1.rds")
                                ))

                                message("predictions made:")
                                message(predictions)
                                
                                peak_table$peak_ID <- paste0("(predicted) ", predictions)
                                writeMonolist(peak_table, "peaks_monolist.csv")
                                
                                # # Render the results in the UI table.
                                # output$ms_results <- DT::renderDataTable({
                                #     ms_results_reactive()
                                # })
                                
                                message("Mass spectral classification completed.")
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

