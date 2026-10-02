#### Shiny port pool for jupyter = TRUE
####
#### nginx on host1 (bustalab.d.umn.edu) proxies 13 named subdomains to fixed
#### localhost ports. A Shiny app started inside a JupyterHub R session can claim
#### one; nginx then exposes it as https://shinyN.bustalab.d.umn.edu over HTTPS.
####
#### Defined here as a top-level helper so this file remains self-contained
#### (per the bustalabfunctions/AGENTS.md self-contained rule). The same
#### definition exists in phylochemistry.R; if you change the pool, keep both
#### in sync.

    .shiny_port_pool <- data.frame(
        port      = c(10101:10110, 10123, 10456, 10789),
        subdomain = c(paste0("shiny", 1:10), "shiny123", "shiny456", "shiny789"),
        stringsAsFactors = FALSE
    )

    .shiny_port_is_free <- function(port) {
        tryCatch({
            s <- httpuv::startServer("127.0.0.1", port,
                                     list(call = function(req) list(status = 200, body = "")))
            s$stop()
            TRUE
        }, error = function(e) FALSE)
    }

    pick_shiny_port <- function() {
        for (i in sample(nrow(.shiny_port_pool))) {
            p <- .shiny_port_pool$port[i]
            if (.shiny_port_is_free(p)) {
                return(list(
                    port = p,
                    url  = paste0("https://", .shiny_port_pool$subdomain[i], ".bustalab.d.umn.edu")
                ))
            }
        }
        NULL
    }


#### analyzeGCMSdata6
####
#### v6 = v5 plus the all-ion chromatogram view.
####
#### The TIC hides co-eluting compounds: two compounds under one hump sum to a
#### single smooth peak, and nothing in the trace says there are two. Drawing
#### every ion as its own thin, semi-transparent line exposes them, because the
#### ions of the two compounds apex at slightly different retention times and
#### carry different peak shapes. A pure peak draws as a tight bundle of lines;
#### an impure one splays. Set all_ion_view = TRUE, or switch "Display" in the
#### sidebar.
####
#### The same ion-selection machinery also feeds integration: `area_peak_ions`
#### sums only the ions that rise across a peak's own window, so a peak is no
#### longer charged for the column bleed and tailing neighbours sitting under
#### it. It is reported ALONGSIDE the v5 TIC `area`, never instead of it, so
#### existing numbers stay comparable.

    analyzeGCMSdata6 <- function(
        CDF_directory_path = getwd(),
        zoom_and_scroll_rate = 100,
        baseline_window = 400,
        samples_per_page = 20,
        x_axis_start_default = NULL,
        x_axis_end_default = NULL,
        path_to_reference_library = busta_spectral_library,
        samples_monolist_subset = NULL,
        ions = 0,
        all_ion_view = FALSE,
        all_ion_threshold = 0.02,
        all_ion_scaling = "sqrt",
        all_ion_max_points = 1200,
        integrate_peak_ions = TRUE,
        ms_classifier_path = "/project_data/shared/mass_spectral_library/ms_classifier.py",
        jupyter = FALSE
    ) {

        ## Detect CDF flavour: "ms" for ANDI-MS (full rt x mz cube),
        ## "chrom" for ANDI-Chrom single-detector (e.g. GC-TCD, GC-FID, GC-ECD),
        ## "unknown" otherwise. Used to fork the ingest path below.
        detect_cdf_type <- function(path) {
            nc <- ncdf4::nc_open(path)
            on.exit(ncdf4::nc_close(nc))
            vars <- names(nc$var)
            if (all(c("mass_values", "intensity_values", "scan_index") %in% vars)) {
                return("ms")
            }
            if ("ordinate_values" %in% vars) {
                return("chrom")
            }
            "unknown"
        }

        ## Read an ANDI-Chrom CDF and return a data frame matching the
        ## framedDataFile shape used downstream (mz, intensity, rt). A
        ## single-detector trace has no mz axis, so mz is a sentinel 0:
        ## the TIC pipeline sums by rt to recover the one signal channel,
        ## and any extracted-ion request returns empty (no spectrum exists).
        read_chrom_cdf_as_framed <- function(path) {
            nc <- ncdf4::nc_open(path)
            on.exit(ncdf4::nc_close(nc))
            ord <- as.numeric(ncdf4::ncvar_get(nc, "ordinate_values"))
            var_names <- names(nc$var)
            delay <- if ("actual_delay_time" %in% var_names) {
                as.numeric(ncdf4::ncvar_get(nc, "actual_delay_time"))
            } else 0
            interval <- if ("actual_sampling_interval" %in% var_names) {
                as.numeric(ncdf4::ncvar_get(nc, "actual_sampling_interval"))
            } else 1
            n <- length(ord)
            rt <- delay + (seq_len(n) - 1) * interval
            data.frame(mz = 0, intensity = ord, rt = rt)
        }


        ## ---- All-ion chromatogram support (new in v6) ----------------------
        ##
        ## Three things make an every-ion overlay affordable at 20 samples per
        ## page, which measurement on a 8434-scan run says it otherwise is not
        ## (~10M line segments, ~25 s per ggplot redraw):
        ##
        ##   * m/z is binned to NOMINAL mass. The extraction path above rounds
        ##     to 0.1, which on unit-resolution quadrupole data yields roughly
        ##     five times more channels than there are real ones.
        ##   * only ions that RISE inside the visible window are drawn. Column
        ##     bleed and persistent background are flat there, so they cost ink
        ##     and CPU without carrying information.
        ##   * retention time is decimated onto a pixel grid, keeping the
        ##     MAXIMUM in each bin so apexes survive. Nothing is gained by
        ##     drawing more x-points than the panel has pixels.
        ##
        ## Together those take the same render to ~1.5M segments and ~5 s, and a
        ## zoomed window to ~1 s.

            .allion_cache <- new.env(parent = emptyenv())
            .allion_cache_order <- character(0)   # least-recently-used first
            .allion_cache_max <- 40L

            ## ---- cache stamping -------------------------------------------------
            ## Every derived file in this folder used to be keyed on its source's
            ## NAME alone, which made replacing an acquisition invisible: drop a new
            ## run in under the same filename and the app kept plotting and
            ## integrating the OLD one, because <file>.CDF.csv already existed so it
            ## was never rebuilt and chromatograms.csv still held its rows. Deleting
            ## the .CDF.csv by hand did not help either, since chromatograms.csv was
            ## the copy actually being read. Each derived file now carries a .stamp
            ## sidecar holding its source's mtime and size.
            file_stamp <- function(path) {
                fi <- file.info(path)
                if (is.na(fi$size)) return(NA_character_)
                paste0(format(as.numeric(fi$mtime), digits = 17), ":", fi$size)
            }
            stamp_path <- function(derived) paste0(derived, ".stamp")
            stamp_is_current <- function(derived, source) {
                if (!file.exists(derived)) return(FALSE)
                want <- file_stamp(source)
                if (is.na(want)) return(TRUE)     # source gone: keep what we have
                sp <- stamp_path(derived)
                if (!file.exists(sp)) return(FALSE)
                identical(readLines(sp, warn = FALSE)[1], want)
            }
            write_stamp <- function(derived, source) {
                want <- file_stamp(source)
                if (!is.na(want)) writeLines(want, stamp_path(derived))
                invisible(NULL)
            }

            ## Write through a temp file in the same directory, then rename. An
            ## interrupted fwrite straight to the destination left a TRUNCATED file
            ## behind, which file.exists() then accepted for ever.
            fwrite_atomic <- function(x, path) {
                tmp <- paste0(path, ".tmp", Sys.getpid())
                data.table::fwrite(x, tmp)
                if (!file.rename(tmp, path)) {
                    file.copy(tmp, path, overwrite = TRUE)
                    unlink(tmp)
                }
                invisible(path)
            }

            allion_store_path <- function(cdf_csv) paste0(cdf_csv, ".allions.csv")

            ## One row per (rt, nominal m/z) with non-zero signal. Built once per
            ## sample from the framed CSV and cached on disk beside it. This is
            ## deliberately NOT folded into chromatograms.csv: that file is
            ## shared by every sample, and appending the full cube to it would
            ## turn a small index into hundreds of megabytes.
            build_allion_store <- function(cdf_csv, force = FALSE) {
                out <- allion_store_path(cdf_csv)
                if (!force && stamp_is_current(out, cdf_csv)) return(invisible(out))
                framed <- data.table::as.data.table(data.table::fread(cdf_csv))
                framed[, mz := round(as.numeric(mz))]
                framed <- framed[!is.na(mz) & !is.na(rt) & !is.na(intensity)]
                agg <- framed[, .(abundance = sum(intensity)), by = .(rt, mz)]
                ## Keep every single-detector scan. Dropping non-positive rows is the
                ## right sparsity win for MS counts (an empty m/z cell carries nothing),
                ## but TCD/FID/ECD ordinate values are SIGNED and routinely sit below
                ## zero between peaks -- pruning them punched holes in the store, so
                ## geom_line drew chords across the gaps and peak_ion_area summed only
                ## the positive part of each peak.
                agg <- agg[mz == 0 | abundance > 0]
                data.table::setorder(agg, rt, mz)
                fwrite_atomic(agg, out)
                write_stamp(out, cdf_csv)
                invisible(out)
            }

            ## In-memory cache: validity is checked against the store's stamp, so a
            ## rebuilt store is not shadowed by the copy already in RAM, and the
            ## cache is bounded rather than growing one m/z cube per sample opened
            ## for the life of the session.
            allion_cache_touch <- function(cdf_csv) {
                .allion_cache_order <<- c(setdiff(.allion_cache_order, cdf_csv), cdf_csv)
                while (length(.allion_cache_order) > .allion_cache_max) {
                    drop <- .allion_cache_order[1]
                    .allion_cache_order <<- .allion_cache_order[-1]
                    if (!is.null(.allion_cache[[drop]])) rm(list = drop, envir = .allion_cache)
                }
                invisible(NULL)
            }

            load_allions <- function(cdf_csv) {
                p <- allion_store_path(cdf_csv)
                fresh <- stamp_is_current(p, cdf_csv)
                if (!fresh) {
                    if (!is.null(.allion_cache[[cdf_csv]])) rm(list = cdf_csv, envir = .allion_cache)
                    build_allion_store(cdf_csv, force = TRUE)
                }
                if (!is.null(.allion_cache[[cdf_csv]])) {
                    allion_cache_touch(cdf_csv)
                    return(.allion_cache[[cdf_csv]])
                }
                d <- data.table::as.data.table(data.table::fread(p))
                data.table::setkey(d, rt)
                assign(cdf_csv, d, envir = .allion_cache)
                allion_cache_touch(cdf_csv)
                d
            }

            ## Ions worth drawing in [xs, xe], with a per-ion baseline removed.
            ## The baseline is that ion's 10th percentile INSIDE the window
            ## rather than a median of the flanks either side: it costs one pass
            ## instead of three, and it degrades gracefully as the window
            ## narrows. `threshold` is a fraction of the largest baseline-
            ## corrected ion in the window, so it means the same thing whether
            ## the window holds the base peak or a trace one.
            ## Returns NULL when nothing clears the threshold.
            ## `rt_shift` is this sample's rt_offset. The decimation grid is built in
            ## OFFSET-CORRECTED space and then shifted back, so every sample on the page
            ## snaps to the SAME centres. Without it the "one x grid" rule held only
            ## within a sample: the union across samples was as dense as the raw data,
            ## which made the ion map's tile width (computed from that union) several
            ## times too narrow and drew the map as thin stripes with white gaps.
            allion_window <- function(cdf_csv, xs, xe, threshold, max_points, rt_shift = 0) {
                d <- load_allions(cdf_csv)
                w <- d[rt >= xs & rt <= xe]
                if (nrow(w) == 0) return(NULL)
                w[, base := stats::quantile(abundance, 0.1), by = mz]
                w[, net := pmax(abundance - base, 0)]
                top <- w[, .(peak = max(net)), by = mz]
                if (nrow(top) == 0 || max(top$peak) <= 0) return(NULL)
                keep <- top[peak > threshold * max(top$peak), mz]
                w <- w[mz %in% keep]
                if (nrow(w) == 0) return(NULL)
                if (is.finite(max_points) && data.table::uniqueN(w$rt) > max_points) {
                    ## Grid spans the REQUESTED window in shifted space, not this
                    ## sample's own min/max, so the centres are identical for every
                    ## sample regardless of where its data starts or what its offset is.
                    br <- seq(xs + rt_shift, xe + rt_shift, length.out = as.integer(max_points) + 1L) - rt_shift
                    if (!all(is.finite(br)) || length(br) < 2L || br[length(br)] <= br[1]) {
                        br <- seq(min(w$rt), max(w$rt), length.out = as.integer(max_points) + 1L)
                    }
                    w[, bin := findInterval(rt, br, rightmost.closed = TRUE)]
                    w <- w[w[, .I[which.max(net)], by = .(mz, bin)]$V1]
                    ## Snap to the bin centre so every ion shares ONE x grid.
                    ## Without this each ion keeps its own within-bin rt, the
                    ## union of rt across ions is as dense as the undecimated
                    ## data (so nothing is saved), and the ion map has no
                    ## regular grid to tile. Apexes move by at most half a bin,
                    ## which at these bin counts is sub-pixel.
                    centres <- (br[-length(br)] + br[-1L]) / 2
                    w[, rt := centres[pmin(pmax(bin, 1L), length(centres))]]
                    w[, bin := NULL]
                }
                w[, base := NULL]
                w[]
            }

            ## y transform for the overlay, plus the per-line opacity that goes
            ## with it. "norm" divides each ion by its own maximum in the window,
            ## which is what actually makes co-elution obvious: every ion becomes
            ## a unit-height peak, so a two-scan apex shift is unmissable. It
            ## also destroys quantitation and would draw a noise-level ion as
            ## tall as the base peak -- hence opacity carrying absolute
            ## abundance in that mode, and only in that mode.
            allion_scale <- function(w, mode) {
                w <- data.table::copy(w)
                if (identical(mode, "norm")) {
                    w[, peak := max(net), by = mz]
                    w[, y := ifelse(peak > 0, net / peak, 0)]
                    top <- max(w$peak)
                    w[, opacity := 0.10 + 0.65 * (peak / top)^0.35]
                    w[, peak := NULL]
                } else {
                    w[, y := switch(mode,
                        raw  = net,
                        sqrt = sqrt(net),
                        log  = log10(pmax(net, 1)),
                        net)]
                    w[, opacity := 0.35]
                }
                w[]
            }

            allion_y_label <- function(mode) switch(mode,
                raw  = "Abundance (counts, baseline removed)",
                sqrt = "sqrt abundance (baseline removed)",
                log  = "log10 abundance (baseline removed)",
                norm = "Abundance, each ion scaled to its own max",
                "Abundance")

            ## Shiny infers the brush's x/y variable names from the FIRST LAYER's
            ## aes, then brushedPoints() looks those names up as columns of the
            ## data frame it is handed -- always chromatograms_updated here.
            ## That holds in the TIC view (x = rt_rt_offset, y = abundance, both
            ## real columns) and breaks in the v6 views, where the first layer
            ## maps y to the transformed `y` (all ions) or to `mz` (ion map).
            ## Neither is a column of chromatograms_updated, so every one of the
            ## ~18 brush-driven actions -- Shift+Q's zoom, Shift+A/G peak add,
            ## Shift+1 MS extraction -- would error the moment the view changed.
            ##
            ## So: delegate to brushedPoints() when the mapping is usable, and
            ## otherwise select on x and the facet alone. The brush's y range in
            ## those modes is in sqrt / normalised / m-z units and means nothing
            ## against raw counts, so ignoring it is not a degradation -- it is
            ## the only correct reading.
            ## Where the plot PANEL sits inside the rendered image, as fractions
            ## of image width.
            ##
            ## Needed because Shiny's fallback coordmap (see brushed_chromatogram
            ## below) reports the brush as a fraction of the whole IMAGE, while
            ## the data only occupies the panel -- the y-axis labels, facet strip
            ## and colour bar all sit outside it. Verified on host1: at 960 px the
            ## panel runs 0.041..0.914, at 1400 px 0.028..0.941, i.e. constant
            ## 39 px and 83 px margins. Ignoring that is a left shift plus a
            ## stretch, which is exactly the "not quite the right region" symptom.
            ##
            ## Absolute gtable widths convert straight to inches; the panel is the
            ## "null" column that soaks up whatever is left. Cached per image size
            ## because building a second grob is not free on a million-row plot.
            .panel_cache <- new.env(parent = emptyenv())

            ## Read a panel_params range whichever way this ggplot2 stores it.
            pp_range <- function(pp, axis) {
                v <- pp[[paste0(axis, ".range")]]
                if (is.null(v)) v <- tryCatch(pp[[axis]]$continuous_range, error = function(e) NULL)
                if (is.null(v)) v <- tryCatch(pp[[axis]]$dimension(), error = function(e) NULL)
                if (length(v) == 2 && all(is.finite(v)) && v[2] > v[1]) v else NULL
            }

            ## view_key and n_facets are part of the CACHE KEY, not decoration: the
            ## margins depend on the plot, not only on the image size. Measured at
            ## 960x600 with this same function -- TIC-like panel x spans
            ## 0.0585..0.9240, all-ions 0.0460..0.9136, all-ions/normalised
            ## 0.0429..0.9136 (the m/z colourbar, and the y-tick label width, both
            ## move with the view). Keying on size alone meant whichever view was
            ## brushed first at a given window size poisoned every other one.
            panel_fraction <- function(p, px_w, px_h, res = 72, view_key = NULL, n_facets = NA) {
                key <- paste0(px_w, "x", px_h, "|", if (is.null(view_key)) "" else view_key, "|", n_facets)
                if (!is.null(.panel_cache[[key]])) return(.panel_cache[[key]])
                tf <- tempfile(fileext = ".png")
                grDevices::png(tf, width = px_w, height = px_h, res = res)
                out <- tryCatch({
                    ## One build, used for both the gtable and the panel ranges.
                    built <- ggplot2::ggplot_build(p)
                    gt <- ggplot2::ggplot_gtable(built)
                    grid::grid.draw(gt)
                    ws <- gt$widths
                    n <- length(ws)
                    isnull <- grid::unitType(ws) == "null"
                    abs_in <- numeric(n)
                    if (any(!isnull)) abs_in[!isnull] <- grid::convertWidth(ws[!isnull], "in", valueOnly = TRUE)
                    total_in <- px_w / res
                    nv <- numeric(n)
                    if (any(isnull)) nv[isnull] <- as.numeric(ws[isnull])
                    rem <- max(total_in - sum(abs_in), 0)
                    win <- abs_in + if (sum(nv) > 0) rem * nv / sum(nv) else 0
                    cum <- cumsum(win)
                    ip <- grepl("^panel", gt$layout$name)
                    lcol <- min(gt$layout$l[ip]); rcol <- max(gt$layout$r[ip])
                    xf <- c(if (lcol > 1) cum[lcol - 1] else 0, cum[rcol]) / total_in

                    ## Same again down the heights. Reported as fractions from the
                    ## TOP of the image, because that is the direction image
                    ## pixels run -- and the y brush has to be flipped through it.
                    hs <- gt$heights
                    nh <- length(hs)
                    hnull <- grid::unitType(hs) == "null"
                    habs <- numeric(nh)
                    if (any(!hnull)) habs[!hnull] <- grid::convertHeight(hs[!hnull], "in", valueOnly = TRUE)
                    total_h_in <- px_h / res
                    hv <- numeric(nh)
                    if (any(hnull)) hv[hnull] <- as.numeric(hs[hnull])
                    hrem <- max(total_h_in - sum(habs), 0)
                    hin <- habs + if (sum(hv) > 0) hrem * hv / sum(hv) else 0
                    hcum <- cumsum(hin)
                    trow <- min(gt$layout$t[ip]); brow <- max(gt$layout$b[ip])
                    yf <- c(if (trow > 1) hcum[trow - 1] else 0, hcum[brow]) / total_h_in

                    ## PER-PANEL extents. yf above is the bounding box of ALL
                    ## panels, which is right for x (one panel column, shared
                    ## axis) and wrong for y: with facet_grid(path ~ .) and
                    ## samples_per_page = 20 it spans all twenty stacked panels
                    ## while the brush was drawn inside exactly one of them.
                    ## Measured on a 20-facet page at 960x2200: the stack spans
                    ## 0.0025..0.9849 and each panel is 0.0468 of the image tall,
                    ## so the same in-facet gesture in panel 1 and in panel 20
                    ## mapped to completely different zooms.
                    lay <- gt$layout[ip, , drop = FALSE]
                    lay <- lay[order(lay$t, lay$l), , drop = FALSE]
                    panels <- data.frame(
                        y_lo = vapply(lay$t, function(k) (if (k > 1) hcum[k - 1] else 0) / total_h_in, numeric(1)),
                        y_hi = vapply(lay$b, function(k) hcum[k] / total_h_in, numeric(1)),
                        x_lo = vapply(lay$l, function(k) (if (k > 1) cum[k - 1] else 0) / total_in, numeric(1)),
                        x_hi = vapply(lay$r, function(k) cum[k] / total_in, numeric(1)),
                        stringsAsFactors = FALSE
                    )

                    ## Which facet each panel holds. The gtable's panels, ordered
                    ## top-to-bottom then left-to-right, line up with the build
                    ## layout ordered by ROW then COL -- that is how the gtable is
                    ## assembled. This is what lets a mapping-less brush recover
                    ## the sample it was drawn on; the fallback's own
                    ## brush$mapping$panelvar1 is always NULL, by definition.
                    bl <- built$layout$layout
                    meta <- c("PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y", "COORD", "AXIS_X", "AXIS_Y")
                    facet_cols <- setdiff(names(bl), meta)
                    bl <- bl[order(bl$ROW, bl$COL), , drop = FALSE]
                    if (length(facet_cols) >= 1 && nrow(bl) == nrow(panels)) {
                        panels$facet <- as.character(bl[[facet_cols[1]]])
                        facet_var <- facet_cols[1]
                    } else {
                        panels$facet <- NA_character_
                        facet_var <- NA_character_
                    }

                    ## Each panel's own drawn data range, in the same display order.
                    ## Every v6 view sets explicit y limits (scale limits, or
                    ## coord_cartesian for the ion map), so scales = "free_y" is
                    ## overridden and these all agree -- but reading them per panel
                    ## means the conversion stays correct if that ever changes.
                    pps <- built$layout$panel_params
                    if (nrow(bl) == nrow(panels) && !is.null(bl$PANEL)) {
                        pidx <- as.integer(bl$PANEL)
                        panels$y_data_lo <- vapply(pidx, function(k) { v <- pp_range(pps[[k]], "y"); if (is.null(v)) NA_real_ else v[1] }, numeric(1))
                        panels$y_data_hi <- vapply(pidx, function(k) { v <- pp_range(pps[[k]], "y"); if (is.null(v)) NA_real_ else v[2] }, numeric(1))
                    } else {
                        panels$y_data_lo <- NA_real_
                        panels$y_data_hi <- NA_real_
                    }

                    ## The panel does NOT span the scale limits: ggplot expands a
                    ## continuous axis by 5% each side, so limits c(200, 3200) draw
                    ## a panel running 50..3350. Read the truth off the built plot
                    ## rather than re-deriving it from the window, or every fraction
                    ## is out by 5% of the window, worst at the panel edges.
                    pp <- built$layout$panel_params[[1]]
                    list(x = xf, y = yf, panels = panels, facet_var = facet_var,
                         x_range = pp_range(pp, "x"), y_range = pp_range(pp, "y"))
                }, error = function(e) NULL)
                grDevices::dev.off(); unlink(tf)
                ok <- function(v) length(v) == 2 && all(is.finite(v)) && v[2] > v[1]
                if (!is.null(out) && ok(out$x) && ok(out$y)) {
                    assign(key, out, envir = .panel_cache)
                    return(out)
                }
                NULL
            }

            ## Which panel (1 = top) a from-the-bottom image fraction falls in.
            ## Returns NULL when the geometry is unknown; a fraction in a gutter or
            ## just off the stack snaps to the nearest panel, because a drag that
            ## overshoots the panel edge is the normal gesture, not an error.
            panel_at_bottom_fraction <- function(pf, frac_from_bottom) {
                if (is.null(pf) || is.null(pf$panels) || nrow(pf$panels) == 0) return(NULL)
                if (!is.finite(frac_from_bottom)) return(NULL)
                lo <- 1 - pf$panels$y_hi          # panels carry fractions from the TOP
                hi <- 1 - pf$panels$y_lo
                hit <- which(frac_from_bottom >= lo & frac_from_bottom <= hi)
                if (length(hit) >= 1) {
                    return(hit[which.min(abs(frac_from_bottom - (lo[hit] + hi[hit]) / 2))])
                }
                near <- which.min(pmin(abs(frac_from_bottom - lo), abs(frac_from_bottom - hi)))
                if (length(near) == 1) near else NULL
            }

            ## Honour a stored y-zoom only in the view it was taken in; otherwise
            ## fall back to the range the data itself occupies.
            y_limits_for <- function(key, fallback) {
                if (!is.null(y_zoom) && length(y_zoom) == 2 &&
                    !is.null(y_zoom_key) && identical(y_zoom_key, key)) {
                    return(y_zoom)
                }
                fallback
            }

            ## The brushed y range, in the data units of whatever was last drawn.
            ##
            ## When the coordmap survives, brush$ymin/ymax are already data values.
            ## When it does not they are fractions of the image measured from the
            ## TOP (the same fallback that makes x an image fraction), so they need
            ## flipping through the panel as well as rescaling: the SMALLER
            ## fraction is nearer the top of the picture and therefore the LARGER
            ## data value.
            ##
            ## Returns NULL when there is nothing dependable to convert, in which
            ## case the y axis is simply left alone.
            brush_y_range <- function(brush) {
                if (is.null(brush) || is.null(brush$ymin) || is.null(brush$ymax)) return(NULL)
                if (!is.finite(brush$ymin) || !is.finite(brush$ymax)) return(NULL)

                yv <- if (!is.null(brush$mapping)) brush$mapping$y else NULL
                if (!is.null(yv)) return(sort(c(brush$ymin, brush$ymax)))   # already data units

                why <- NULL
                if (length(last_brush_yrange) != 2 ||
                    !all(is.finite(last_brush_yrange)) || last_brush_yrange[2] <= last_brush_yrange[1]) {
                    why <- "the y range of the last render is not known"
                } else if (is.null(last_brush_plot)) {
                    why <- "no plot recorded from the last render"
                } else if (length(last_brush_size) != 2 ||
                           !all(is.finite(last_brush_size)) || !all(last_brush_size > 0)) {
                    why <- "the rendered plot size is not known"
                } else if (!is.finite(brush$ymin) || !is.finite(brush$ymax) ||
                           brush$ymin < -0.5 || brush$ymax > 1.5) {
                    why <- "the y coordinates are not usable fractions"
                }
                if (!is.null(why)) {
                    cat(paste0("  y-zoom skipped (", why, "); raw y ",
                               signif(brush$ymin, 4), " to ", signif(brush$ymax, 4), ".\n"))
                    return(NULL)
                }

                pf <- panel_fraction(last_brush_plot, last_brush_size[1], last_brush_size[2],
                                     view_key = last_brush_key, n_facets = last_brush_facets)
                if (is.null(pf)) {
                    cat("  y-zoom skipped (panel bounds unknown).\n")
                    return(NULL)
                }

                ## Dragging a box around a peak naturally overshoots the panel top
                ## and bottom, and those coordinates come back outside 0..1.
                ## Clamping is the right reading -- the user meant "everything
                ## vertically" -- where rejecting outright silently dropped every
                ## other y-zoom in the field (2026-09-14).
                ##
                ## The y fractions are measured from the BOTTOM of the image, not
                ## the top -- the opposite of the x convention and of what the
                ## pixel direction suggests. panel_fraction() reports the panel
                ## from the top, so flip it to match. Established from the field
                ## on 2026-09-14: selecting 600-400 zoomed to 200-0 and selecting
                ## 200-0 zoomed to 600-400, a clean reflection about the axis
                ## midpoint. A selection near the middle looked almost right under
                ## either reading, which is why it first appeared to "mostly" work.
                ## ...and it has to be rescaled inside the ONE panel the gesture
                ## was made in, not across the whole stack. With a single facet the
                ## two are identical, which is the case the field evidence above
                ## validated; with twenty, the stack is ~21x the panel, so an
                ## in-facet drag in panel 1 and the same drag in panel 20 used to
                ## produce completely different zooms.
                centre <- mean(c(brush$ymin, brush$ymax))
                idx <- if (!is.null(pf$panels) && nrow(pf$panels) > 1) panel_at_bottom_fraction(pf, centre) else NULL
                if (!is.null(idx)) {
                    panel_lo <- 1 - pf$panels$y_hi[idx]
                    panel_hi <- 1 - pf$panels$y_lo[idx]
                    where <- paste0("panel ", idx, "/", nrow(pf$panels),
                                    if (!is.na(pf$panels$facet[idx])) paste0(" ", pf$panels$facet[idx]) else "")
                } else {
                    panel_lo <- 1 - pf$y[2]        # panel bottom, from the image bottom
                    panel_hi <- 1 - pf$y[1]        # panel top
                    where <- "whole panel stack"
                }
                fr <- pmin(pmax((c(brush$ymin, brush$ymax) - panel_lo) / (panel_hi - panel_lo), 0), 1)

                ## The panel spans its DRAWN range, which carries ggplot's 5% axis
                ## expansion -- limits c(0, 1000) draw a panel running -50..1050.
                ## Read it off the built plot; fall back to the limits recorded at
                ## render time only when the build did not report one.
                drawn <- if (!is.null(idx) && !is.na(pf$panels$y_data_lo[idx]) && !is.na(pf$panels$y_data_hi[idx])) {
                    c(pf$panels$y_data_lo[idx], pf$panels$y_data_hi[idx])
                } else if (!is.null(pf$y_range)) {
                    pf$y_range
                } else {
                    last_brush_yrange
                }
                lo <- drawn[1]; hi <- drawn[2]
                vals <- sort(lo + fr * (hi - lo))
                cat(paste0("  y-zoom: raw ", signif(brush$ymin, 4), "-", signif(brush$ymax, 4),
                           " in ", where,
                           ", panel ", signif(panel_lo, 4), "-", signif(panel_hi, 4), " (from bottom)",
                           ", against drawn y ", signif(lo, 6), "-", signif(hi, 6),
                           " -> ", signif(vals[1], 6), " to ", signif(vals[2], 6), "\n"))
                if (!all(is.finite(vals)) || vals[2] <= vals[1]) return(NULL)
                vals
            }

            ## A peak add or a spectrum extraction belongs to ONE sample. When the
            ## facet cannot be recovered the selection spans every panel on the page,
            ## and taking $path_to_cdf_csv[1] silently attributes it to the first
            ## sample in table order rather than the brushed one -- a wrong peak
            ## written against a wrong sample, with an area summed over all twenty.
            ## Refuse and say so; the diagnostics above name what was resolved.
            ## (Shift+G is exempt by design: it writes the same window to every
            ## sample on purpose.)
            single_sample <- function(pp, what) {
                if (is.null(pp) || nrow(pp) == 0) return(NULL)
                s <- unique(as.character(pp$path_to_cdf_csv))
                if (length(s) == 1) return(s)
                cat(paste0(what, " needs a brush inside ONE sample's panel, but this selection covers ",
                           length(s), " samples (", paste(utils::head(s, 3), collapse = ", "),
                           if (length(s) > 3) ", ..." else "", "). Nothing written.\n"))
                NULL
            }

            ## The mass-spectrum plots have the chromatogram's coordmap problem too,
            ## and a worse symptom: with an empty mapping a bare brushedPoints()
            ## ERRORS ("not able to automatically infer `xvar`") rather than
            ## misreading, which took the whole Shift+1/2/3 block down with it --
            ## including the re-plot the user was trying to get back to. Same
            ## conversion as brushed_chromatogram, one axis, no facets.
            brushed_ms <- function(df, brush, xcol = "mz") {
                empty <- if (is.null(df)) NULL else df[0, , drop = FALSE]
                if (is.null(brush) || is.null(df) || nrow(df) == 0 || !(xcol %in% names(df))) return(empty)

                xv <- if (!is.null(brush$mapping)) brush$mapping$x else NULL
                yv <- if (!is.null(brush$mapping)) brush$mapping$y else NULL
                if (!is.null(xv) && !is.null(yv) && xv %in% names(df) && yv %in% names(df)) {
                    return(tryCatch(shiny::brushedPoints(df, brush),
                                    error = function(e) { cat(paste0("MS brush could not be read (", conditionMessage(e), ").\n")); empty }))
                }
                if (is.null(brush$xmin) || is.null(brush$xmax) ||
                    !is.finite(brush$xmin) || !is.finite(brush$xmax)) return(empty)

                rng <- range(df[[xcol]], na.rm = TRUE)
                bxmin <- brush$xmin; bxmax <- brush$xmax
                if (bxmin >= 0 && bxmax <= 1 && rng[2] > 1) {
                    fr <- c(bxmin, bxmax)
                    pfm <- NULL
                    if (!is.null(last_ms_plot) && length(last_ms_size) == 2 &&
                        all(is.finite(last_ms_size)) && all(last_ms_size > 0)) {
                        pfm <- panel_fraction(last_ms_plot, last_ms_size[1], last_ms_size[2],
                                              view_key = "massSpectra_1", n_facets = 1)
                    }
                    if (!is.null(pfm)) {
                        fr <- pmin(pmax((fr - pfm$x[1]) / (pfm$x[2] - pfm$x[1]), 0), 1)
                        lo <- if (!is.null(pfm$x_range)) pfm$x_range[1] else rng[1]
                        hi <- if (!is.null(pfm$x_range)) pfm$x_range[2] else rng[2]
                    } else {
                        cat("  (mass-spectrum panel bounds unknown - using the uncorrected image fraction, expect an offset)\n")
                        lo <- rng[1]; hi <- rng[2]
                    }
                    bxmin <- lo + fr[1] * (hi - lo)
                    bxmax <- lo + fr[2] * (hi - lo)
                    cat(paste0("MS brush read as ", signif(bxmin, 6), " to ", signif(bxmax, 6), " ", xcol, ".\n"))
                }
                if (!is.finite(bxmin) || !is.finite(bxmax) || bxmax <= bxmin) return(empty)
                df[!is.na(df[[xcol]]) & df[[xcol]] >= bxmin & df[[xcol]] <= bxmax, , drop = FALSE]
            }

            brushed_chromatogram <- function(df, brush) {
                empty <- df[0, , drop = FALSE]
                if (is.null(brush)) return(empty)

                xv <- if (!is.null(brush$mapping)) brush$mapping$x else NULL
                yv <- if (!is.null(brush$mapping)) brush$mapping$y else NULL

                ## Best case: the mapping names real columns, so brushedPoints()
                ## can do its own job (log scales, discrete axes, panels).
                if (!is.null(xv) && !is.null(yv) && xv %in% names(df) && yv %in% names(df)) {
                    return(shiny::brushedPoints(df, brush))
                }

                ## Otherwise fall back to x alone. Two reasons to land here:
                ## the y mapping names a column that does not exist in
                ## chromatograms_updated (the v6 views map y to the transformed
                ## `y`, or to `mz`), or the mapping is missing altogether.
                if (is.null(brush$mapping) || is.null(xv)) {
                    cat(paste0("Brush has no plot mapping. Fields present: ",
                               paste(names(brush), collapse = ", "), "\n"))
                }
                xv_mapped <- if (!is.null(xv) && xv %in% names(df)) xv else NULL
                if (is.null(xv) || !(xv %in% names(df))) xv <- "rt_rt_offset"
                if (!(xv %in% names(df))) return(empty)

                ## xmin/xmax are already in DATA space when Shiny had a coordinate
                ## map to convert with, so they stay usable even when the mapping
                ## itself did not survive. Sanity-check them against the data's own
                ## x range before trusting them: if the plot had no coordinate map
                ## they are pixels, and pixels would silently select the wrong
                ## slice of the run rather than failing loudly.
                if (is.null(brush$xmin) || is.null(brush$xmax) ||
                    !is.finite(brush$xmin) || !is.finite(brush$xmax)) {
                    cat("Brush carries no usable x range - ignoring it.\n")
                    return(empty)
                }
                rng <- range(df[[xv]], na.rm = TRUE)
                bxmin <- brush$xmin; bxmax <- brush$xmax
                pfull <- NULL        # panel geometry, filled in by the fallback below

                ## Field evidence (host1, 2026-09-14): when Shiny's ggplot coordmap
                ## extraction fails it still sends a brush, but with an EMPTY
                ## mapping and xmin/xmax expressed as FRACTIONS OF THE PANEL
                ## (0..1) against a degenerate 0..1 domain. Two successive drags
                ## came back 0.259-0.411 and 0.347-0.498: identical widths, the
                ## second shifted right, exactly as two equal drags at different
                ## positions would look.
                ##
                ## Those are recoverable. The panel spans the current x-axis
                ## window, so a fraction maps straight onto it. Guarded on the
                ## data's own x range sitting clear of 0..1 -- true for retention
                ## times in seconds, and the check means real data coordinates are
                ## never reinterpreted by mistake.
                if (is.null(xv_mapped) && bxmin >= 0 && bxmax <= 1 && rng[1] > 1) {
                    win_lo <- if (!is.null(x_axis_start) && length(x_axis_start) == 1 && is.finite(x_axis_start)) x_axis_start else rng[1]
                    win_hi <- if (!is.null(x_axis_end)   && length(x_axis_end)   == 1 && is.finite(x_axis_end))   x_axis_end   else rng[2]
                    if (win_hi > win_lo) {
                        fr <- c(brush$xmin, brush$xmax)

                        ## Re-express image fractions as PANEL fractions. Skipped
                        ## (with a warning) if the last plot or its rendered size
                        ## is unknown, in which case the old uncorrected reading
                        ## is used -- offset, but better than nothing.
                        pf <- NULL
                        if (!is.null(last_brush_plot) && length(last_brush_size) == 2 &&
                            all(is.finite(last_brush_size)) && all(last_brush_size > 0)) {
                            pfull <- panel_fraction(last_brush_plot, last_brush_size[1], last_brush_size[2],
                                                    view_key = last_brush_key, n_facets = last_brush_facets)
                            if (!is.null(pfull)) pf <- pfull$x
                        }
                        if (!is.null(pf)) {
                            fr <- pmin(pmax((fr - pf[1]) / (pf[2] - pf[1]), 0), 1)
                        } else {
                            cat("  (panel bounds unknown - using the uncorrected image fraction, expect an offset)\n")
                        }

                        ## The panel spans its DRAWN x range, not the axis window:
                        ## scale_x_continuous(limits = c(200, 3200)) draws a panel
                        ## running 50..3350, because ggplot expands a continuous
                        ## axis by 5% each side. Fraction 0 is rt 50, and reading it
                        ## as 200 is a 5%-of-window error, worst at the panel edges
                        ## and zero in the middle -- the same signature as the
                        ## y-flip bug, which is why it could hide behind it.
                        if (!is.null(pfull) && !is.null(pfull$x_range)) {
                            win_lo <- pfull$x_range[1]; win_hi <- pfull$x_range[2]
                        }

                        bxmin <- win_lo + fr[1] * (win_hi - win_lo)
                        bxmax <- win_lo + fr[2] * (win_hi - win_lo)
                        cat(paste0("Brush came back as image fractions (", signif(brush$xmin, 4), "-",
                                   signif(brush$xmax, 4), ")",
                                   if (!is.null(pf)) paste0(", panel spans ", signif(pf[1], 4), "-", signif(pf[2], 4)) else "",
                                   ", drawn x ", signif(win_lo, 6), "-", signif(win_hi, 6),
                                   "; reading as ", signif(bxmin, 8), " to ", signif(bxmax, 8), ".\n"))
                    }
                }

                if (bxmax < rng[1] || bxmin > rng[2]) {
                    cat(paste0("Brush x range (", signif(bxmin, 6), " to ", signif(bxmax, 6),
                               ") lies outside the data (", signif(rng[1], 6), " to ", signif(rng[2], 6),
                               "). Ignoring it.\n"))
                    return(empty)
                }

                keep <- !is.na(df[[xv]]) & df[[xv]] >= bxmin & df[[xv]] <= bxmax

                ## WHICH FACET was brushed, in descending order of trust:
                ##   1. brush$mapping$panelvar1 + brush$panelvar1 - the coordmap survived;
                ##   2. brush$panelvar1 on its own - Shiny sometimes sends the panel
                ##      value at top level even when the mapping came back empty;
                ##   3. the panel the brush's y centre falls in, from the per-panel
                ##      extents measured off the built plot.
                ## Route 1 was the only one there before, and it read the variable
                ## NAME out of brush$mapping -- whose emptiness is precisely what puts
                ## us in this fallback. So it was dead code: the selection ran across
                ## every facet on the page, and Shift+A then took $path_to_cdf_csv[1],
                ## the first sample in table order rather than the brushed one.
                pv <- if (!is.null(brush$mapping)) brush$mapping$panelvar1 else NULL
                if (is.null(pv) && !is.null(pfull) && !is.na(pfull$facet_var)) pv <- pfull$facet_var
                panel_value <- if (!is.null(brush$panelvar1)) as.character(brush$panelvar1) else NULL
                panel_source <- if (!is.null(panel_value)) "brush$panelvar1" else NULL
                if (is.null(panel_value) && !is.null(pfull) && !is.null(pfull$panels) &&
                    nrow(pfull$panels) > 1 && !is.null(brush$ymin) && !is.null(brush$ymax) &&
                    is.finite(brush$ymin) && is.finite(brush$ymax) &&
                    brush$ymin > -0.5 && brush$ymax < 1.5) {
                    pidx <- panel_at_bottom_fraction(pfull, mean(c(brush$ymin, brush$ymax)))
                    if (!is.null(pidx) && !is.na(pfull$panels$facet[pidx])) {
                        panel_value <- pfull$panels$facet[pidx]
                        panel_source <- paste0("panel ", pidx, "/", nrow(pfull$panels), " by y position")
                    }
                }
                if (!is.null(pv) && pv %in% names(df) && !is.null(panel_value)) {
                    keep <- keep & as.character(df[[pv]]) == panel_value
                    cat(paste0("Brush resolved to facet ", panel_value, " (", panel_source, ").\n"))
                } else if (!is.null(pfull) && !is.null(pfull$panels) && nrow(pfull$panels) > 1) {
                    cat("Brush facet could not be resolved - the selection spans every panel on the page.\n")
                }
                out <- df[which(keep), , drop = FALSE]
                if (nrow(out) == 0) {
                    cat(paste0("Brush (", xv, " ", signif(bxmin, 6), " to ", signif(bxmax, 6),
                               if (!is.null(brush$panelvar1)) paste0(", panel ", brush$panelvar1) else "",
                               ") matched no rows.\n"))
                } else {
                    cat(paste0("Brush selected ", nrow(out), " points, ", xv, " ",
                               signif(min(out[[xv]]), 6), " to ", signif(max(out[[xv]]), 6), ".\n"))
                }
                out
            }

            ## Area over the ions a peak actually generates. v5 integrates the
            ## TIC, which charges the peak for whatever co-elutes beneath it.
            ## Restricting the sum to ions that genuinely belong to the peak
            ## removes that contribution. No decimation here -- this is a number,
            ## not a picture, so every scan counts.
            ##
            ## Rewritten 2026-09-30. The first version reused allion_window(),
            ## which was built for DRAWING, and it was wrong here in two ways:
            ##
            ##   * Its per-ion baseline is the 10th percentile INSIDE the
            ##     window. For a correctly-bounded peak the peak fills its own
            ##     window, so that percentile sits well up the flank and the
            ##     peak eats its own baseline. Measured on a clean Gaussian:
            ##     bounds at +/-3 sigma recovered 94% of the true area, +/-2
            ##     sigma 64%, +/-1.5 sigma 39%. The number therefore swung ~2x
            ##     with nothing but how tightly the bounds were drawn.
            ##     => the baseline is now taken from the FLANKS, OUTSIDE
            ##        [start, end], and interpolated across the peak, which is
            ##        what the TIC path already does.
            ##
            ##   * It selected ions by in-window AMPLITUDE. A co-eluting second
            ##     compound's ions are tall inside that window, so they cleared
            ##     the cut and were summed straight back in -- i.e. the column
            ##     did not do the one thing it exists to do.
            ##     => selection is now a SHAPE test: an ion is kept only if it
            ##        apexes with the peak and tracks its profile.
            ##
            ## Returns NA (not 0) when there is nothing dependable to integrate,
            ## so a failure reads as "no number" rather than "no signal".
            peak_ion_area <- function(cdf_csv, start_rt, end_rt, threshold,
                                      flank_frac = 1.0, apex_tol_sigmas = 0.5,
                                      min_cor = 0.8) {

                na_out <- list(area = NA_real_, ions = NA_character_)
                if (!is.finite(start_rt) || !is.finite(end_rt) || end_rt <= start_rt) return(na_out)

                d <- try(load_allions(cdf_csv), silent = TRUE)
                if (inherits(d, "try-error") || is.null(d) || nrow(d) == 0) return(na_out)

                ## Single-detector data (GC-TCD/FID/ECD) carries the mz = 0
                ## sentinel and has no spectral dimension, so an ion-restricted
                ## area is meaningless. Short-circuit, as Shift+1 and Shift+4 do.
                if (all(d$mz == 0)) return(na_out)

                width  <- end_rt - start_rt
                pad_lo <- start_rt - flank_frac * width
                pad_hi <- end_rt   + flank_frac * width
                padded <- d[rt >= pad_lo & rt <= pad_hi]
                if (nrow(padded) == 0) return(na_out)

                inside <- padded[rt >= start_rt & rt <= end_rt]
                if (nrow(inside) == 0 || data.table::uniqueN(inside$rt) < 3) return(na_out)

                ## Per-ion baseline from the flanks, interpolated across the peak.
                ## Median rather than mean so a neighbouring peak intruding into
                ## a flank does not drag the level up. Falls back to whichever
                ## flank exists, then to the in-window minimum when the peak sits
                ## against the start or end of the run.
                left  <- padded[rt <  start_rt]
                right <- padded[rt >  end_rt]
                lsum <- if (nrow(left))  left[,  .(bl = stats::median(abundance), rtl = stats::median(rt)), by = mz] else NULL
                rsum <- if (nrow(right)) right[, .(br = stats::median(abundance), rtr = stats::median(rt)), by = mz] else NULL

                w <- inside
                if (!is.null(lsum)) w <- merge(w, lsum, by = "mz", all.x = TRUE)
                if (!is.null(rsum)) w <- merge(w, rsum, by = "mz", all.x = TRUE)
                if (is.null(lsum)) { w[, bl := NA_real_]; w[, rtl := NA_real_] }
                if (is.null(rsum)) { w[, br := NA_real_]; w[, rtr := NA_real_] }

                w[, floor_mz := min(abundance), by = mz]
                w[, base := data.table::fifelse(
                        !is.na(bl) & !is.na(br) & is.finite(rtr - rtl) & (rtr - rtl) != 0,
                        bl + (br - bl) * (rt - rtl) / (rtr - rtl),
                    data.table::fifelse(!is.na(bl), bl,
                    data.table::fifelse(!is.na(br), br, floor_mz)))]
                w[, net := pmax(abundance - base, 0)]

                ## Reference = the BASE ION, not the summed profile. The sum is
                ## a blend when something co-elutes, and its apex sits between
                ## the two compounds, which is exactly the case we are trying to
                ## resolve. The single strongest ion belongs unambiguously to one
                ## compound, so it is the honest model peak.
                tops <- w[, .(top = max(net), ion_apex = rt[which.max(net)]), by = mz]
                if (nrow(tops) == 0 || max(tops$top) <= 0) return(na_out)
                base_mz <- tops$mz[which.max(tops$top)]

                ref <- w[mz == base_mz, .(rt, ref = net)]
                data.table::setorder(ref, rt)
                if (nrow(ref) < 3 || max(ref$ref) <= 0) return(na_out)
                apex_rt <- ref$rt[which.max(ref$ref)]

                ## Peak width from the base ion's full width at half maximum.
                ## The apex tolerance has to be in units of the PEAK, not of the
                ## integration window: a window twice as wide does not make two
                ## compounds twice as resolvable. Falls back to a quarter of the
                ## window when the peak is too sparsely sampled to measure.
                half  <- ref$rt[ref$ref >= max(ref$ref) / 2]
                sigma_est <- if (length(half) >= 2) (max(half) - min(half)) / 2.355 else width / 4
                if (!is.finite(sigma_est) || sigma_est <= 0) sigma_est <- width / 4

                w <- merge(w, ref, by = "rt")
                prof <- w[, .(
                        top      = max(net),
                        ion_apex = rt[which.max(net)],
                        r        = suppressWarnings(stats::cor(net, ref))
                    ), by = mz]

                keep <- prof[
                    top > threshold * max(prof$top) &                              # amplitude floor
                    abs(ion_apex - apex_rt) <= apex_tol_sigmas * sigma_est &       # co-apexes with the base ion
                    !is.na(r) & r >= min_cor,                                      # and tracks its profile
                    mz]

                ## A real peak's base ion defines `ref`, so it passes by
                ## construction; if nothing does, the window is not a peak.
                ## Keep the strongest ion rather than silently reporting 0.
                if (length(keep) == 0) keep <- prof$mz[which.max(prof$top)]

                kept <- w[mz %in% keep]
                list(area = sum(kept$net), ions = paste(sort(unique(kept$mz)), collapse = ";"))
            }

        setwd(CDF_directory_path)

        ## Chemstation exports detector channels as SIGNAL01.CDF, SIGNAL02.CDF, etc.
        ## Rename to <sample_name>.CDF using the file's sample_name global attribute
        ## so downstream labels are meaningful instead of channel numbers.
        signal_files <- dir()[grep("^SIGNAL.*\\.cdf$", dir(), ignore.case = TRUE)]
        for (sf in signal_files) {
            sample_name <- tryCatch({
                nc <- ncdf4::nc_open(sf)
                att <- ncdf4::ncatt_get(nc, 0, "sample_name")
                ncdf4::nc_close(nc)
                if (isTRUE(att$hasatt) && is.character(att$value)) att$value else ""
            }, error = function(e) "")
            if (!nzchar(sample_name)) next
            safe_name <- gsub("[^A-Za-z0-9._-]+", "_", sample_name)
            ext <- sub("^.*(\\.[Cc][Dd][Ff])$", "\\1", sf)
            target <- paste0(safe_name, ext)
            if (target == sf) next
            if (file.exists(target)) {
                i <- 2L
                while (file.exists(paste0(safe_name, "_", i, ext))) i <- i + 1L
                target <- paste0(safe_name, "_", i, ext)
            }
            file.rename(sf, target)
            if (file.exists(paste0(sf, ".csv"))) {
                file.rename(paste0(sf, ".csv"), paste0(target, ".csv"))
            }
            cat(paste0("Renamed ", sf, " -> ", target, "\n"))
        }

        paths_to_cdfs <- dir()[grep("\\.CDF$", dir(), ignore.case = TRUE)]
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

                for (file in seq_along(paths_to_cdfs)) {

                    ## If the cdf.csv doesn't exist for this cdf -- or is stale
                    ## against it -- create it. Stale means the .CDF has a different
                    ## mtime/size from the one the .csv was built from, i.e. the
                    ## acquisition was replaced under the same filename. Before the
                    ## stamp, that case silently kept serving the previous run.

                        cdf_csv_stale <- !stamp_is_current(paths_to_cdf_csvs[file], paths_to_cdfs[file])

                        if ( cdf_csv_stale ) {

                            if ( file.exists(paths_to_cdf_csvs[file]) ) {
                                cat(paste0(paths_to_cdfs[file], " has changed since its CSV was built - rebuilding it and everything derived from it.\n"))
                            }

                            cdf_type_this <- detect_cdf_type(paths_to_cdfs[file])

                            if (cdf_type_this == "chrom") {

                                ## ANDI-Chrom single-detector file (e.g. GC-TCD).
                                ## One row per timepoint with mz=0 sentinel.
                                cat(paste("CDF (single-detector / chromatography format) to CSV conversion. Reading data file ", paths_to_cdfs[file], "\n", sep = ""))
                                    framedDataFile <- read_chrom_cdf_as_framed(paths_to_cdfs[file])

                                cat("   Writing out data file as CSV... \n")
                                    data.table::fwrite(framedDataFile, file = paste(paths_to_cdfs[file], ".csv", sep = ""), col.names = TRUE, row.names = FALSE)

                            } else {

                                cat(paste("CDF to CSV conversion. Reading data file ", paths_to_cdfs[file], "\n", sep = ""))
                                    rawDataFile <- xcms::loadRaw(xcms::xcmsSource(paths_to_cdfs[file]))

                                cat("   Framing data file ... \n")
                                    rt <- rawDataFile$rt
                                    scanindex <- rawDataFile$scanindex

                                    filteredRawDataFile <- list()
                                    for ( i in seq_len(max(0L, length(rt) - 1L)) ) {
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

                            ## Stamp the new CSV, then drop everything derived from
                            ## the old one: the all-ion store, its in-memory copy,
                            ## and this sample's rows in chromatograms.csv (which is
                            ## the table actually plotted, so leaving them there is
                            ## exactly how a replaced acquisition stayed invisible).
                            write_stamp(paths_to_cdf_csvs[file], paths_to_cdfs[file])
                            unlink(c(allion_store_path(paths_to_cdf_csvs[file]),
                                     stamp_path(allion_store_path(paths_to_cdf_csvs[file]))))
                            if (!is.null(.allion_cache[[paths_to_cdf_csvs[file]]])) {
                                rm(list = paths_to_cdf_csvs[file], envir = .allion_cache)
                            }
                            if (is.data.frame(chromatograms) && nrow(chromatograms) > 0) {
                                chromatograms <- chromatograms[chromatograms$path_to_cdf_csv != paths_to_cdf_csvs[file], , drop = FALSE]
                            }
                            if (file.exists("chromatograms.csv")) {
                                stale_rows <- readMonolist("chromatograms.csv")
                                writeMonolist(
                                    monolist = stale_rows[stale_rows$path_to_cdf_csv != paths_to_cdf_csvs[file], , drop = FALSE],
                                    monolist_out_path = "chromatograms.csv"
                                )
                            }

                        }

                    ## Build the nominal-mass all-ion store for this sample.
                    ## Done here rather than lazily at first render so the cost
                    ## lands during ingest, where there is already a progress
                    ## log, instead of freezing the first redraw.

                        if ( isTRUE(all_ion_view) && !stamp_is_current(allion_store_path(paths_to_cdf_csvs[file]), paths_to_cdf_csvs[file]) ) {
                            cat(paste("   Building all-ion store for ", paths_to_cdf_csvs[file], "\n", sep = ""))
                            build_allion_store(paths_to_cdf_csvs[file])
                        }

                    ## If any chromatograms (tic and ion) are not present for this csv, extract them

                        ## Which ions this sample already has, read from the IN-MEMORY
                        ## tables rather than by re-reading chromatograms.csv from disk
                        ## once per sample. That read (plus the matching write below) is
                        ## what made folder loading O(N^2) in bytes: at 152 samples the
                        ## whole table was read and rewritten 152 times.
                        ## Both sides are compared as character, because the on-disk
                        ## column mixes "0" with "baseline" while `ions` is numeric.
                        ions_for_this_cdf_csv <- character(0)
                        if (is.data.frame(chromatograms) && nrow(chromatograms) > 0) {
                            ions_for_this_cdf_csv <- c(ions_for_this_cdf_csv,
                                as.character(chromatograms$ion[chromatograms$path_to_cdf_csv == paths_to_cdf_csvs[file]]))
                        }
                        if (is.data.frame(chromatograms_to_add) && nrow(chromatograms_to_add) > 0) {
                            ions_for_this_cdf_csv <- c(ions_for_this_cdf_csv,
                                as.character(chromatograms_to_add$ion[chromatograms_to_add$path_to_cdf_csv == paths_to_cdf_csvs[file]]))
                        }
                        ions_for_this_cdf_csv <- unique(ions_for_this_cdf_csv)

                        if ( length(ions_for_this_cdf_csv) > 0 ) {

                            missing_ions <- as.numeric(as.character(ions[!as.character(ions) %in% ions_for_this_cdf_csv]))
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
                                    for ( ion in seq_along(numeric_ions) ){
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

                }   ## end for each file

                ## ONE write, after the loop. This used to sit INSIDE it, so the whole
                ## table was rebuilt and rewritten once per sample -- and
                ## unconditionally, so even a fully-cached relaunch with nothing new to
                ## extract rewrote it N times for no reason. Now a cached relaunch
                ## writes nothing at all unless there are stale rows to prune.
                ##
                ## Prune on the way out: rows for CDFs no longer in the folder were
                ## dropped from the in-memory copy but written back regardless, so a
                ## deleted acquisition stayed in chromatograms.csv for ever.
                    added_any <- is.data.frame(chromatograms_to_add) && nrow(chromatograms_to_add) > 0
                    had_any   <- is.data.frame(chromatograms) && nrow(chromatograms) > 0

                    if ( added_any ) {

                        to_write <- if (had_any) rbind(chromatograms, chromatograms_to_add) else chromatograms_to_add
                        to_write <- to_write[to_write$path_to_cdf_csv %in% paths_to_cdf_csvs, , drop = FALSE]
                        writeMonolist(monolist = to_write, monolist_out_path = "chromatograms.csv")

                    } else if ( had_any ) {

                        pruned <- chromatograms[chromatograms$path_to_cdf_csv %in% paths_to_cdf_csvs, , drop = FALSE]
                        if (nrow(pruned) != nrow(chromatograms)) {
                            cat(paste0("Dropping ", nrow(chromatograms) - nrow(pruned),
                                       " chromatogram rows for CDFs no longer in this folder.\n"))
                            writeMonolist(monolist = pruned, monolist_out_path = "chromatograms.csv")
                        }

                    } else if ( !file.exists("chromatograms.csv") ) {

                        writeMonolist(chromatograms_to_add, "chromatograms.csv")

                    }

                print("done")
            }

            ## Read in chromatograms

                chromatograms <- readMonolist("chromatograms.csv")
                        
            ## Set up new samples monolist

                ## ONE derivation for Sample_ID, used by both the create and the
                ## append path below. They used to disagree: the first write stored
                ## the full `foo.CDF.csv`, while later appends ran
                ## gsub("\\..*$", "", ...) -- truncating at the FIRST dot, so
                ## `WT.rep2.CDF` and `WT.rep3.CDF` both collapsed to `WT`. Strip only
                ## the .CDF.csv suffix (and any directory), as the facet labels do.
                sample_id_from_path <- function(x) {
                    gsub("\\.CDF\\.csv$", "", gsub(".*/", "", x), ignore.case = TRUE)
                }

                ## If it doesn't exist, create it
                
                    if ( !file.exists("samples_monolist.csv") ) {

                        samples_monolist <- data.frame(
                            Sample_ID = sample_id_from_path(unique(chromatograms$path_to_cdf_csv)),
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
                                Sample_ID = sample_id_from_path(missing_from_samples_monolist),
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

                            ## Re-read: the rows above went to the FILE, and the
                            ## in-memory copy was left one acquisition behind.
                            samples_monolist <- readMonolist("samples_monolist.csv")
                        }
                    }

            ## Filter chromatograms so only the CDFs in this folder are included

                chromatograms <- chromatograms[chromatograms$path_to_cdf_csv %in% dir()[grep("\\.CDF\\.csv$", dir(), ignore.case = TRUE)],]

            ## Set up several variables, plot_height, and x_axis limits if not specified in function call
                
                peak_data <- NULL
                peak_points <- NULL

                ## y-zoom state. y_zoom holds the brushed y range and y_zoom_key
                ## names the view it was taken in, so it is only ever re-applied
                ## to the axis it actually refers to.
                y_zoom <- NULL
                y_zoom_key <- NULL
                last_brush_plot <- NULL
                last_brush_size <- c(NA_real_, NA_real_)
                last_brush_key <- NULL
                last_brush_yrange <- c(NA_real_, NA_real_)
                ## How many facets the last render actually drew. Part of the
                ## panel-geometry cache key, because panel heights move with it.
                last_brush_facets <- NA_integer_
                ## Same pair for the mass-spectrum plot, used by brushed_ms().
                last_ms_plot <- NULL
                last_ms_size <- c(NA_real_, NA_real_)

                ## Everything the handlers reach for with `<<-` is declared HERE, so
                ## the assignment binds in this function's frame. Without a local
                ## binding, `<<-` walks all the way out and writes to globalenv() --
                ## and nothing clears globalenv() when the app quits. Open folder A,
                ## quit, open folder B in the same R session and the panel showed
                ## folder A's traces, while Shift+6's guard passed on the stale
                ## object and wrote folder A's peaks into folder B's
                ## peaks_monolist.csv. Every guard on these names must therefore be
                ## is.null(), NOT exists(): once declared they always exist.
                chromatograms_updated <- NULL
                x_axis_start <- NULL
                x_axis_end <- NULL
                y_axis_start <- NULL
                y_axis_end <- NULL
                MS_out_1 <- NULL
                MS_ret_start_line <- NULL
                MS_ret_end_line <- NULL
                framedDataFile_to_subtract <- NULL
                ## What the editable peak table was last rendered FROM. Shift+Z writes
                ## the browser's copy straight back, so if the file moved on since, that
                ## write silently reverts it. See render_peak_table() / Shift+Z.
                peak_table_snapshot <- NULL
                ms_wide <- NULL
                predictions <- NULL
                ## Bounded: never draw more than one page of facets at once, so the
                ## render can't blow up the graphics device on large sample sets.
                plot_height <- 200 + 100*min(samples_per_page, length(unique(chromatograms$path_to_cdf_csv)))

                ## x-axis defaults are resolved HERE, not in the Shift+Q handler, so
                ## they are never NULL by the time a key handler can fire. The pan/zoom
                ## keys (F/D/V/C) do `x_axis_start_default <<- x_axis_start_default + rate`;
                ## NULL + rate is numeric(0), which is NOT NULL, so the Shift+Q
                ## initialiser's is.null() check would skip it and the renderer's
                ## dplyr::filter() would then die on a length-0 predicate.
                if ( is.null(x_axis_start_default) || length(x_axis_start_default) == 0 || !is.finite(x_axis_start_default) ) {
                    x_axis_start_default <- if (nrow(chromatograms) > 0) min(chromatograms$rt, na.rm = TRUE) else 0
                }
                if ( is.null(x_axis_end_default) || length(x_axis_end_default) == 0 || !is.finite(x_axis_end_default) ) {
                    x_axis_end_default <- if (nrow(chromatograms) > 0) max(chromatograms$rt, na.rm = TRUE) else 1
                }

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

                ## Draggable sidebar. sidebarLayout() emits a fixed Bootstrap
                ## col-sm-4 / col-sm-8 pair, so the split cannot be changed from
                ## R. Convert that row to flexbox and put a grab handle between
                ## the two columns. On release we fire a window resize, which is
                ## what makes Shiny re-measure the plot outputs -- that also
                ## refreshes last_brush_size, which the brush maths depends on.
                tags$head(
                    tags$style(HTML("
                        .gcms-flex-row { display: flex !important; align-items: stretch; }
                        .gcms-sidebar  { flex: 0 0 auto !important; float: none !important;
                                         max-width: none !important; }
                        .gcms-main     { flex: 1 1 auto !important; float: none !important;
                                         width: auto !important; min-width: 0; }
                        .gcms-drag     { flex: 0 0 6px; cursor: col-resize; background: #dcdcdc;
                                         border-radius: 3px; margin: 0 6px; }
                        .gcms-drag:hover { background: #9e9e9e; }
                        body.gcms-dragging { cursor: col-resize; user-select: none; }
                    ")),
                    tags$script(HTML("
                        $(function() {
                            var well = $('.well').first();
                            if (!well.length) return;
                            var side = well.parent();
                            var row  = side.parent();
                            var main = side.next();
                            row.addClass('gcms-flex-row');
                            side.addClass('gcms-sidebar').css('width', '340px');
                            main.addClass('gcms-main');
                            var handle = $('<div></div>').addClass('gcms-drag')
                                          .attr('title', 'Drag to resize the sidebar');
                            side.after(handle);
                            var dragging = false;
                            handle.on('mousedown', function(e) {
                                dragging = true;
                                $('body').addClass('gcms-dragging');
                                e.preventDefault();
                            });
                            $(document).on('mousemove', function(e) {
                                if (!dragging) return;
                                var w = e.pageX - side.offset().left;
                                var maxw = $(window).width() - 300;
                                if (w < 160) w = 160;
                                if (w > maxw) w = maxw;
                                side.css('width', w + 'px');
                            });
                            $(document).on('mouseup', function() {
                                if (!dragging) return;
                                dragging = false;
                                $('body').removeClass('gcms-dragging');
                                $(window).trigger('resize');
                            });
                        });
                    "))
                ),

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

                        tags$hr(),
                        strong("Chromatogram display"),
                        selectInput(
                            "display_mode", NULL,
                            choices = c("Total ion current" = "tic",
                                        "All ions (thin lines)" = "all_ions",
                                        "Ion map (heatmap)" = "ion_map"),
                            selected = if (isTRUE(all_ion_view)) "all_ions" else "tic"
                        ),
                        selectInput(
                            "allion_scaling", "Y scaling",
                            choices = c("Raw counts" = "raw",
                                        "Square root" = "sqrt",
                                        "Log10" = "log",
                                        "Each ion to its own max" = "norm"),
                            selected = all_ion_scaling
                        ),
                        sliderInput(
                            "allion_threshold", "Ion inclusion (fraction of the window's largest ion)",
                            min = 0.001, max = 0.2, value = all_ion_threshold, step = 0.001
                        ),
                        sliderInput(
                            "allion_max_points", "Max x-points per ion (render cost)",
                            min = 200, max = 4000, value = all_ion_max_points, step = 100
                        ),
                        helpText("A pure peak draws as a tight bundle; a splayed bundle means two compounds. Display and Y scaling apply immediately; the two sliders apply on the next Shift+Q."),

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

                                    ## Pagination controls: draw only samples_per_page
                                    ## samples at a time and step through pages.
                                    fluidRow(
                                        column(12,
                                            actionButton("prev_page", paste0("\u25c0 Prev ", samples_per_page)),
                                            actionButton("next_page", paste0("Next ", samples_per_page, " \u25b6")),
                                            tags$span(
                                                style = "margin-left:15px; font-weight:bold;",
                                                textOutput("page_indicator", inline = TRUE)
                                            )
                                        )
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

                ## Pagination state. current_page = which page of samples is shown;
                ## redraw_trigger is bumped by Shift+Q so the (reactive) chromatogram
                ## plot recomputes only when the data actually changes, not on every brush.
                    current_page   <- reactiveVal(1)
                    redraw_trigger <- reactiveVal(0)

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

                ## ONE way to render the editable peak table. Shift+G used to render it
                ## as DT::renderDataTable while every other handler used rhandsontable --
                ## and `input$peak_table` from a DT output carries nothing hot_to_r() can
                ## read, so Shift+Z silently did nothing for the rest of the session.
                ## Each render also records what it rendered FROM, so Shift+Z can tell
                ## an edit from a revert.
                    render_peak_table <- function() {
                        current <- if (file.exists("peaks_monolist.csv")) read.csv("peaks_monolist.csv") else NULL
                        peak_table_snapshot <<- current
                        output$peak_table <- rhandsontable::renderRHandsontable(rhandsontable::rhandsontable({
                            current
                        }))
                        invisible(NULL)
                    }

                ## Save manual changes to table on "Z" (90) keystroke

                    observeEvent(input$keypress, {
                        if (input$keypress == 90 ) {
                            ## Write out any modifications to peak table (i.e. sample IDs)
                                hot = isolate(input$peak_table)
                                if (is.null(hot)) {
                                    cat("Shift+Z: no editable peak table in the browser yet - press Shift+Q first.\n")
                                    return()
                                }
                                edited <- rhandsontable::hot_to_r(hot)

                                ## The browser holds a SNAPSHOT from the last render. If the file
                                ## has moved on since (a Shift+A, a Shift+6, another Shift+Q), then
                                ## writing that snapshot back reverts everything added in between.
                                ## Refuse rather than guess: a cell-level merge would need a stable
                                ## row key, and this table has none that the user cannot also edit.
                                on_disk <- if (file.exists("peaks_monolist.csv")) read.csv("peaks_monolist.csv") else NULL
                                moved_on <- !is.null(peak_table_snapshot) && !is.null(on_disk) &&
                                            !identical(dim(peak_table_snapshot), dim(on_disk))
                                if (moved_on) {
                                    cat(paste0("Shift+Z: peaks_monolist.csv has changed since the table was drawn (",
                                               nrow(peak_table_snapshot), " rows then, ", nrow(on_disk),
                                               " now). Saving would discard the difference - press Shift+Q to refresh the table, then re-make your edits.\n"))
                                    return()
                                }

                                writeMonolist(edited, "peaks_monolist.csv")
                                peak_table_snapshot <<- edited
                                cat("Peak list saved!\n")
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
                        if (is.null(chromatograms_updated)) {
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
                              # bind_rows, not rbind: once a Shift+Q has run, the on-disk table carries extra
                              # columns (peak_number_within_sample, rt_offset, area_peak_ions, ...) that
                              # peak_data_new does not, and rbind errors on the mismatch. The new rows
                              # get NA in those columns; the next Shift+Q recomputes them.
                              peak_data_combined <- dplyr::bind_rows(peak_data_existing, peak_data_new)
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
                            if (is.null(chromatograms_updated)) {
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

                                for ( chrom in seq_along(unique(chromatograms_updated$path_to_cdf_csv)) ) {
                          
                                    chromatogram <- dplyr::filter(chromatograms_updated, path_to_cdf_csv == unique(chromatograms_updated$path_to_cdf_csv)[chrom])
                                    tic <- filter(chromatogram, ion == 0)

                                    prelim_baseline_window <- samples_monolist$baseline_window[match(chromatogram$path_to_cdf_csv[1], samples_monolist$path_to_cdf_csv)]

                                    ## A sample with fewer scans than baseline_window (a short TCD run,
                                    ## a narrow SIM method) gives floor(...) == 0; `1:0` then negative-indexes
                                    ## to numeric(0), min() returns Inf and data.frame() throws
                                    ## "arguments imply differing number of rows: 0, 1" with no hint that
                                    ## baseline_window is the knob. Clamp to at least one window, iterate with
                                    ## seq_len(), and clamp the slice so the last window cannot over-run.
                                    n_scans_in_tic <- length(tic$rt)
                                    n_prelim_baseline_windows <- max(1L, floor(n_scans_in_tic/prelim_baseline_window))
                                    if ( n_scans_in_tic < prelim_baseline_window ) {
                                        cat(paste0("Sample ", chromatogram$path_to_cdf_csv[1], " has ", n_scans_in_tic,
                                                   " scans, fewer than baseline_window (", prelim_baseline_window,
                                                   "); using one window over the whole trace. Lower baseline_window in samples_monolist.csv for a finer baseline.\n"))
                                    }
                                    prelim_baseline <- list()
                                    for ( i in seq_len(n_prelim_baseline_windows) ) {
                                        window_lo <- (prelim_baseline_window*(i-1))+1
                                        window_hi <- min(prelim_baseline_window*i, n_scans_in_tic)
                                        abundances_in_window <- tic$abundance[window_lo:window_hi]
                                        prelim_baseline[[i]] <- data.frame(
                                            rt = tic$rt[(which.min(abundances_in_window)+((i-1)*prelim_baseline_window))],
                                            min = min(abundances_in_window)
                                        )
                                    }
                                    prelim_baseline <- do.call(rbind, prelim_baseline)

                                    ## bound baseline to first and last point in chromatogram
                                    ## (index by position, not by rt value — original code accidentally
                                    ##  used the rt value as an integer index, which only happened to
                                    ##  land in-range when rt values were large; on TCD data where rt
                                    ##  starts near 0, tic$abundance[0] is length-0 and breaks data.frame())
                                    prelim_baseline <- rbind(
                                        data.frame(rt = min(tic$rt), min = tic$abundance[which.min(tic$rt)]),
                                        prelim_baseline,
                                        data.frame(rt = max(tic$rt), min = tic$abundance[which.max(tic$rt)])
                                    )

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
                                ## rbind onto chromatograms_updated (the SUBSET), not the unfiltered
                                ## `chromatograms` — rebuilding from the latter re-admitted every excluded
                                ## sample with no ion == "baseline" rows and rt_offset = NA, which killed
                                ## the render as soon as one of them landed on the visible page with a peak.
                                chromatograms_updated <- rbind(chromatograms_updated, baselined_chromatograms)

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
                                        y_axis_start <<- 0
                                        y_axis_end <<- max(chromatograms$abundance)
                                        ## No brush means "reset", and that has to include y --
                                        ## otherwise a stale y window silently clips the full view.
                                        y_zoom <<- NULL
                                        y_zoom_key <<- NULL
                                    }

                                ## If brush is not null, assign brush values to start and end

                                    if ( !is.null(input$chromatogram_brush) ) {
                                        peak_points <<- isolate(brushed_chromatogram(chromatograms_updated, input$chromatogram_brush))

                                        ## An empty selection used to propagate straight through:
                                        ## min()/max() of nothing return Inf/-Inf, the axis limits
                                        ## become Inf..-Inf, the next render filters every row away,
                                        ## and ggplot dies with "Faceting variables must have at
                                        ## least one value". Because the limits are globals the app
                                        ## then stays wedged until it is restarted -- there is no
                                        ## brush you can draw to recover. Keep a usable view instead
                                        ## and say what happened.
                                        ## rt_rt_offset, not rt: the renderer filters the page on
                                        ## rt_rt_offset, so setting the limits from native rt is wrong
                                        ## by exactly rt_offset as soon as retention-time alignment is
                                        ## used -- and if the offset exceeds the brush width the window
                                        ## contains no points at all, so a brush drawn on a visible peak
                                        ## yields the "no points in range" placeholder.
                                        brush_x_col <- if ("rt_rt_offset" %in% names(peak_points)) "rt_rt_offset" else "rt"
                                        brush_ok <- nrow(peak_points) > 0 &&
                                            is.finite(min(peak_points[[brush_x_col]])) && is.finite(max(peak_points[[brush_x_col]]))

                                        if (brush_ok) {
                                            x_axis_start <<- min(peak_points[[brush_x_col]])
                                            x_axis_end <<- max(peak_points[[brush_x_col]])

                                            ## y follows the brush box itself, not the x-slice's
                                            ## full extent. Stored with a key naming the view it
                                            ## was taken in (display mode + y scaling), because
                                            ## the same number means different things in raw
                                            ## counts, square root, log and per-ion-normalised
                                            ## space -- applying it to the wrong one would clip
                                            ## the trace to nothing.
                                            ysel <- brush_y_range(input$chromatogram_brush)
                                            if (!is.null(ysel) && ysel[2] > ysel[1]) {
                                                y_zoom <<- ysel
                                                y_zoom_key <<- last_brush_key
                                                ## y_axis_start/end are the UNKEYED pair, and the renderer
                                                ## hands them to y_limits_for() as its FALLBACK. Writing a
                                                ## sqrt / log / per-ion-normalised / m-z range into them
                                                ## leaked that zoom into every other view: the key check
                                                ## correctly rejected the foreign y_zoom and then fell
                                                ## straight back to the same numbers, so a brush taken in
                                                ## all-ions/sqrt flattened the TIC through oob = squish,
                                                ## and one taken on the ion map gave the TIC a y-axis in
                                                ## m/z. They now only ever hold raw counts.
                                                if (is.null(y_zoom_key) || grepl("^tic/", y_zoom_key)) {
                                                    y_axis_start <<- ysel[1]
                                                    y_axis_end   <<- ysel[2]
                                                } else {
                                                    cat(paste0("  y-zoom stored for ", y_zoom_key,
                                                               " only; the raw-count axis is left as it was.\n"))
                                                }
                                            } else {
                                                y_zoom <<- NULL
                                                y_zoom_key <<- NULL
                                                y_axis_start <<- min(peak_points$abundance)
                                                y_axis_end <<- max(peak_points$abundance)
                                            }
                                        } else {
                                            cat("Brush selected no chromatogram points - keeping the current view.\n")
                                            cat("  Re-brush on a drawn trace and press Shift+Q again, or press Shift+Q with no brush to reset.\n")
                                            logMessage("Brush selected no points; view left unchanged.")
                                            if (is.null(x_axis_start) || length(x_axis_start) == 0 || !is.finite(x_axis_start)) x_axis_start <<- x_axis_start_default
                                            if (is.null(x_axis_end)   || length(x_axis_end)   == 0 || !is.finite(x_axis_end))   x_axis_end   <<- x_axis_end_default
                                            if (is.null(y_axis_start) || length(y_axis_start) == 0 || !is.finite(y_axis_start)) y_axis_start <<- 0
                                            if (is.null(y_axis_end)   || length(y_axis_end)   == 0 || !is.finite(y_axis_end))   y_axis_end   <<- max(chromatograms$abundance)
                                        }
                                    }
                                
                                ## One line that distinguishes every way the zoom can fail: no
                                ## brush registered at all, a brush that produced no selection,
                                ## or limits that were set correctly and then ignored downstream.

                                    cat(paste0(
                                        "Shift+Q: brush ",
                                        if (is.null(input$chromatogram_brush)) "ABSENT" else "present",
                                        "; x-axis now ", signif(x_axis_start, 8), " to ", signif(x_axis_end, 8),
                                        "; y-axis ",
                                        if (is.null(y_zoom)) "auto" else paste0(signif(y_zoom[1], 6), " to ", signif(y_zoom[2], 6)),
                                        "\n"))

                            ## No plot is assembled here. Shift+Q computes and writes; the
                            ## paginated output$chromatograms below draws. v5 moved the drawing
                            ## out but left v4's assembly in place, where it went on building a
                            ## chromatogram_plot that was immediately discarded -- two gsub passes
                            ## and a named-vector build over ONE element per row (~1-2M rows at
                            ## 152 samples), plus a full copy of the table attached to the dead
                            ## plot, on every single Shift+Q. Deleted 2026-10-01; nothing
                            ## downstream read chromatograms_updated_filtered or the plot.

                            ## Add peaks, if any
                                
                                peak_table <- read.csv("peaks_monolist.csv")
                        
                                if (dim(peak_table)[1] > 0) {

                                    ## Filter out duplicate peaks and NA peaks
                                        
                                        ## Dedupe on the (start, end) PAIR. It used to run two
                                        ## independent passes -- one on peak_start, one on peak_end --
                                        ## so two genuinely distinct peaks that happened to share an
                                        ## end bound had the second silently dropped and committed to
                                        ## disk. Sorted by peak_start first, so "first wins" means the
                                        ## earlier peak rather than whichever row was appended first.
                                        peak_table <- peak_table[!is.na(peak_table$peak_start),]
                                        peak_table <- peak_table[order(peak_table$path_to_cdf_csv, peak_table$peak_start, peak_table$peak_end),]
                                        peak_table <- peak_table %>% group_by(path_to_cdf_csv) %>%
                                            mutate(duplicated = duplicated(paste(peak_start, peak_end, sep = "_")))
                                        peak_table <- as.data.frame(peak_table)
                                        peak_table <- dplyr::filter(peak_table, duplicated == FALSE)
                                        peak_table <- peak_table[,!colnames(peak_table) == "duplicated"]

                                    ## Update with peak_number_within_sample
                                        
                                        peak_table_updated <- list()
                                        samples_with_peaks <- unique(as.character(peak_table$path_to_cdf_csv))
                                        for (sample_number in seq_along(samples_with_peaks)) {
                                          peaks_in_this_sample <- peak_table[peak_table$path_to_cdf_csv == samples_with_peaks[sample_number],]
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
                                        
                                        if (isTRUE(integrate_peak_ions)) {
                                            if (is.null(peak_table$area_peak_ions)) peak_table$area_peak_ions <- NA_real_
                                            if (is.null(peak_table$peak_ions))      peak_table$peak_ions      <- NA_character_
                                        }

                                        ## Iterate over the UNIQUE sample list, and index into that
                                        ## same vector. The bound used to be length(unique(...)) while
                                        ## the index went into the raw column, so one duplicated row in
                                        ## samples_monolist.csv made the loop visit an early sample
                                        ## twice and never reach the last one -- whose areas then kept
                                        ## their stale values with no error anywhere.
                                        samples_for_area <- unique(as.character(samples_monolist$path_to_cdf_csv))
                                        for (sample_number in seq_along(samples_for_area)) {

                                          this_sample <- samples_for_area[sample_number]
                                          peaks_in_this_sample <- peak_table[peak_table$path_to_cdf_csv == this_sample,]
                                          
                                          ## Area = TIC sum over the peak window, minus the
                                          ## interpolated baseline over the same window.
                                          ##
                                          ## Fixed 2026-09-30. The baseline lives as ROWS tagged
                                          ## ion == "baseline" carrying their value in `abundance`
                                          ## (long format). It is NOT a `baseline` COLUMN -- that
                                          ## is the WIDE layout of modules/gcms.R, from which this
                                          ## block was copied verbatim. v4 changed the table to
                                          ## long form and adapted only the first sum; the second
                                          ## kept asking for $baseline, which returns NULL, and
                                          ## sum(NULL) is 0. So from v4 until now the subtraction
                                          ## was a silent no-op and every `area` ever written was
                                          ## a raw, un-baseline-corrected TIC sum. The old form
                                          ## also had no ion == 0 restriction on the baseline sum,
                                          ## so it would have summed across the TIC, the baseline
                                          ## and every extracted-ion row.
                                          ##
                                          ## seq_len, not 1:length(): a sample with no peaks yet
                                          ## gave 1:0 == c(1, 0), and the zero index handed
                                          ## dplyr::filter() a length-0 predicate, which errors.
                                          ## That is the normal state of a part-annotated folder,
                                          ## so Shift+Q died on the first un-annotated sample.
                                          areas <- numeric(nrow(peaks_in_this_sample))
                                          for (peak in seq_len(nrow(peaks_in_this_sample))) {
                                            rows_for_peak <- chromatograms_updated[
                                                chromatograms_updated$path_to_cdf_csv ==
                                                    as.character(peaks_in_this_sample$path_to_cdf_csv[peak]), ]
                                            in_window <- dplyr::filter(
                                                rows_for_peak,
                                                rt >= peaks_in_this_sample$peak_start[peak],
                                                rt <= peaks_in_this_sample$peak_end[peak]
                                            )
                                            areas[peak] <-
                                                sum(in_window$abundance[in_window$ion == 0]) -
                                                sum(in_window$abundance[in_window$ion == "baseline"])
                                          }

                                          peak_table$area[peak_table$path_to_cdf_csv == this_sample] <- areas

                                          ## v6: the same area restricted to the ions the peak
                                          ## actually generates. The TIC area above charges a peak
                                          ## for everything under it -- column bleed, a tailing
                                          ## neighbour, a co-eluting second compound. Summing only
                                          ## the ions that rise across the peak window drops that.
                                          ## Reported alongside `area`, not instead of it, so
                                          ## numbers from v5 stay comparable.
                                          if (isTRUE(integrate_peak_ions) && nrow(peaks_in_this_sample) > 0) {
                                            ion_areas <- numeric(nrow(peaks_in_this_sample))
                                            ion_lists <- character(nrow(peaks_in_this_sample))
                                            for (peak in seq_len(nrow(peaks_in_this_sample))) {
                                              res <- peak_ion_area(
                                                as.character(peaks_in_this_sample$path_to_cdf_csv[peak]),
                                                peaks_in_this_sample$peak_start[peak],
                                                peaks_in_this_sample$peak_end[peak],
                                                ## The LIVE control, not the startup argument. These
                                                ## two diverge the moment the slider is touched, so
                                                ## the ion set you tune on screen was not the one
                                                ## that produced area_peak_ions.
                                                if (is.null(input$allion_threshold)) all_ion_threshold else input$allion_threshold
                                              )
                                              ion_areas[peak] <- res$area
                                              ion_lists[peak] <- res$ions
                                            }
                                            rows <- peak_table$path_to_cdf_csv == this_sample
                                            peak_table$area_peak_ions[rows] <- ion_areas
                                            peak_table$peak_ions[rows]      <- ion_lists
                                          }

                                          # peaks_in_this_sample$area <- areas
                                          # peak_table_updated[[sample_number]] <- peaks_in_this_sample
                                        }
                                        # peak_table_updated <- do.call(rbind, peak_table_updated)
                                        # peak_table <- peak_table_updated
                          
                                    ## Write out peaks now with assigned peak_number_within_sample and RT offset, update the peak_table in ui
                            
                                        write.table(peak_table, file = "peaks_monolist.csv", col.names = TRUE, sep = ",", row.names = FALSE)

                                        render_peak_table()

                                }

                            ## The x-window globals must be usable before the paginated renderer
                            ## filters on them. Hoisted out of the peak block 2026-10-01: it used
                            ## to sit inside `if (dim(peak_table)[1] > 0)`, so on a folder with no
                            ## peaks yet the safety net did not run at all.
                                if (length(x_axis_start) == 0) {x_axis_start <<- min(chromatograms$rt)}
                                if (length(x_axis_end) == 0) {x_axis_end <<- max(chromatograms$rt)}

                                ## Drawing moved out of Shift+Q: the chromatogram is now drawn by
                                ## the reactive output$chromatograms defined just below this handler,
                                ## which renders only the current PAGE of samples (see pagination
                                ## block). Here we only signal that fresh data is ready.
                                redraw_trigger(isolate(redraw_trigger()) + 1)

                                cat("Chromatogram updated.\n")
                        }
                    })

                ## ------------------------------------------------------------------
                ## Paginated chromatogram renderer + page controls
                ##
                ## A single faceted plot (one panel per sample) is what crashed v4 at
                ## 100-200 samples: render cost scales with the number of panels, and a
                ## ~20,000px image exhausts the graphics device (plot greys out, app
                ## goes unresponsive). Here we draw only samples_per_page samples at a
                ## time. Shift+Q recomputes the data and bumps redraw_trigger; Prev/Next
                ## just change current_page. Paging does NOT re-run the (expensive)
                ## baseline / area computation, so moving between views stays snappy.
                ## ------------------------------------------------------------------

                    ## Which samples exist, and the current page's slice
                        page_bounds <- function() {
                            if (is.null(chromatograms_updated)) return(NULL)
                            all_samples <- sort(unique(as.character(chromatograms_updated$path_to_cdf_csv)))
                            n <- length(all_samples)
                            ## No samples at all: return an EMPTY page rather than
                            ## all_samples[1:0], which gives NA_character_ -- a phantom sample
                            ## that the renderer then tried to facet on, under an indicator
                            ## reading "samples 1-0 of 0".
                            if (n == 0) {
                                return(list(all_samples = character(0), n = 0L, n_pages = 1L,
                                            pg = 1L, lo = 0L, hi = 0L, page_samples = character(0)))
                            }
                            n_pages <- max(1, ceiling(n / samples_per_page))
                            pg <- min(max(1, current_page()), n_pages)
                            lo <- (pg - 1) * samples_per_page + 1
                            hi <- min(pg * samples_per_page, n)
                            list(all_samples = all_samples, n = n, n_pages = n_pages,
                                 pg = pg, lo = lo, hi = hi,
                                 page_samples = all_samples[seq(lo, hi)])
                        }

                    ## Page navigation buttons
                        observeEvent(input$prev_page, {
                            current_page(max(1, current_page() - 1))
                        })
                        observeEvent(input$next_page, {
                            b <- page_bounds()
                            if (is.null(b)) { cat("Press Shift+Q to load chromatograms first.\n"); return() }
                            current_page(min(b$n_pages, current_page() + 1))
                        })

                    ## Page indicator text (e.g. "Page 2 / 8   (samples 21-40 of 152)")
                        output$page_indicator <- renderText({
                            redraw_trigger(); current_page()
                            b <- page_bounds()
                            if (is.null(b)) return("Press Shift+Q to load chromatograms")
                            if (b$n == 0) return("No samples loaded")
                            paste0("Page ", b$pg, " / ", b$n_pages,
                                   "   (samples ", b$lo, "–", b$hi, " of ", b$n, ")")
                        })

                    ## The chromatogram plot itself - reactive on redraw_trigger (Shift+Q)
                    ## and current_page (Prev/Next). Axis limits come from the globals that
                    ## Shift+Q set from the brush/defaults, so brushing does not itself force
                    ## a redraw (matches v4: press Q to refresh the view).
                        output$chromatograms <- renderPlot({

                            redraw_trigger()
                            b <- page_bounds()

                            ## Returning NULL here draws a plot with NO coordinate map, and the
                            ## browser then holds a degenerate 0..1 map for this output until the
                            ## next successful render -- one route to the fraction-coordinate
                            ## brushes seen on host1. Hand back a real ggplot with a real x
                            ## mapping instead, so the output always has a usable coordmap.
                            if (is.null(b)) {
                                return(
                                    ggplot(data.frame(rt_rt_offset = range(chromatograms$rt), abundance = c(0, 0)),
                                           aes(x = rt_rt_offset, y = abundance)) +
                                        geom_blank() +
                                        annotate("text", x = mean(range(chromatograms$rt)), y = 0, size = 5,
                                                 label = "Press Shift+Q to load chromatograms") +
                                        theme_classic()
                                )
                            }
                            page_samples <- b$page_samples

                            ## Axis limits set by Shift+Q (globals); guard against unset
                                xs <- x_axis_start; xe <- x_axis_end
                                if (length(xs) == 0 || is.null(xs)) xs <- min(chromatograms$rt)
                                if (length(xe) == 0 || is.null(xe)) xe <- max(chromatograms$rt)
                                ys <- if (!is.null(y_axis_start) && length(y_axis_start) > 0) y_axis_start else 0
                                ye <- if (!is.null(y_axis_end)   && length(y_axis_end)   > 0) y_axis_end   else max(chromatograms$abundance)

                            ## Display settings (v6). The controls may not have rendered on the
                            ## first pass, so each falls back to the function argument.
                                mode_display <- if (is.null(input$display_mode)) {
                                    if (isTRUE(all_ion_view)) "all_ions" else "tic"
                                } else input$display_mode
                                mode_scale <- if (is.null(input$allion_scaling))    all_ion_scaling    else input$allion_scaling
                                ## The two EXPENSIVE controls are read through isolate(), so this
                                ## renderer does not take a reactive dependency on them. Dragging
                                ## either slider used to queue one multi-million-segment redraw per
                                ## intermediate value, each recomputing the ion window for the whole
                                ## page. They now apply on the next Shift+Q, which is what the help
                                ## text has always said. Display mode and y scaling stay live: both
                                ## are cheap and users expect them to be instant.
                                thr        <- isolate(if (is.null(input$allion_threshold))  all_ion_threshold  else input$allion_threshold)
                                npx        <- isolate(if (is.null(input$allion_max_points)) all_ion_max_points else input$allion_max_points)

                            ## Base data for this page only.
                            ## INCLUSIVE bounds, matching integration. Integration uses >= / <=
                            ## while the overlay and this view filter used > / <, so a peak sitting
                            ## flush against the zoom edge was integrated but vanished from the
                            ## plot, and a ribbon excluded the two boundary scans its own area
                            ## included.
                                cuf <- dplyr::filter(
                                    chromatograms_updated,
                                    path_to_cdf_csv %in% page_samples,
                                    rt_rt_offset >= xs & rt_rt_offset <= xe
                                )

                                ## EVERY sample on this page gets a panel, drawn or not. The facet
                                ## column is a factor over page_samples and the facets are built
                                ## with drop = FALSE, so a sample with nothing over the ion
                                ## threshold shows as an empty panel rather than vanishing. It used
                                ## to disappear entirely: a page returned 14 panels while the
                                ## indicator still read "samples 21-40 of 152", which reads as
                                ## "those six samples failed".
                                as_page_facet <- function(x) factor(as.character(x), levels = page_samples)
                                cuf$path_to_cdf_csv <- as_page_facet(cuf$path_to_cdf_csv)

                                ## Per SAMPLE, not per ROW. labeller() only ever looks up one
                                ## entry per facet, but this built a named vector as long as the
                                ## table -- two gsub passes and a names<- over ~1-2M elements on
                                ## every redraw, including every page flip.
                                facet_labels <- sample_id_from_path(page_samples)
                                names(facet_labels) <- page_samples

                            ## Gather the all-ion data for the page, one sample at a time.
                            ## Each sample carries its own rt offset and its own baseline, so
                            ## the window is queried in that sample's native rt and the offset
                            ## is added back afterwards. A sample that yields nothing (no ion
                            ## clears the threshold in this window) is simply skipped.
                                allion <- NULL
                                if (mode_display %in% c("all_ions", "ion_map")) {
                                    npx_used <- if (mode_display == "ion_map") min(npx, 600) else npx
                                    ## Read the offsets FRESH. Shift+Q re-reads
                                    ## samples_monolist into a handler-LOCAL variable, so this
                                    ## renderer used to resolve the name to the startup copy:
                                    ## editing an rt_offset moved the TIC and the peak marks
                                    ## but not the ion lines, which in all-ions mode manufactures
                                    ## exactly the apex shift the view exists to detect. One
                                    ## small CSV per render is nothing beside the draw.
                                    offsets_now <- if (file.exists("samples_monolist.csv")) {
                                        try(read.csv("samples_monolist.csv", stringsAsFactors = FALSE), silent = TRUE)
                                    } else NULL
                                    if (inherits(offsets_now, "try-error") || is.null(offsets_now) ||
                                        !all(c("rt_offset", "path_to_cdf_csv") %in% names(offsets_now))) {
                                        offsets_now <- samples_monolist
                                    }
                                    pieces <- list()
                                    for (s in page_samples) {
                                        off <- offsets_now$rt_offset[match(s, offsets_now$path_to_cdf_csv)]
                                        if (length(off) == 0 || is.na(off)) off <- 0
                                        w <- try(allion_window(s, xs - off, xe - off, thr, npx_used, rt_shift = off), silent = TRUE)
                                        if (inherits(w, "try-error") || is.null(w) || nrow(w) == 0) next
                                        w <- allion_scale(w, mode_scale)
                                        ## Round the reconstructed grid. allion_window() snaps every
                                        ## sample to one shared set of centres in offset-corrected
                                        ## space, but it hands them back shifted into native rt, and
                                        ## (c - off) + off is not bit-identical to c. Those ~1e-15
                                        ## differences make unique() see one grid per sample again,
                                        ## which is what made the ion map's tiles too narrow. Six
                                        ## decimals is orders of magnitude finer than any real
                                        ## retention-time resolution.
                                        w$rt_rt_offset  <- round(w$rt + off, 6)
                                        w$path_to_cdf_csv <- s
                                        pieces[[length(pieces) + 1L]] <- w
                                    }
                                    if (length(pieces) > 0) allion <- as.data.frame(data.table::rbindlist(pieces))
                                    if (!is.null(allion)) allion$path_to_cdf_csv <- as_page_facet(allion$path_to_cdf_csv)

                                    ## Single-detector data (GC-TCD/FID/ECD) carries the mz = 0
                                    ## sentinel, and the all-ion views colour and position by m/z
                                    ## on a LOG scale: log10(0) is -Inf, so every line came back
                                    ## grey with an empty legend and the ion map's y axis collapsed
                                    ## to c(0, 0). Fall back to the TIC, as Shift+1 and Shift+4
                                    ## already do for the same reason.
                                    if (!is.null(allion) && all(allion$mz == 0, na.rm = TRUE)) {
                                        cat("All-ion view is not available for single-detector data (no m/z dimension) - showing the TIC.\n")
                                        mode_display <- "tic"
                                        allion <- NULL
                                    }
                                }

                                ## Nothing to draw in the requested mode -> fall back to the TIC
                                ## rather than showing an empty panel.
                                if (mode_display %in% c("all_ions", "ion_map") && is.null(allion)) mode_display <- "tic"

                                ## Nothing in range at all -- e.g. axis limits left over from a
                                ## failed brush. facet_grid() errors on an empty frame ("Faceting
                                ## variables must have at least one value"), which in a renderPlot
                                ## leaves a broken panel and no way back, so say what happened.
                                if (nrow(cuf) == 0 && is.null(allion)) {
                                    ## Carry a REAL x mapping over the current window. The old
                                    ## placeholder was annotate(x = 0, y = 0) + theme_void(): a
                                    ## degenerate 0..0 domain with no panel, which is exactly the
                                    ## "never hand back a plot with no coordinate map" trap -- and
                                    ## it appeared precisely when the user was trying to brush
                                    ## their way out of a bad zoom, leaving the next brush with
                                    ## nothing sane to convert against.
                                    mid <- if (is.finite(xs) && is.finite(xe)) (xs + xe) / 2 else 0
                                    return(
                                        ggplot(data.frame(rt_rt_offset = c(xs, xe), abundance = c(0, 1)),
                                               aes(x = rt_rt_offset, y = abundance)) +
                                            geom_blank() +
                                            annotate("text", x = mid, y = 0.5, size = 5, lineheight = 1.2,
                                                     label = paste0(
                                                         "No chromatogram points between x = ",
                                                         signif(xs, 6), " and ", signif(xe, 6), ".\n",
                                                         "Press Shift+Q with no brush to reset the view.")) +
                                            scale_x_continuous(limits = c(xs, xe), name = "Retention (Scan number)") +
                                            scale_y_continuous(name = "Abundance (counts)") +
                                            theme_classic()
                                    )
                                }

                            ## A y-zoom only means something in the view it was taken in --
                            ## raw counts, square root, log, per-ion-normalised and m/z are five
                            ## different axes. Tag it with the view, and ignore it elsewhere.
                            ## drawn_yrange is what the panel spans when no zoom applies, and is
                            ## also what a later brush's y fractions get measured against.
                                view_key <- paste(mode_display, mode_scale, sep = "/")
                                drawn_yrange <- if (mode_display == "ion_map" && !is.null(allion)) {
                                    range(allion$mz, na.rm = TRUE)
                                } else if (mode_display == "all_ions" && !is.null(allion)) {
                                    range(allion$y, na.rm = TRUE)
                                } else {
                                    c(ys, ye)
                                }

                            ## Base plot, per display mode
                                if (mode_display == "all_ions") {

                                    chromatogram_plot <- ggplot() +
                                        geom_line(
                                            data = allion,
                                            mapping = aes(x = rt_rt_offset, y = y, colour = mz, alpha = opacity,
                                                          group = interaction(path_to_cdf_csv, mz)),
                                            linewidth = 0.18
                                        ) +
                                        scale_alpha_identity() +
                                        ## log colour: nearly every ion in a GC-MS run sits
                                        ## between 40 and 150, so a linear m/z ramp renders the
                                        ## whole plot one shade of purple and wastes the top
                                        ## two-thirds of the scale on a handful of heavy ions.
                                        scale_colour_viridis_c(name = "m/z", trans = "log10") +
                                        scale_x_continuous(limits = c(xs, xe), name = "Retention (Scan number)") +
                                        scale_y_continuous(name = allion_y_label(mode_scale),
                                                           limits = y_limits_for(view_key, drawn_yrange),
                                                           oob = scales::squish) +
                                        facet_grid(path_to_cdf_csv~., scales = "free_y", drop = FALSE, labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                        theme_classic() +
                                        guides(fill = "none")

                                } else if (mode_display == "ion_map") {

                                    ## geom_tile rather than geom_raster: decimation can leave
                                    ## the rt grid with gaps, which geom_raster rejects.
                                    ## Tile width from the ACTUAL spacing of the drawn grid.
                                    ## Dividing the window by the COUNT of unique x values was
                                    ## several times too narrow whenever the union across samples
                                    ## was denser than any one sample's grid -- the map rendered as
                                    ## thin stripes with white gaps, hidden on a single-sample page.
                                    ## It is also wrong whenever decimation does not fire at all.
                                    ## Samples now share one grid (see allion_window's rt_shift),
                                    ## so the median gap is the true tile width.
                                    tile_x <- sort(unique(allion$rt_rt_offset))
                                    tile_w <- if (length(tile_x) > 1) stats::median(diff(tile_x)) else (xe - xs)
                                    if (!is.finite(tile_w) || tile_w <= 0) tile_w <- (xe - xs) / max(1, length(tile_x))
                                    chromatogram_plot <- ggplot() +
                                        geom_tile(
                                            data = allion,
                                            mapping = aes(x = rt_rt_offset, y = mz, fill = y),
                                            width = tile_w, height = 1
                                        ) +
                                        scale_fill_viridis_c(name = allion_y_label(mode_scale)) +
                                        ## coord_cartesian, not scale limits: a tile straddling
                                        ## the window edge should be clipped, not dropped.
                                        ## Both limits go through coord_cartesian, not the scales:
                                        ## tiles have width and height, so a tile straddling either
                                        ## edge should be clipped, not dropped.
                                        coord_cartesian(xlim = c(xs, xe),
                                                        ylim = y_limits_for(view_key, drawn_yrange)) +
                                        scale_x_continuous(name = "Retention (Scan number)") +
                                        scale_y_continuous(name = "m/z") +
                                        facet_grid(path_to_cdf_csv~., scales = "free_y", drop = FALSE, labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                        theme_classic()

                                } else {

                                    chromatogram_plot <- ggplot() +
                                        geom_line(
                                            data = filter(cuf, ion == "baseline"),
                                            mapping = aes(x = rt_rt_offset, y = abundance), color = "grey"
                                        ) +
                                        scale_x_continuous(limits = c(xs, xe), name = "Retention (Scan number)") +
                                        scale_y_continuous(limits = y_limits_for(view_key, c(ys, ye)),
                                                           name = "Abundance (counts)", oob = scales::squish) +
                                        facet_grid(path_to_cdf_csv~., scales = "free_y", drop = FALSE, labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                        theme_classic() +
                                        guides(fill = "none") +
                                        scale_fill_continuous(type = "viridis") +
                                        scale_color_manual(values = discrete_palette)

                                }

                            ## Peak overlay (read the table Shift+Q wrote), page + x-axis filtered.
                            ## The shaded ribbons are drawn in raw counts, so they only make sense
                            ## against the TIC's raw y-axis; in the other modes the peak bounds and
                            ## numbers are still drawn, but the shading is dropped.
                                if (file.exists("peaks_monolist.csv")) {
                                    peak_table <- read.csv("peaks_monolist.csv")
                                    needed <- c("peak_start_rt_offset", "peak_end_rt_offset", "peak_number_within_sample")
                                    if (nrow(peak_table) > 0 && all(needed %in% colnames(peak_table))) {
                                        peak_table <- dplyr::filter(
                                            peak_table,
                                            path_to_cdf_csv %in% page_samples,
                                            peak_start_rt_offset >= xs & peak_end_rt_offset <= xe
                                        )
                                        all_ribbons <- list(); all_vlines <- list(); all_labels <- list()
                                        if (nrow(peak_table) > 0) {
                                            for (peak in seq_len(nrow(peak_table))) {
                                                signal_for_this_peak <- dplyr::filter(
                                                    chromatograms_updated,
                                                    path_to_cdf_csv == peak_table[peak, ]$path_to_cdf_csv,
                                                    rt_rt_offset >= peak_table[peak, ]$peak_start_rt_offset,
                                                    rt_rt_offset <= peak_table[peak, ]$peak_end_rt_offset
                                                )
                                                if (nrow(signal_for_this_peak) > 0) {
                                                    signal_for_this_peak$peak_number_within_sample <- peak_table$peak_number_within_sample[peak]
                                                    ribbon <- dplyr::filter(signal_for_this_peak, ion == 0)
                                                    ## Length-mismatch guard: if this sample carries no ion == "baseline"
                                                    ## rows, the filter returns a length-0 vector and assigning it into an
                                                    ## n-row frame kills the entire render with no way back. Skip this
                                                    ## peak's shading with a message instead.
                                                    baseline_for_this_peak <- dplyr::filter(signal_for_this_peak, ion == "baseline")$abundance
                                                    if (nrow(ribbon) == 0 || length(baseline_for_this_peak) != nrow(ribbon)) {
                                                        cat(paste0("Skipping peak overlay for ", peak_table[peak, ]$path_to_cdf_csv,
                                                                   ", peak ", peak_table$peak_number_within_sample[peak], ": ",
                                                                   nrow(ribbon), " TIC rows vs ", length(baseline_for_this_peak),
                                                                   " baseline rows. Press Shift+Q to rebuild the baseline.\n"))
                                                        next
                                                    }
                                                    ribbon$baseline <- baseline_for_this_peak
                                                    all_ribbons[[peak]] <- ribbon
                                                    all_vlines[[peak]]  <- signal_for_this_peak[1, ]
                                                    all_labels[[peak]]  <- dplyr::filter(signal_for_this_peak, ion == 0) %>%
                                                        dplyr::summarize(
                                                            peak_number_within_sample = peak_number_within_sample[1],
                                                            x = median(rt_rt_offset),
                                                            y = max(abundance),
                                                            path_to_cdf_csv = path_to_cdf_csv[1]
                                                        )
                                                }
                                            }
                                            all_ribbons <- dplyr::bind_rows(all_ribbons)
                                            all_vlines  <- dplyr::bind_rows(all_vlines)
                                            all_labels  <- dplyr::bind_rows(all_labels)
                                            ## Same facet factor as the base layer, or these land
                                            ## in panels of their own.
                                            if (nrow(all_ribbons) > 0) all_ribbons$path_to_cdf_csv <- as_page_facet(all_ribbons$path_to_cdf_csv)
                                            if (nrow(all_vlines)  > 0) all_vlines$path_to_cdf_csv  <- as_page_facet(all_vlines$path_to_cdf_csv)
                                            if (nrow(all_labels)  > 0) all_labels$path_to_cdf_csv  <- as_page_facet(all_labels$path_to_cdf_csv)
                                            if (nrow(all_ribbons) > 0) {
                                                if (mode_display == "tic") {
                                                    chromatogram_plot <- chromatogram_plot +
                                                        geom_vline(data = all_vlines, mapping = aes(xintercept = rt_rt_offset), alpha = 0.3) +
                                                        geom_ribbon(
                                                            data = all_ribbons,
                                                            mapping = aes(x = rt_rt_offset, ymax = abundance, ymin = baseline,
                                                                          fill = peak_number_within_sample, group = peak_number_within_sample),
                                                            alpha = 0.8
                                                        ) +
                                                        geom_text(data = all_labels, mapping = aes(label = peak_number_within_sample, x = x, y = y), color = "black")
                                                } else {
                                                    chromatogram_plot <- chromatogram_plot +
                                                        geom_vline(data = all_vlines, mapping = aes(xintercept = rt_rt_offset), alpha = 0.3)
                                                }
                                            }
                                        }
                                    }
                                }

                            ## TIC trace on top (page samples only). In the all-ion modes the
                            ## overlay IS the chromatogram, and a raw-counts TIC would sit far
                            ## off the transformed y-scale.
                                if (mode_display == "tic") {
                                    chromatogram_plot <- chromatogram_plot +
                                        geom_line(
                                            data = filter(cuf, ion != "baseline"),
                                            mapping = aes(x = rt_rt_offset, y = abundance, color = ion),
                                            alpha = 0.8
                                        )
                                }

                            ## Kept so a mapping-less brush can be re-expressed against the
                            ## panel's real bounds. clientData carries the size the browser
                            ## actually rendered at, which is what the fractions refer to.
                                last_brush_plot   <<- chromatogram_plot
                                last_brush_size   <<- suppressWarnings(as.numeric(c(
                                    session$clientData$output_chromatograms_width,
                                    session$clientData$output_chromatograms_height)))
                                last_brush_key    <<- view_key
                                last_brush_facets <<- if (mode_display %in% c("all_ions", "ion_map") && !is.null(allion)) {
                                    length(unique(as.character(allion$path_to_cdf_csv)))
                                } else {
                                    length(unique(as.character(cuf$path_to_cdf_csv)))
                                }
                                last_brush_yrange <<- y_limits_for(view_key, drawn_yrange)

                            chromatogram_plot
                        })

                ## Transfer chromatogram_brush info to selected_peak table

                    output$selected_peak <- DT::renderDataTable(DT::datatable({

                        if ( !is.null(input$chromatogram_brush )) {
                            peak_points <- brushed_chromatogram(chromatograms_updated, input$chromatogram_brush)
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

                            peak_points <- brushed_chromatogram(chromatograms_updated, input$chromatogram_brush)
                            if (nrow(peak_points) == 0) {
                                cat("Shift+E: the brush selected no points; nothing excised.\n")
                                return()
                            }
                            path_to_cdf_csv <- single_sample(peak_points, "Shift+E")
                            if (is.null(path_to_cdf_csv)) return()
                            selection_start = min(peak_points$rt)
                            selection_end   = max(peak_points$rt)

                            peak_table <- read.csv("peaks_monolist.csv")
                            if (nrow(peak_table) == 0) {
                                cat("Shift+E: no peaks on file.\n")
                                return()
                            }

                            # Remove any peak that fully resides within [selection_start, selection_end]
                            # for that single cdf. INCLUSIVE bounds: with `>` and `<` a peak could not
                            # be excised by the very brush that created it, because Shift+A set its
                            # bounds to exactly this selection's min and max.
                            peak_table <- peak_table[!
                              apply(cbind(
                                peak_table$peak_start >= selection_start,
                                peak_table$peak_end <= selection_end,
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

                            peak_points <- brushed_chromatogram(chromatograms_updated, input$chromatogram_brush)
                            if (nrow(peak_points) == 0) {
                                cat("Shift+R: the brush selected no points; nothing removed.\n")
                                return()
                            }

                            # read peak table
                            peak_table <- read.csv("peaks_monolist.csv")
                            if (nrow(peak_table) == 0) {
                                cat("Shift+R: no peaks on file.\n")
                                return()
                            }

                            # A GLOBAL removal has to compare like with like. The brush is in the
                            # brushed sample's NATIVE rt, while every other sample's peak bounds are
                            # in its own native rt -- so comparing the two across samples is wrong by
                            # each sample's rt_offset. Use the offset-corrected pair on both sides
                            # when the table carries it (it does after any Shift+Q), and fall back to
                            # native rt only on a table that has never been through one.
                            aligned <- all(c("peak_start_rt_offset", "peak_end_rt_offset") %in% names(peak_table)) &&
                                       "rt_rt_offset" %in% names(peak_points) &&
                                       !any(is.na(peak_table$peak_start_rt_offset))
                            if (aligned) {
                                selection_start <- min(peak_points$rt_rt_offset)
                                selection_end   <- max(peak_points$rt_rt_offset)
                                starts <- peak_table$peak_start_rt_offset
                                ends   <- peak_table$peak_end_rt_offset
                            } else {
                                selection_start <- min(peak_points$rt)
                                selection_end   <- max(peak_points$rt)
                                starts <- peak_table$peak_start
                                ends   <- peak_table$peak_end
                                cat("  (no rt_offset columns yet - comparing native rt; press Shift+Q first if samples are aligned)\n")
                            }

                            # remove ANY peak in ANY file that fully resides in [start, end];
                            # inclusive, for the same reason as Shift+E
                            peak_table <- peak_table[!
                              apply(cbind(
                                starts >= selection_start,
                                ends <= selection_end
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
                            
                                ## Evaluate the brush ONCE. It was called four times here,
                                ## each re-running the whole coordinate conversion and printing
                                ## its own diagnostic line.
                                pp <- brushed_chromatogram(chromatograms_updated, input$chromatogram_brush)

                                ## An empty selection must not be written. brushed_chromatogram()
                                ## returns a 0-row frame on several paths (brush outside the data,
                                ## unusable x range, nothing matched), and min/max of nothing give
                                ## Inf/-Inf -- which used to be appended to peaks_monolist.csv as
                                ## a permanent poison row that the Shift+Q dedupe does not drop.
                                if (nrow(pp) == 0) {
                                    cat("Shift+A: the brush selected no points; no peak added.\n")
                                    return()
                                }
                                pp_sample <- single_sample(pp, "Shift+A")
                                if (is.null(pp_sample)) return()

                                ## area = TIC minus baseline, both read out of the LONG format
                                ## (rows tagged ion == 0 / ion == "baseline", value in `abundance`).
                                ## Fixed 2026-09-30: this asked for $tic and $baseline, which are
                                ## columns of modules/gcms.R's WIDE table, not of this one. Both
                                ## were NULL, sum(NULL) is 0, so every hand-added peak was written
                                ## with area = 0.
                                write.table(
                                    x = data.frame(
                                            peak_start = min(pp$rt),
                                            peak_end = max(pp$rt),
                                            peak_ID = "unknown",
                                            path_to_cdf_csv = pp_sample,
                                            area = sum(pp$abundance[pp$ion == 0]) - sum(pp$abundance[pp$ion == "baseline"])
                                        ),
                                    file = "peaks_monolist.csv",
                                    append = TRUE,
                                    row.names = FALSE,
                                    col.names = FALSE,
                                    sep = ","
                                )

                                render_peak_table()
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
                            
                                ## Evaluate the brush once (see Shift+A) and refuse an empty one.
                                pp <- brushed_chromatogram(chromatograms_updated, input$chromatogram_brush)
                                if (nrow(pp) == 0) {
                                    cat("Shift+G: the brush selected no points; no peaks added.\n")
                                    return()
                                }

                                x_peaks <-  data.frame(
                                                peak_start = min(pp$rt_rt_offset),
                                                peak_end = max(pp$rt_rt_offset),
                                                peak_ID = "unknown",
                                                path_to_cdf_csv = unique(chromatograms_updated$path_to_cdf_csv),
                                                area = NA_real_
                                            )

                                ## Shift+G writes the SAME window to every sample, so the area has
                                ## to be integrated per sample rather than shared. The old code
                                ## computed one scalar from the brushed sample and gave it to all
                                ## of them -- and that scalar was 0 anyway, because $tic/$baseline
                                ## are columns of modules/gcms.R's wide table, not of this one.
                                ## Fixed 2026-09-30.
                                for (row_i in seq_len(nrow(x_peaks))) {
                                    win <- dplyr::filter(
                                        chromatograms_updated,
                                        path_to_cdf_csv == x_peaks$path_to_cdf_csv[row_i],
                                        rt_rt_offset >= x_peaks$peak_start[row_i],
                                        rt_rt_offset <= x_peaks$peak_end[row_i]
                                    )
                                    x_peaks$area[row_i] <-
                                        sum(win$abundance[win$ion == 0]) -
                                        sum(win$abundance[win$ion == "baseline"])
                                }

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

                                render_peak_table()
                                cat("Added global peak.\n")
                            }
                    })

                ## [MS1 extract ("shift+1"), update ("shift+2"), subtract ("shift+3"), library search ("shift+4"), save ("shift+5")]
                    
                    observeEvent(input$keypress, {

                        ## If "shift+1", MS from chromatogram brush -> MS_out_1
                            
                            if( input$keypress == 33 ) {

                                ## Evaluate the brush once (it was converted three times here),
                                ## refuse an empty one, and refuse a selection that spans more than
                                ## one sample rather than extracting the spectrum from whichever
                                ## sample happens to sort first.
                                ms_pp <- brushed_chromatogram(chromatograms_updated, input$chromatogram_brush)
                                if (is.null(ms_pp) || nrow(ms_pp) == 0) {
                                    cat("Shift+1: the brush selected no points; no spectrum extracted.\n")
                                    return()
                                }
                                sample_name_MS <- single_sample(ms_pp, "Shift+1")
                                if (is.null(sample_name_MS)) return()
                                ret_start_MS <- min(ms_pp$rt)
                                ret_end_MS <- max(ms_pp$rt)

                                chromatogram_updated_MS <- filter(chromatograms_updated, path_to_cdf_csv == sample_name_MS)

                                MS_ret_start_line <<- chromatogram_updated_MS$rt_first_row_in_raw[which.min(abs(
                                    chromatogram_updated_MS$rt - ret_start_MS
                                ))] + 1

                                MS_ret_end_line <<- chromatogram_updated_MS$rt_last_row_in_raw[which.min(abs(
                                    chromatogram_updated_MS$rt - ret_end_MS
                                ))] + 1

                                if (.Platform$OS.type == "unix") {
                                    
                                    ## shQuote every path: CDF_directory_path is caller-supplied and
                                    ## a space in it used to split the command, so `head`/`sed` wrote
                                    ## temp_MS.csv somewhere else entirely and the read below picked
                                    ## up the PREVIOUS peak's spectrum with no error.
                                    ms_src <- shQuote(file.path(CDF_directory_path, sample_name_MS))
                                    ms_tmp <- shQuote(file.path(CDF_directory_path, "temp_MS.csv"))
                                    system(paste0("head -1 ", ms_src, " > ", ms_tmp))
                                    system(paste0("sed -n ", MS_ret_start_line, ",", MS_ret_end_line, "p ", shQuote(sample_name_MS), " >> ", ms_tmp))
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

                                if (is.null(MS_out_1)) {
                                    cat("No mass spectrum extracted yet.\n")
                                    return()
                                } else {
                                
                                    ## Same again: one conversion, and the subtraction has to come
                                    ## from a single sample's panel.
                                    sub_pp <- isolate(brushed_chromatogram(chromatograms_updated, input$chromatogram_brush))
                                    if (is.null(sub_pp) || nrow(sub_pp) == 0) {
                                        cat("Shift+3: the brush selected no points; nothing subtracted.\n")
                                        return()
                                    }
                                    sub_sample <- single_sample(sub_pp, "Shift+3")
                                    if (is.null(sub_sample)) return()
                                    framedDataFile_to_subtract <- isolate(as.data.frame(
                                                        data.table::fread(sub_sample)
                                    ))
                                    framedDataFile_to_subtract <- isolate(dplyr::filter(
                                                            framedDataFile_to_subtract,
                                                            rt > min(sub_pp$rt),
                                                            rt < max(sub_pp$rt)
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

                                if (is.null(MS_out_1)) {
                                    cat("No mass spectrum extracted yet.\n")
                                    return()
                                } else if (nrow(MS_out_1) == 0 || (nrow(MS_out_1) == 1 && all(MS_out_1$mz == 0))) {
                                    # Single-detector data (e.g. GC-TCD): no mass spectrum exists.
                                    output$massSpectra_1 <- renderPlot({
                                        ggplot() +
                                            annotate("text", x = 0.5, y = 0.5,
                                                     label = "No mass spectrum available\n(non-MS detector - e.g. TCD, FID, ECD)",
                                                     size = 6, colour = "grey40") +
                                            theme_void() +
                                            xlim(0, 1) + ylim(0, 1)
                                    })
                                } else {

                                    # Normalize to max abu 100
                                        MS_out_1$intensity <- MS_out_1$intensity*100/max(MS_out_1$intensity)

                                    # Set ranges on mass spec (allow zooming in with selection), Sometimes brush returns Inf or -Inf if the range gets too small, so if that happens, just go up to the main view
                                
                                        if (isolate(is.null(input$massSpectra_1_brush))) {
                                            MS1_low_x_limit <- 0; MS1_high_x_limit <- 1200; MS1_high_y_limit <- 110
                                        } else {
                                            ## One conversion, not two, and through brushed_ms() so an
                                            ## empty coordmap cannot error the whole block away.
                                            ms1_sel <- isolate(brushed_ms(MS_out_1, input$massSpectra_1_brush))
                                            if (is.null(ms1_sel) || nrow(ms1_sel) == 0) {
                                                cat("MS brush selected no peaks - showing the full spectrum.\n")
                                                MS1_low_x_limit <- 0; MS1_high_x_limit <- 1200; MS1_high_y_limit <- 110
                                            } else {
                                                MS1_low_x_limit <- min(ms1_sel$mz)
                                                MS1_high_x_limit <- max(ms1_sel$mz)
                                                in_win <- dplyr::filter(MS_out_1, mz >= MS1_low_x_limit & mz <= MS1_high_x_limit)$intensity
                                                MS1_high_y_limit <- if (length(in_win) > 0) max(in_win) + 8 else 110
                                            }
                                        }
                                        if (MS1_low_x_limit %in% c(Inf, -Inf) | MS1_high_x_limit %in% c(Inf, -Inf) | MS1_high_y_limit %in% c(Inf, -Inf)) {
                                            MS1_low_x_limit <- 0; MS1_high_x_limit <- 1200; MS1_high_y_limit <- 110
                                        }

                                    # Make plot

                                        output$massSpectra_1 <- renderPlot({

                                            ms1_plot <- ggplot() +
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

                                            ## Kept so a mapping-less brush on this plot can be
                                            ## re-expressed against the panel's real bounds, exactly
                                            ## as the chromatogram's is.
                                            last_ms_plot <<- ms1_plot
                                            last_ms_size <<- suppressWarnings(as.numeric(c(
                                                session$clientData$output_massSpectra_1_width,
                                                session$clientData$output_massSpectra_1_height)))
                                            ms1_plot
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
                                message("UPDATED3!")
                                # Read in the peak table (ensure your file is in the correct format)
                                    peak_table <<- read.csv("peaks_monolist.csv", stringsAsFactors = FALSE)
                                    if(nrow(peak_table) == 0){
                                        cat("No peaks available.\n")
                                        return()
                                    }

                                ### somethign to do with the duplicate names
                                ### also all predicted peaks ID should be prefixed with "prediction"
                                ### perhaps it should only output predictions on things labelled "unknown"?

                                # Loop over each detected peak.
                                #
                                # The framed CSV is read once per RUN of peaks from the same sample,
                                # not once per peak. It used to be fread() inside this loop, so a
                                # sample with 40 peaks read its whole m/z cube 40 times -- the single
                                # dominant cost of Shift+4. Shift+Q now writes peaks_monolist.csv
                                # sorted by (sample, peak_start), so in practice that is one read per
                                # sample. Deliberately NOT re-sorted here: `predictions` is still
                                # matched to peak_table POSITIONALLY (finding A4, deferred), so
                                # changing the visiting order would change which peak gets which
                                # label. One cube is held at a time either way.
                                    ms_list <- list()
                                    loaded_file <- NULL; df <- NULL
                                    for(i in seq_len(nrow(peak_table))){
                                        # message(paste0("\n*getting spectrum for peak ", i))
                                        peak <- peak_table[i, ]
                                        # Get the corresponding .CDF.csv file
                                        cdf_file <- as.character(peak$path_to_cdf_csv)
                                        if(file.exists(cdf_file)){
                                            if (is.null(loaded_file) || !identical(loaded_file, cdf_file)) {
                                                df <- data.table::fread(cdf_file)
                                                loaded_file <- cdf_file
                                            }
                                            # Filter rows that fall within the peak’s retention time window
                                            df_subset <- df[df$rt >= peak$peak_start & df$rt <= peak$peak_end, ]
                                            if(nrow(df_subset) > 0){
                                                # Round m/z values and average intensities by m/z
                                                df_agg <- df_subset %>%
                                                    dplyr::mutate(mz = round(mz, 0)) %>%
                                                    dplyr::group_by(mz) %>%
                                                    dplyr::summarize(intensity = sum(intensity, na.rm = TRUE))
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

                                # Guard: single-detector data (e.g. GC-TCD) has no
                                # mass spectrum to classify, so stop here cleanly.
                                    if (all(ms_data$mz == 0)) {
                                        message("Mass spectral classification skipped - single-detector data (e.g. TCD) has no spectrum to classify.")
                                        return()
                                    }

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
                                    pivot_wider(names_from = mz, values_from = intensity, values_fill = 0) -> ms_wide_raw
                                ms_wide_pre <- ms_wide_raw[ms_wide_raw$peak_unique_id != 0, ]
                                ms_wide_norm <- as_tibble(t(apply(ms_wide_pre[,-1], 1, function(x) (x / max(x)) * 100))) # normalize to 100
                                ms_wide <<- ms_wide_norm

                                message("ms_wide made and saved for inspection and plotting.")

                                # Classify peaks with the Python classifier (XGBoost, or sklearn
                                # HistGradientBoosting fallback). The trained model pkl must sit beside
                                # ms_classifier.py, which resolves it relative to its own __file__.
                                if (is.null(ms_classifier_path) || !file.exists(ms_classifier_path)) {
                                    message(paste0("Classification skipped - ms_classifier.py not found at: ",
                                                   ifelse(is.null(ms_classifier_path), "NULL", ms_classifier_path)))
                                    return()
                                }

                                # Export ms_wide (with peak_unique_id) for the Python classifier
                                ms_wide_export <- cbind(peak_unique_id = ms_wide_pre$peak_unique_id, ms_wide_norm)
                                write.csv(ms_wide_export, "ms_wide_for_prediction.csv", row.names = FALSE)

                                exit_code <- system(paste0("python3 ", shQuote(ms_classifier_path),
                                                           " predict ms_wide_for_prediction.csv predictions.csv"))
                                if (exit_code != 0 || !file.exists("predictions.csv")) {
                                    message("ERROR: Python classifier failed. Check that ms_classifier.py is trained and its model pkl is present.")
                                    return()
                                }

                                # Read predictions back (columns: predicted_compound, confidence)
                                pred_df <- read.csv("predictions.csv", stringsAsFactors = FALSE)
                                predictions <<- pred_df$predicted_compound

                                message("predictions made:")

                                peak_table$peak_ID <- paste0("(predicted) ", predictions,
                                                             " [", round(pred_df$confidence * 100, 0), "%]")
                                writeMonolist(peak_table, "peaks_monolist.csv")

                                message("Mass spectral classification completed.")
                            }

                    })

            }

        ## Call the app
            
            if ( jupyter == TRUE) {

                slot <- pick_shiny_port()
                if (is.null(slot)) {
                    cat("All available ports are in use. Please try again later.")
                } else {
                    cat(paste0("Connect at: ", slot$url))
                    runApp(shinyApp(ui = ui, server = server), host = "127.0.0.1", port = slot$port)
                }
            
            }

            if (jupyter == FALSE) {
                shinyApp(ui = ui, server = server)
            }
        

    }

