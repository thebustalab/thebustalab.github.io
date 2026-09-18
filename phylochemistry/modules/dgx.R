#################################################
## DGX ENDPOINTS (research students, JupyterHub) ##
#################################################

# Split out of modules/language_model_analysis.R on 2026-09-08.
#
# WHY THIS IS ITS OWN MODULE: language_model_analysis.R is sourced ONLY on the lab
# build (bustalab = TRUE). A research student sourcing phylochemistry.R therefore
# never got dgxGenerate/runModelGrid/dgxHealth, and dgx_example_R.ipynb failed on
# its first cell with `could not find function "dgxHealth"`. This file is sourced
# unconditionally from phylochemistry.R, and it deliberately needs only packages
# that are already in the CORE scope (httr, jsonlite, readr, dplyr, tibble) so it
# adds no install burden for a chapters-1-11 student.

message("Loading dgx module...")

        # These reach the lab's dgx through SSH tunnels that already exist on host1:
        #   9000 = ESM2 protein embeddings, 9001 = ESMC, 9002 = Gemma text generation.
        # They are NOT reachable from a laptop — only from inside JupyterHub.
        #
        # THE ONE THING TO UNDERSTAND: the dgx holds ONE model at a time. Asking for a
        # different model (or a different quantisation) evicts the old one and loads the new,
        # which takes about 3 minutes. Answering a prompt then takes under a second. So never
        # loop model-inside-prompt; use runModelGrid(), which groups by model.

        dgxGenerate <- function(
          prompt,
          model = NULL,          # e.g. "google/gemma-4-12B-it"; NULL = server default
          quant = NULL,          # 4 = nf4, 8 = int8, 16 = bf16 (full precision)
          max_tokens = 256,
          system = NULL,
          server_url = Sys.getenv("DGX_GENERATE_URL", "http://127.0.0.1:9002"),
          timeout_sec = 900,     # a cold model load is ~180s; 60 would fail on every first call
          # If the dgx is too full to load the model right now (e.g. a big model while a video
          # render is running), it answers 503 "memory_admission". That means "not yet", so
          # WAIT and retry, for up to this long. 0 = fail straight away. Added 2026-09-18.
          max_wait_sec = as.numeric(Sys.getenv("DGX_MEMORY_WAIT_SEC", 4 * 3600)),
          on_wait = NULL         # internal: runModelGrid uses it to keep its lock fresh
        ) {
          suppressPackageStartupMessages({ library(httr); library(jsonlite) })
          payload <- list(prompt = prompt, max_tokens = max_tokens)
          if (!is.null(model))  payload$model  <- model
          if (!is.null(quant))  payload$quant  <- as.integer(quant)
          if (!is.null(system)) payload$system <- system
          waited <- 0
          repeat {
            resp <- httr::POST(
              url = paste0(sub("/+$", "", server_url), "/generate"),
              httr::add_headers(`Content-Type` = "application/json",
                                `User-Agent` = "Mozilla/5.0 (phylochemistry dgxGenerate)"),
              body = jsonlite::toJSON(payload, auto_unbox = TRUE),
              encode = "raw", httr::timeout(timeout_sec)
            )
            txt <- httr::content(resp, "text", encoding = "UTF-8")
            if (httr::status_code(resp) != 503 || !grepl("memory_admission", txt, fixed = TRUE)) break
            wait <- suppressWarnings(as.numeric(httr::headers(resp)[["retry-after"]]))
            if (length(wait) != 1 || is.na(wait)) wait <- 300
            why <- tryCatch(jsonlite::fromJSON(txt)$detail, error = function(e) txt)
            why <- sub("memory_admission: ", "", why, fixed = TRUE)
            if (waited + wait > max_wait_sec)
              stop("Gave up after waiting ", round(waited / 60), " min for memory: ", why, call. = FALSE)
            if (waited == 0) {
              message("The dgx is too full to load ", model, " @ ", quant, "-bit right now:\n  ", why,
                      "\nWaiting and retrying every ", round(wait / 60), " min (up to ",
                      round(max_wait_sec / 3600), " h). Nothing is lost.")
            } else {
              message("  still waiting for memory (", round(waited / 60), " min so far)...")
            }
            if (is.function(on_wait)) on_wait()
            Sys.sleep(wait)
            waited <- waited + wait
            if (is.function(on_wait)) on_wait()
          }
          if (httr::status_code(resp) >= 300) {
            stop("dgx said ", httr::status_code(resp), ": ",
                 substr(httr::content(resp, "text", encoding = "UTF-8"), 1, 300))
          }
          jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))
        }

        # ---- SHARING: the dgx has ONE model slot for the whole lab ------------------------
        # Two grids at once do not queue politely — they evict each other's model on every
        # request, so both pay a ~3-minute reload per prompt and neither finishes. On
        # 2026-09-16 that churn exhausted the dgx's memory, the kernel killed the model
        # server, and the 45,358 prompts queued behind it failed in under a minute into a
        # 114 MB CSV with nothing usable in it. These guards exist so that cannot repeat.
        # They mirror dgx_client.py exactly — keep the two sides in step.

        dgxGridLock      <- function() Sys.getenv("DGX_GRID_LOCK",
          "/project_data/shared/general_lab_resources/remote_lm_server/.grid_lock.json")
        dgxLockStaleSec  <- function() as.numeric(Sys.getenv("DGX_GRID_LOCK_STALE_SEC", "1800"))
        dgxBusyIdleSec   <- function() as.numeric(Sys.getenv("DGX_GRID_BUSY_IDLE_SEC", "300"))
        dgxMaxGridCells  <- function() as.numeric(Sys.getenv("DGX_MAX_GRID_CELLS", "5000"))

        .dgxMe <- function() list(host = as.character(Sys.info()[["nodename"]]),
                                  pid  = as.integer(Sys.getpid()))

        .dgxReadLock <- function() {
          f <- dgxGridLock()
          if (!file.exists(f)) return(NULL)
          tryCatch(jsonlite::fromJSON(f), error = function(e) NULL)
        }

        .dgxWriteLock <- function(info) {
          f <- dgxGridLock(); tmp <- paste0(f, ".tmp", Sys.getpid())
          ok <- tryCatch({
            writeLines(jsonlite::toJSON(info, auto_unbox = TRUE), tmp)
            file.rename(tmp, f)
            try(Sys.chmod(f, "666"), silent = TRUE)
            TRUE
          }, error = function(e) FALSE, warning = function(w) FALSE)
          if (!ok) {
            try(unlink(tmp), silent = TRUE)
            message("[warn] could not write the grid lock at ", f,
                    ". Carrying on, but nobody else can see that you are running.")
          }
          invisible(ok)
        }

        # The live lock, if somebody OTHER than us holds it. NULL if free, stale, or ours.
        .dgxLockHolder <- function() {
          lock <- .dgxReadLock(); if (is.null(lock)) return(NULL)
          me <- .dgxMe()
          if (identical(as.character(lock$host), me$host) &&
              identical(as.integer(lock$pid), me$pid)) return(NULL)
          hb <- suppressWarnings(as.numeric(lock$heartbeat))
          if (is.na(hb) || (as.numeric(Sys.time()) - hb) > dgxLockStaleSec()) return(NULL)
          lock
        }

        .dgxReleaseLock <- function() {
          lock <- .dgxReadLock(); if (is.null(lock)) return(invisible(NULL))
          me <- .dgxMe()
          if (identical(as.character(lock$host), me$host) &&
              identical(as.integer(lock$pid), me$pid)) try(unlink(dgxGridLock()), silent = TRUE)
          invisible(NULL)
        }

        # Seconds since the generate server last answered anyone; NA if we cannot tell.
        .dgxServerIdle <- function(
          server_url = Sys.getenv("DGX_GENERATE_URL", "http://127.0.0.1:9002")
        ) {
          tryCatch({
            resp <- httr::GET(paste0(sub("/+$", "", server_url), "/health"), httr::timeout(10))
            if (httr::status_code(resp) != 200) return(NA_real_)
            idle <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))$idle_sec
            if (is.null(idle)) NA_real_ else as.numeric(idle)
          }, error = function(e) NA_real_)
        }

        # Is it safe to start a grid right now? Run this before a long job, and again if
        # runModelGrid() refuses to start.
        dgxGridStatus <- function() {
          holder <- .dgxLockHolder()
          if (!is.null(holder)) {
            mins <- (as.numeric(Sys.time()) - as.numeric(holder$started)) / 60
            message("BUSY — ", holder$user, " on ", holder$host, " started a grid ",
                    round(mins), " min ago (", holder$total, " prompts -> ", holder$out_csv, ").")
          } else message("Grid lock: free.")
          idle <- .dgxServerIdle()
          if (is.na(idle)) message("Server: could not reach the generate endpoint.")
          else if (idle < dgxBusyIdleSec())
            message("Server: in use — it answered someone ", round(idle), " s ago.")
          else message("Server: quiet — idle for ", round(idle / 60), " min.")
          invisible(NULL)
        }

        # Refuse to start a grid that would collide with someone else or monopolise the box.
        .dgxPreflight <- function(total, force, verbose) {
          problems <- character(0)
          if (total > dgxMaxGridCells()) problems <- c(problems, paste0(
            "This grid is ", format(total, big.mark = ","), " generations. That is days of GPU\n",
            "    time for the whole lab. Pilot on 20-50 prompts across the full grid first,\n",
            "    decide which (model, quant) actually works, then run everything through that\n",
            "    ONE cell."))
          holder <- .dgxLockHolder()
          if (!is.null(holder)) {
            mins <- (as.numeric(Sys.time()) - as.numeric(holder$started)) / 60
            problems <- c(problems, paste0(
              "Another grid is running: ", holder$user, " on ", holder$host, ", started ",
              round(mins), " min ago\n    (", holder$total, " prompts -> ", holder$out_csv,
              "). Starting now would make both runs evict\n",
              "    each other's model on every prompt."))
          } else {
            idle <- .dgxServerIdle()
            if (is.na(idle)) problems <- c(problems, paste0(
              "Cannot reach the generate server. Run dgxHealth() first — starting a grid\n",
              "    against a dead server just fills the CSV with errors."))
            else if (idle < dgxBusyIdleSec()) problems <- c(problems, paste0(
              "The server answered somebody ", round(idle), " s ago, so it is probably in use.\n",
              "    Wait for them to finish, or check with dgxGridStatus()."))
          }
          if (!length(problems)) return(invisible(NULL))
          bar <- strrep("-", 90)
          body <- paste0("  * ", problems, collapse = "\n\n")
          if (force) {
            if (verbose) message("\n", bar, "\nrunModelGrid WARNINGS (overridden by force = TRUE):\n\n",
                                 body, "\n", bar)
            return(invisible(NULL))
          }
          stop("\n", bar, "\nrunModelGrid refused to start:\n\n", body,
               "\n\nIf you are sure, call runModelGrid(..., force = TRUE).\n", bar, call. = FALSE)
        }

        # Run every prompt against every (model, quantisation) and write a tidy CSV.
        # Batched BY CELL so the ~3-minute load is paid once per cell, not once per prompt.
        # RESUMABLE: appends after every prompt and skips rows that SUCCEEDED, so if the
        # kernel dies four hours in, re-running the same cell continues where it stopped.
        # Rows that FAILED are retried — a failure is not a result. A retried prompt appends a
        # second row; the later row is the one that counts.
        # IT WILL REFUSE TO START if another grid holds the lock, if the server is answering
        # somebody else, or if the grid is larger than DGX_MAX_GRID_CELLS. Call dgxGridStatus()
        # to see why, and pass force = TRUE if you are sure.
        runModelGrid <- function(
          prompts,
          models = c("google/gemma-4-E2B-it", "google/gemma-4-E4B-it",
                     "google/gemma-4-12B-it", "google/gemma-4-26B-A4B-it",
                     "google/gemma-4-31B-it"),
          quants = c(4),
          out_csv = "grid_results.csv",
          max_tokens = 16,
          system = NULL,
          labels = NULL,
          verbose = TRUE,
          force = FALSE,
          max_consecutive_errors = as.numeric(Sys.getenv("DGX_MAX_CONSECUTIVE_ERRORS", "10"))
        ) {
          suppressPackageStartupMessages({ library(readr); library(dplyr); library(tibble) })
          if (is.null(labels)) labels <- as.character(seq_along(prompts))
          if (length(labels) != length(prompts)) stop("labels must match prompts in length")

          total <- length(models) * length(quants) * length(prompts)
          .dgxPreflight(total, force, verbose)

          # A row counts as done only if it has no error. The LAST row for a key wins, so a
          # successful retry supersedes the failure that came before it.
          done_keys <- character(0); retrying <- 0
          if (file.exists(out_csv)) {
            prev <- suppressMessages(readr::read_csv(out_csv, show_col_types = FALSE))
            if (nrow(prev) > 0) {
              prev <- prev %>%
                mutate(.key = paste(model, quant, as.character(prompt_id)),
                       .ok  = is.na(error) | !nzchar(trimws(as.character(error))),
                       .row = row_number()) %>%
                group_by(.key) %>% filter(.row == max(.row)) %>% ungroup()
              done_keys <- prev$.key[prev$.ok]
              retrying  <- sum(!prev$.ok)
              if (verbose) {
                message("Resuming: ", length(done_keys), " successful results already in ", out_csv)
                if (retrying > 0) message("          ", retrying, " failed rows will be RETRIED.")
              }
            }
          }

          n <- length(done_keys); consecutive_errors <- 0
          me <- .dgxMe()
          lock <- list(host = me$host, pid = me$pid,
                       user = as.character(Sys.info()[["user"]]),
                       out_csv = normalizePath(out_csv, mustWork = FALSE),
                       total = total, started = as.numeric(Sys.time()),
                       heartbeat = as.numeric(Sys.time()))
          .dgxWriteLock(lock)
          on.exit(.dgxReleaseLock(), add = TRUE)

          for (m in models) for (q in quants) {
            keep <- !(paste(m, q, labels) %in% done_keys)
            if (!any(keep)) { if (verbose) message("[skip] ", m, " @ ", q, "-bit"); next }
            if (verbose) message("\n[cell] ", m, " @ ", q, "-bit — ", sum(keep),
                                 " prompts\n       loading the model, ~3 minutes...")
            for (i in which(keep)) {
              failed <- FALSE
              row <- tryCatch({
                r <- dgxGenerate(prompts[i], model = m, quant = q,
                                 max_tokens = max_tokens, system = system,
                                 on_wait = function() {   # keep the lock fresh while waiting
                                   lock$heartbeat <- as.numeric(Sys.time()); .dgxWriteLock(lock)
                                 })
                tibble(model = m, quant = q, prompt_id = labels[i], prompt = prompts[i],
                       completion = r$completion, duration_sec = r$duration_sec,
                       error = NA_character_, timestamp = format(Sys.time()))
              }, error = function(e) {    # one bad prompt must not kill a long run
                failed <<- TRUE
                tibble(model = m, quant = q, prompt_id = labels[i], prompt = prompts[i],
                       completion = NA_character_, duration_sec = NA_real_,
                       error = substr(conditionMessage(e), 1, 300), timestamp = format(Sys.time()))
              })
              readr::write_csv(row, out_csv, append = file.exists(out_csv))
              n <- n + 1
              lock$heartbeat <- as.numeric(Sys.time()); .dgxWriteLock(lock)
              consecutive_errors <- if (failed) consecutive_errors + 1 else 0

              if (consecutive_errors >= max_consecutive_errors) {
                message("\nSTOPPED after ", consecutive_errors, " failures in a row. ",
                        "The last one was:\n    ", row$error[1],
                        "\nNothing is being generated, so there is no point continuing — the ",
                        "rest of the grid\nwould just fill ", out_csv, " with the same error.",
                        "\nCheck the server with dgxHealth(), then re-run this same cell: the ",
                        "failed rows\nwill be retried automatically.")
                stop(consecutive_errors, " consecutive failures — run stopped. Last error: ",
                     row$error[1], call. = FALSE)
              }
              if (verbose && (n %% 10 == 0)) message("       ", n, "/", total, " overall")
            }
          }
          if (verbose) message("\nDone. Results in ", out_csv)
          invisible(out_csv)
        }

        dgxHealth <- function() {
          suppressPackageStartupMessages(library(httr))
          urls <- c(
            generate = Sys.getenv("DGX_GENERATE_URL", "http://127.0.0.1:9002"),
            esm2     = Sys.getenv("DGX_ESM2_URL",     "http://127.0.0.1:9000"),
            esmc     = Sys.getenv("DGX_ESMC_URL",     "http://127.0.0.1:9001")
          )
          # Loop over NAMES, not values: `for (nm in urls)` binds nm to the URL and
          # loses the label, so every line printed the address instead of the service.
          for (nm in names(urls)) {
            ok <- tryCatch(httr::status_code(httr::GET(paste0(sub("/+$", "", urls[[nm]]), "/health"),
                                                       httr::timeout(10))) == 200,
                           error = function(e) FALSE)
            message(nm, " (", urls[[nm]], "): ",
                    if (ok) "ok" else "UNREACHABLE (are you on JupyterHub?)")
          }
          invisible(NULL)
        }

message("Done with dgx loading!")
