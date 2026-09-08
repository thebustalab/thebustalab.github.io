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
          timeout_sec = 900      # a cold model load is ~180s; 60 would fail on every first call
        ) {
          suppressPackageStartupMessages({ library(httr); library(jsonlite) })
          payload <- list(prompt = prompt, max_tokens = max_tokens)
          if (!is.null(model))  payload$model  <- model
          if (!is.null(quant))  payload$quant  <- as.integer(quant)
          if (!is.null(system)) payload$system <- system
          resp <- httr::POST(
            url = paste0(sub("/+$", "", server_url), "/generate"),
            httr::add_headers(`Content-Type` = "application/json",
                              `User-Agent` = "Mozilla/5.0 (phylochemistry dgxGenerate)"),
            body = jsonlite::toJSON(payload, auto_unbox = TRUE),
            encode = "raw", httr::timeout(timeout_sec)
          )
          if (httr::status_code(resp) >= 300) {
            stop("dgx said ", httr::status_code(resp), ": ",
                 substr(httr::content(resp, "text", encoding = "UTF-8"), 1, 300))
          }
          jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))
        }

        # Run every prompt against every (model, quantisation) and write a tidy CSV.
        # Batched BY CELL so the ~3-minute load is paid once per cell, not once per prompt.
        # RESUMABLE: appends after every prompt and skips rows already present, so if the
        # kernel dies four hours in, re-running the same cell continues where it stopped.
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
          verbose = TRUE
        ) {
          suppressPackageStartupMessages({ library(readr); library(dplyr); library(tibble) })
          if (is.null(labels)) labels <- as.character(seq_along(prompts))
          if (length(labels) != length(prompts)) stop("labels must match prompts in length")

          done <- tibble(model = character(), quant = integer(), prompt_id = character())
          if (file.exists(out_csv)) {
            prev <- suppressMessages(readr::read_csv(out_csv, show_col_types = FALSE))
            if (nrow(prev) > 0) {
              done <- prev %>% select(model, quant, prompt_id) %>%
                        mutate(prompt_id = as.character(prompt_id))
              if (verbose) message("Resuming: ", nrow(done), " results already in ", out_csv)
            }
          }

          total <- length(models) * length(quants) * length(prompts); n <- nrow(done)
          for (m in models) for (q in quants) {
            keep <- !(paste(m, q, labels) %in% paste(done$model, done$quant, done$prompt_id))
            if (!any(keep)) { if (verbose) message("[skip] ", m, " @ ", q, "-bit"); next }
            if (verbose) message("\n[cell] ", m, " @ ", q, "-bit — ", sum(keep),
                                 " prompts\n       loading the model, ~3 minutes...")
            for (i in which(keep)) {
              row <- tryCatch({
                r <- dgxGenerate(prompts[i], model = m, quant = q,
                                 max_tokens = max_tokens, system = system)
                tibble(model = m, quant = q, prompt_id = labels[i], prompt = prompts[i],
                       completion = r$completion, duration_sec = r$duration_sec,
                       error = NA_character_, timestamp = format(Sys.time()))
              }, error = function(e)      # one bad prompt must not kill a long run
                tibble(model = m, quant = q, prompt_id = labels[i], prompt = prompts[i],
                       completion = NA_character_, duration_sec = NA_real_,
                       error = substr(conditionMessage(e), 1, 300), timestamp = format(Sys.time())))
              readr::write_csv(row, out_csv, append = file.exists(out_csv))
              n <- n + 1
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
