# decode_codes.R — turn submitted escape-room codes back into answer paths.
#
# This mirrors ../shared/codec.js byte-for-byte. If you change the JS
# arithmetic, change it here too (and re-run the round-trip test at the bottom).
#
# Typical use after downloading a Canvas assignment export:
#   source("decode_codes.R")
#   roster <- readr::read_csv("canvas_export.csv")   # needs x500 + code columns
#   scored <- grade_submissions(roster, ALASKA_KEY,
#                               id_col = "x500", code_col = "code")
#   readr::write_csv(scored, "graded.csv")

SECRET <- "chem5725-noatak-2026"   # must match escape-engine.js
TWO32  <- 2^32
CROCKFORD <- strsplit("0123456789ABCDEFGHJKMNPQRSTVWXYZ", "")[[1]]

# ---- primitives (mirror codec.js) ----

hash32 <- function(str) {
  h <- 2166136261
  for (ch in utf8ToInt(str)) h <- (h * 31 + ch) %% TWO32
  h
}

hash_bytes <- function(bytes) {
  h <- 2166136261
  for (b in bytes) h <- (h * 31 + b) %% TWO32
  h
}

keystream <- function(secret, student_id, n) {
  seed <- hash32(paste0(secret, "|", tolower(trimws(student_id))))
  state <- seed
  out <- integer(n)
  for (i in seq_len(n)) {
    state <- (1664525 * state + 1013904223) %% TWO32
    out[i] <- floor(state / 65536) %% 256
  }
  out
}

crockford_encode <- function(bytes) {
  bits <- 0; value <- 0; out <- ""
  for (b in bytes) {
    value <- value * 256 + b
    bits <- bits + 8
    while (bits >= 5) {
      bits <- bits - 5
      idx <- floor(value / 2^bits) %% 32
      out <- paste0(out, CROCKFORD[idx + 1])
      value <- value %% 2^bits            # drop the consumed high bits
    }
  }
  if (bits > 0) {
    idx <- (value * 2^(5 - bits)) %% 32
    out <- paste0(out, CROCKFORD[idx + 1])
  }
  out
}

crockford_decode <- function(str) {
  str <- gsub("-", "", toupper(trimws(str)))
  chars <- strsplit(str, "")[[1]]
  vals <- match(chars, CROCKFORD) - 1
  if (any(is.na(vals))) stop("Invalid character in code")
  bits <- 0; value <- 0; out <- integer(0)
  for (v in vals) {
    value <- value * 32 + v
    bits <- bits + 5
    if (bits >= 8) {
      bits <- bits - 8
      out <- c(out, floor(value / 2^bits) %% 256)
      value <- value %% 2^bits            # drop the consumed high bits
    }
  }
  out
}

bitxor_byte <- function(a, b) bitwXor(as.integer(a), as.integer(b))

# ---- encode (for testing / answer keys) ----

# Arg order mirrors codec.js encode(): version first, then scenario_id.
encode_code <- function(version, scenario_id, steps, student_id, secret = SECRET) {
  # steps: list of list(answer=, attempts=)
  header <- bitwOr(bitwShiftL(bitwAnd(version, 15L), 4L), bitwAnd(scenario_id, 15L))
  payload <- header
  for (s in steps) {
    ans <- bitwAnd(as.integer(s$answer), 31L)
    # 0 = "not attempted" (a skipped hub-and-spoke node under an N-of-M gate);
    # resolved nodes carry attempts >= 1. Mirrors codec.js.
    att <- bitwAnd(as.integer(min(max(s$attempts, 0), 7)), 7L)
    payload <- c(payload, bitwOr(bitwShiftL(att, 5L), ans))
  }
  chk <- hash_bytes(payload) %% 256
  full <- c(payload, chk)
  ks <- keystream(secret, student_id, length(full))
  scrambled <- mapply(bitxor_byte, full, ks)
  raw <- crockford_encode(scrambled)
  gsub("(.{4})(?=.)", "\\1-", raw, perl = TRUE)
}

# ---- decode ----

decode_code <- function(code, student_id, secret = SECRET) {
  scrambled <- crockford_decode(code)
  ks <- keystream(secret, student_id, length(scrambled))
  full <- mapply(bitxor_byte, scrambled, ks)
  n <- length(full)
  payload <- full[1:(n - 1)]
  chk <- full[n]
  ok <- (hash_bytes(payload) %% 256) == chk
  header <- payload[1]
  version <- bitwShiftR(header, 4L)
  scenario_id <- bitwAnd(header, 15L)
  step_bytes <- payload[-1]
  answers  <- sapply(step_bytes, function(b) bitwAnd(b, 31L))
  attempts <- sapply(step_bytes, function(b) bitwShiftR(b, 5L))
  list(valid = ok, version = version, scenario_id = scenario_id,
       answers = as.integer(answers), attempts = as.integer(attempts))
}

# ---- grading ----

# A scenario key: correct answer index per step, and a scoring function.
ALASKA_KEY <- list(
  scenario_id = 1,
  correct = c(18, 3, 1),          # 0-based indices, in step order
  # points per step: full if right, scaled down by attempts, 0 if never right
  score_step = function(correct, answer, attempts) {
    if (answer != correct) return(0)
    if (attempts <= 1) return(10)
    if (attempts == 2) return(7)
    if (attempts == 3) return(5)
    3
  }
)

grade_one <- function(code, student_id, key, secret = SECRET) {
  d <- tryCatch(decode_code(code, student_id, secret), error = function(e) NULL)
  if (is.null(d) || !d$valid || d$scenario_id != key$scenario_id) {
    return(list(valid = FALSE, points = NA_integer_, detail = "invalid/mismatched code"))
  }
  n <- length(key$correct)
  pts <- 0
  detail <- character(0)
  for (i in seq_len(n)) {
    a  <- if (i <= length(d$answers)) d$answers[i] else -1
    at <- if (i <= length(d$attempts)) d$attempts[i] else 0
    p  <- key$score_step(key$correct[i], a, at)
    pts <- pts + p
    detail <- c(detail, sprintf("Q%d: ans=%d att=%d -> %dpt", i, a, at, p))
  }
  list(valid = TRUE, points = pts, detail = paste(detail, collapse = "; "))
}

# ---- graph-mode grading (hub-and-spoke) ----
#
# A graph scenario's code carries, in canonical NODE order:
#   - one byte per spoke (answer + attempts; attempts 0 = the student skipped it
#     under the N-of-M gate), then
#   - one trailing boss byte whose answer bit is 1 if the student REACHED the
#     boss and produced a figure (the figure itself is graded by hand — the code
#     only records that they got there, for the integrity/watermark link).
#
# A graph key: `n_spokes`, the per-spoke `correct` indices (node order), and a
# per-spoke `score_step`. The boss is not auto-scored here.
grade_graph <- function(code, student_id, key, secret = SECRET) {
  d <- tryCatch(decode_code(code, student_id, secret), error = function(e) NULL)
  if (is.null(d) || !d$valid || d$scenario_id != key$scenario_id) {
    return(list(valid = FALSE, points = NA_integer_, boss_reached = NA,
                detail = "invalid/mismatched code"))
  }
  n <- key$n_spokes
  pts <- 0
  detail <- character(0)
  for (i in seq_len(n)) {
    a  <- if (i <= length(d$answers)) d$answers[i] else -1L
    at <- if (i <= length(d$attempts)) d$attempts[i] else 0L
    p  <- if (at == 0) 0 else key$score_step(key$correct[i], a, at)  # 0 attempts = skipped
    pts <- pts + p
    detail <- c(detail, sprintf("S%d: ans=%d att=%d -> %dpt", i, a, at, p))
  }
  boss_reached <- length(d$answers) > n && d$answers[n + 1] == 1
  detail <- c(detail, sprintf("boss_reached=%s (figure graded by hand)", boss_reached))
  list(valid = TRUE, points = pts, boss_reached = boss_reached,
       detail = paste(detail, collapse = "; "))
}

# Demo/trial hub-and-spoke room (scenario id 2): 3 spokes reusing the Alaska
# questions + a boss figure task. Spoke `correct` indices match demo_hub/scenario.js.
DEMO_KEY <- list(
  scenario_id = 2,
  n_spokes = 3,
  correct = c(18, 3, 1),
  score_step = function(correct, answer, attempts) {
    if (answer != correct) return(0)
    if (attempts <= 1) return(10)
    if (attempts == 2) return(7)
    if (attempts == 3) return(5)
    3
  }
)

# Journey/chain room (scenario id 3): 2 case-rooms + a boss figure. Same grading
# shape as a graph key (spokes in room order, trailing boss byte). Spoke `correct`
# indices match datavis1/scenario.js (Alaska 18, Hawai‘i 2).
DATAVIS1_KEY <- list(
  scenario_id = 3,
  n_spokes = 2,
  correct = c(18, 2),
  score_step = function(correct, answer, attempts) {
    if (answer != correct) return(0)
    if (attempts <= 1) return(10)
    if (attempts == 2) return(7)
    if (attempts == 3) return(5)
    3
  }
)

# Explore/pano room (scenario id 4): 3 rooms, one MC each + a reserved trailing
# byte (no boss figure in the prototype). Spoke `correct` = Alaska 18/3/1.
ALASKA_STATION_KEY <- list(
  scenario_id = 4,
  n_spokes = 3,
  correct = c(18, 3, 1),
  score_step = function(correct, answer, attempts) {
    if (answer != correct) return(0)
    if (attempts <= 1) return(10)
    if (attempts == 2) return(7)
    if (attempts == 3) return(5)
    3
  }
)

# Panorama explore room (scenario id 5): 2 rooms, one MC each (Alaska 18 / 1).
ALASKA_PANO_KEY <- list(
  scenario_id = 5,
  n_spokes = 2,
  correct = c(18, 1),
  score_step = function(correct, answer, attempts) {
    if (answer != correct) return(0)
    if (attempts <= 1) return(10)
    if (attempts == 2) return(7)
    if (attempts == 3) return(5)
    3
  }
)

# Vectorised over a data frame of submissions.
grade_submissions <- function(df, key, id_col = "x500", code_col = "code",
                              secret = SECRET) {
  res <- lapply(seq_len(nrow(df)), function(i) {
    g <- grade_one(df[[code_col]][i], df[[id_col]][i], key, secret)
    data.frame(
      id = df[[id_col]][i], code = df[[code_col]][i],
      valid = g$valid, points = g$points, detail = g$detail,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

# ---- self-test: round-trip a known path ----
if (identical(environment(), globalenv()) && sys.nframe() == 0) {
  steps <- list(list(answer = 18, attempts = 1),
                list(answer = 3,  attempts = 2),
                list(answer = 1,  attempts = 1))
  code <- encode_code(scenario_id = 1, version = 1, steps = steps,
                      student_id = "bust0037")
  cat("Encoded code:", code, "\n")
  dec <- decode_code(code, "bust0037")
  cat("Valid:", dec$valid, "| answers:", paste(dec$answers, collapse = ","),
      "| attempts:", paste(dec$attempts, collapse = ","), "\n")
  wrong <- decode_code(code, "someone_else")   # wrong id -> should be invalid
  cat("Wrong-id decode valid (should be FALSE):", wrong$valid, "\n")
  g <- grade_one(code, "bust0037", ALASKA_KEY)
  cat("Points:", g$points, "|", g$detail, "\n")

  # Regression: long (10-step) code round-trip.
  # FAILURE MODE this guards against: the base32 accumulator must be reduced
  # (value %% 2^bits) after each byte is emitted. Without it, codes longer than
  # ~6 bytes let `value` grow past 2^53, at which point the browser's JS doubles
  # lose precision and silently disagree with this R decoder — long/large
  # scenarios would produce codes that decode to the wrong answers. This
  # asserts a 10-step path survives the round-trip exactly.
  long_steps <- lapply(0:9, function(i) list(answer = (i * 3) %% 20,
                                             attempts = (i %% 4) + 1))
  long_code <- encode_code(version = 1, scenario_id = 2, steps = long_steps,
                           student_id = "test_student")
  ld <- decode_code(long_code, "test_student")
  exp_ans <- vapply(long_steps, function(s) as.integer(s$answer), integer(1))
  exp_att <- vapply(long_steps, function(s) as.integer(s$attempts), integer(1))
  ok <- ld$valid && identical(ld$answers, exp_ans) && identical(ld$attempts, exp_att)
  cat("Long-code (10-step) round-trip OK (should be TRUE):", ok, "\n")
  if (!ok) stop("REGRESSION: long-code round-trip failed — check base32 accumulator reduction")

  # Regression: graph-mode round-trip — a SKIPPED spoke (attempts 0) + boss byte.
  # FAILURE MODE this guards against: the attempts=0 sentinel for an un-attempted
  # hub-and-spoke node must survive encode/decode. If the old `max(attempts, 1)`
  # floor crept back into either codec.js or encode_code(), a skipped node would
  # decode as attempts=1 and be mis-scored as an attempt. Also asserts the boss
  # byte (reached=1) round-trips.
  gsteps <- list(list(answer = 18, attempts = 1),   # spoke 1: solved first try
                 list(answer = 0,  attempts = 0),    # spoke 2: SKIPPED (not attempted)
                 list(answer = 1,  attempts = 2),    # spoke 3: solved, 2 tries
                 list(answer = 1,  attempts = 0))    # boss byte: figure produced
  gcode <- encode_code(version = 1, scenario_id = 2, steps = gsteps,
                       student_id = "grid_test")
  gd <- decode_code(gcode, "grid_test")
  gok <- gd$valid &&
    identical(gd$answers,  c(18L, 0L, 1L, 1L)) &&
    identical(gd$attempts, c(1L, 0L, 2L, 0L))
  cat("Graph-mode round-trip OK (should be TRUE):", gok, "\n")
  if (!gok) stop("REGRESSION: graph-mode round-trip failed — check attempts=0 sentinel")
  gg <- grade_graph(gcode, "grid_test", DEMO_KEY)
  cat("Graph grade — points:", gg$points, "| boss_reached:", gg$boss_reached, "\n")
  cat("  detail:", gg$detail, "\n")
  # spoke1 correct(18) 1 try = 10; spoke2 skipped = 0; spoke3 correct(1) 2 tries = 7 -> 17
  if (!isTRUE(gg$points == 17 && gg$boss_reached)) {
    stop("REGRESSION: graph grade wrong — expected 17 pts + boss_reached=TRUE")
  }

  # Regression: journey/chain round-trip — 2 rooms solved in order + boss byte.
  jsteps <- list(list(answer = 18, attempts = 1),   # room 1 (Alaska) solved
                 list(answer = 2,  attempts = 2),    # room 2 (Hawai‘i) solved, 2 tries
                 list(answer = 1,  attempts = 0))    # boss byte: figure produced
  jcode <- encode_code(version = 1, scenario_id = 3, steps = jsteps,
                       student_id = "journey_test")
  jg <- grade_graph(jcode, "journey_test", DATAVIS1_KEY)
  cat("Journey grade — points:", jg$points, "| boss_reached:", jg$boss_reached, "\n")
  # room1 10 + room2 (2 tries) 7 = 17
  if (!isTRUE(jg$points == 17 && jg$boss_reached)) {
    stop("REGRESSION: journey grade wrong — expected 17 pts + boss_reached=TRUE")
  }
}
