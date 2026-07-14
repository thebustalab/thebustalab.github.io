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
    att <- bitwAnd(as.integer(min(max(s$attempts, 1), 7)), 7L)
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
}
