# tools/verify/verify-data-check-conditions.R
# Condition digest for data_check(). data_check() returns nothing a fit reads --
# it only emits messages, warnings and errors -- so its entire observable
# behaviour IS that condition set, and `/golden-check` cannot see any of it: the
# four reference models reach only a handful of its branches, and a condition
# that stops being raised moves no objective.
#
# This captures every condition data_check() emits over every bundled dataset,
# crossed with switch settings chosen to reach branches a default build does
# not, and digests it. It is the gate for decomposing data_check() -- capture on
# the clean baseline, require an identical digest afterwards.
#
# DETERMINISM is the whole point, so it is enforced rather than assumed:
#   * conditions are sorted before digesting, so dispatch order cannot matter;
#   * absolute paths, times and the session's temp directory are scrubbed from
#     every message, because several carry a file path;
#   * numbers are kept as emitted -- a message whose number changes IS a
#     behaviour change and must show up;
#   * the R and package versions are recorded alongside, so a locale- or
#     version-driven difference is attributable rather than mysterious.
# Run it twice on an unmodified tree and diff the two digests before trusting
# it; `selftest` below does exactly that in one invocation.
#
# WHAT THIS DOES NOT DO, measured rather than guessed. R/1-data_check.R has 36
# stop/warning/message sites; the cases below reach 12 distinct conditions. So
# an identical digest means the reached branches are unchanged -- it does NOT
# mean a condition elsewhere was not dropped. Demonstrated: silencing the first
# warning() in the file leaves the digest IDENTICAL, because no case reaches it.
# Before this can gate a decomposition of a 2,212-line function it needs either
# cases for the unreached sites or a run over the 375-workbook release corpus
# (183 of which carry a fleet_control sheet). Treat the present digest as a
# regression net for the common paths, not as equivalence.
#
# Output defaults into dev/ (gitignored scratch) so a digest survives between
# the two runs you need to compare.
#
# Usage:
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript tools/verify/verify-data-check-conditions.R dev/dc-before.rds
#   NOT_CRAN=true Rscript tools/verify/verify-data-check-conditions.R dev/dc-after.rds compare dev/dc-before.rds
#   NOT_CRAN=true Rscript tools/verify/verify-data-check-conditions.R selftest

args         <- commandArgs(trailingOnly = TRUE)
dir.create("dev", showWarnings = FALSE)
selftest     <- length(args) >= 1 && args[1] == "selftest"
out_path     <- if (!selftest && length(args) >= 1) args[1] else "dev/dc-conditions.rds"
do_compare   <- length(args) >= 2 && args[2] == "compare"
compare_path <- if (length(args) >= 3) args[3] else NULL

suppressMessages(pkgload::load_all(".", quiet = TRUE, compile = FALSE))

# data_check() is internal, so reach it through the namespace. Calling
# Rceattle::data_check() instead raises "not an exported object" -- which an
# earlier version of this harness captured 17 times and called a pass.
dc <- get("data_check", envir = asNamespace("Rceattle"))
stopifnot(is.function(dc))

# Scrub what varies between two runs of the same code. A message carrying a
# tempdir, an absolute path or a timestamp would otherwise make every digest
# unique and the harness worthless.
scrub <- function(x) {
  x <- gsub(tempdir(), "<TMP>", x, fixed = TRUE)
  x <- gsub(getwd(), "<WD>", x, fixed = TRUE)
  x <- gsub("/[^ '\"]*/Rceattle[^ '\"]*", "<PATH>", x)
  x <- gsub("[0-9]{4}-[0-9]{2}-[0-9]{2}[ T][0-9:]+", "<TIME>", x)
  trimws(x)
}

# Every condition data_check() raises, as (class, message) pairs. withCallingHandlers
# lets execution continue past a warning or message so one dataset yields its
# WHOLE condition set rather than stopping at the first.
capture_conditions <- function(d) {
  got <- character(0)
  record <- function(cls, msg) got <<- c(got, paste0(cls, ": ", scrub(msg)))
  res <- withCallingHandlers(
    tryCatch({ dc(d); "ok" },
             error = function(e) { record("error", conditionMessage(e)); "error" }),
    warning = function(w) {
      record("warning", conditionMessage(w)); invokeRestart("muffleWarning")
    },
    message = function(m) {
      record("message", conditionMessage(m)); invokeRestart("muffleMessage")
    })
  # Sorted, so the order data_check() happens to emit in is not part of the
  # digest. A decomposition that reorders its sections is not a behaviour
  # change; one that drops a condition is.
  list(outcome = res, conditions = sort(unique(got)))
}

# The datasets, and switch settings that reach branches a default build misses.
# Each case states what it is for, so a future reader can tell whether a new
# branch needs a new case.
bundled <- c("Atka2022", "BS2017MS", "BS2017SS", "GeorgesBank3spp", "GOA2018SS",
             "GOAatf", "GOAatf2023", "GOAcod", "GOApollock",
             "NorthernRockfish2022", "whamGrowthData")

# data_check() runs on a switch_check()'d list, which is what fit_mod() hands
# it. On raw bundled data it reaches only the unresolved-switch refusal --
# "Invalid 'Time_varying_q' specified for fleets" -- which an earlier version
# of this harness captured for nearly every dataset and mistook for a
# condition set. `raw = TRUE` keeps two cases on that path deliberately.
cases <- list()
for (nm in bundled) cases[[nm]] <- list(data = nm, mutate = NULL)
cases[["BS2017SS_raw"]]  <- list(data = "BS2017SS",  raw = TRUE)
cases[["GOA2018SS_raw"]] <- list(data = "GOA2018SS", raw = TRUE)
# Branch-reaching mutations. These deliberately produce conditions.
cases[["BS2017SS_msm"]]      <- list(data = "BS2017SS",
  mutate = function(d) { d$msmMode <- 1; d })
cases[["BS2017SS_kinzey"]]   <- list(data = "BS2017SS",
  mutate = function(d) { d$msmMode <- 5; d })        # refused: Kinzey & Punt
cases[["BS2017SS_estdyn"]]   <- list(data = "BS2017SS",
  mutate = function(d) { d$estDynamics <- rep(2, d$nspp); d })
cases[["GOA2018SS_m1est"]]   <- list(data = "GOA2018SS",
  mutate = function(d) { d$M1_model <- rep(3, d$nspp); d })
cases[["GOA2018SS_badsex"]]  <- list(data = "GOA2018SS",
  mutate = function(d) { d$comp_data$Sex[1] <- 3L; d })   # sex on a 1-sex species
cases[["BS2017SS_nodiet"]]   <- list(data = "BS2017SS",
  mutate = function(d) { d$msmMode <- 1; d$diet_data <- d$diet_data[0, ]; d })

# One builder for both capture passes. Two copies of this differing only in
# indentation is what let the switch_check omission survive a patch.
build_case <- function(spec) {
  e <- new.env()
  utils::data(list = spec$data, package = "Rceattle", envir = e)
  d <- get(spec$data, envir = e)
  if (!isTRUE(spec$raw)) d <- suppressMessages(Rceattle::switch_check(d))
  # Mutations go on AFTER switch_check, or it canonicalises them away.
  if (!is.null(spec$mutate)) d <- spec$mutate(d)
  d
}

digest <- list()
for (case in names(cases)) {
  digest[[case]] <- capture_conditions(build_case(cases[[case]]))
}

# COVERAGE, reported because it is the harness's real limitation. data_check()
# emits from many call sites and these cases reach a minority of them, so an
# identical digest proves the reached branches are unchanged and says nothing
# about the rest. Counting the sites by parsing keeps the denominator honest.
emit_sites <- function(path) {
  pd <- utils::getParseData(parse(path, keep.source = TRUE))
  calls <- pd$text[pd$token == "SYMBOL_FUNCTION_CALL"]
  sum(calls %in% c("stop", "warning", "message"))
}
n_sites <- tryCatch(emit_sites("R/1-data_check.R"), error = function(e) NA_integer_)
distinct <- unique(unlist(lapply(digest, function(x) x$conditions)))

meta <- list(
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  pkg_version = as.character(utils::packageVersion("Rceattle")),
  n_cases = length(digest),
  n_conditions = sum(vapply(digest, function(x) length(x$conditions), integer(1))),
  n_distinct = length(distinct),
  n_emit_sites = n_sites
)

report <- function(d, m) {
  cat(sprintf("cases: %d   conditions: %d   distinct: %d   (R %s, Rceattle %s)\n",
              m$n_cases, m$n_conditions, m$n_distinct, m$r_version, m$pkg_version))
  if (!is.na(m$n_emit_sites)) {
    cat(sprintf("COVERAGE: %d distinct conditions against %d stop/warning/message sites in R/1-data_check.R.\n",
                m$n_distinct, m$n_emit_sites))
    cat("  An identical digest proves the REACHED branches are unchanged, and\n")
    cat("  nothing about the rest. Not sufficient on its own to gate a\n")
    cat("  decomposition -- see the header.\n")
  }
  for (case in names(d)) {
    cat(sprintf("  %-24s %-6s %3d condition(s)\n",
                case, d[[case]]$outcome, length(d[[case]]$conditions)))
  }
}

if (selftest) {
  # The pre-flight the plan asks for: run the capture twice in one process and
  # require the two digests identical. A harness that cannot reproduce itself
  # cannot gate anything.
  second <- list()
  for (case in names(cases)) {
    second[[case]] <- capture_conditions(build_case(cases[[case]]))
  }
  report(digest, meta)
  uniq <- unique(unlist(lapply(digest, function(x) x$conditions)))
  resolved <- intersect(names(digest), bundled)
  if (!any(vapply(digest[resolved], function(x) x$outcome == "ok", logical(1)))) {
    stop("SELFTEST FAILED: no switch-resolved bundled dataset passed ",
         "data_check(). These datasets all fit, so the capture is reaching a ",
         "refusal rather than the real condition set.")
  }
  if (length(uniq) <= 1L) {
    stop("SELFTEST FAILED: every case produced the same ", length(uniq),
         " condition(s). That is a harness defect, not a result -- the most ",
         "likely cause is data_check() failing identically for all of them.")
  }
  if (!identical(digest, second)) {
    differing <- names(digest)[!vapply(names(digest),
      function(n) identical(digest[[n]], second[[n]]), logical(1))]
    stop("SELFTEST FAILED: the capture is not reproducible within one process. ",
         "Differing cases: ", paste(differing, collapse = ", "),
         ". The digest cannot gate a refactor until this is fixed.")
  }
  cat("\nSELFTEST PASSED: two captures in one process are identical.\n")
  cat("Now run it twice as separate processes and diff the .rds files.\n")
  quit(status = 0)
}

saveRDS(list(digest = digest, meta = meta), out_path)
report(digest, meta)
cat("\nwrote", out_path, "\n")

if (do_compare) {
  if (is.null(compare_path) || !file.exists(compare_path)) {
    stop("compare needs a baseline .rds: ... compare dev/dc-before.rds")
  }
  before <- readRDS(compare_path)
  if (identical(before$digest, digest)) {
    cat("\nIDENTICAL to", compare_path, "-- data_check()'s conditions are unchanged.\n")
  } else {
    cat("\nDIFFERENT from", compare_path, "\n")
    for (case in union(names(before$digest), names(digest))) {
      a <- before$digest[[case]]$conditions
      b <- digest[[case]]$conditions
      if (identical(a, b)) next
      cat("\n--", case, "\n")
      for (x in setdiff(a, b)) cat("  LOST  ", x, "\n")
      for (x in setdiff(b, a)) cat("  GAINED", x, "\n")
      oa <- before$digest[[case]]$outcome; ob <- digest[[case]]$outcome
      if (!identical(oa, ob)) cat("  OUTCOME", oa, "->", ob, "\n")
    }
    quit(status = 1)
  }
}
