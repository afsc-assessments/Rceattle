# tools/verify/verify-safebounds.R
# Array-bounds check on the TMB model.
#
# The ordinary build sets `safebounds = FALSE` (src/TMB/compile.R), so an
# out-of-range access is a silent read or write into adjacent memory rather than
# an error. That is not hypothetical here: 5.15.0 fixed exactly such a write,
# where `age_hat` / `age_obs_hat` were sized from `comp_obs` and overran it
# whenever `nlengths < nages` -- "Eigen does not bounds-check in a release
# build, so that was a silent write into adjacent memory, not a crash".
#
# It is also the standing explanation for the one CI failure this net exists to
# catch: R CMD check on windows-latest died with exit code -1073741819
# (0xC0000005, an access violation) inside test-selectivity-catchability.R. An
# access violation is memory corruption, not a bad optimum -- a fit landing in a
# different basin gives a huge gradient and a failed sdreport, not a fault.
#
# This harness rebuilds with bounds checking ON and drives the model over the
# configurations most likely to trip it: ragged composition data (a model whose
# length bins and age bins differ, in both directions), joint-sex rows, tail
# accumulation, and the fleet set the crashing test fits. A violation surfaces
# as an R error naming the object, which is the whole point -- attributable
# instead of a mystery segfault.
#
# Slow by construction: every array access is range-tested, and the model does
# millions per objective evaluation. Minutes, not seconds.
#
#   export PATH=/usr/bin:$PATH
#   RCEATTLE_SAFEBOUNDS=true NOT_CRAN=true Rscript tools/verify/verify-safebounds.R

if (!identical(tolower(Sys.getenv("RCEATTLE_SAFEBOUNDS")), "true")) {
  stop("Run with RCEATTLE_SAFEBOUNDS=true, or this rebuilds without bounds\n",
       "checking and verifies nothing. See src/TMB/compile.R.", call. = FALSE)
}

# ---- positive control -------------------------------------------------------
# Prove the build is actually bounds-checked before trusting a clean result.
# Without this the harness passes just as happily against a stale unchecked .so:
# pkgload only recompiles when sources changed, and this script changes none, so
# setting the flag alone guarantees nothing. Force the rebuild by removing the
# artifacts, then assert -DTMB_SAFEBOUNDS in the compile line.
obj   <- file.path("src", "TMB", "ceattle.o")
so    <- file.path("src", paste0("ceattle", .Platform$dynlib.ext))
shlib <- file.path("src", paste0("Rceattle", .Platform$dynlib.ext))
unlink(c(obj, so, shlib))

# Set it here as well as reading it: the package SHLIB is built by a child
# process further down, and it reads this from the environment.
Sys.setenv(RCEATTLE_SAFEBOUNDS = "true")

# compile.R resolves "ceattle.cpp" against the WORKING DIRECTORY, so it must be
# run from src/TMB. Run from anywhere else and its file.exists() guard is false:
# it prints the safebounds message, compiles nothing, and exits 0 -- which is
# exactly the silent no-op the assertion below exists to catch.
cat("Rebuilding with bounds checking (slow)...\n")
build_log <- local({
  old <- setwd(file.path("src", "TMB"))
  on.exit(setwd(old), add = TRUE)
  suppressWarnings(system2("Rscript", "compile.R",
                           stdout = TRUE, stderr = TRUE,
                           env = "RCEATTLE_SAFEBOUNDS=true"))
})

if (!any(grepl("-DTMB_SAFEBOUNDS", build_log, fixed = TRUE))) {
  cat(tail(build_log, 20), sep = "\n")
  stop("the rebuild did NOT define TMB_SAFEBOUNDS -- this run would have\n",
       "checked nothing. Fix the build before reading any result below.",
       call. = FALSE)
}
cat("  confirmed: -DTMB_SAFEBOUNDS in the compile line\n")

# Put the tree back to a NORMAL build on the way out, rather than just deleting
# the checked artifacts. Two reasons, both learned the hard way:
#
#   * Left alone, the next ordinary load_all() would NOT rebuild -- pkgload only
#     recompiles when sources changed, and this script changes none -- so every
#     later fit would silently use the slow bounds-checked model.
#   * Deleting only the TMB artifacts is worse still. src/Makevars builds the
#     TMB library as a prerequisite of $(SHLIB), so with src/Rceattle.so still
#     present make considers the target current, never re-runs `tmblib`, and
#     load_all() then fails outright on useDynLib(ceattle) -- a broken tree that
#     a plain rebuild does not repair.
#
# So drop the package SHLIB too and rebuild unchecked, leaving things as found.
#
# Not on.exit(): registered at a script's top level it has no function frame to
# attach to and never fires, which is how this script first left a checked .so
# behind. reg.finalizer() on the global environment runs at session end, on a
# normal exit and on an error alike.
restore_normal_build <- function(...) {
  unlink(c(obj, so, shlib))
  message("verify-safebounds: restoring the unchecked build...")
  Sys.setenv(RCEATTLE_SAFEBOUNDS = "false")
  # Build through make, not compile.R alone: that rebuilds the TMB library AND
  # relinks the package SHLIB. Rebuilding only the former leaves src/Rceattle.so
  # missing, and the next plain load_all() fails on useDynLib(ceattle) -- the
  # tree is left unusable by the very step meant to put it back.
  if (requireNamespace("pkgbuild", quietly = TRUE)) {
    try(pkgbuild::compile_dll(".", quiet = TRUE), silent = TRUE)
  } else {
    local({
      old <- setwd(file.path("src", "TMB"))
      on.exit(setwd(old), add = TRUE)
      system2("Rscript", "compile.R", stdout = FALSE, stderr = FALSE)
    })
  }
}
invisible(reg.finalizer(globalenv(), restore_normal_build, onexit = TRUE))

# NAMESPACE loads TWO libraries: `ceattle`, the TMB model built above, and
# `Rceattle`, the package SHLIB. Only the first was built here, and on a fresh
# checkout the second has never existed -- so `compile = FALSE` loaded neither
# and every case below died with "'ceattle' was not found in the list of loaded
# DLLs". It passed on a developer's machine only because a src/Rceattle.so was
# already lying around.
#
# `compile = TRUE` builds it. The TMB library is not rebuilt unchecked by that:
# src/Makevars makes `tmblib` a prerequisite of $(SHLIB), so compile.R runs
# again with RCEATTLE_SAFEBOUNDS still true, and TMB::compile() no-ops because
# the artifacts above are current.
cat("Building the package library...\n")
suppressMessages(pkgload::load_all(".", quiet = TRUE, compile = TRUE))

# The whole run is worthless if the model did not load, and five identical
# "not found in the list of loaded DLLs" errors below would be reported as five
# bounds violations. Fail here instead, where the cause is legible.
loaded <- vapply(getLoadedDLLs(), function(d) d[["name"]], character(1))
missing_dll <- setdiff(c("ceattle", "Rceattle"), loaded)
if (length(missing_dll)) {
  stop("did not load: ", paste(missing_dll, collapse = ", "),
       ". Nothing was bounds-checked; this is a build failure, not a result.",
       call. = FALSE)
}

cat("Rceattle", as.character(utils::packageVersion("Rceattle")),
    "built with safebounds = TRUE\n\n")

results <- list()

run_case <- function(label, expr) {
  t0 <- Sys.time()
  out <- tryCatch({ force(expr); "ok" },
                  error = function(e) paste0("ERROR: ", conditionMessage(e)))
  el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  results[[label]] <<- out
  cat(sprintf("  %-52s %-6s %s\n", label, paste0(el, "s"),
              if (identical(out, "ok")) "ok" else out))
  invisible(out)
}

fit3 <- function(d, ...) suppressWarnings(suppressMessages(
  fit_mod(data_list = d, inits = NULL, file = NULL, estimateMode = 3,
          random_rec = FALSE, msmMode = 0,
          fit_control = fit_control(verbose = 0), ...)))

cat("Bundled data sets, built and evaluated once each:\n")

# GOA2018SS is the discriminating one: nages = (10, 21, 12) against
# nlengths = (7, 26, 117), so it is ragged in BOTH directions -- species 1 has
# fewer length bins than ages, species 3 has ten times as many -- and it carries
# joint-sex (Sex = 3) composition rows, which double the age width the
# composition slot writes. It is also the data the crashing test fits.
run_case("GOA2018SS (ragged comps, joint sex)", {
  data("GOA2018SS", envir = environment())
  f <- fit3(GOA2018SS)
  stopifnot(is.finite(f$obj$fn(f$obj$par)))
})

run_case("BS2017SS", {
  data("BS2017SS", envir = environment())
  f <- fit3(BS2017SS)
  stopifnot(is.finite(f$obj$fn(f$obj$par)))
})

# Multispecies: exercises the 5-D suitability and predation arrays, which are
# the other place a stray index would land somewhere plausible-looking.
run_case("BS2017MS (predation arrays)", {
  data("BS2017MS", envir = environment())
  f <- suppressWarnings(suppressMessages(
    fit_mod(data_list = BS2017MS, inits = NULL, file = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 1, suitMode = 0,
            fit_control = fit_control(verbose = 0))))
  stopifnot(is.finite(f$obj$fn(f$obj$par)))
})

# The configuration the crashing CI job was running: catchability estimated
# across the GOA fleet set, which is where the worker died.
run_case("GOA2018SS, estimated catchability (the CI crash config)", {
  data("GOA2018SS", envir = environment())
  d <- GOA2018SS
  d$fleet_control$Catchability <- 1
  f <- fit3(d)
  stopifnot(is.finite(f$obj$fn(f$obj$par)))
})

# Simulation drives every SIMULATE block, including the composition and CAAL
# draws that write into comp_sim / caal_sim at the raw bin width -- the ragged
# tail those blocks clear by hand.
run_case("sim_mod() draws on ragged comps", {
  data("GOA2018SS", envir = environment())
  f <- fit3(GOA2018SS)
  invisible(suppressWarnings(suppressMessages(sim_mod(f, simulate = TRUE))))
})

# ---- 5.46.0's new array dimensions ------------------------------------------
# Until 5.46.0 every length-indexed array was sized on nlengths. The SS3 options
# added a SECOND length grid: growth_matrix, sel_at_length, non_par_sel and
# log_non_par_sel are sized on max_nlengths_pop, ceattle.cpp walks each data
# bin's run of population bins through pop_bin_lo / pop_bin_hi, and age_error
# gained a third extent (n_ae) indexed by a user column rather than by species.
# That is the shape of the bug this harness was built for -- an array sized from
# one dimension and indexed by another -- so both are driven here.
sb_env <- new.env(parent = asNamespace("Rceattle"))
for (h in c("helpers.R", "helpers-make-msm-data.R")) {
  f <- file.path("tests", "testthat", h)
  if (file.exists(f)) sys.source(f, envir = sb_env)
}

# A population grid four times finer than the data bins, so every population
# array is wider than the data array beside it and a stray index lands inside
# the object rather than past it -- silent, and exactly what bounds checking
# turns into an attributable error.
run_case("finer population length grid (pop_lengths)", {
  set.seed(11)
  d <- sb_env$make_msm_test_data()$data_list
  edges <- sort(unique(d$caal_data$Length[d$caal_data$Species == 1]))
  fine  <- sort(unique(c(edges, seq(min(edges), max(edges), by = 1))))
  d$fleet_control$Selectivity_dimension <- "Length"
  f <- suppressWarnings(suppressMessages(fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3, msmMode = 0,
    random_rec = FALSE,
    growthFun = build_growth(fun = "vonBertalanffy", pop_lengths = fine),
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  stopifnot(is.finite(f$obj$fn(f$obj$par)))
  # The premise: the two grids really do differ, or this case checks nothing.
  stopifnot(length(fine) > d$nlengths[1])
})

# A second ageing-error matrix for species 1, so n_ae exceeds nspp and a fleet
# reads age_error on an index that is NOT its species number.
run_case("second ageing-error matrix (n_ae > nspp)", {
  set.seed(11)
  d <- sb_env$make_msm_test_data()$data_list
  ae <- as.data.frame(d$age_error)
  extra <- ae[ae$Species == 1, , drop = FALSE]
  ae$Ageing_error_index <- ae$Species
  extra$Ageing_error_index <- d$nspp + 1L
  d$age_error <- rbind(ae, extra)
  # Point species 1's fleets at the new matrix, the deepest slice of the array.
  sp1 <- which(d$fleet_control$Species == 1)
  d$fleet_control$Ageing_error_index <- d$fleet_control$Species
  d$fleet_control$Ageing_error_index[sp1] <- d$nspp + 1L
  f <- fit3(d)
  stopifnot(is.finite(f$obj$fn(f$obj$par)))
})

# The file the Windows worker died in, driven here rather than as a second CI
# step. A separate step has to re-establish the bounds-checked build, and it
# cannot: this script restores the unchecked one on exit, and TMB::compile() is
# incremental, so `load_all(compile = TRUE)` afterwards no-ops and the test runs
# unchecked. Measured 2026-08-27 -- ceattle.so's mtime did not move. Run it
# inside the one build instead, where it is checked by construction.
run_case("test-selectivity-catchability.R (the file that crashed CI)", {
  res <- testthat::test_file(
    "tests/testthat/test-selectivity-catchability.R", reporter = "silent")
  df <- as.data.frame(res)
  # A skipped file reports as a pass, which is the failure this whole job
  # exists to defeat.
  if (sum(df$passed) == 0) {
    stop("ran no assertions -- skipped, so nothing was bounds-checked")
  }
  if (sum(df$failed) + sum(df$error) > 0) {
    stop(sum(df$failed) + sum(df$error), " assertion(s) failed under bounds checking")
  }
})

cat("\n")

# "Could not run" is not "found a violation". Reporting them alike is how five
# DLL-load failures were printed as five bounds violations -- a reader would
# conclude the model is riddled with them when in fact nothing was measured.
# A case that could not build the model was never bounds-checked at all.
not_run <- vapply(results, function(x) {
  is.character(x) && grepl("was not found in the list of loaded DLLs|could not find function|there is no package",
                           x)
}, logical(1))
bad <- !not_run & !vapply(results, identical, logical(1), "ok")

if (any(not_run)) {
  cat(sum(not_run), "of", length(results),
      "case(s) COULD NOT RUN, so they were not bounds-checked:\n")
  for (nm in names(results)[not_run]) cat("  -", nm, "\n     ", results[[nm]], "\n")
}
if (any(bad)) {
  cat(sum(bad), "case(s) failed under bounds checking:\n")
  for (nm in names(results)[bad]) cat("  -", nm, "\n     ", results[[nm]], "\n")
}
if (any(not_run) || any(bad)) quit(status = 1)

cat("No bounds violation in", length(results), "configurations.\n")
cat("NOTE: absence here is not proof -- only the paths driven above were\n",
    "checked, and the CI fault is intermittent. Widen the cases rather than\n",
    "concluding the model is clean.\n")
