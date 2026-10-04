# tools/verify/verify-refit-like.R
# Before/after equivalence harness for the .refit_like() refit path. Runs six of
# the eight .refit_like() callers -- retrospective / jitter / self_test / profile
# (9-retro_and_jitter.R), run_mse (10-run_mse.R), remove_F (10-project-no-F.R) --
# on one tiny seeded model and writes a numeric digest. Golden-check does NOT
# cover these paths, so this is the regression net: capture on the clean
# baseline, then require bit-identical after the change.
#
# NOT covered here: sample_rec(update_model = TRUE) (8-sim_mod.R) and
# reweight_comps() (6-reweight.R). They also refit through .refit_like(), so a
# change to that helper needs them checked separately. run_mse() calls
# sample_rec() with update_model = FALSE, so it does not exercise that path.
#
# Output defaults into dev/ (gitignored scratch, created if absent) so the digest
# survives between the two runs you need to compare.
#
# Usage:
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript tools/verify/verify-refit-like.R dev/refit-before.rds
#   NOT_CRAN=true Rscript tools/verify/verify-refit-like.R dev/refit-after.rds compare dev/refit-before.rds

args <- commandArgs(trailingOnly = TRUE)
dir.create("dev", showWarnings = FALSE)
out_path     <- if (length(args) >= 1) args[1] else "dev/refit-digest.rds"
do_compare   <- length(args) >= 2 && args[2] == "compare"
compare_path <- if (length(args) >= 3) args[3] else NULL

suppressMessages(pkgload::load_all(".", quiet = TRUE))

# Recursively walk any return structure; at every fit node (has $quantities --
# true for both optimized fits and estimateMode = 3 build-only objects like
# remove_F's) capture a numeric fingerprint and stop descending. Deterministic
# (list order preserved). obj is NA for build-only objects (no $opt); ssb/R
# still fingerprint them.
collect <- function(x, acc = list()) {
  if (is.list(x) && !is.null(x[["quantities"]])) {
    acc[[length(acc) + 1L]] <- list(
      obj = tryCatch(as.numeric(x$opt$objective), error = function(e) NA_real_),
      ssb = tryCatch(as.numeric(x$quantities$ssb), error = function(e) NA_real_),
      R   = tryCatch(as.numeric(x$quantities$R),   error = function(e) NA_real_))
    return(acc)
  }
  if (is.list(x)) for (el in x) acc <- collect(el, acc)
  acc
}
# Run an entry function; store its fingerprint, or its error message (so a path
# that errors in baseline compares equal only if it errors identically after).
section <- function(expr) {
  r <- tryCatch(expr, error = function(e) structure(list(.err = conditionMessage(e)),
                                                    class = "harness_err"))
  if (inherits(r, "harness_err")) return(r$.err)
  collect(r)
}

digest <- list()

data(BS2017SS)
BS2017SS$projyr <- 2020
BS2017SS$fleet_control$Proj_F_proportion <- rep(1, 7)

base <- Rceattle::fit_mod(
  data_list = BS2017SS, inits = NULL, file = NULL,
  estimateMode = 1, random_rec = FALSE, msmMode = 0,
  fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0))
digest$base <- collect(base)

# sites 1,2 -- retrospective
digest$retro     <- section(Rceattle::retrospective(base, peels = 1, getsd = FALSE))
# ... and the same two sites reached with a SUBSET of depths, which a scalar
# `peels` cannot express: it always dispatches 1:n. This fingerprints the vector
# path and two deeper peels, and nothing more -- on BS2017SS it does NOT reach
# the fixed-effect covariate branch (no `*_linkages`), the once-before-dispatch
# warnings (random_vars is character(0), forecast_rec defaults to "mean") or
# .parallel_lapply() (cores = 1). Those need their own fixture; see
# test-functions-retrospective.R, which covers them.
#
# NOTE: adding a section name makes "bit-identical across all sections"
# unreachable against a baseline captured before it existed -- the compare
# reports `DIFF: retro_sub`. Re-baseline on the commit that adds this line.
digest$retro_sub <- section(Rceattle::retrospective(base, peels = 2:3, getsd = FALSE,
                                                    cores = 1))
# site 3 -- jitter
digest$jitter    <- section(Rceattle::jitter(base, njitter = 1, sd = 0.1, phase = FALSE, seed = 1))
# site 4 -- self_test
digest$self_test <- section(Rceattle::self_test(base, nsim = 1, simulate = FALSE, seed = 1))
# site 5 -- profile
digest$profile   <- section(stats::profile(base, param = "M1", slots = list(c(1, 1, 1)),
                                            values = list(c(0.2, 0.25))))
# site 11 -- remove_F
digest$remove_F  <- section(Rceattle::remove_F(base))

# sites 6/9/10 -- run_mse (normal, Tier-3 EM)
em_tier3 <- suppressWarnings(Rceattle::fit_mod(
  data_list = BS2017SS, inits = base$estimated_params, estimateMode = 2,
  HCR = build_hcr(HCR = 5, Ftarget = 0.4, Flimit = 0.35, Plimit = 0.2, Alpha = 0.05),
  msmMode = 0, fit_control = fit_control(getsd = FALSE, verbose = 0)))
digest$mse_normal <- section(Rceattle::run_mse(
  om = base, em = em_tier3, nsim = 1, assessment_period = 1, sampling_period = 1,
  simulate_data = FALSE, sample_rec = FALSE, seed = 42))

# sites 7/8 -- run_mse with regenerate_past + input-F EM (HCR == 2). Built inside
# the guard so a construction failure records as an error string rather than
# aborting before saveRDS. Ftarget is per-species (matches log_Ftarget).
digest$mse_regen <- section({
  em_inputF <- suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017SS, inits = base$estimated_params, estimateMode = 2,
    HCR = build_hcr(HCR = 2, Ftarget = rep(0.2, base$data_list$nspp)),
    msmMode = 0, fit_control = fit_control(getsd = FALSE, verbose = 0)))
  Rceattle::run_mse(
    om = base, em = em_inputF, nsim = 1, assessment_period = 1, sampling_period = 1,
    simulate_data = FALSE, regenerate_past = TRUE, sample_rec = FALSE, seed = 7)
})

# site 9 (OM refit), MULTISPECIES -- the important case: the OM's suitability
# window (suit_styr/suit_endyr) is pinned to the pristine om$ so predation
# suitability stays constant through the MSE. A single-species MSE cannot
# exercise that pin (suitMode = 0 is inert), so run a Bering Sea multispecies
# MSE warm-started from the SS MLEs (golden recipe) to fingerprint it.
digest$mse_ms <- section({
  data(BS2017MS)
  BS2017MS$projyr <- 2020
  BS2017MS$fleet_control$Proj_F_proportion <- rep(1, nrow(BS2017MS$fleet_control))
  ms_om <- suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017MS, inits = base$estimated_params, file = NULL,
    estimateMode = 1, random_rec = FALSE, msmMode = 1, suitMode = 0, niter = 5,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0)))
  # Multispecies mode supports only HCRs 0/1/2/3/6; use ConstantF (2), per-species.
  ms_em <- suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017MS, inits = ms_om$estimated_params, estimateMode = 2,
    HCR = build_hcr(HCR = 2, Ftarget = rep(0.15, ms_om$data_list$nspp)),
    msmMode = 1, suitMode = 0, niter = 5,
    fit_control = fit_control(getsd = FALSE, verbose = 0)))
  Rceattle::run_mse(
    om = ms_om, em = ms_em, nsim = 1, assessment_period = 1, sampling_period = 1,
    simulate_data = FALSE, sample_rec = FALSE, seed = 11)
})

saveRDS(digest, out_path)
cat("wrote digest ->", out_path, "\n")
for (nm in names(digest)) {
  v <- digest[[nm]]
  n <- if (is.character(v)) paste0("ERR(", substr(v, 1, 40), ")") else paste0(length(v), " fit(s)")
  cat(sprintf("  %-12s %s\n", nm, n))
}

if (do_compare && !is.null(compare_path)) {
  before <- readRDS(compare_path)
  cat("\n=== COMPARE vs", compare_path, "===\n")
  if (identical(before, digest)) {
    cat("PASS: bit-identical across all sections.\n")
  } else {
    cat("FAIL: sections differ:\n")
    for (nm in union(names(before), names(digest))) {
      if (!identical(before[[nm]], digest[[nm]])) cat("  DIFF:", nm, "\n")
    }
  }
}
