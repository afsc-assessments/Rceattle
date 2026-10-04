# Simulation recovery of the initial recruitment level, srr_linkages `R_init`
# (5.48.0).
#
# The level multiplies the initial age-structure only, so it is informed by the
# first year's composition and index observations and by nothing else -- R0 is
# identified by the whole recruitment series, and the per-age departures by
# init_dev. That makes recovery the check that separates a correct offset from a
# plausible one: a level that is not identified would come back at its starting
# value with a flat gradient, and the unit tests, which drive the tape directly,
# could not tell the difference.
#
# This harness imposes a known level on make_test_data()'s single species,
# simulates every observation with sim_mod(), refits with the level free, and
# reports the mean estimate against the truth. A mean more than two standard
# errors of the mean away is the number to investigate.
#
# Usage: Rscript tools/verify/verify-sim-recovery-r-init.R [n_reps] [true_multiplier] [seed_offset]
# Run time: a few seconds per replicate.

args <- commandArgs(trailingOnly = TRUE)
n_reps  <- if (length(args) >= 1) as.integer(args[1]) else 20L
true_mult <- if (length(args) >= 2) as.numeric(args[2]) else 0.5
seed0   <- if (length(args) >= 3) as.integer(args[3]) else 2000L
true_lvl <- log(true_mult)

suppressMessages(pkgload::load_all(".", compile = FALSE, quiet = TRUE))

# 40 years so the first-year observations are a small share of the data: the
# level has to be identified, not absorbed by the rest of the series.
d <- make_test_data(nyrs = 40)
d$env_data <- data.frame(Year = d$styr:d$projyr, lvl = 1)
spec <- build_srr(linkages = list(R_init = linkage_spec(~ 1)))

qfit <- function(dat, inits = NULL, mode = "Hindcast", sd = FALSE) {
  suppressMessages(suppressWarnings(fit_mod(
    data_list = dat, inits = inits, msmMode = 0, estimateMode = mode,
    recFun = spec, initMode = "NonEquilibrium", random_rec = FALSE,
    fit_control = fit_control(phase = FALSE, getsd = sd, verbose = 0))))
}

# Operating model: the level-free MLE with the level imposed.
base <- qfit(d)
ini  <- base$estimated_params
ini$beta_linkage[1] <- true_lvl
om   <- qfit(d, inits = ini, mode = "DebugBuild")

# No reported SE. Not because the fit is poor -- it reaches max|grad| ~1e-05 --
# but because this fixture's Hessian is rank-deficient in four directions that
# load on log_sel_slp / sel_inf for fleets with no composition data, so
# sdreport() returns nothing at all. The level has ZERO loading on those four
# and a conditional curvature of ~460 (conditional SD ~0.05), so it is itself
# identified; the empirical spread over replicates is the honest measure here.
#
# The replicate statistics are NOT reproducible to two significant figures:
# `base` below is an optimization, so `ini`, the OM and every draw inherit a
# machine-dependent starting vector. The runaway COUNT in particular moves
# between machines (3 and 6 of 60 both observed on the same arguments). Read
# the rate as "a few percent", and re-measure rather than quoting a past run.
res <- data.frame()
for (s in seq_len(n_reps)) {
  set.seed(seed0 + s)
  sim <- sim_mod(om, simulate = TRUE)
  em  <- qfit(sim, inits = ini)
  res <- rbind(res, data.frame(rep = s, est = em$estimated_params$beta_linkage[1],
                               maxgr = max(abs(em$obj$gr()))))
  cat(sprintf("rep %2d  est %.4f (x%.4g)  max|grad| %.1e\n",
              s, res$est[s], exp(res$est[s]), res$maxgr[s]))
}

# A replicate that walks into the flat tail below reports a CLEAN gradient at a
# level the data cannot distinguish from an annihilated initial cohort, so a
# gradient filter alone does not find it. Separate those out rather than
# averaging them in, and report how often it happened -- that rate is the
# number a real assessment cares about.
ok      <- res[res$maxgr < 1e-2, ]
runaway <- ok[abs(ok$est - true_lvl) > 3, ]
good    <- ok[abs(ok$est - true_lvl) <= 3, ]
cat(sprintf("\n%d of %d converged; %d of those ran away past 3 log units (x%.3g and below)\n",
            nrow(ok), n_reps, nrow(runaway), if (nrow(runaway)) max(exp(runaway$est)) else 0))
cat(sprintf("true %.4f (x%.3f)   mean %.4f   empirical sd %.4f   n = %d\n",
            true_lvl, true_mult, mean(good$est), sd(good$est), nrow(good)))
cat(sprintf("z of the mean: %.2f  (share below truth %.2f)\n",
            (mean(good$est) - true_lvl) / (sd(good$est) / sqrt(nrow(good))),
            mean(good$est < true_lvl)))
if (nrow(runaway)) {
  cat(sprintf("\nRUNAWAY REPLICATES (%s): the level is weakly identified at the low end.\n",
              paste(runaway$rep, collapse = ", ")))
  cat("Bound it on the natural scale to keep a fit out of that tail, e.g.\n")
  cat("  linkage_spec(~ 1, bounds = list(`(Intercept)` = c(0.05, 20)))\n")
}

# The level must stay finite and differentiable across its whole plausible
# range: a collapsed initial state (x1e-13) and a 1e13-fold one both have to
# give a finite objective, or the optimiser can walk off the tape.
dbg <- qfit(d, inits = ini, mode = "DebugBuild")
p <- dbg$obj$env$last.par.best
i <- grep("beta_linkage", names(p))
cat("\nlimiting cases (objective and gradient must stay finite):\n")
for (b in c(-30, -10, 0, 10, 30)) {
  p[i] <- b
  cat(sprintf("  beta %6.1f  fn %.6g  gr %.3e\n",
              b, dbg$obj$fn(p), dbg$obj$gr(p)[i]))
}
