# Shared fixture for the diet-simulation tests (test-functions-sim-mod-diet*.R).
#
# These live in a helper rather than in one test file because that file was 66%
# of the whole coverage run: ~1200 sim_mod() calls at nrep = 150, each
# evaluating the compiled model through covr-instrumented R, which took ~78
# minutes of a 119-minute suite and set the wall-clock floor no matter how many
# shards it was split across. The assertions are unchanged; they are spread over
# several files so the shard planner can put them on different runners.
#
# `.diet_sim_env` caches the base fit per PROCESS, so files that land on the
# same shard still fit it once between them.

# Per stomach, in diet_data row order: the prey proportions the density predicts
# (offset and renormalized), and the same transform applied to a simulated row.
.diet_transform <- function(mod, values, offset = 1e-5) {
  re  <- Rceattle::rearrange_data(mod$data_list)
  sid <- re$stomach_id
  out <- rep(NA_real_, length(values))
  for (i in unique(sid)) {
    idx <- which(sid == i)
    p <- c(values[idx], max(1 - sum(values[idx]), 0))
    p <- (p + offset) / sum(p + offset)
    out[idx] <- p[seq_along(idx)]
  }
  out
}

# What the density predicts for each prey bin.
.diet_pred <- function(mod, offset = 1e-5) {
  .diet_transform(mod, mod$quantities$diet_hat[, 2], offset = offset)
}

# What the refit would see, given a simulated table.
.diet_seen <- function(mod, sim, offset = 1e-5) {
  .diet_transform(mod, sim, offset = offset)
}

# Fit BS2017MS once with an estimated suitability, so the diet likelihood is
# active and diet_hat is concentrated enough to measure a spread on. Then hold
# every parameter there and vary only the composition family and its
# concentration, so those are the single thing that differs between cases.
# NOTE the parameter `diet_comp_weights` means different things by family: for
# the multinomial it is used raw, as a likelihood weight / effective-sample-size
# multiplier (1 = neutral), and for the Dirichlet-multinomial it is log theta
# (exp() is taken). Passing log(1) = 0 for a multinomial sets the weight to
# ZERO, which drops the data from the likelihood and from the draw.
.diet_sim_env <- new.env(parent = emptyenv())
.diet_sim_fit <- function(family, weight_par) {
  if (is.null(.diet_sim_env$base)) {
    utils::data("BS2017SS", package = "Rceattle", envir = environment())
    utils::data("BS2017MS", package = "Rceattle", envir = environment())
    ss <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
      BS2017SS, inits = NULL, msmMode = 0, estimateMode = 1,
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
    .diet_sim_env$dat <- BS2017MS
    .diet_sim_env$base <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
      BS2017MS, inits = ss$estimated_params, msmMode = 1, suitMode = 4,
      estimateMode = 1, niter = 3,
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  }
  d <- .diet_sim_env$dat
  d$Diet_distribution <- rep(family, d$nspp)
  inits <- .diet_sim_env$base$estimated_params
  inits$diet_comp_weights[] <- weight_par
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, inits = inits, msmMode = 1, suitMode = 4, estimateMode = 3, niter = 3,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
}

.diet_reps <- function(mod, nrep, seed = 7) {
  set.seed(seed)
  replicate(nrep, suppressWarnings(
    Rceattle::sim_mod(mod, simulate = TRUE)$diet_data$Stomach_proportion_by_weight))
}

# The round-trip check, for one composition family.
#
# Split by family because the two-family x two-offset loop it came from was 600
# of the file's 1200 sim_mod() calls -- half the cost of the whole diet suite in
# a single test_that(). One family per file halves the floor again.
#
# Run at a raised comp_offset as well as the default. Drawing from the
# offset-and-renormalized prediction rather than the raw one applies the offset
# twice over the round trip, which is 0.02% at the default and over 20% at
# 1e-2 -- invisible at the default, unmissable here.
.diet_roundtrip_check <- function(fam) {
  mod <- .diet_sim_fit(fam, if (fam == "Multinomial") 1 else log(1))
  # One seeded draw serves both offsets: `off` enters only the transform below.
  reps <- .diet_reps(mod, nrep = 150)
  for (off in c(1e-5, 1e-2)) {
    pred <- .diet_pred(mod, offset = off)
    seen <- apply(reps, 2, function(x) .diet_seen(mod, x, offset = off))
    keep <- which(pred > 1e-4)
    testthat::expect_gt(length(keep), 50)
    # 0.04 sits well above the Monte Carlo noise of this aggregate and well
    # below the +23% the pre-fix draw carried at comp_offset = 1e-2.
    testthat::expect_equal(sum(rowMeans(seen)[keep]) / sum(pred[keep]), 1,
                           tolerance = 0.04,
                           label = paste(fam, "comp_offset", off))
  }
}
