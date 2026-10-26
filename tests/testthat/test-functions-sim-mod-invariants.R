# Invariants of the TMB simulation path that no golden model reaches.
#
# The four reference models carry no natural-scale index, no AMAK/Ianelli
# stock-recruit penalty, and no zero-sample-size composition row, so
# /golden-check is silent on every case here. Each test pins one property that,
# if it broke, would not error -- it would quietly hand a self_test() or an MSE
# data it was never told were simulated.

testthat::skip_on_cran()

# ---------------------------------------------------------------------------
# Simulated data must not come back under the observed data's name.
#
# TMB never clears its report environment, so a draw written into a DATA_ object
# and REPORTed under that name stays readable as the observed data for the rest
# of the object's life. Every observation type therefore draws into a copy and
# reports it as *_obs_sim.
# ---------------------------------------------------------------------------
testthat::test_that("obj$simulate() leaves the observed-data reports untouched", {
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data(nyrs = 12, nages = 5, seed = 7)
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

  obs <- fit$obj$report()          # no draw taken
  set.seed(3)
  sim <- fit$obj$simulate()        # draw taken

  # catch_obs and index_obs are not REPORTed at all, so their in-place draw
  # cannot leak; comp_obs, caal_obs and diet_obs are, and each has to be drawn
  # into a copy. Assert which is which, so adding a REPORT() for one of the first
  # two without also giving it a copy fails here rather than silently handing a
  # replicate back under the observed name.
  reported <- c("comp_obs", "caal_obs", "diet_obs")
  unreported <- c("catch_obs", "index_obs")
  testthat::expect_true(all(vapply(reported, function(n) !is.null(obs[[n]]), logical(1))))
  testthat::expect_true(all(vapply(unreported, function(n) is.null(obs[[n]]), logical(1))))

  for (nm in reported) {
    testthat::expect_equal(sim[[nm]], obs[[nm]], tolerance = 1e-12,
                           info = paste(nm, "was overwritten by the draw"))
  }

  # ...and the draw really did happen, or the loop above is vacuous.
  testthat::expect_false(is.null(sim$comp_obs_sim))
  testthat::expect_false(isTRUE(all.equal(as.numeric(sim$comp_obs_sim),
                                          as.numeric(obs$comp_obs))))
})

# ---------------------------------------------------------------------------
# A composition row with no sample size comes back EMPTY, not filled in from the
# prediction.
#
# comp_hat rows are normalized to sum to one, so falling back to the prediction
# would return noise-free proportions summing to one -- a row that was never
# sampled, indistinguishable from a real full-weight composition. run_mse()
# decides "no catch, so no composition" by testing whether the row sums above
# zero, and it zeroes Sample_size for empty rows and feeds them back in, so this
# state is reachable on the next assessment rather than hypothetical.
# ---------------------------------------------------------------------------
testthat::test_that("a zero-Sample_size composition row simulates to zero", {
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data(nyrs = 12, nages = 5, seed = 7)
  # The fixture carries exactly one composition row per fleet, so zeroing that
  # row would leave the fleet with no composition at all and data_check() would
  # refuse the model before the draw is ever reached. Give fleet 1 a second row
  # in another year and zero THAT one, so the fleet keeps a usable row and the
  # zero-sample-size case is the only thing under test.
  extra <- d$comp_data[1, , drop = FALSE]
  extra$Year <- d$comp_data$Year[1] + 1L
  testthat::skip_if(extra$Year > d$endyr, "no spare year in the fixture")
  d$comp_data <- rbind(d$comp_data, extra)
  zero_row <- nrow(d$comp_data)
  d$comp_data$Sample_size[zero_row] <- 0

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

  set.seed(5)
  sim <- fit$obj$simulate()
  got <- as.numeric(sim$comp_obs_sim[zero_row, ])

  testthat::expect_equal(sum(got), 0, tolerance = 1e-12)
  # The sampled row of the same fleet still draws, or the assertion above would
  # pass on a simulator that had stopped working altogether.
  testthat::expect_gt(sum(as.numeric(sim$comp_obs_sim[1, ])), 0)
})

# ---------------------------------------------------------------------------
# Recruitment scored by two densities is not redrawn.
#
# Under srr_fun = 0 with srr_pred_fun > 0 (AMAK/Ianelli) the stock-recruit curve
# is fitted as a penalty on the recruitment deviations, so rec_dev is scored by
# JNLL_REC_DEV and again by JNLL_SRR_PENALTY. Two penalties on one latent do not
# compose into a distribution to draw from: drawing at R_sd from the first alone
# ignores the second, and the second couples the deviations to the stock-recruit
# parameters rather than being a density on rec_dev alone, so there is no
# closed-form conditional to use instead.
# ---------------------------------------------------------------------------
testthat::test_that("the AMAK/Ianelli stock-recruit penalty suppresses the recruitment draw", {
  testthat::skip_if_not_installed("TMB")

  mk <- function(srr_pred) {
    d <- make_test_data(nyrs = 15, nages = 5, seed = 21)
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
      recFun = Rceattle::build_srr(srr_fun = 0, srr_pred_fun = srr_pred),
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
  }

  # srr_pred_fun = 0: one density on rec_dev, so recruitment is drawn.
  plain <- mk(0)
  testthat::expect_true(as.logical(plain$obj$report()$rec_srr_single_density))
  set.seed(2)
  sim_plain <- Rceattle::sim_mod(plain, simulate = TRUE, process = "recruitment")
  truth <- attr(sim_plain, "process_sim")
  testthat::expect_false(is.null(truth$rec_dev))
  testthat::expect_false(isTRUE(all.equal(as.numeric(truth$rec_dev),
                                          as.numeric(plain$estimated_params$rec_dev))))

  # srr_pred_fun > 0 with srr_fun = 0: two densities, so nothing is drawn, the
  # fitted deviations stand, and sim_mod() says why.
  pen <- mk(2)
  testthat::expect_false(as.logical(pen$obj$report()$rec_srr_single_density))
  set.seed(2)
  w <- testthat::capture_warnings(
    sim_pen <- Rceattle::sim_mod(pen, simulate = TRUE, process = "recruitment"))
  testthat::expect_true(any(grepl("stock-recruit curve as a penalty", w)))
  # No truth is returned, because none was generated -- handing back the fitted
  # deviations would look like a truth to recover.
  testthat::expect_null(attr(sim_pen, "process_sim")$rec_dev)
})

# ---------------------------------------------------------------------------
# A 2DAR1 selectivity fleet asked for process = "selectivity" is not silently
# frozen.
#
# The template draws only the linkage-grammar random effects; sel_coff_dev
# carries a density (IID, 2D-AR1 or 3D-AR1 GMRF) but no draw.
# .sim_warn_process_unsimulated() keys off Time_varying_sel, which a 2DAR1 fleet
# does not set -- so the warning has to reach it some other way, or the caller
# gets fitted values dressed as a simulation.
# ---------------------------------------------------------------------------
testthat::test_that("a 2DAR1 selectivity fleet reports that its deviations were not drawn", {
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data(nyrs = 12, nages = 5, seed = 7)
  d$fleet_control$Selectivity[1] <- "2DAR1"
  d$fleet_control$N_sel_bins[1]  <- 4

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

  set.seed(9)
  w <- testthat::capture_warnings(
    sim <- Rceattle::sim_mod(fit, simulate = TRUE, process = "selectivity"))

  testthat::expect_true(any(grepl("selectivity", w, ignore.case = TRUE)),
                        info = "a 2DAR1 fleet was frozen without a word")
  # Nothing is returned as truth for a process that was not drawn.
  testthat::expect_null(attr(sim, "process_sim")$sel_coff_dev)
})

# ---------------------------------------------------------------------------
# A missing hindcast survey year stays missing; a projection row is still drawn.
#
# The index likelihood skips a hindcast row whose Observation is NA, but the
# draw used to fill every row, so a self_test() refit fitted survey years the
# real data lacked and recovery looked better than it was. run_mse() still needs
# the negative-Year rows drawn, since it reveals them as each assessment's data.
# ---------------------------------------------------------------------------
testthat::test_that("a missing hindcast survey year simulates to NA", {
  testthat::skip_if_not_installed("TMB")

  for (dist in c("Lognormal", "Normal", "TruncatedNormal")) {
    d <- make_test_data(nyrs = 12, nages = 5, seed = 7)
    srv <- which(d$index_data$Fleet_name == "Survey")
    d$fleet_control$Index_distribution[d$fleet_control$Fleet_name == "Survey"] <- dist
    if (dist != "Lognormal") d$index_data$Log_sd[srv] <- 10   # absolute sd
    missing_row <- srv[5]
    d$index_data$Observation[missing_row] <- NA
    proj <- d$index_data[utils::tail(srv, 1), , drop = FALSE]
    proj$Year <- -(d$endyr + 1)
    proj$Observation <- NA
    d$index_data <- rbind(d$index_data, proj)
    proj_row <- nrow(d$index_data)

    fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
      d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

    set.seed(11)
    sim <- testthat::expect_no_warning(Rceattle::sim_mod(fit, simulate = TRUE))
    obs <- sim$index_data$Observation
    testthat::expect_true(is.na(obs[missing_row]), info = dist)
    testthat::expect_gt(obs[proj_row], 0)
    testthat::expect_true(all(obs[setdiff(srv, missing_row)] > 0), info = dist)
  }
})

# ---------------------------------------------------------------------------
# The same, for `simulate = FALSE`. This is the path `self_test(simulate =
# FALSE)` and `run_mse(regenerate_past = TRUE)` take, and it wrote `index_hat`
# over every row -- so a survey year the operating model never produced came
# back as the model's PREDICTION for it, which a refit then fits. That is the
# same defect as drawing one, through a different door, and it is why the
# predicate is shared between the two branches rather than written twice.
# ---------------------------------------------------------------------------
testthat::test_that("a missing hindcast survey year stays NA without simulation", {
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data(nyrs = 12, nages = 5, seed = 7)
  srv <- which(d$index_data$Fleet_name == "Survey")
  missing_row <- srv[5]
  d$index_data$Observation[missing_row] <- NA
  proj <- d$index_data[utils::tail(srv, 1), , drop = FALSE]
  proj$Year <- -(d$endyr + 1)
  proj$Observation <- NA
  d$index_data <- rbind(d$index_data, proj)
  proj_row <- nrow(d$index_data)

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

  sim <- suppressMessages(suppressWarnings(Rceattle::sim_mod(fit, simulate = FALSE)))
  obs <- sim$index_data$Observation
  # The unfitted hindcast year keeps its NA ...
  testthat::expect_true(is.na(obs[missing_row]))
  # ... while the projection row DOES take the expected value, which is what
  # run_mse() splices in, and every fitted year is filled as before.
  testthat::expect_gt(obs[proj_row], 0)
  testthat::expect_true(all(obs[setdiff(srv, missing_row)] > 0))
  # Not vacuous: the model does predict a value for that year, so the old
  # behaviour had something to write and this is not an empty assertion.
  testthat::expect_gt(as.numeric(fit$quantities$index_hat)[missing_row], 0)
})
