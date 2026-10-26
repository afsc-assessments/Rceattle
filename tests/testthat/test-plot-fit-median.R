# The predicted line on an index or catch fit figure is the MEDIAN of the fitted
# observation distribution, not the reported mean: the lognormal likelihood fits
# log(y) to log(hat) - b*sigma^2/2, so hat is the mean while the bars drawn
# beside it are quantiles bracketing the median.
#
# The expected values are computed by hand from that definition -- exp(-0.125)
# at sigma 0.5 -- not read back from the plotter, so dropping, doubling or
# sign-flipping the offset all fail. sigma is injected, so nothing has to
# converge for the arithmetic to be exact.
testthat::skip_on_cran()

# Built, not optimized: this is about a transform, not about the optimum.
.fitm_fit <- function() {
  data("BS2017SS", package = "Rceattle", envir = environment())
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017SS, inits = NULL, estimateMode = 3, msmMode = 0,
    fit_control = Rceattle::fit_control(verbose = 0))))
}

# Known sd and flag, so the expected median is exactly hat * exp(-b*sd^2/2).
.fitm_set <- function(fit, sd, ba_obs) {
  fit$quantities$index_sd[] <- sd
  fit$quantities$catch_sd[] <- sd
  fit$data_list$bias_adjust_obs <- ba_obs
  fit
}


testthat::test_that("the index line is the median at bias_adjust_obs = 1", {
  testthat::skip_if_not_installed("TMB")
  fit <- .fitm_set(.fitm_fit(), sd = 0.5, ba_obs = 1)

  p  <- suppressWarnings(Rceattle::plot_index(fit))
  df <- p$data
  testthat::expect_gt(nrow(df), 0L)

  # Premise: every drawn row here is lognormal, so every one is offset.
  nat <- Rceattle:::.index_rows_natural_scale(fit$data_list)
  testthat::expect_false(any(nat))

  idx  <- fit$data_list$index_data
  keep <- idx$Year > 0 & idx$Year <= fit$data_list$endyr
  want <- as.numeric(fit$quantities$index_hat)[keep] * exp(-0.125)
  testthat::expect_equal(df$Predicted, want, tolerance = 1e-10)
  # The median is below the mean.
  testthat::expect_true(all(df$Predicted <
                              as.numeric(fit$quantities$index_hat)[keep]))
  # The gap this fixes: 13.3% at sigma 0.5.
  testthat::expect_equal(exp(0.125) - 1, 0.1331485, tolerance = 1e-6)
})


testthat::test_that("the index line is the reported hat at bias_adjust_obs = 0", {
  testthat::skip_if_not_installed("TMB")
  fit <- .fitm_set(.fitm_fit(), sd = 0.5, ba_obs = 0)

  df   <- suppressWarnings(Rceattle::plot_index(fit))$data
  idx  <- fit$data_list$index_data
  keep <- idx$Year > 0 & idx$Year <= fit$data_list$endyr
  testthat::expect_equal(df$Predicted,
                         as.numeric(fit$quantities$index_hat)[keep],
                         tolerance = 0)
})


testthat::test_that("the observation bars are untouched by the median fix", {
  testthat::skip_if_not_installed("TMB")
  # The bars are observation quantiles, so the flag must move only the line.
  on_  <- suppressWarnings(Rceattle::plot_index(.fitm_set(.fitm_fit(), 0.5, 1)))$data
  off_ <- suppressWarnings(Rceattle::plot_index(.fitm_set(.fitm_fit(), 0.5, 0)))$data

  testthat::expect_equal(on_$Lower95, off_$Lower95, tolerance = 0)
  testthat::expect_equal(on_$Upper95, off_$Upper95, tolerance = 0)
  testthat::expect_equal(on_$Observation, off_$Observation, tolerance = 0)
  # Not vacuous: the line DID move.
  testthat::expect_false(isTRUE(all.equal(on_$Predicted, off_$Predicted)))
})


testthat::test_that("a natural-scale index fleet is unchanged by either flag", {
  testthat::skip_if_not_installed("TMB")
  # No offset in its likelihood, and its sigma is an ABSOLUTE sd.
  data("BS2017SS", package = "Rceattle", envir = environment())
  d <- BS2017SS
  # Fleet_type is NUMERIC before switch_check() (1 fishery, 2 survey) and a
  # string after, and Index_distribution does not exist until then -- so set it
  # at full length here, by the numeric code.
  srv <- which(d$fleet_control$Fleet_type == 2)[1]
  d$fleet_control$Index_distribution <- rep("Lognormal", nrow(d$fleet_control))
  d$fleet_control$Index_distribution[srv] <- "Normal"
  base <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 0,
    fit_control = Rceattle::fit_control(verbose = 0))))

  nat_rows <- Rceattle:::.index_rows_natural_scale(base$data_list)
  testthat::expect_true(any(nat_rows))   # the fixture really has one

  code <- base$data_list$fleet_control$Fleet_code[srv]
  got <- lapply(c(1, 0), function(b) {
    f <- .fitm_set(base, sd = 0.5, ba_obs = b)
    df <- suppressWarnings(Rceattle::plot_index(f))$data
    df$Predicted[df$Fleet == as.character(
      base$data_list$fleet_control$Fleet_name[srv])]
  })
  testthat::expect_equal(got[[1]], got[[2]], tolerance = 0)
})


testthat::test_that("the catch line is the median, on every fishery", {
  testthat::skip_if_not_installed("TMB")
  # No family switch for catch: lognormal throughout, so always offset.
  fit <- .fitm_set(.fitm_fit(), sd = 0.5, ba_obs = 1)
  df  <- suppressWarnings(Rceattle::plot_catch(fit))$data
  testthat::expect_gt(nrow(df), 0L)

  cat_d <- fit$data_list$catch_data
  fc    <- fit$data_list$fleet_control
  keep  <- cat_d$Fleet_code %in% fc$Fleet_code[fc$Fleet_type == "Fishery"] &
    cat_d$Year > 0 & cat_d$Year <= fit$data_list$endyr
  want  <- as.numeric(fit$quantities$catch_hat)[keep] * exp(-0.125)
  testthat::expect_equal(df$Predicted, want, tolerance = 1e-10)
})


testthat::test_that("plot_index_residual() pearson matches residuals.Rceattle()", {
  testthat::skip_if_not_installed("TMB")
  # The plotter used to compute a bare log(obs) - log(hat) and disagree.
  fit <- .fitm_set(.fitm_fit(), sd = 0.4, ba_obs = 1)

  pl <- suppressWarnings(
    Rceattle::plot_index_residual(fit, residual_type = "pearson"))$data
  rs <- stats::residuals(fit, source = "index", type = "pearson")
  rs <- rs[!is.na(rs$Residual), ]

  # Join on (fleet, year): the two filter rows differently.
  key_pl <- paste(pl$Fleet, pl$Year)
  key_rs <- paste(rs$Fleet_name, rs$Year)
  common <- intersect(key_pl, key_rs)
  testthat::expect_gt(length(common), 0L)
  testthat::expect_equal(pl$Residual[match(common, key_pl)],
                         rs$Residual[match(common, key_rs)],
                         tolerance = 1e-10)

  # Not vacuous: at sigma 0.4 the offset is 0.2 of a standardized residual.
  testthat::expect_equal(0.4^2 / 2 / 0.4, 0.2, tolerance = 1e-12)
})


testthat::test_that("the offset uses each row's own sd, not a shared one", {
  testthat::skip_if_not_installed("TMB")
  fit <- .fitm_fit()
  fit$data_list$bias_adjust_obs <- 1
  # A distinct sd on every row: a constant sd would agree under ANY permutation
  # of the row subset, so the blocks above cannot see a misalignment.
  n <- length(fit$quantities$index_sd)
  fit$quantities$index_sd[] <- seq(0.1, 0.6, length.out = n)

  df   <- suppressWarnings(Rceattle::plot_index(fit))$data
  idx  <- fit$data_list$index_data
  keep <- idx$Year > 0 & idx$Year <= fit$data_list$endyr
  sd_k <- as.numeric(fit$quantities$index_sd)[keep]
  want <- as.numeric(fit$quantities$index_hat)[keep] * exp(-sd_k^2 / 2)
  testthat::expect_equal(df$Predicted, want, tolerance = 1e-10)
  # Not vacuous: the sds genuinely differ row to row.
  testthat::expect_gt(stats::sd(sd_k), 0)
})

testthat::test_that("a projection row carries no bias offset", {
  testthat::skip_if_not_installed("TMB")
  # A projection prediction is the removal taken under the HCR, not an
  # observation drawn with error, so the flag must not move it. Compared
  # between flag settings rather than against catch_hat, which needs no join.
  on_  <- suppressWarnings(Rceattle::plot_catch(
    .fitm_set(.fitm_fit(), sd = 0.5, ba_obs = 1), incl_proj = TRUE))$data
  off_ <- suppressWarnings(Rceattle::plot_catch(
    .fitm_set(.fitm_fit(), sd = 0.5, ba_obs = 0), incl_proj = TRUE))$data
  endyr <- .fitm_fit()$data_list$endyr

  proj <- on_$Year > endyr
  testthat::skip_if(!any(proj), "no projection catch rows drawn")
  testthat::expect_equal(on_$Predicted[proj], off_$Predicted[proj],
                         tolerance = 1e-12)
  # ...while the hindcast rows in the same figure DO move with the flag.
  testthat::expect_false(
    isTRUE(all.equal(on_$Predicted[!proj], off_$Predicted[!proj])))
})


testthat::test_that("an MSE operating model takes no bias offset at all", {
  testthat::skip_if_not_installed("TMB")
  # run_mse() advances the OM's endyr to the last assessment year
  # (10-run_mse.R), so Year <= endyr is TRUE for every row it draws and no year
  # test can separate hindcast from projection there. Its realized catch is
  # catch_hat by construction, so converting would draw the line below an
  # observation it equals exactly. Simulated by advancing endyr, which is the
  # one property of the OM that matters here.
  adv <- function(ba) {
    f <- .fitm_set(.fitm_fit(), sd = 0.5, ba_obs = ba)
    f$data_list$endyr <- f$data_list$projyr
    f
  }
  on_  <- Rceattle:::.fleet_fit_df(list(adv(1)), kind = "catch",
                                   incl_proj = TRUE, mse = TRUE)
  hat  <- as.numeric(adv(1)$quantities$catch_hat)

  # Premise: without the mse flag every one of these rows WOULD convert, so
  # this is not vacuous -- the advanced endyr defeats the year test.
  naive <- Rceattle:::.fleet_fit_df(list(adv(1)), kind = "catch",
                                    incl_proj = TRUE)
  testthat::expect_false(isTRUE(all.equal(naive$Predicted, on_$Predicted)))

  # ...and under the flag the prediction is the reported value, so it cannot
  # move with bias_adjust_obs.
  off_ <- Rceattle:::.fleet_fit_df(list(adv(0)), kind = "catch",
                                   incl_proj = TRUE, mse = TRUE)
  testthat::expect_equal(on_$Predicted, off_$Predicted, tolerance = 1e-12)
  testthat::expect_equal(on_$Predicted, hat[seq_len(nrow(on_))],
                         tolerance = 1e-10)
})
