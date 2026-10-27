# hindcast_skill(): forecast skill across retrospective peels.
#
# This answers a different question from Mohn's rho. Rho measures how the
# ESTIMATE of a year moves as data accumulate -- estimation consistency. This
# measures how well a peel PROJECTED years it could not see, given the catch
# that was actually taken, which is what a recruitment-projection assumption
# (proj_mean_rec TRUE/FALSE, or a DSEM) should be judged on.

# hindcast_skill() does not return the peel models, so recompute the one we need
# rather than asserting against numbers the function itself produced.
retro_peel_for <- function(hs, fit, peel) {
  r <- suppressWarnings(suppressMessages(
    Rceattle::retrospective(fit, peels = peel, cores = 1, getsd = FALSE,
                            forecast_rec = "model")))
  r$Rceattle_list[[1]]
}

testthat::test_that("hindcast_skill scores peels against the full model", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = 0, random_rec = TRUE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
  hs <- suppressWarnings(suppressMessages(
    Rceattle::hindcast_skill(fit, peels = 2, quantity = "ssb",
                             reference = "model", cores = 1, getsd = FALSE)))

  testthat::expect_true(all(c("mase", "by_year") %in% names(hs)))
  # Keyed by `Forecast year`, not peel: the scoring averages ACROSS peels at a
  # fixed horizon (Kell et al. 2021 eq. 5), so `N` is the count that matters.
  testthat::expect_true(all(c("Forecast year", "species", "Object", "N",
                              "mae_forecast", "mae_naive", "mase") %in%
                              names(hs$mase)))
  testthat::expect_false("peel" %in% names(hs$mase))
  testthat::expect_gt(nrow(hs$mase), 0)

  # NOTE on what is worth asserting here. `years_ahead == year - (endyr - peel)`,
  # `year > endyr - peel` and `mase == mae_forecast/mae_naive` all restate the
  # lines that computed them and cannot fail; `all(is.finite(x) | is.na(x))`
  # passes for any numeric vector. They are omitted deliberately. What CAN fail
  # is whether the three series are the right quantities from the right objects,
  # which is what follows.
  endyr <- fit$data_list$endyr
  yr1 <- hs$by_year[hs$by_year$peel == 1, ]
  sp1 <- yr1[yr1$species == fit$data_list$spnames[1], ]
  col <- sp1$year - fit$data_list$styr + 1L

  # reference: the FULL model's estimate at those years.
  testthat::expect_equal(sp1$reference,
                         as.numeric(fit$quantities$ssb[1, col]),
                         tolerance = 1e-8)

  # forecast: the PEEL's own estimate at those years, and NOT the full model's.
  # Nothing tested this, so a forecast column silently taken from the wrong
  # object would have passed everything above.
  peel1 <- retro_peel_for(hs, fit, peel = 1)
  testthat::expect_equal(sp1$forecast,
                         as.numeric(peel1$quantities$ssb[1, col]),
                         tolerance = 1e-8)
  testthat::expect_false(isTRUE(all.equal(sp1$forecast, sp1$reference,
                                          tolerance = 1e-6)))

  # naive: one value held flat, and specifically the peel's TERMINAL estimate.
  testthat::expect_equal(length(unique(sp1$naive)), 1L)
  testthat::expect_equal(unique(sp1$naive),
                         as.numeric(peel1$quantities$ssb[1, (endyr - 1) -
                                                           fit$data_list$styr + 1L]),
                         tolerance = 1e-8)
})

testthat::test_that("hindcast_skill can score the held-out index instead", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = 0, random_rec = TRUE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
  # Record each peel's reported mean as it reaches the median conversion; the
  # real conversion still runs, so the scored result is unchanged.
  means <- list()
  real_median <- Rceattle:::.hat_as_median
  testthat::local_mocked_bindings(
    .hat_as_median = function(hat, sigma, ba_obs, log_scale) {
      means[[length(means) + 1L]] <<- hat
      real_median(hat, sigma, ba_obs, log_scale)
    }, .package = "Rceattle")
  hs <- suppressWarnings(suppressMessages(
    Rceattle::hindcast_skill(fit, peels = 2, reference = "observed",
                             cores = 1, getsd = FALSE)))

  # [[ not $, and the column is `Object`: a $ read of a column that has been
  # renamed away returns NULL, grepl(NULL) returns logical(0), and the subset is
  # then 0 rows -- which turns every all() assertion below into a vacuous TRUE
  # rather than a failure. That is exactly what happened when `quantity` became
  # `Object`: this block reported PASS while checking nothing.
  idx <- hs$mase[grepl("^index_fleet", hs$mase[["Object"]]), ]
  testthat::expect_gt(nrow(idx), 0)

  # The reference here must be the OBSERVED index, not a model quantity -- every
  # scored value has to appear in index_data at that fleet and year.
  by <- hs$by_year[grepl("^index_fleet", hs$by_year[["Object"]]), ]
  obs <- fit$data_list$index_data
  # Non-vacuity first. all(logical(0)) is TRUE, so without this an empty `by`
  # passes the fleet-matching assertion below instead of failing it.
  testthat::expect_gt(nrow(by), 0)
  # Match on FLEET and year, not year alone: with several fleets reporting in
  # the same year, a year-only match is satisfied by any of them and would pass
  # even if the reference were taken from the wrong fleet entirely.
  hit <- vapply(seq_len(min(nrow(by), 20L)), function(i) {
    flt <- as.integer(sub("^index_fleet_", "", by[["Object"]][i]))
    rows <- obs$Year == by$year[i] & obs$Fleet_code == flt
    any(rows) && any(abs(obs$Observation[rows] - by$reference[i]) < 1e-8)
  }, logical(1))
  testthat::expect_length(hit, min(nrow(by), 20L))
  testthat::expect_true(all(hit))

  # And the peel must be predicting years it did not fit.
  testthat::expect_true(all(by[["Forecast year"]] >= 1))

  # The forecast is the predicted MEDIAN: BS2017SS's surveys are lognormal with
  # fixed sds, so each is its peel's mean times exp(-Log_sd^2 / 2) (bias_adjust_obs
  # = 1), with Log_sd read from the data rather than from the function.
  testthat::expect_length(means, 2L)
  r <- mapply(function(o, y) which(obs$Fleet_code == o & obs$Year == y),
              as.integer(sub("^index_fleet_", "", by[["Object"]])), by$year)
  med <- exp(-obs$Log_sd[r]^2 / 2)
  for (p in unique(by$peel)) {
    i <- by$peel == p
    hit <- vapply(means, function(m) isTRUE(all.equal(by$forecast[i], m[r[i]] * med[i],
                                                      tolerance = 1e-12)), logical(1))
    testthat::expect_true(any(hit), info = paste("peel", p))
    testthat::expect_true(all(by$forecast[i] < means[[which(hit)[1]]][r[i]]))
  }
})

testthat::test_that("hindcast_skill takes the median from the fit's bias_adjust_obs", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  # Under fit_control(bias_adjust_obs = FALSE) the likelihood is centred on
  # log(index_hat), so the reported value already is the median and is scored as is.
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = 0, random_rec = TRUE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0,
                                        bias_adjust_obs = FALSE))))
  means <- list()
  real_median <- Rceattle:::.hat_as_median
  testthat::local_mocked_bindings(
    .hat_as_median = function(hat, sigma, ba_obs, log_scale) {
      means[[length(means) + 1L]] <<- hat
      real_median(hat, sigma, ba_obs, log_scale)
    }, .package = "Rceattle")
  hs <- suppressWarnings(suppressMessages(
    Rceattle::hindcast_skill(fit, peels = 2, reference = "observed",
                             cores = 1, getsd = FALSE)))
  by <- hs$by_year[grepl("^index_fleet", hs$by_year[["Object"]]), ]
  obs <- fit$data_list$index_data
  testthat::expect_gt(nrow(by), 0)
  testthat::expect_length(means, 2L)
  r <- mapply(function(o, y) which(obs$Fleet_code == o & obs$Year == y),
              as.integer(sub("^index_fleet_", "", by[["Object"]])), by$year)
  for (p in unique(by$peel)) {
    i <- by$peel == p
    hit <- vapply(means, function(m) identical(by$forecast[i], m[r[i]]), logical(1))
    testthat::expect_true(any(hit), info = paste("peel", p))
  }
})


# MASE averages ACROSS PEELS at a fixed horizon -- Kell et al. (2021) eq. 5.
#
# Regression. It was grouped the other way round: by peel, averaging over
# horizons, which is the transpose of the published definition. Two consequences
# made it look like a property of the statistic rather than a defect. It reported
# one MASE per peel where Kell's Table 2 reports one per steps-ahead; and it made
# the h = 1 denominator a SINGLE |naive - reference| instead of a mean over
# peels, which is where the old "a one-year MASE is dominated by its own
# denominator, compare the deeper peels" advice came from. That advice was
# compensating for the grouping.
#
# Pinned by recomputing eq. 5 from $by_year rather than by restating the line
# that produced $mase: the point is the GROUPING, so the check has to do its own
# grouping.
testthat::test_that("MASE is computed per horizon, averaging over peels", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data()
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, file = NULL, estimateMode = 1,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  hs <- suppressMessages(suppressWarnings(
    Rceattle::hindcast_skill(fit, peels = 2:4, quantity = "R",
                             reference = "model", cores = 1, getsd = FALSE)))

  by <- hs$by_year
  testthat::expect_gt(length(unique(by$peel)), 1L)   # or there is nothing to average

  # One row per (Forecast year, species, Object) -- NOT per peel.
  key <- paste(by[["Forecast year"]], by$species, by$Object)
  testthat::expect_equal(nrow(hs$mase), length(unique(key)))

  # `N` is the number of peels reaching that horizon, and it FALLS as the
  # horizon grows, because a horizon h needs a peel at least h deep.
  for (i in seq_len(nrow(hs$mase))) {
    z <- by[by[["Forecast year"]] == hs$mase[["Forecast year"]][i] &
            by$species              == hs$mase$species[i] &
            by$Object               == hs$mase$Object[i], ]
    # PEELS, not rows: a peel contributing several rows at one horizon -- which
    # the observed path allows, two index rows in a year at different months --
    # is averaged within the peel first, so `N` != nrow(z) there.
    testthat::expect_equal(hs$mase$N[i], length(unique(z$peel)))
    # eq. 5, recomputed with its own within-peel average
    af <- tapply(abs(z$forecast - z$reference), z$peel, mean)
    an <- tapply(abs(z$naive    - z$reference), z$peel, mean)
    testthat::expect_equal(hs$mase$mase[i], mean(af) / mean(an))
  }
  m1 <- hs$mase[hs$mase[["Forecast year"]] == min(hs$mase[["Forecast year"]]), ]
  mN <- hs$mase[hs$mase[["Forecast year"]] == max(hs$mase[["Forecast year"]]), ]
  testthat::expect_gt(m1$N[1], mN$N[1])
})


# A peel contributing TWO rows at one horizon is weighted once, not twice.
#
# Kell et al. (2021) eq. 5 sums one term per t, each t a different peel's
# terminal year. On reference = "model" that is automatic -- one forecast per
# peel per horizon. On reference = "observed" it is not: a fleet may legally
# carry two index rows in the same Year at different Months, because
# data_check() tests duplicates on (Fleet_code, Year, Month), and both land in
# the same `Forecast year` from the same peel.
#
# $mase therefore averages within a peel first and then across peels. This is
# the ONLY branch where that differs from a pooled-row mean, and no bundled
# dataset reaches it -- BS2017SS has zero (Fleet_code, Year) duplicates and
# every Month is 6 -- so the fixture is constructed here. Without it, `N`
# could silently revert to a row count and nothing in the suite would notice.
testthat::test_that("a peel with two rows at one horizon is weighted once", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  data(BS2017SS, package = "Rceattle")
  d <- BS2017SS
  # Duplicate fleet 4's recent rows at a second month, perturbed so a
  # row-weighted mean and a peel-weighted mean cannot coincide by accident.
  extra <- d$index_data[d$index_data$Fleet_code == 4 &
                        d$index_data$Year >= 2014, ]
  testthat::skip_if(nrow(extra) == 0L, "fixture fleet has no recent rows")
  extra$Month <- 9
  extra$Observation <- extra$Observation * 1.35
  d$index_data <- rbind(d$index_data, extra)
  d$index_data <- d$index_data[order(d$index_data$Fleet_code,
                                     d$index_data$Year), ]
  # The shape must be legal, or this is testing a model that cannot be built.
  testthat::expect_equal(sum(duplicated(
    d$index_data[c("Fleet_code", "Year", "Month")])), 0L)
  testthat::expect_gt(sum(duplicated(d$index_data[c("Fleet_code", "Year")])), 0L)

  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, file = NULL, estimateMode = 1,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  hs <- suppressMessages(suppressWarnings(
    Rceattle::hindcast_skill(fit, peels = 2:3, reference = "observed",
                             cores = 1, getsd = FALSE)))

  by <- hs$by_year
  # [[ and `Object`: with a $ read of a renamed column this was NULL == "..." ->
  # logical(0) -> 0 rows, and the test SKIPPED ITSELF with a reason that was not
  # true. The header above says why this test exists; a silent skip retires it.
  dup <- by[by[["Object"]] == "index_fleet_4", ]
  testthat::skip_if(nrow(dup) == 0L, "fleet 4 contributed no scored rows")

  for (i in seq_len(nrow(hs$mase))) {
    z <- by[by[["Forecast year"]] == hs$mase[["Forecast year"]][i] &
            by$species              == hs$mase$species[i] &
            by$Object               == hs$mase$Object[i], ]
    af <- tapply(abs(z$forecast - z$reference), z$peel, mean)
    an <- tapply(abs(z$naive    - z$reference), z$peel, mean)
    # A peel needs BOTH a forecast and a baseline to be scored: on the observed
    # path the naive value is the observation at the peel's terminal year, and a
    # fleet with no observation there has no scale for the ratio, so that peel
    # drops out even though its forecast is fine. `N` counts the peels that
    # actually contributed, which is what the column claims to mean.
    ok <- !is.na(af) & !is.na(an)
    # PEELS, not rows -- this is the assertion the duplicate shape exists for.
    testthat::expect_equal(hs$mase$N[i], sum(ok))
    if (any(ok)) {
      testthat::expect_equal(hs$mase$mae_forecast[i], mean(af[ok]))
      testthat::expect_equal(hs$mase$mae_naive[i],    mean(an[ok]))
    } else {
      # A mean of nothing is NA, not 0 -- a 0 numerator would read as perfect
      # skill and a 0 denominator as infinitely bad.
      testthat::expect_true(is.na(hs$mase$mae_forecast[i]))
      testthat::expect_true(is.na(hs$mase$mae_naive[i]))
      testthat::expect_true(is.na(hs$mase$mase[i]))
    }
  }

  # And at least one group really does have more rows than peels, or the
  # fixture has drifted and the test is passing on the easy case.
  hit <- vapply(seq_len(nrow(hs$mase)), function(i) {
    z <- by[by[["Forecast year"]] == hs$mase[["Forecast year"]][i] &
            by$species              == hs$mase$species[i] &
            by$Object               == hs$mase$Object[i], ]
    nrow(z) > length(unique(z$peel))
  }, logical(1))
  testthat::expect_true(any(hit))
})

testthat::test_that("an analytical index sd comes from the peel, not the rebuild", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  # Walters-Ludwig solves the sd from every row a model fits; the rebuild that
  # predicts the held-out years fits them too, so the median must use the peel's.
  d <- Rceattle::BS2017SS
  srv <- which(d$fleet_control$Fleet_name == "BT_Pollock")
  d$fleet_control$Estimate_index_sd[srv] <- 2
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = TRUE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  sds <- list()
  real_median <- Rceattle:::.hat_as_median
  testthat::local_mocked_bindings(
    .hat_as_median = function(hat, sigma, ba_obs, log_scale) {
      sds[[length(sds) + 1L]] <<- sigma
      real_median(hat, sigma, ba_obs, log_scale)
    }, .package = "Rceattle")
  hs <- suppressWarnings(suppressMessages(
    Rceattle::hindcast_skill(fit, peels = 1, reference = "observed",
                             cores = 1, getsd = FALSE)))
  testthat::expect_length(sds, 1L)

  peel <- retro_peel_for(hs, fit, peel = 1)
  code <- d$fleet_control$Fleet_code[srv]
  peel_sd <- peel$quantities$index_sd[match(code, peel$data_list$index_data$Fleet_code)]
  full_sd <- fit$quantities$index_sd[match(code, fit$data_list$index_data$Fleet_code)]
  # Non-vacuity: the analytical sd does depend on which years were fitted.
  testthat::expect_false(isTRUE(all.equal(peel_sd, full_sd, tolerance = 1e-8)))
  rows <- fit$data_list$index_data$Fleet_code == code
  testthat::expect_true(all(abs(sds[[1]][rows] - peel_sd) < 1e-12))
})

testthat::test_that("hindcast_skill scores only rows the likelihood fits", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  # Rows on an Off fleet and NA observations are not fitted, so neither is a
  # forecast to score.
  d <- make_test_data(nyrs = 20, seed = 1)
  d$fleet_control$Catchability[1] <- "Estimated"
  off <- d$fleet_control[1, ]
  off$Fleet_name <- "Survey_off"; off$Fleet_code <- 3L; off$Fleet_type <- "Off"
  off$Selectivity_index <- 3L; off$Catchability_index <- 3L
  d$fleet_control <- rbind(d$fleet_control, off)
  i3 <- d$index_data[d$index_data$Fleet_code == 1, ]
  i3$Fleet_code <- 3L; i3$Fleet_name <- "Survey_off"
  d$index_data <- rbind(d$index_data, i3)
  # Unobserved at endyr - 1 and endyr - 2: each peel's terminal year, so neither
  # has a persistence baseline, and endyr - 1 is also held out from peel 2.
  srv <- d$index_data$Fleet_code == 1
  d$index_data$Observation[srv & d$index_data$Year %in% (d$endyr - 1:2)] <- NA
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  ws <- character(0)
  hs <- withCallingHandlers(suppressMessages(
    Rceattle::hindcast_skill(fit, peels = 2, reference = "observed",
                             cores = 1, getsd = FALSE)),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  by <- hs$by_year
  testthat::expect_gt(nrow(by), 0)
  # The unobserved terminal year is reported as a missing baseline.
  testthat::expect_true(any(grepl("observation in the peel's terminal year", ws, fixed = TRUE)))
  testthat::expect_identical(unique(by[["Object"]]), "index_fleet_1")
  testthat::expect_false(any(by$year == d$endyr - 1))
  testthat::expect_false(anyNA(by$reference))
  testthat::expect_gt(sum(by$peel == 2), 0)
  testthat::expect_true(all(is.na(by$naive[by$peel == 2])))
})

# Catchability = "Environmental": log q_y = q_mu + beta * temp_y. A peel is fitted
# with post-peel temp blanked (mean-filled), so its forecast of a held-out year
# must not move with that year's true temp. Through 5.66.1.9004 the forecasts
# tracked it, 82 to 144 against observations of 82 to 144.
testthat::test_that("the held-out prediction does not see post-peel covariates", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 20, seed = 1)
  yrs <- d$styr:d$endyr
  set.seed(3)
  d$env_data <- data.frame(Year = yrs, temp = stats::rnorm(length(yrs)))
  srv <- d$index_data$Fleet_code == 1
  d$index_data$Observation[srv] <- 100 *
    exp(0.3 * d$env_data$temp[match(d$index_data$Year[srv], yrs)])
  d$fleet_control$Catchability[1] <- "Environmental"
  d$fleet_control$Time_varying_q[1] <- 1
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  testthat::expect_gt(abs(fit$estimated_params$index_q_beta[1]), 0.1)   # temp matters
  hs <- suppressWarnings(suppressMessages(
    Rceattle::hindcast_skill(fit, peels = 3, reference = "observed",
                             cores = 1, getsd = FALSE)))
  by <- hs$by_year[hs$by_year$peel == 3, ]
  testthat::expect_equal(nrow(by), 3L)
  # At the peel's own covariate the forecasts differ only by population dynamics.
  testthat::expect_lt(max(by$forecast) / min(by$forecast) - 1, 0.01)
})

# The observed score predicts held-out survey years from a rebuild of each peel,
# so the rebuild must be the peel. Through 5.66.1.9004 it took the Ianelli
# stock-recruit windows and (under predation) post-peel diet from the full
# model: held-out biomass up to 2.8% and 2.1% off the peel's own.
testthat::test_that("the held-out rebuild reproduces each peel's own projection", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  srr <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0, fit_control = fc,
    recFun = Rceattle::build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt"))))
  ss <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017MS, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0, fit_control = fc)))
  msm <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017MS, inits = ss$estimated_params, file = NULL,
    estimateMode = 0, random_rec = FALSE, msmMode = 1, niter = 3, fit_control = fc)))

  rets <- lapply(list(srr, msm), function(fit) suppressWarnings(suppressMessages(
    Rceattle::retrospective(fit, peels = 2, cores = 1, getsd = FALSE,
                            forecast_rec = "model"))))
  # Record every rebuild hindcast_skill() makes; the real one still runs.
  rebuilds <- list()
  real_refit <- Rceattle:::.refit_like
  testthat::local_mocked_bindings(.refit_like = function(...) {
    out <- real_refit(...)
    rebuilds[[length(rebuilds) + 1L]] <<- out
    out
  }, .package = "Rceattle")
  for (k in 1:2) {
    fit <- list(srr, msm)[[k]]; r <- rets[[k]]
    rebuilds <- list()
    suppressWarnings(suppressMessages(
      Rceattle::hindcast_skill(fit, retro = r, reference = "observed")))
    peels <- r$Rceattle_list[-length(r$Rceattle_list)]
    testthat::expect_length(rebuilds, length(peels))
    for (p in rebuilds) {
      ep <- p$data_list$endyr_peel
      peel <- peels[[which(vapply(peels, function(m) m$data_list$endyr_peel, 0) == ep)]]
      cols <- (ep + 1):fit$data_list$endyr - fit$data_list$styr + 1
      testthat::expect_equal(p$quantities$biomass[, cols], peel$quantities$biomass[, cols],
                             tolerance = 1e-10)
    }
  }
})
