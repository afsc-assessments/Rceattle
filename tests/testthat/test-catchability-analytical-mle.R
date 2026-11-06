# Catchability = "Analytical" must be the q that maximizes the index likelihood
# the model evaluates, on every family it accepts (formulas in
# vignette("model-parameterizations")); TruncatedNormal has no closed-form q and
# is refused. Up to 5.66.11 every family took the geometric mean of obs / pred:
# on this fixture (complete series) the lognormal objective sat 0.065 (sd by
# year) to 0.24 (sd estimated or analytical) above its minimum over q, with q
# 2-8% low, and q was 51% high under Normal and 88% under MVN. No golden model
# uses analytical q, so this is the regression net.
#
# The check needs no fit: at the starting values, the analytical-q objective
# must equal the minimum over log q (and log sd) of the same model with those
# estimated. The missing year checks that q uses only fitted rows.

.aq_build <- function(q, sd_form, bias = TRUE, family = "Lognormal") {
  d <- make_test_data(nyrs = 12)
  d$index_data$Observation <- 100 * exp(c(0.3, -0.2, 0.1, 0.4, -0.3, 0, 0.2,
                                          -0.1, 0.25, -0.15, 0.05, -0.35))
  d$index_data$Log_sd <- seq(0.1, 0.5, length.out = 12)  # sd varies by year
  d$fleet_control$Index_sd[1] <- 0.3
  if (family %in% c("Normal", "TruncatedNormal")) {       # an absolute sd
    d$index_data$Log_sd <- seq(10, 30, length.out = 12)
    d$fleet_control$Index_sd[1] <- 20
  }
  if (family %in% c("MVN", "MVNORM")) {                   # heteroscedastic, rho = 0.3
    sds <- seq(10, 30, length.out = 12)
    rho <- matrix(0.3, 12, 12); diag(rho) <- 1
    d$index_cov <- list(Survey = diag(sds) %*% rho %*% diag(sds))
  } else {
    # A year without a survey. Not on MVN, whose R-side covariance alignment
    # does not yet accept one (Rceattle-dev CLEANUP_BACKLOG.md).
    d$index_data$Observation[5] <- NA
  }
  d$fleet_control$Catchability_init[1] <- 0.8
  d$fleet_control$Catchability[1] <- q
  d$fleet_control$Estimate_index_sd[1] <- sd_form
  d$fleet_control$Index_distribution <- "Lognormal"
  d$fleet_control$Index_distribution[1] <- family         # the survey
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(bias_adjust_obs = bias, verbose = 0))))$obj
}

# Minimum of obj$fn over the named parameters, the rest held at their inits.
.aq_min <- function(obj, free) {
  p0 <- obj$par
  if (!length(free)) return(list(objective = obj$fn(p0), par = p0))
  idx <- which(names(p0) %in% free)
  testthat::expect_length(idx, length(free))
  f <- function(x) { p <- p0; p[idx] <- x; obj$fn(p) }
  g <- function(x) { p <- p0; p[idx] <- x; obj$gr(p)[idx] }
  o <- nlminb(p0[idx], f, g, control = list(rel.tol = 1e-14, x.tol = 1e-14))
  p0[idx] <- o$par
  list(objective = o$objective, par = p0)
}

# Analytical q against q estimated, at the same sd form.
.aq_expect_mle <- function(sd_form, bias, family, ref = NULL) {
  info <- paste0(family, ", bias_adjust_obs = ", bias, ", Estimate_index_sd = ", sd_form)
  sd_free <- if (sd_form == "Estimated") "index_log_sd" else character(0)
  oA <- .aq_build("Analytical", sd_form, bias, family)
  mA <- .aq_min(oA, sd_free)
  mE <- .aq_min(.aq_build("Estimated", sd_form, bias, family), c("index_log_q", sd_free))
  testthat::expect_true(is.finite(mA$objective), label = info)
  # Absolute: the objective is ~1e5 here, so a relative tolerance is loose.
  testthat::expect_lt(abs(mA$objective - mE$objective), 1e-7, label = info)
  testthat::expect_equal(oA$report(mA$par)$index_q[1, 1],
                         exp(mE$par[["index_log_q"]]), tolerance = 1e-4, info = info)
  if (!is.null(ref)) testthat::expect_lt(abs(mA$objective - ref), 1e-7, label = info)
}

test_that("lognormal analytical q reaches the likelihood's optimum for every sd form", {
  testthat::skip_if_not_installed("TMB")
  for (bias in c(TRUE, FALSE)) {
    # Both quantities free: the reference objective when sd is not input.
    ref <- .aq_min(.aq_build("Estimated", "Estimated", bias),
                   c("index_log_q", "index_log_sd"))$objective
    .aq_expect_mle("Fixed", bias, "Lognormal")
    .aq_expect_mle("Estimated", bias, "Lognormal", ref)
    .aq_expect_mle("Analytical", bias, "Lognormal", ref)
  }
})

test_that("Normal analytical q reaches the likelihood's optimum", {
  testthat::skip_if_not_installed("TMB")
  .aq_expect_mle("Fixed", TRUE, "Normal")
  .aq_expect_mle("Estimated", TRUE, "Normal")
})

test_that("MVN and MVNORM analytical q reach the likelihood's optimum", {
  testthat::skip_if_not_installed("TMB")
  .aq_expect_mle("Fixed", TRUE, "MVN")
  .aq_expect_mle("Fixed", TRUE, "MVNORM")
})

test_that("analytical q uses only the fitted rows of its own fleet (fishery CPUE, fleet 2)", {
  testthat::skip_if_not_installed("TMB")
  # The survey's q is fixed; the analytical fleet is the fishery, whose index
  # also carries projection and negative-year rows with extreme values.
  cpue <- function(q) {
    d <- make_test_data(nyrs = 12)
    fi <- d$index_data
    fi$Fleet_name <- "Fishery"; fi$Fleet_code <- 2
    fi$Observation <- 30 * exp(c(-0.35, 0.05, -0.15, 0.25, -0.1, 0.2, 0, -0.3,
                                 0.4, 0.1, -0.2, 0.3))
    fi$Observation[7] <- NA
    fi$Log_sd <- seq(0.1, 0.5, length.out = 12)
    ex <- fi[1:4, ]
    ex$Year <- c(13, 14, -3, -4); ex$Observation <- c(1e4, 1e-3, 5e3, 2e-3)
    d$index_data <- rbind(d$index_data, fi, ex)
    fc <- d$fleet_control
    fc$Catchability_index[2] <- 2; fc$Catchability[2] <- q
    fc$Catchability_init[2] <- 0.5; fc$Catchability_prior_sd[2] <- 0.2
    fc$Time_varying_q[2] <- 0; fc$Time_varying_q_sd[2] <- 1
    fc$Estimate_index_sd[2] <- "Analytical"; fc$Index_sd[2] <- 0.3
    d$fleet_control <- fc
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, estimateMode = 3, msmMode = 0, random_rec = FALSE,
      fit_control = Rceattle::fit_control(verbose = 0))))$obj
  }
  oA <- cpue("Analytical")
  mE <- .aq_min(cpue("Estimated"), "index_log_q")
  testthat::expect_lt(abs(oA$fn(oA$par) - mE$objective), 1e-7)
  testthat::expect_equal(oA$report(oA$par)$index_q[2, 1],
                         exp(mE$par[["index_log_q"]]), tolerance = 1e-4)
})

test_that("analytical q is refused on a TruncatedNormal index, unless the fleet is Off", {
  testthat::skip_if_not_installed("TMB")
  testthat::expect_error(.aq_build("Analytical", "Fixed", TRUE, "TruncatedNormal"),
                         "TruncatedNormal")
  d <- make_test_data(nyrs = 12)
  d$fleet_control$Catchability[1] <- 3                    # integer codes
  d$fleet_control$Index_distribution <- c(4, 0)
  testthat::expect_error(Rceattle:::data_check(d), "TruncatedNormal")
  d$fleet_control$Fleet_type[1] <- "Off"
  testthat::expect_no_error(suppressWarnings(Rceattle:::.check_analytical_q_family(d)))
  testthat::expect_length(Rceattle:::.check_analytical_q_family(d), 0)
})

test_that("AnalyticalArith skips a missing year", {
  testthat::skip_if_not_installed("TMB")
  o <- .aq_build("AnalyticalArith", "Fixed")
  r <- o$report(o$par)
  obs <- o$env$data$index_obs[, 1]
  fitted <- which(obs > 0 & o$env$data$index_ctl[, 1] == 1)
  testthat::expect_true(is.finite(o$fn(o$par)))
  # q = sum(obs) / sum(pred) over the observed years, so the scaled sums agree.
  testthat::expect_equal(sum(r$index_hat[fitted]), sum(obs[fitted]), tolerance = 1e-12)
})
