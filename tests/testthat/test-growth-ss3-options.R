# build_growth()'s SS3 growth options -- the population length grid
# (pop_lengths), CV-form growth variability (sd_form), the plus-group mean
# length (plus_group_length) -- and maturity-at-length (the L50_mat_len /
# slope_mat_len control columns). Each is checked against an independent R
# calculation from the reported mean length-at-age, following the SS3 3.30
# source (SS_ALK.tpl, SS_biofxn.tpl): normal length-at-age on the population
# bins with the first bin a minus group below the second edge and the last a
# plus group above the last edge, weight at population-bin midpoints, and
# fecundity = sum over length of P(l | a) x weight(l) x maturity(l).

testthat::skip_on_cran()

ss3_opt_build <- function(d, ...) {
  suppressWarnings(suppressMessages(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    growthFun = build_growth(fun = "vonBertalanffy", ...),
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
}

# Age-length key on lower edges `edges`, SS3's edge convention.
ss3_alk <- function(mu, sd, edges) {
  n <- length(edges)
  cdf <- stats::pnorm((edges - mu) / sd)
  p <- numeric(n)
  p[1] <- cdf[2]
  p[n] <- 1 - cdf[n]
  p[2:(n - 1)] <- cdf[3:n] - cdf[2:(n - 1)]
  p
}
ss3_mid <- function(edges) c((edges[-1] + edges[-length(edges)]) / 2,
                             edges[length(edges)] + diff(utils::tail(edges, 2)) / 2)

# SD of length-at-age for species 1, sex 1, from the fitted parameters.
ss3_sd <- function(m, age_idx, len, plus, form = "SD", fracyr = 0) {
  lsd  <- m$estimated_params$growth_log_sd[1, 1, ]
  gp   <- m$quantities$growth_parameters[1, 1, 1, ]
  cur  <- age_idx - 1 + m$data_list$minage[1] + fracyr
  v <- if (cur <= m$data_list$growth_age_L1[1]) exp(lsd[1]) else if (plus) exp(lsd[2]) else
    exp(lsd[1]) + (exp(lsd[2]) - exp(lsd[1])) / (gp[3] - gp[2]) * (len - gp[2])
  if (form == "CV") v * len else v
}

set.seed(11)
d <- make_msm_test_data()$data_list
edges <- sort(unique(d$caal_data$Length[d$caal_data$Species == 1]))
fine  <- sort(unique(c(edges, seq(min(edges), max(edges), by = 1))))
alpha <- d$alpha_wt_len[1]; beta <- d$beta_wt_len[1]

testthat::test_that("a population grid equal to the data bins changes nothing", {
  m0 <- ss3_opt_build(d)
  m1 <- ss3_opt_build(d, pop_lengths = list(edges))
  q0 <- m0$quantities; q1 <- m1$quantities
  testthat::expect_equal(q1$length_hat, q0$length_hat, tolerance = 1e-12)
  testthat::expect_equal(q1$weight_hat, q0$weight_hat, tolerance = 1e-12)
  testthat::expect_equal(q1$growth_matrix, q0$growth_matrix, tolerance = 1e-12)
  testthat::expect_equal(q1$ssb, q0$ssb, tolerance = 1e-12)
})

testthat::test_that("the default spawning output is spawning weight x maturity, and SSB sums it", {
  m <- ss3_opt_build(d)
  q <- m$quantities
  na <- d$nages[1]
  sm <- m$data_list$spawn_month[1]
  ssb1 <- sum(q$N_at_age[1, 1, 1:na, 1] * exp(-q$Z_at_age[1, 1, 1:na, 1] * sm / 12) * q$spawn_output[1, 1:na, 1])
  testthat::expect_equal(unname(q$ssb[1, 1]), ssb1, tolerance = 1e-10)
  mat <- as.numeric(d$maturity[d$maturity$Species == 1, grep("^Age", names(d$maturity))])[1:na]
  sr  <- as.numeric(d$sex_ratio[d$sex_ratio$Species == 1, grep("^Age", names(d$sex_ratio))])[1:na]
  testthat::expect_equal(unname(q$spawn_output[1, 1:na, 1]),
                         unname(q$weight_hat[2, 1, 1:na, 1] * mat * sr), tolerance = 1e-12)
})

testthat::test_that("the age-length key and weight-at-age integrate over the population grid", {
  m <- ss3_opt_build(d, pop_lengths = list(fine))
  q <- m$quantities
  na <- d$nages[1]
  for (a in seq_len(na)) {
    mu <- q$length_hat[1, 1, a, 1]
    p  <- ss3_alk(mu, ss3_sd(m, a, mu, a == na), fine)
    # growth_matrix is stored at POPULATION resolution: SS3 forms sel(L) * P(L|age)
    # on the population grid and bins the result, so the key itself is not
    # pre-aggregated. The summing into data bins happens later, in pred_CAAL,
    # which is what the length compositions are built from.
    testthat::expect_equal(unname(q$growth_matrix[1, 1, a, seq_along(fine), 1]),
                           as.numeric(p), tolerance = 1e-10)
    # ... and that aggregation is still exactly the data-bin key.
    to_data <- findInterval(fine + 1e-8, edges)
    testthat::expect_equal(
      as.numeric(tapply(unname(q$growth_matrix[1, 1, a, seq_along(fine), 1]),
                        factor(to_data, levels = seq_along(edges)), sum)),
      as.numeric(tapply(p, factor(to_data, levels = seq_along(edges)), sum)),
      tolerance = 1e-10)
    testthat::expect_equal(unname(q$weight_hat[1, 1, a, 1]), sum(p * alpha * ss3_mid(fine)^beta),
                           tolerance = 1e-10)
  }
  testthat::expect_equal(sum(q$growth_matrix[1, 1, 1, , 1]), 1, tolerance = 1e-12)
})

testthat::test_that("sd_form = 'CV' makes the SD a CV times mean length", {
  m <- ss3_opt_build(d, sd_form = "CV")
  q <- m$quantities
  na <- d$nages[1]
  testthat::expect_equal(unname(m$estimated_params$growth_log_sd[1, 1, ]), rep(log(0.1), 2))
  for (a in seq_len(na)) {
    mu <- q$length_hat[1, 1, a, 1]
    p  <- ss3_alk(mu, ss3_sd(m, a, mu, a == na, form = "CV"), edges)
    testthat::expect_equal(unname(q$weight_hat[1, 1, a, 1]), sum(p * alpha * ss3_mid(edges)^beta),
                           tolerance = 1e-10)
  }
})

testthat::test_that("plus_group_length follows SS3's Linf_decay forms", {
  na <- d$nages[1]
  A  <- na - 1 + d$minage[1]                                 # oldest age
  base <- function(m) {
    gp <- m$quantities$growth_parameters[1, 1, 1, ]           # K, L1, Linf, m
    list(gp = gp, LA = unname(gp[3] + (gp[2] - gp[3]) * exp(-gp[1] * (A - m$data_list$growth_age_L1[1]))))
  }

  m <- ss3_opt_build(d, plus_group_length = "none")
  b <- base(m)
  testthat::expect_equal(unname(m$quantities$length_hat[1, 1, na, 1]), b$LA, tolerance = 1e-10)

  m <- ss3_opt_build(d, plus_group_length = "SS3.24")
  b <- base(m); a <- 0:A; w <- exp(-0.2 * a)
  testthat::expect_equal(unname(m$quantities$length_hat[1, 1, na, 1]),
                         sum(w * (b$LA + (a / A) * (b$gp[3] - b$LA))) / sum(w), tolerance = 1e-10)

  m <- ss3_opt_build(d, plus_group_length = "decay", plus_group_decay = 0.4)
  b <- base(m)
  size <- b$LA; num <- b$LA; den <- 1; wt <- 1
  for (i in seq_len(2 * A)) {
    wt <- wt * exp(-0.4); size <- size + (b$gp[3] - size) * (1 - exp(-b$gp[1]))
    num <- num + wt * size; den <- den + wt
  }
  testthat::expect_equal(unname(m$quantities$length_hat[1, 1, na, 1]), unname(num / den), tolerance = 1e-10)

  # Under an SS3 form the plus group grows within the year; under "M1" it does not.
  d3 <- d
  fl <- which(d3$fleet_control$Species == 1)[1]
  d3$fleet_control$Month[fl] <- 6
  slot <- 2 * d3$nspp + fl
  m <- ss3_opt_build(d3, plus_group_length = "none")
  testthat::expect_gt(m$quantities$length_hat[slot, 1, na, 1], m$quantities$length_hat[1, 1, na, 1])
  m1 <- ss3_opt_build(d3)
  testthat::expect_equal(m1$quantities$length_hat[slot, 1, na, 1], m1$quantities$length_hat[1, 1, na, 1])
})

testthat::test_that("maturity-at-length integrates maturity x weight over length at spawning", {
  d2 <- d
  d2$L50_mat_len   <- rep(NA_real_, d$nspp); d2$L50_mat_len[1]   <- 40
  d2$slope_mat_len <- rep(NA_real_, d$nspp); d2$slope_mat_len[1] <- 0.2
  m <- ss3_opt_build(d2, pop_lengths = list(fine))
  q <- m$quantities
  na <- d$nages[1]
  frac <- m$data_list$spawn_month[1] / 12
  sr <- as.numeric(d$sex_ratio[d$sex_ratio$Species == 1, grep("^Age", names(d$sex_ratio))])[1:na]
  mid <- ss3_mid(fine)
  for (a in seq_len(na)) {
    mu <- q$length_hat[2, 1, a, 1]                          # length at spawning
    p  <- ss3_alk(mu, ss3_sd(m, a, mu, a == na, fracyr = frac), fine)
    testthat::expect_equal(unname(q$spawn_output[1, a, 1]),
                           sum(p * alpha * mid^beta / (1 + exp(-0.2 * (mid - 40)))) * sr[a],
                           tolerance = 1e-10)
  }
  # Species without the columns keep the age-based maturity sheet
  m0 <- ss3_opt_build(d, pop_lengths = list(fine))
  testthat::expect_equal(q$spawn_output[2, , ], m0$quantities$spawn_output[2, , ], tolerance = 1e-12)
})

testthat::test_that("maturity-at-length and the new build_growth() arguments are validated", {
  d2 <- d
  d2$L50_mat_len <- rep(40, d$nspp)
  testthat::expect_error(data_check(d2), "Set both L50_mat_len and slope_mat_len")
  d2$slope_mat_len <- rep(-0.2, d$nspp)
  testthat::expect_error(data_check(d2), "must be positive")
  testthat::expect_error(build_growth("vonBertalanffy", plus_group_length = "decay"),
                         "positive `plus_group_decay`")
  testthat::expect_error(build_growth("vonBertalanffy", sd_form = "var"))
  testthat::expect_error(build_growth("vonBertalanffy", pop_lengths = c(10, 5)), "strictly")
  testthat::expect_error(ss3_opt_build(d, pop_lengths = list(seq(1, 99, by = 2))),
                         "must include every data length-bin lower edge")
})

testthat::test_that("maturity-at-length round-trips through write_data() / read_data()", {
  d2 <- d
  d2$L50_mat_len   <- c(40, rep(NA_real_, d$nspp - 1))
  d2$slope_mat_len <- c(0.2, rep(NA_real_, d$nspp - 1))
  f <- tempfile(fileext = ".xlsx")
  on.exit(unlink(f))
  suppressMessages(write_data(d2, f))
  back <- suppressMessages(read_data(f))
  testthat::expect_identical(as.numeric(back$L50_mat_len), d2$L50_mat_len)
  testthat::expect_identical(as.numeric(back$slope_mat_len), d2$slope_mat_len)
})
