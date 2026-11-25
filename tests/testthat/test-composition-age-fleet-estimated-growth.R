# Under estimated growth, a fleet with selectivity at age (Selectivity_dimension
# "Age", or Selectivity "Fixed" on either dimension) has no selectivity-at-length,
# and its age and length comps were predicted as all zero: only a
# Dirichlet-multinomial weight could move the comp likelihood. Found tracing the
# CAAL test failure of 2026-10-30; with a length-based survey beside it no warning
# named the comps.
#
# The expected comps are computed by hand from reported quantities, so a formula
# error cannot pass by agreeing with itself. GOAatf is the two-sex vehicle because
# a joint-sex comp is built in two halves normalized by a shared sum.

.age_fleet_fit <- function(d, inits = NULL, pop_lengths = NULL) {
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    d, inits = inits, file = NULL, estimateMode = "DebugBuild", msmMode = 0,
    random_rec = FALSE,
    growthFun = Rceattle::build_growth(fun = "vonBertalanffy", pop_lengths = pop_lengths),
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
}

# Numbers at true age caught (fishery) or sampled (survey) by a fleet in one year.
.catch_at_age <- function(q, d, flt, sex, yr) {
  sp <- d$fleet_control$Species[flt]
  a  <- seq_len(d$nages[sp])
  Z  <- q$Z_at_age[sp, sex, a, yr]
  N  <- q$N_at_age[sp, sex, a, yr]
  if (as.character(d$fleet_control$Fleet_type[flt]) %in% c("1", "Fishery")) {
    return(q$F_flt_age[flt, sex, a, yr] / Z * (1 - exp(-Z)) * N)
  }
  N * q$sel_at_age[flt, sex, a, yr] * q$index_q[flt, yr] *
    exp(-d$fleet_control$Month[flt] / 12 * Z)
}

# Expected comp row: one block per sex for a joint comp, otherwise pooled.
.hand_comp <- function(q, d, i) {
  row <- d$comp_data[i, ]
  flt <- row$Fleet_code
  sp  <- row$Species
  yr  <- row$Year - d$styr + 1
  by_sex <- lapply(seq_len(d$nsex[sp]), function(sex) {
    x <- .catch_at_age(q, d, flt, sex, yr)
    if (row$Age0_Length1 == 0) return(x)
    # P(length | age) from the fleet's own slot, the species' two slots coming first.
    key <- q$growth_matrix[d$nspp * 2 + flt, sex, seq_len(d$nages[sp]),
                           seq_len(d$nlengths[sp]), yr]
    as.vector(x %*% key)
  })
  out <- switch(as.character(row$Sex),
                "0" = Reduce(`+`, by_sex),
                "1" = by_sex[[1]],
                "2" = by_sex[[2]],
                "3" = unlist(by_sex))
  unname(out / sum(out))
}

.expect_row_matches <- function(q, d, i) {
  want <- .hand_comp(q, d, i)
  got  <- q$comp_hat[i, seq_along(want)]
  testthat::expect_true(all(is.finite(want)) && sum(want) > 0)
  testthat::expect_equal(as.numeric(got), want, tolerance = 1e-10)
}

testthat::test_that("an age-selectivity fleet's age comps follow sel(a) under estimated growth", {
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1, growth = "vonBertalanffy")
  testthat::expect_true(all(d$fleet_control$Selectivity_dimension == "Age"))
  testthat::expect_true(all(d$comp_data$Age0_Length1 == 0))

  q <- .age_fleet_fit(d)$quantities
  for (i in seq_len(nrow(d$comp_data))) .expect_row_matches(q, d, i)
})

testthat::test_that("an age-selectivity fleet's length comps use the estimated key, not age_trans_matrix", {
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1, growth = "vonBertalanffy")
  d$comp_data$Age0_Length1 <- 1

  q <- .age_fleet_fit(d)$quantities
  for (i in seq_len(nrow(d$comp_data))) .expect_row_matches(q, d, i)

  # The input key belongs to empirical growth, so reversing it must change nothing.
  alt <- d
  bins <- grep("^Length_", names(alt$age_trans_matrix))
  alt$age_trans_matrix[, bins] <- alt$age_trans_matrix[, rev(bins)]
  testthat::expect_equal(.age_fleet_fit(alt)$quantities$comp_hat, q$comp_hat, tolerance = 0)
})

testthat::test_that("an age-selectivity fleet's comps respond to its selectivity", {
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1, growth = "vonBertalanffy")
  base <- .age_fleet_fit(d)
  p <- base$estimated_params
  p$sel_inf[] <- p$sel_inf + 0.5
  p$log_sel_slp[] <- p$log_sel_slp - 0.3
  moved <- .age_fleet_fit(d, inits = p)

  comp_row <- function(fit) fit$quantities$jnll_comp["Composition data", 1:2]
  testthat::expect_true(all(is.finite(comp_row(base)), is.finite(comp_row(moved))))
  testthat::expect_true(all(abs(comp_row(moved) - comp_row(base)) > 1e-4))
})

testthat::test_that("both halves of a joint comp and every sex code follow sel(a)", {
  testthat::skip_if_not_installed("TMB")
  d <- Rceattle::GOAatf
  d$fleet_control$Selectivity_dimension <- "Age"
  # Logistic on every fleet, so sel_inf below can give each sex its own curve.
  d$fleet_control$Selectivity <- 1
  # Identity ageing error, so the hand calculation reads true ages for both comp types.
  obs  <- grep("^Obs_age", names(d$age_error))
  nage <- min(nrow(d$age_error), length(obs))
  d$age_error[, obs[seq_len(nage)]] <- diag(nage)

  # The first row of each fleet (survey age, survey length, fishery length), a fishery
  # age row, then those rows under each other sex code. Only a joint row has a second half.
  nl   <- d$nlengths[1]
  live <- which(d$comp_data$Year > 0 & d$comp_data$Sample_size > 0)
  rows <- d$comp_data[sapply(1:3, function(f) live[d$comp_data$Fleet_code[live] == f][1]), ]
  fish_age <- rows[3, ]
  fish_age$Age0_Length1 <- 0
  rows <- rbind(rows, fish_age)
  testthat::expect_equal(as.numeric(rows$Age0_Length1), c(0, 1, 1, 0))
  testthat::expect_true(all(rows$Sex == 3))
  comp_cols <- grep("^Comp_", names(rows))
  cut_to <- function(r, half) {
    r[comp_cols[seq_along(comp_cols) > half]] <- NA
    r
  }
  rows[4, ] <- cut_to(rows[4, ], 2 * d$nages[1])
  variants <- do.call(rbind, lapply(0:2, function(s) {
    r <- rows
    r$Sex <- s
    for (j in seq_len(nrow(r))) {
      r[j, ] <- cut_to(r[j, ], if (r$Age0_Length1[j] == 0) d$nages[1] else nl)
    }
    r
  }))
  n0 <- nrow(d$comp_data)
  d$comp_data <- rbind(d$comp_data, rows[4, ], variants)
  check <- c(match(rownames(rows)[1:3], rownames(d$comp_data)), (n0 + 1):nrow(d$comp_data))

  # Growth estimation requires a CAAL row for every length bin.
  d$caal_data <- cbind(
    data.frame(Fleet_name = rows$Fleet_name[1], Fleet_code = rows$Fleet_code[1],
               Species = rows$Species[1], Sex = 0L, Year = rows$Year[1],
               Length = seq_len(nl), Sample_size = 1L),
    as.data.frame(matrix(1 / d$nages[1], nl, d$nages[1],
                         dimnames = list(NULL, paste0("CAAL_", seq_len(d$nages[1]))))))

  # Different selectivity by sex, or a wrong-sex read would pass.
  p <- .age_fleet_fit(d)$estimated_params
  p$sel_inf[1, , 2] <- p$sel_inf[1, , 2] + 1
  q <- .age_fleet_fit(d, inits = p)$quantities
  testthat::expect_gt(max(abs(q$sel_at_age[1, 1, , 1] - q$sel_at_age[1, 2, , 1])), 0.01)
  testthat::expect_gt(max(abs(q$sel_at_age[3, 1, , 1] - q$sel_at_age[3, 2, , 1])), 0.01)
  for (i in check) .expect_row_matches(q, d, i)
})

testthat::test_that("an age-selectivity fleet's comps read ageing error only on age rows", {
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1, growth = "vonBertalanffy")
  d$comp_data <- rbind(d$comp_data, d$comp_data)
  d$comp_data$Age0_Length1 <- c(0, 0, 1, 1)
  base <- .age_fleet_fit(d)$quantities$comp_hat

  # 60% on the true age, 40% read one year older, the oldest age absorbing its own.
  obs <- grep("^Obs_age", names(d$age_error))
  a <- diag(5)
  for (i in 1:4) {
    a[i, i] <- 0.6
    a[i, i + 1] <- 0.4
  }
  d$age_error[, obs[1:5]] <- a
  q <- .age_fleet_fit(d)$quantities

  testthat::expect_true(all(is.finite(q$comp_hat)))
  testthat::expect_equal(q$comp_hat[3:4, ], base[3:4, ], tolerance = 0)
  for (i in 1:2) {
    x <- .hand_comp(q, d, i) %*% a
    testthat::expect_equal(as.numeric(q$comp_hat[i, 1:5]), as.numeric(x / sum(x)), tolerance = 1e-10)
  }
})

testthat::test_that("an age comp on an age-selectivity fleet ignores a finer length grid", {
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1, growth = "vonBertalanffy")
  q <- .age_fleet_fit(d, pop_lengths = seq(1, 5, by = 0.25))$quantities
  testthat::expect_gt(dim(q$growth_matrix)[4], d$nlengths)
  for (i in seq_len(nrow(d$comp_data))) .expect_row_matches(q, d, i)
})

testthat::test_that("a Fixed fleet labelled Length is predicted from sel(a)", {
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1, growth = "vonBertalanffy")
  d$fleet_control$Selectivity_dimension <- c("Length", "Length")
  d$fleet_control$Selectivity[2] <- "Fixed"
  yrs <- d$styr:d$endyr
  sel <- c(0.1, 0.4, 0.8, 1, 1)
  d$emp_sel <- do.call(rbind, lapply(yrs, function(y) {
    cbind(data.frame(Fleet_name = "Fishery", Fleet_code = 2L, Species = 1L, Sex = 0L, Year = y),
          stats::setNames(as.list(sel), paste0("Comp_", seq_along(sel))))
  }))

  q <- .age_fleet_fit(d)$quantities
  testthat::expect_equal(as.numeric(q$sel_at_age[2, 1, , 1]), sel)
  .expect_row_matches(q, d, which(d$comp_data$Fleet_code == 2)[1])
})

testthat::test_that("CAAL on an age-selectivity fleet stays unpredicted", {
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1, growth = "vonBertalanffy")
  base <- .age_fleet_fit(d)
  testthat::expect_gt(sum(base$quantities$comp_hat), 0)
  testthat::expect_equal(sum(abs(base$quantities$pred_CAAL)), 0)

  p <- base$estimated_params
  p$sel_inf[] <- p$sel_inf + 0.5
  moved <- .age_fleet_fit(d, inits = p)
  testthat::expect_identical(moved$quantities$caal_hat, base$quantities$caal_hat)
})

testthat::test_that("a length-based fleet's comps still come from pred_CAAL", {
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1, growth = "vonBertalanffy")
  d$fleet_control$Selectivity_dimension <- "Length"
  q <- .age_fleet_fit(d)$quantities

  for (i in seq_len(nrow(d$comp_data))) {
    flt <- d$comp_data$Fleet_code[i]
    yr  <- d$comp_data$Year[i] - d$styr + 1
    x   <- rowSums(q$pred_CAAL[flt, 1, 1:5, 1:5, yr])
    testthat::expect_gt(sum(x), 0)
    testthat::expect_equal(as.numeric(q$comp_hat[i, 1:5]), x / sum(x), tolerance = 1e-10)
  }
})

testthat::test_that("simulated comps on an age-selectivity fleet centre on comp_hat", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1, growth = "vonBertalanffy")
  d$comp_data$Sample_size <- 200
  fit <- .age_fleet_fit(d)

  draws <- .comp_draws(fit, nrep = 60)
  for (i in seq_len(nrow(d$comp_data))) {
    p <- as.numeric(fit$quantities$comp_hat[i, 1:5])
    testthat::expect_gt(sum(p), 0)
    x <- apply(draws[i, 1:5, ], 2, function(v) v / sum(v))
    # Multinomial standard error of a mean over 60 draws of N = 200.
    se <- sqrt(p * (1 - p) / (200 * 60))
    testthat::expect_true(all(abs(rowMeans(x) - p) < 4 * se))
  }
})
