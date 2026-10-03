# A length composition must not carry ageing error. The age-length key is
# P(length | TRUE age) and a measured length involves no otolith reading, so the
# ageing-error matrix has no part in predicting one. SS3 applies its `age_age`
# matrix only to age and CAAL data, WHAM only to its age and CAAL predictions,
# and SPoRC passes `AgeingError = NA` for a length composition.
#
# GOAatf is the vehicle because arrowtooth flounder carries the configuration
# that matters: a strongly non-identity ageing-error matrix (max|A - I| = 0.77)
# and length compositions on a JOINT-SEX fleet, which the template builds in two
# halves -- females over the first `nlengths` bins, males over the next
# `nlengths`. Both halves must read the true-age composition, and the row is
# then normalized by a shared sum, so a half corrected on its own contaminates
# the other. Every length-comp row in GOAatf, GOAatf2023 and GOA2018SS species 2
# is joint-sex, so the male half is the only one those fits exercise.

# `comp_ctl` is `comp_data`'s own rows in order (R/5-rearrange_data.R), so row i
# of `comp_hat` is row i of `comp_data`.
ae_scrambled <- function(d) {
  # One row per TRUE age, and only that many Obs_age columns carry the matrix --
  # the sheet is padded to the widest species, so the rest stay NA.
  obs <- grep("^Obs_age", names(d$age_error))
  nage <- nrow(d$age_error)
  # Row-stochastic and markedly different: 60% on the true age, 40% read one
  # year older, the oldest age absorbing its own.
  a <- diag(nage)
  for (i in seq_len(nage - 1)) {
    a[i, i] <- 0.6
    a[i, i + 1] <- 0.4
  }
  d$age_error[, obs[seq_len(nage)]] <- a
  d
}

comp_hat_of <- function(d) {
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, file = NULL, estimateMode = "DebugBuild", msmMode = 0,
    random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
  fit$quantities$comp_hat
}

testthat::test_that("a length composition ignores the ageing-error matrix", {
  testthat::skip_if_not_installed("TMB")

  d <- Rceattle::GOAatf
  is_len <- d$comp_data$Age0_Length1 == 1
  testthat::expect_true(any(is_len))
  # The male half of a joint-sex length comp is the site this pins.
  testthat::expect_true(all(d$comp_data$Sex[is_len] == 3))

  base <- comp_hat_of(d)
  alt  <- comp_hat_of(ae_scrambled(d))

  # Both fits must be finite first: a degenerate alternative would differ from
  # the base everywhere and pass the power check below while proving nothing.
  testthat::expect_true(all(is.finite(base)))
  testthat::expect_true(all(is.finite(alt)))

  testthat::expect_equal(base[is_len, , drop = FALSE],
                         alt[is_len, , drop = FALSE], tolerance = 0)

  # Power check: without it, a matrix that never reached the model would make
  # the invariance above pass for free. Age comps DO read the matrix.
  testthat::expect_false(isTRUE(all.equal(base[!is_len, , drop = FALSE],
                                          alt[!is_len, , drop = FALSE])))
})
