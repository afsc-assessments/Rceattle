# A species with estDynamics > 0 has its numbers-at-age supplied rather than
# estimated, and build_map_fixed_natage() maps its whole M1 random-effect block
# out -- log_M1_dev, M1_dev_log_sd and M1_rho alike -- after build_map_m1() has
# run. So nothing is estimated for it.
#
# But the template's density is gated on M1_re(sp) alone
# (src/TMB/ceattle.cpp, the JNLL_M_RE blocks; no estDynamics test anywhere in
# them), so it still evaluated: N(0, exp(0)) over a field of constant zeros,
# contributing n * log(2*pi)/2 with n the element count. That is a constant --
# zero gradient, no parameter bias -- but it made the objective incomparable
# with an M1_re = 0 fit of the same model, and broke the invariant that the
# "M random effects" jnll row is exactly 0 when nothing is estimated.
#
# Measured end to end on BS2017SS (nyrs_hind 39) with species 1's
# numbers-at-age fixed, M1_model = 1, M1_re = 2:
#
#   before   jnll row "M random effects"  35.838603, 35.838603, 35.838603
#   after                                 0,         35.838603, 35.838603
#   objective 1333557.488469 -> 1333521.649866, a difference of 35.838603
#
# which is exactly 39 * log(2*pi)/2. Species 2 and 3 keep their 35.838603:
# their deviations ARE estimated, and that is the density at the starting
# values, not a degenerate constant.
#
# fit_mod() now treats M1_re as 0 for such a species and says so, in the same
# style as it already forces M1_use_prior to 0 when M1_model is 0.
#
# No bundled dataset ships a fixed-numbers species, so the fixture is built
# here from a normal fit's N_at_age.

fixed_natage_data <- function() {
  base <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = "DebugBuild", msmMode = 0, niter = 3, random_rec = FALSE,
    fit_control = Rceattle::fit_control(verbose = 0))))
  d <- suppressMessages(Rceattle::switch_check(Rceattle::BS2017SS))
  N <- base$quantities$N_at_age
  yrs <- d$styr:d$projyr
  nag <- d$nages[1]
  rows <- do.call(rbind, lapply(seq_along(yrs), function(k) {
    v <- rep(NA_real_, 21)
    v[seq_len(nag)] <- N[1, 1, seq_len(nag), k]
    data.frame(Species_name = d$spnames[1], Species = 1L, Sex = 0L,
               Year = yrs[k],
               as.list(stats::setNames(v, paste0("Age", 1:21))),
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
  dl <- Rceattle::BS2017SS
  dl$NByageFixed <- rows
  dl$estDynamics <- c(1, 0, 0)
  dl
}

m_re_row <- function(fit) {
  i <- which(rownames(fit$quantities$jnll_comp) == "M random effects")
  fit$quantities$jnll_comp[i, ]
}


testthat::test_that("a fixed-numbers species scores no M random effect", {
  testthat::skip_on_cran()
  dl <- fixed_natage_data()

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = dl, inits = NULL, file = NULL, estimateMode = "DebugBuild",
    msmMode = 0, niter = 3, random_rec = FALSE,
    M1Fun = Rceattle::build_M1(M1_model = 1, M1_re = 2),
    fit_control = Rceattle::fit_control(verbose = 0))))

  row <- m_re_row(fit)
  # Species 1's numbers are input, so nothing of its M is estimated.
  testthat::expect_equal(unname(row[1]), 0)
  # Species 2 and 3 estimate theirs, and 39 * log(2*pi)/2 is the density at
  # the starting values -- deviations 0, sigma exp(0) = 1.
  testthat::expect_equal(unname(row[2]), 39 * log(2 * pi) / 2,
                         tolerance = 1e-6)
  testthat::expect_equal(unname(row[3]), 39 * log(2 * pi) / 2,
                         tolerance = 1e-6)
})


testthat::test_that("it warns, and only for the fixed-numbers species", {
  testthat::skip_on_cran()
  dl <- fixed_natage_data()

  w <- character()
  withCallingHandlers(
    suppressMessages(Rceattle::fit_mod(
      data_list = dl, inits = NULL, file = NULL, estimateMode = "DebugBuild",
      msmMode = 0, niter = 3, random_rec = FALSE,
      M1Fun = Rceattle::build_M1(M1_model = 1, M1_re = 2),
      fit_control = Rceattle::fit_control(verbose = 0))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  hit <- grep("treated as 0 there", w, value = TRUE)
  testthat::expect_length(hit, 1L)
  # Names the species whose numbers are input, and no other.
  testthat::expect_match(hit, "species 1", fixed = TRUE)

  # M1_re = 0 has nothing to normalise, so no warning.
  w2 <- character()
  withCallingHandlers(
    suppressMessages(Rceattle::fit_mod(
      data_list = dl, inits = NULL, file = NULL, estimateMode = "DebugBuild",
      msmMode = 0, niter = 3, random_rec = FALSE,
      M1Fun = Rceattle::build_M1(M1_model = 1, M1_re = 0),
      fit_control = Rceattle::fit_control(verbose = 0))),
    warning = function(x) {
      w2 <<- c(w2, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  testthat::expect_length(grep("treated as 0 there", w2), 0L)
})
