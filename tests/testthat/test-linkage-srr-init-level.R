# A recruitment `init` linkage is an unpenalised log-scale multiplier on the
# initial age-structure.
#
# Provenance: SS3 carries the level of an initial population sitting away from R0
# as a free, unpenalised parameter (an SR_regime block on the year before the
# hindcast) and penalises only the age-specific departures. Rceattle had nowhere
# to put the level, so a bridge had to fold it into init_dev and pay the
# recruitment-deviate penalty for it. Measured on GOA Pacific cod 2024: moving a
# -1.3879 level out of init_dev dropped the objective 49.91 nats, all of it in
# the initial-abundance-deviate row, with index, catch, length-composition and
# CAAL unchanged to 1.9e-04.

# `fit_mod()` overwrites data_list$srr_linkages from recFun, so the spec has to
# travel through build_srr().
init_fit <- function(initMode = "NonEquilibrium", link = "log",
                     formula = ~ 0 + lvl) {
  d <- make_test_data()
  d$env_data <- data.frame(Year = d$styr:d$projyr, lvl = 1)
  fit_mod(d,
          recFun = build_srr(linkages = list(
            init = linkage_spec(formula = formula, link = link))),
          estimateMode = "DebugBuild", initMode = initMode,
          msmMode = 0, random_rec = FALSE, verbose = 0)
}

test_that("an init linkage multiplies the initial age-structure, unpenalised", {
  skip_on_cran()
  beta_init <- -1.3879     # log scale: a quarter of R0, SS3's GOA cod value
  fit <- init_fit()

  # Drive beta through the tape rather than its starting value of 0, where the
  # level is inert and the test would pass on any implementation.
  p <- fit$obj$env$last.par.best
  i <- grep("beta_linkage", names(p))
  expect_length(i, 1L)
  base <- fit$obj$report(p)
  p[i] <- beta_init
  lvl <- fit$obj$report(p)

  sp <- 1L
  ages <- 2:fit$data_list$nages[sp]     # age 0 comes from R, not R_init
  expect_equal(lvl$N_at_age[sp, 1, ages, 1],
               base$N_at_age[sp, 1, ages, 1] * exp(beta_init),
               tolerance = 1e-10)

  # Not an additive offset on the natural scale.
  expect_false(isTRUE(all.equal(lvl$N_at_age[sp, 1, ages, 1],
                                base$N_at_age[sp, 1, ages, 1] + beta_init)))

  # The level is NOT charged the init_dev penalty: that row is a function of
  # init_dev alone, which the linkage does not touch. Row read through the
  # registry, never a bare integer (CLAUDE.md, the JnllRow enum).
  irow <- match("Initial abundance deviates", names(Rceattle:::.JNLL_ROW_AXIS))
  expect_false(is.na(irow))
  expect_equal(lvl$jnll_comp[irow, sp], base$jnll_comp[irow, sp],
               tolerance = 1e-12)

  # R0 is the unfished level and must not move with the initial state.
  expect_equal(lvl$R0, base$R0, tolerance = 1e-12)
})

test_that("an intercept-only init linkage carries its own level", {
  skip_on_cran()
  # `~ 1` is the natural spelling for a constant level. It must estimate a
  # coefficient of its own: `init` has no base parameter in rec_pars to
  # re-target, so a folded intercept would leave nothing to fit.
  fit <- init_fit(formula = ~ 1)
  p <- fit$obj$env$last.par.best
  i <- grep("beta_linkage", names(p))
  expect_length(i, 1L)

  base <- fit$obj$report(p)
  p[i] <- -0.5
  lvl <- fit$obj$report(p)
  ages <- 2:fit$data_list$nages[1]
  expect_equal(lvl$N_at_age[1, 1, ages, 1],
               base$N_at_age[1, 1, ages, 1] * exp(-0.5),
               tolerance = 1e-10)
})

test_that("an init linkage is refused where the initMode cannot carry it", {
  skip_on_cran()
  # FreeParams never reads R_init, so the level would be estimated and inert.
  expect_error(init_fit(initMode = "FreeParams"), "FreeParams")
  # OffsetEquilibrium already scales the same ages by rec_dev[, 1].
  expect_error(init_fit(initMode = "OffsetEquilibrium"), "OffsetEquilibrium")
})

test_that("an init linkage is refused on a link that is not log", {
  skip_on_cran()
  expect_error(init_fit(link = "identity"), "link = \"log\"")
})

test_that("the init param code is in lockstep between R and linkage.hpp", {
  # Rule 12: the codes must match, so read both rather than trusting one.
  hpp_path <- testthat::test_path("..", "..", "src", "TMB", "linkage.hpp")
  skip_if_not(file.exists(hpp_path), "linkage.hpp not in this build")
  hpp <- readLines(hpp_path)
  code <- sub(".*RCEATTLE_REC_INIT +([0-9]+).*", "\\1",
              grep("define RCEATTLE_REC_INIT", hpp, value = TRUE)[1])
  expect_equal(as.integer(code),
               unname(Rceattle:::LINKAGE_PARAM_CODES$recruitment[["init"]]))
  # The offset tensor must have a slot for it.
  n <- sub(".*RCEATTLE_N_REC_PARAMS +([0-9]+).*", "\\1",
           grep("define RCEATTLE_N_REC_PARAMS", hpp, value = TRUE)[1])
  expect_gt(as.integer(n),
            unname(Rceattle:::LINKAGE_PARAM_CODES$recruitment[["init"]]))
})

test_that("an init linkage is refused when more than one column is given", {
  skip_on_cran()
  d <- make_test_data()
  d$env_data <- data.frame(Year = d$styr:d$projyr, lvl = 1,
                           temp = seq(0, 1, length.out = d$projyr - d$styr + 1))
  # Only year 0 is read, so an intercept plus a slope share one number.
  expect_error(
    fit_mod(d, recFun = build_srr(linkages = list(
              init = linkage_spec(formula = ~ temp))),
            estimateMode = "DebugBuild", initMode = "NonEquilibrium",
            msmMode = 0, random_rec = FALSE),
    "one design column")
  # A per-year random effect estimates deviates no year but the first reads.
  expect_error(
    fit_mod(d, recFun = build_srr(linkages = list(
              init = linkage_spec(formula = ~ 0 + (1 | Year)))),
            estimateMode = "DebugBuild", initMode = "NonEquilibrium",
            msmMode = 0, random_rec = FALSE),
    "one design column")
})
