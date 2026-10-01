# A recruitment `R_init` linkage is an unpenalised log-scale multiplier on the
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
#
# `R_init` is the only linkage parameter with no base parameter in `rec_pars`, so
# its `(Intercept)` is the one that stays estimable. That makes it the one place
# where the shared intercept machinery -- the map, the starting value, the bound
# and the C++ prior re-target -- has to take the other branch, and all four are
# covered below.

# `fit_mod()` overwrites data_list$srr_linkages from recFun, so the spec has to
# travel through build_srr().
#
# The mode argument is `mode`, not `initMode`: R partial-matches a supplied
# `init = ` onto a formal named `initMode` when there is no exact match, which
# silently sends the starting value to the wrong argument.
r_init_fit <- function(mode = "NonEquilibrium", link = "log",
                       formula = ~ 0 + lvl, init = NULL, bounds = NULL,
                       priors = NULL) {
  d <- make_test_data()
  d$env_data <- data.frame(Year = d$styr:d$projyr, lvl = 1)
  fit_mod(d,
          recFun = build_srr(linkages = list(
            R_init = linkage_spec(formula = formula, link = link,
                                  init = init, bounds = bounds,
                                  priors = priors))),
          estimateMode = "DebugBuild", initMode = mode,
          msmMode = 0, random_rec = FALSE,
          fit_control = fit_control(verbose = 0))
}

test_that("an R_init linkage multiplies the initial age-structure, unpenalised", {
  skip_on_cran()
  beta_init <- -1.3879     # log scale: a quarter of R0, SS3's GOA cod value
  fit <- r_init_fit()

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

test_that("an intercept-only R_init linkage carries its own level", {
  skip_on_cran()
  # `~ 1` is the natural spelling for a constant level. It must estimate a
  # coefficient of its own: `R_init` has no base parameter in rec_pars to
  # re-target, so a folded intercept would leave nothing to fit.
  fit <- r_init_fit(formula = ~ 1)
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

test_that("an R_init linkage is refused where the initMode cannot carry it", {
  skip_on_cran()
  # FreeParams never reads R_init, so the level would be estimated and inert.
  expect_error(r_init_fit(mode = "FreeParams"), "FreeParams")
  # OffsetEquilibrium already scales the same ages by rec_dev[, 1].
  expect_error(r_init_fit(mode = "OffsetEquilibrium"), "OffsetEquilibrium")
})

test_that("an R_init linkage is refused on a link that is not log", {
  skip_on_cran()
  expect_error(r_init_fit(link = "identity"), "link = \"log\"")
})

test_that("the R_init param code is in lockstep between R and linkage.hpp", {
  # Rule 12: the codes must match, so read both rather than trusting one.
  hpp_path <- testthat::test_path("..", "..", "src", "TMB", "linkage.hpp")
  skip_if_not(file.exists(hpp_path), "linkage.hpp not in this build")
  hpp <- readLines(hpp_path)
  code <- sub(".*RCEATTLE_REC_R_INIT +([0-9]+).*", "\\1",
              grep("define RCEATTLE_REC_R_INIT", hpp, value = TRUE)[1])
  expect_equal(as.integer(code),
               unname(Rceattle:::LINKAGE_PARAM_CODES$recruitment[["R_init"]]))
  # The offset tensor must have a slot for it.
  n <- sub(".*RCEATTLE_N_REC_PARAMS +([0-9]+).*", "\\1",
           grep("define RCEATTLE_N_REC_PARAMS", hpp, value = TRUE)[1])
  expect_gt(as.integer(n),
            unname(Rceattle:::LINKAGE_PARAM_CODES$recruitment[["R_init"]]))
  # R_init has no rec_pars column, so it must stay out of that registry: the
  # builders key the whole natural-scale/level contract off the NA.
  expect_true(is.na(Rceattle:::.REC_PARAM_TO_INDEX["R_init"]))
})

test_that("an R_init linkage is refused when more than one column is given", {
  skip_on_cran()
  d <- make_test_data()
  d$env_data <- data.frame(Year = d$styr:d$projyr, lvl = 1,
                           temp = seq(0, 1, length.out = d$projyr - d$styr + 1))
  # Only year 0 is read, so an intercept plus a slope share one number.
  expect_error(
    fit_mod(d, recFun = build_srr(linkages = list(
              R_init = linkage_spec(formula = ~ temp))),
            estimateMode = "DebugBuild", initMode = "NonEquilibrium",
            msmMode = 0, random_rec = FALSE),
    "one design column")
  # A per-year random effect estimates deviates no year but the first reads.
  expect_error(
    fit_mod(d, recFun = build_srr(linkages = list(
              R_init = linkage_spec(formula = ~ 0 + (1 | Year)))),
            estimateMode = "DebugBuild", initMode = "NonEquilibrium",
            msmMode = 0, random_rec = FALSE),
    "one design column")
})


# -- The natural-scale contract on the one estimable intercept ---------------
#
# `init` and `bounds` are documented as natural-scale on an (Intercept) and
# link-scale on a slope, and for every other parameter the natural-scale value
# is logged onto the BASE parameter. R_init has no base parameter, so the same
# value has to be logged onto the coefficient instead. Read raw it would invert
# the user's meaning: init = 0.5, "half of R0", would start the level at
# exp(0.5) = 1.65x R0.

test_that("a natural-scale `init` lands logged on the R_init coefficient", {
  skip_on_cran()
  fit <- r_init_fit(formula = ~ 1, init = list(`(Intercept)` = 0.25))
  p <- fit$obj$env$last.par.best
  i <- grep("beta_linkage", names(p))
  expect_equal(unname(p[i]), log(0.25), tolerance = 1e-12)

  # The table's unset default is 0, which is no multiplier at all, so it has to
  # read as 1: no shift off R0, i.e. a coefficient of 0.
  p0 <- r_init_fit(formula = ~ 1)$obj$env$last.par.best
  expect_equal(unname(p0[grep("beta_linkage", names(p0))]), 0,
               tolerance = 1e-12)

  # A multiplier cannot be negative, and log() would hand back NaN in silence.
  expect_error(r_init_fit(formula = ~ 1, init = list(`(Intercept)` = -0.5)),
               "cannot be negative")
})

test_that("a natural-scale `bounds` lands logged on the R_init coefficient", {
  skip_on_cran()
  fit <- r_init_fit(formula = ~ 1,
                    bounds = list(`(Intercept)` = c(0.05, 2)))
  tbl    <- fit$data_list$linkage_table
  is_lvl <- tbl$param == "R_init" & tbl$design_col == "(Intercept)"
  expect_equal(sum(is_lvl), 1L)
  expect_equal(fit$bounds$lower$beta_linkage[is_lvl], log(0.05),
               tolerance = 1e-12)
  expect_equal(fit$bounds$upper$beta_linkage[is_lvl], log(2),
               tolerance = 1e-12)

  # An absent bound stays unbounded rather than becoming log(-Inf) = NaN.
  f2 <- r_init_fit(formula = ~ 1)
  t2 <- f2$data_list$linkage_table
  l2 <- t2$param == "R_init" & t2$design_col == "(Intercept)"
  expect_false(any(is.nan(f2$bounds$lower$beta_linkage[l2])))
  expect_false(any(is.nan(f2$bounds$upper$beta_linkage[l2])))
})


# -- The C++ prior re-target ------------------------------------------------
#
# The linkage-prior block re-points every recruitment (Intercept) prior at
# rec_pars(sp, param), because for R0/alpha/beta the base parameter carries the
# level. rec_pars is nspp x 3 and R_init is code 3, so that read runs one
# element past a PARAMETER_MATRIX. The prior belongs on beta_linkage(i), which
# IS the level -- the same treatment log_sel_apical gets, where a lognormal
# centred on 1 means "no offset".

test_that("a prior on an R_init intercept is read on its own coefficient", {
  skip_on_cran()
  p1 <- 0; p2 <- 0.5
  fit <- r_init_fit(formula = ~ 1,
                    priors = list(`(Intercept)` = prior_lognormal(p1, p2)))

  prow <- match("Linkage-table priors", names(Rceattle:::.JNLL_ROW_AXIS))
  expect_false(is.na(prow))

  p <- fit$obj$env$last.par.best
  i <- grep("beta_linkage", names(p))
  beta <- -0.8
  p[i] <- beta
  rep <- fit$obj$report(p)

  # lognormal on the multiplier exp(beta): the density is evaluated on
  # log(b_nat) = beta itself, against a bias-corrected mean.
  bias  <- as.numeric(fit$data_list$bias_adjust_proc)
  mu    <- p1 - bias * p2^2 / 2
  want  <- -stats::dnorm(beta, mu, p2, log = TRUE)
  expect_equal(sum(rep$jnll_comp[prow, ]), want, tolerance = 1e-8)

  # And the objective stays finite: an out-of-bounds rec_pars read gives
  # whatever was adjacent in memory, which is not reliably finite or stable.
  expect_true(is.finite(rep$jnll))
  expect_equal(fit$obj$report(p)$jnll, rep$jnll, tolerance = 0)
})


# -- The equilibrium the level has to be consistent with --------------------
#
# N_eq is the equilibrium age-structure, and equil_catch_hat integrates Baranov
# over every age of it, age 0 included. If the level scaled ages 1+ but not
# age 0 the "equilibrium" would sit at two recruitment levels at once, so the
# predicted equilibrium catch would not be proportional to the level. initMode 6
# is the only mode that reads an equilibrium catch, and it is the mode the GOA
# Pacific cod bridge this feature was built for uses.

test_that("the equilibrium catch scales with the level at every age", {
  skip_on_cran()
  skip_if_not_installed("TMB")
  skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS
  yrs <- d$styr:d$projyr
  d$env_data <- data.frame(Year = yrs, lvl = 1)

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, file = NULL, inits = NULL,
    recFun = Rceattle::build_srr(linkages = list(
      R_init = Rceattle::linkage_spec(formula = ~ 1))),
    estimateMode = "DebugBuild", initMode = 6, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))

  # The dataset really does carry the equilibrium observation this reads.
  expect_gt(length(fit$quantities$equil_catch_hat), 0)

  p <- fit$obj$env$last.par.best
  i <- grep("beta_linkage", names(p))
  base <- fit$obj$report(p)
  beta <- -0.3
  p[i] <- beta                      # same level for every species
  lvl  <- fit$obj$report(p)

  # Every age of N_eq carries the level, so the Baranov integral over it is
  # exactly proportional. An unscaled age 0 would break the proportionality by
  # that age's share of the equilibrium catch.
  expect_equal(lvl$equil_catch_hat,
               base$equil_catch_hat * exp(beta),
               tolerance = 1e-10)
})


# -- One column PER SPECIES, not one column in the table --------------------
#
# build_srr() takes a list of specs under one parameter key, so a multispecies
# model can give each species its own covariate. Counting unique design columns
# across the pooled table refuses that, although each species has exactly the
# one column the initial state can identify.

test_that("per-species R_init specs may name different covariates", {
  skip_on_cran()
  skip_if_not_installed("TMB")
  sim <- make_msm_test_data(nspp = 2, years = 1:20,
                            log_phi = matrix(-Inf, 2, 2))
  d   <- sim$data_list
  yrs <- d$styr:d$projyr
  d$env_data <- data.frame(
    Year  = yrs,
    temp  = seq(-1, 1, length.out = length(yrs)),
    depth = seq(2, 3, length.out = length(yrs))
  )

  specs <- list(
    Rceattle::linkage_spec(formula = ~ 0 + temp,  species = 1),
    Rceattle::linkage_spec(formula = ~ 0 + depth, species = 2)
  )
  expect_no_error(suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, file = NULL, inits = NULL,
    recFun = Rceattle::build_srr(linkages = list(R_init = specs)),
    estimateMode = "DebugBuild", initMode = "NonEquilibrium",
    msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0)))))
})
