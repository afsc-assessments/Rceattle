# report_tables() collects the quantities a SAFE chapter reports. The risk it
# carries is not a crash: it is a table that looks right and is keyed wrong --
# a likelihood row attributed to the wrong fleet, two species' biomass sharing a
# year, or a per-recruit reference point read as an estimate when the model
# never defined one. These tests pin the keys.

make_fit <- function(msmMode = 0) {
  set.seed(123)
  sim <- make_msm_test_data(years = 1:30)
  suppressMessages(fit_mod(
    data_list = sim$data_list, estimateMode = 3, msmMode = msmMode,
    growthFun = build_growth(fun = "vonBertalanffy"),
    fit_control = fit_control(getsd = FALSE, verbose = 0)))
}


testthat::test_that("report_tables() returns the documented sections", {
  testthat::skip_on_cran()
  tabs <- report_tables(make_fit())

  testthat::expect_s3_class(tabs, "rceattle_report")
  testthat::expect_true(all(c("model", "likelihood", "timeseries",
                              "reference_points", "fits") %in% names(tabs)))
  # Diagnostics were not supplied, so their sections are absent rather than
  # present and empty.
  testthat::expect_false(any(c("retrospective", "jitter", "osa") %in% names(tabs)))
  # Every section is keyed by model.
  for (s in names(tabs)) testthat::expect_true("model" %in% names(tabs[[s]]))
  testthat::expect_equal(nrow(tabs$model), 1L)
})


testthat::test_that("species is keyed by NAME in every section that has one", {
  # residuals() keys species by integer code; the rest of the package keys it by
  # name. A mismatch here makes the sections unjoinable and silently empties a
  # species filter -- which is how the `fits` section went missing from
  # standard_output().
  testthat::skip_on_cran()
  fit <- make_fit()
  tabs <- report_tables(fit)
  sp <- fit$data_list$spnames

  for (s in c("timeseries", "reference_points", "fits")) {
    testthat::expect_true(all(tabs[[s]]$species %in% sp),
                          info = paste(s, "has a species not in spnames"))
  }
})


testthat::test_that("several models stack with distinct model keys", {
  testthat::skip_on_cran()
  fit <- make_fit()
  tabs <- report_tables(list(base = fit, alt = fit))
  testthat::expect_equal(nrow(tabs$model), 2L)
  testthat::expect_equal(sort(unique(tabs$timeseries$model)), c("alt", "base"))
  # Unnamed models are numbered, not left blank.
  t2 <- report_tables(list(fit, fit))
  testthat::expect_equal(sort(unique(t2$model$model)), c("Model_1", "Model_2"))
  # And an explicit name wins.
  t3 <- report_tables(list(fit), model_names = "m24.0")
  testthat::expect_equal(t3$model$model, "m24.0")
})


testthat::test_that("the likelihood table keys each row on its own axis", {
  # jnll_comp columns count FLEETS on rows 1-8 and SPECIES on rows 9-20. Reading
  # every column as a fleet would report a species' recruitment penalty against
  # a survey.
  testthat::skip_on_cran()
  fit <- make_fit()
  tabs <- report_tables(fit)
  lk <- tabs$likelihood
  d <- fit$data_list

  fleet_rows <- lk[!is.na(lk$axis) & lk$axis == "fleet", ]
  sp_rows    <- lk[!is.na(lk$axis) & lk$axis == "species", ]
  testthat::expect_true(all(fleet_rows$unit %in% d$fleet_control$Fleet_name))
  testthat::expect_true(all(sp_rows$unit %in% d$spnames))

  # The total is carried, and equals the sum of the matrix it came from.
  tot <- lk$weighted[lk$component == "Total"]
  testthat::expect_equal(tot, sum(fit$quantities$jnll_comp, na.rm = TRUE))
})


testthat::test_that("era splits hindcast from projection at endyr", {
  testthat::skip_on_cran()
  fit <- make_fit()
  tabs <- report_tables(fit)
  ts <- tabs$timeseries
  endyr <- fit$data_list$endyr

  testthat::expect_setequal(unique(ts$era), c("time", "fore"))
  testthat::expect_true(all(ts$year[ts$era == "time"] <= endyr))
  testthat::expect_true(all(ts$year[ts$era == "fore"] > endyr))
})


testthat::test_that("per-recruit reference points are NA, not zero, under predation", {
  # The SPR block runs only under msmMode = 0; the template leaves these at zero
  # otherwise, and a zero there would be read as an estimate of zero rather than
  # as "not defined for this model".
  testthat::skip_on_cran()
  ms <- report_tables(make_fit(msmMode = 1))$reference_points
  spr <- ms$estimate[ms$quantity %in% c("SPR0", "SPRlimit", "SPRtarget")]
  testthat::expect_true(all(is.na(spr)))

  # Under single species SPR0 is defined -- it is spawning output per recruit at
  # F = 0, which needs no harvest control rule. SPRtarget / SPRlimit still are
  # not, because this fixture's HCR does not estimate the F they are evaluated
  # at; they are blanked for that reason and not this one.
  ss <- report_tables(make_fit(msmMode = 0))$reference_points
  spr0 <- ss[ss$quantity == "SPR0", ]
  testthat::expect_true(all(is.finite(spr0$estimate)))
  testthat::expect_true(all(spr0$basis == "estimated"))

  ms_spr0 <- ms[ms$quantity == "SPR0", ]
  testthat::expect_true(all(grepl("multispecies", ms_spr0$basis)))
})


testthat::test_that("biomass proxies use the model's own Ptarget / Plimit", {
  # B40% is Ptarget * SB0, not a hardcoded 0.40 -- a model configured to another
  # fraction must report the fraction it actually used.
  testthat::skip_on_cran()
  fit <- make_fit()
  rp <- report_tables(fit)$reference_points
  d <- fit$data_list
  sp1 <- d$spnames[1]

  sb0 <- rp$estimate[rp$species == sp1 & rp$quantity == "SB0"]
  btg <- rp$estimate[rp$species == sp1 & rp$quantity == "B_target"]
  testthat::expect_equal(btg, as.numeric(d$Ptarget)[1] * sb0)
})


testthat::test_that("report_tables() rejects things that are not fits", {
  testthat::expect_error(report_tables(list()), "non-empty list")
  testthat::expect_error(report_tables(list(1, 2)), "must be an Rceattle fit")
})


testthat::test_that("report_tables() never runs a diagnostic itself", {
  # A retrospective or a jitter is tens to hundreds of refits. Passing a
  # non-diagnostic object must error rather than quietly refitting.
  testthat::skip_on_cran()
  fit <- make_fit()
  testthat::expect_error(report_tables(fit, retro = list(1, 2)),
                         "one Rceattle_retro object per model")
})


testthat::test_that("standard_output() emits the NOAA schema in its own order", {
  testthat::skip_on_cran()
  fit <- make_fit()
  std <- standard_output(fit, species = 1)

  testthat::expect_true(all(Rceattle:::.RCE_STANDARD_COLS %in% names(std)))
  # The standard's columns come first, in the standard's order.
  testthat::expect_equal(names(std)[seq_along(Rceattle:::.RCE_STANDARD_COLS)],
                         Rceattle:::.RCE_STANDARD_COLS)
  testthat::expect_true(all(c("species", "model") %in% names(std)))
  testthat::expect_setequal(unique(std$module_name),
                            c("timeseries", "reference_points", "likelihood",
                              "fits"))
})


testthat::test_that("the emitted schema is the one stockplotr reads", {
  # Checking .RCE_STANDARD_COLS against itself proves nothing: a misspelled
  # column is missing from the consumer's filters, which is an error there
  # rather than an empty column. `len_bins` for `length_bins` was exactly that.
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("stockplotr")
  utils::data("example_data", package = "stockplotr", envir = environment())
  std <- standard_output(make_fit(), species = 1)
  testthat::expect_equal(setdiff(names(example_data), names(std)), character(0))
})


testthat::test_that("standard_output() translates names through the dictionary", {
  testthat::skip_on_cran()
  std <- standard_output(make_fit(), species = 1)
  # The quantities the standard has a word for use it...
  testthat::expect_true("spawning_biomass" %in% std$label)
  testthat::expect_true("recruitment" %in% std$label)
  testthat::expect_true("fishing_mortality" %in% std$label)
  testthat::expect_false("ssb" %in% std$label)
  # ...and one it has no word for keeps the Rceattle name rather than vanishing.
  testthat::expect_true("Ftarget" %in% std$label)
})


testthat::test_that("standard_output() refuses to merge species", {
  # The standard describes ONE stock. Returning a frame in which two species'
  # biomass share a year would be silently wrong, so this errors instead.
  testthat::skip_on_cran()
  fit <- make_fit()
  testthat::expect_error(standard_output(fit),
                         "has no species dimension")
  testthat::expect_error(standard_output(fit, species = "Nothing"),
                         "`species` must be one of")

  std <- standard_output(fit, species = fit$data_list$spnames[2])
  testthat::expect_equal(unique(std$species), fit$data_list$spnames[2])
})


testthat::test_that("no standard_label is used by two different quantities", {
  # Three per-recruit quantities once all mapped to "spr", which collapsed them
  # to one key in the standardized output.
  d <- quantity_dictionary()
  lab <- d$standard_label[!is.na(d$standard_label)]
  testthat::expect_false(any(duplicated(lab)))
})


testthat::test_that("a reference point the fit never estimated is NA with a basis", {
  # The single worst failure mode this table has: Ftarget / Flimit sit at their
  # initial exp(0) = 1 when the HCR does not estimate them, and 1.0/yr reported
  # as a target F in an executive summary is catastrophic. The gating comes from
  # build_hcr_map(), not from a second reading of the HCR switch.
  testthat::skip_on_cran()
  fit <- make_fit()          # HCR defaults to NoFishing, so neither is estimated
  rp <- report_tables(fit)$reference_points

  testthat::expect_true("basis" %in% names(rp))
  f <- rp[rp$quantity %in% c("Ftarget", "Flimit"), ]
  testthat::expect_true(all(is.na(f$estimate)))
  testthat::expect_true(all(grepl("not estimated", f$basis)))
  # And the raw array really does hold 1, so this is not a vacuous test.
  testthat::expect_equal(unname(fit$quantities$Ftarget[1]), 1)

  # A quantity that IS defined keeps its value and says so.
  spr0 <- rp[rp$quantity == "SPR0", ]
  testthat::expect_true(all(spr0$basis == "estimated"))
})


testthat::test_that("the gating reads an HCR held as an integer code", {
  # build_hcr_map() matches the harvest control rule by name. A data_list
  # carrying the integer code -- an older saved fit, or one that went through
  # convert_switches() -- matched no rule and blanked an F40% that WAS
  # estimated, with a basis saying it was not.
  testthat::skip_on_cran()
  fit <- make_fit()
  fit$data_list$HCR <- "NPFMC"
  named <- Rceattle:::.rce_refpoint_availability(fit)

  fit$data_list$HCR <- unname(hcr_map["NPFMC"])
  coded <- Rceattle:::.rce_refpoint_availability(fit)

  testthat::expect_true(any(named$ftarget))    # not a vacuously equal pair
  testthat::expect_equal(coded$ftarget, named$ftarget)
  testthat::expect_equal(coded$flimit, named$flimit)
})


testthat::test_that("convergence is reported from the hindcast fit", {
  # The projection re-optimizes and overwrites $opt and $sdrep, so a gradient
  # read off those describes the projection, not whether the assessment
  # converged. fit_mod() snapshots the hindcast into .conv_hindcast.
  testthat::skip_on_cran()
  fit <- make_fit()
  fit$.conv_hindcast <- list(max_gradient = 1.5e-6, pdHess = TRUE)
  fit$opt$max_gradient <- 42        # what a projection refit would leave behind
  fit$sdrep <- list(pdHess = FALSE)

  m <- report_tables(fit)$model
  testthat::expect_equal(m$max_gradient, 1.5e-6)
  testthat::expect_true(m$pdHess)
})


testthat::test_that("an underived unfished reference blanks SB0 and its proxies", {
  # Under msmMode > 0 the template overwrites SB0 / B0 with the MSSB0 / MSB0
  # inputs, which stand at the 999 mt placeholder until fit_mod() derives them.
  # B_target = Ptarget * SB0 would then be 399.6 mt against a real scale of
  # 1e5-1e6 mt. MSSB0_derived is the flag that distinguishes them.
  testthat::skip_on_cran()
  fit <- make_fit(msmMode = 1)
  testthat::expect_false(any(as.logical(fit$data_list$MSSB0_derived)))

  rp <- report_tables(fit)$reference_points
  blanked <- rp[rp$quantity %in% c("SB0", "B0", "B_target", "B_limit"), ]
  testthat::expect_true(all(is.na(blanked$estimate)))
  testthat::expect_true(all(grepl("MSSB0 placeholder", blanked$basis)))
})


testthat::test_that("depletion survives multispecies when the template defines it", {
  # Section 12.1 of the template does NOT divide by SB0 under
  # HCR = NoFishing & msmMode > 0: it divides by biomass in the last projection
  # year, the equilibrated unfished reference. Blanking depletion there because
  # SB0 is a placeholder would destroy a valid series on exactly the fits this
  # package exists for.
  testthat::skip_on_cran()
  fit <- make_fit(msmMode = 1)
  ts <- report_tables(fit)$timeseries
  depl <- ts[ts$quantity %in% c("ssb_depletion", "biomass_depletion"), ]
  testthat::expect_true(any(!is.na(depl$estimate)))

  rp <- report_tables(fit)$reference_points
  td <- rp[rp$quantity == "terminal_depletion", ]
  testthat::expect_true(all(!is.na(td$estimate)))
})


testthat::test_that("standard_output() keeps the FLEET likelihood rows", {
  # jnll_comp is keyed on two axes. Filtering `unit` as though it were always a
  # species dropped every fleet row -- 31 of 38 on a 16-fleet assessment -- and
  # the fleet decomposition is most of what a model comparison reads.
  testthat::skip_on_cran()
  fit <- make_fit()
  tabs <- report_tables(fit)
  n_fleet <- sum(!is.na(tabs$likelihood$axis) & tabs$likelihood$axis == "fleet")
  testthat::expect_gt(n_fleet, 0)

  std <- standard_output(tabs, species = fit$data_list$spnames[1])
  lk <- std[std$module_name == "likelihood", ]
  testthat::expect_equal(sum(!is.na(lk$fleet)), n_fleet)
  # Fleet rows carry the fleet in the standard's own `fleet` column.
  testthat::expect_true(all(lk$fleet[!is.na(lk$fleet)] %in%
                              fit$data_list$fleet_control$Fleet_name))
})


testthat::test_that("a species name containing ', ' still round-trips", {
  # standard_output() used to recover the species list by splitting the model
  # row's comma-joined summary, so "Pollock, GOA" split into two and became
  # unselectable.
  testthat::skip_on_cran()
  fit <- make_fit()
  fit$data_list$spnames <- c("Pollock, GOA", "Cod")
  tabs <- report_tables(fit)
  std <- standard_output(tabs, species = "Pollock, GOA")
  testthat::expect_equal(unique(std$species), "Pollock, GOA")
})


testthat::test_that("a diagnostics list is matched by name, never silently by order", {
  # The realistic mistake is passing ONE model's osa_residuals() stored as a
  # list of parts (index/catch/comp), which would otherwise be attributed one
  # part per model.
  testthat::skip_on_cran()
  fit <- make_fit()
  models <- list(a = fit, b = fit)

  # A stub rather than osa_residuals(): the fixture is a mode-3 build, which
  # oneStepPredict cannot differentiate, and this exercises the pairing rule
  # rather than the residuals themselves.
  osa1 <- structure(
    data.frame(source = "index", fleet = 1L, residual = stats::rnorm(20),
               stringsAsFactors = FALSE),
    class = c("rceattle_osa", "data.frame"))
  wrong <- list(index = osa1, catch = osa1)
  testthat::expect_error(report_tables(models, osa = wrong),
                         "named for model\\(s\\) not in `object`")

  short <- list(a = osa1)
  testthat::expect_error(report_tables(models, osa = short),
                         "one rceattle_osa object per model")

  # Named for the models: matched on the name, so order does not matter.
  ok <- list(b = osa1, a = osa1)
  testthat::expect_silent(res <- report_tables(models, osa = ok))
  testthat::expect_setequal(unique(res$osa$model), c("a", "b"))

  # Unnamed: paired positionally, but it says so.
  testthat::expect_message(report_tables(models, osa = list(osa1, osa1)),
                           "paired with the models in order")
})


testthat::test_that("the two negative log-likelihoods are reported and not conflated", {
  # `marginal_nll` is what nlminb minimized (random effects integrated out);
  # `joint_nll` is what jnll_comp decomposes. Reporting only one made the model
  # and likelihood tables look inconsistent on a random-effects fit, where they
  # differ by the Laplace correction.
  testthat::skip_on_cran()
  fit <- make_fit()
  tabs <- report_tables(fit)
  m <- tabs$model

  testthat::expect_true(all(c("n_random", "marginal_nll",
                              "joint_nll") %in% names(m)))
  # Named for what they are: an NLL, not a likelihood. Nothing should suggest
  # a quantity where larger is better.
  testthat::expect_false(any(grepl("^(marginal|joint)_likelihood$", names(m))))
  # The likelihood table sums to the JOINT nll, by construction.
  total <- tabs$likelihood$weighted[tabs$likelihood$component == "Total"]
  testthat::expect_equal(m$joint_nll, total)
  # No random effects here, so the two agree.
  testthat::expect_equal(m$n_random, 0L)
})


testthat::test_that("a build with no optimization has no parameter table", {
  # estimateMode = 3 never optimizes, so there are no estimated parameters to
  # report and the section is absent rather than empty.
  testthat::skip_on_cran()
  testthat::expect_false("parameters" %in% names(report_tables(make_fit())))
})


testthat::test_that("the parameter table names sigma_R and its process", {
  # Kalei's question: where are sigma_R and M? The parameter table answers it by
  # joining coef()/vcov() to the parameter dictionary.
  testthat::skip_on_cran()
  # random_rec = TRUE, because sigma_R is only ESTIMATED when recruitment is a
  # random effect. Fixed, it never enters coef() and is read off quantities$R_sd.
  utils::data("GOAatf", package = "Rceattle", envir = environment())
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = GOAatf, estimateMode = 1, msmMode = 0, phase = FALSE,
    random_rec = TRUE, fit_control = fit_control(getsd = TRUE, verbose = 0))))

  p <- report_tables(fit)$parameters
  testthat::expect_false(is.null(p))
  testthat::expect_named(p, c("model", "parameter", "index", "natural",
                              "process", "estimate", "std_error"))
  testthat::expect_equal(nrow(p), length(stats::coef(fit)))

  sr <- p[p$parameter == "R_log_sd", ]
  testthat::expect_equal(nrow(sr), 1L)
  testthat::expect_equal(sr$natural, "sigma_R")
  testthat::expect_equal(sr$process, "recruitment")

  # Only a repeated block is numbered, so a row identifies which element it is.
  testthat::expect_true(all(table(p$parameter[!is.na(p$index)]) > 1))
  testthat::expect_true(all(table(p$parameter[is.na(p$index)]) == 1))

  # With recruitment fixed, sigma_R is not estimated and so is not in the table.
  fixed <- suppressMessages(suppressWarnings(fit_mod(
    data_list = GOAatf, estimateMode = 1, msmMode = 0, phase = FALSE,
    random_rec = FALSE, fit_control = fit_control(getsd = TRUE, verbose = 0))))
  testthat::expect_equal(
    nrow(report_tables(fixed)$parameters[
      report_tables(fixed)$parameters$parameter == "R_log_sd", ]), 0L)
  testthat::expect_true(all(is.finite(fixed$quantities$R_sd)))
})


# An indefinite Hessian inverts without being positive definite, so a variance
# can come back below zero. Both reporters took sqrt(abs(.)) of it, which prints
# a meaningless number as an uncertainty in the table that becomes the executive
# summary. Found reviewing the 5.34.0-5.43.0 release (#158).
testthat::test_that("a negative variance reports no standard error", {
  nm  <- c("ln_mean_rec", "sel_inf")
  cov <- diag(c(4, -1)); dimnames(cov) <- list(nm, nm)
  fit <- structure(list(opt = list(par = stats::setNames(c(1, 2), nm)),
                        sdrep = list(cov.fixed = cov, pdHess = FALSE)),
                   class = "Rceattle")

  # The positive variance still reports sqrt(4) = 2, so this is not "NA for
  # everything once one entry is bad".
  testthat::expect_no_warning(p <- Rceattle:::.rce_tab_parameters(fit, "m"))
  testthat::expect_equal(p$std_error, c(2, NA_real_))

  testthat::expect_no_warning(s <- summary(fit))
  testthat::expect_equal(s$coefficients$std_error, c(2, NA_real_))
})


# print(summary(fit)) explains WHY a standard error column is empty, and there are
# two reasons. Neither arm had a test. Added reviewing #158.
testthat::test_that("print(summary()) says which reason the standard errors are missing", {
  nm <- c("ln_mean_rec", "sel_inf")
  mk <- function(sdrep) structure(
    list(opt = list(par = stats::setNames(c(1, 2), nm)), sdrep = sdrep),
    class = "Rceattle")

  # No sdreport: getsd = FALSE.
  none <- utils::capture.output(print(summary(mk(NULL))))
  testthat::expect_true(any(grepl("no sdreport", none)))

  # An sdreport whose every variance is negative: the Hessian is not positive
  # definite, which is a different thing to tell the user.
  cov <- diag(c(-4, -1)); dimnames(cov) <- list(nm, nm)
  bad <- utils::capture.output(print(summary(mk(list(cov.fixed = cov, pdHess = FALSE)))))
  testthat::expect_true(any(grepl("pdHess", bad)))
  testthat::expect_false(any(grepl("no sdreport", bad)))
})
