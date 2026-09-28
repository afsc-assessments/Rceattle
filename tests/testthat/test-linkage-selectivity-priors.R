# =============================================================================
# Priors on the base selectivity parameters via the linkage grammar.
#
# A selectivity linkage with an intercept-only formula (`~ 1`) and a `priors`
# entry places a prior on the base selectivity parameter that carries the level
# (log_sel_slp for the slopes, sel_inf for the inflections) -- exactly like the
# prior-only `build_composition()` path. The (Intercept) beta_linkage row is
# pinned at 0 (adds no offset to selectivity), and the C++ prior loop
# re-targets the prior to the base parameter:
#   * slp_asc / slp_desc  -> log_sel_slp (LOG scale)      -> use lognormal()
#   * inf_asc / inf_desc  -> sel_inf     (NATURAL scale)  -> use normal()
# This reproduces the AFSC GOA pollock selectivity priors (Monnahan).
# =============================================================================

.fit3 <- function(sel) {
  e <- new.env(parent = asNamespace("Rceattle"))
  for (f in list.files(testthat::test_path(), "^helper", full.names = TRUE)) sys.source(f, e)
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = e$make_test_data(), selFun = sel, estimateMode = 3, msmMode = 0,
    random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, phase = FALSE, verbose = 0))))
}

testthat::test_that("no selectivity prior => linkage-prior jnll row is 0", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  fit <- .fit3(Rceattle::build_selectivity())
  testthat::expect_equal(sum(fit$quantities$jnll_comp["Linkage-table priors", ]), 0)
})

testthat::test_that("normal prior on inf_asc re-targets the NATURAL-scale sel_inf", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  MU <- 4.0; SD <- 0.5
  # A `~1` selectivity prior defaults to by = ~ fleet (one coefficient per fleet);
  # target fleet 1 so the prior is a single term matching sel_inf[1, 1, 1].
  fit <- .fit3(Rceattle::build_selectivity(linkages = list(
    inf_asc = Rceattle::linkage_spec(~1, fleet = 1L,
                                     priors = list(`(Intercept)` = normal(MU, SD))))))
  inf_val <- fit$estimated_params$sel_inf[1, 1, 1]
  # Prior is read on the natural inflection (NOT exp(0) = 1).
  testthat::expect_equal(sum(fit$quantities$jnll_comp["Linkage-table priors", ]),
                         -dnorm(inf_val, MU, SD, log = TRUE), tolerance = 1e-8)
  # The base parameter is still estimated (a real prior on the mean).
  mp <- Rceattle:::build_map(fit$data_list, fit$estimated_params)$mapList
  testthat::expect_true(any(!is.na(mp$sel_inf)))
})

testthat::test_that("lognormal prior on slp_asc re-targets the LOG-scale log_sel_slp", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  MU <- -1.0; SD <- 1.5
  # A `~1` selectivity prior defaults to by = ~ fleet; target fleet 1 so the
  # prior is a single term matching log_sel_slp[1, 1, 1].
  fit <- .fit3(Rceattle::build_selectivity(linkages = list(
    slp_asc = Rceattle::linkage_spec(~1, fleet = 1L,
                                     priors = list(`(Intercept)` = lognormal(MU, SD))))))
  lslp <- fit$estimated_params$log_sel_slp[1, 1, 1]
  # Mean-centred under the default bias_adjust_proc = TRUE.
  testthat::expect_equal(sum(fit$quantities$jnll_comp["Linkage-table priors", ]),
                         -dnorm(lslp, MU - SD^2 / 2, SD, log = TRUE), tolerance = 1e-8)
})

# --- Guards on sel-prior re-targeting (density-review follow-ups) ------------

testthat::test_that("prior on DoubleNormal right_floor/inf_desc is rejected (logit scale)", {
  lt <- data.frame(process = "sel", param = "right_floor", fleet = 1L,
                   prior_family = "normal", stringsAsFactors = FALSE)
  fc <- data.frame(Fleet_name = "F1", Selectivity = "DoubleNormal",
                   Selectivity_index = 1L, Fleet_code = 1L, stringsAsFactors = FALSE)
  testthat::expect_error(Rceattle:::.check_sel_linkage_support(lt, fc), "logit")
  # The same inflection prior on a DoubleLogistic fleet (natural scale) is allowed.
  fc$Selectivity <- "DoubleLogistic"; lt$param <- "inf_desc"
  testthat::expect_silent(Rceattle:::.check_sel_linkage_support(lt, fc))
})

testthat::test_that("prior on a mirror fleet (shared Selectivity_index) is rejected", {
  fc <- data.frame(Fleet_name = c("F1", "F2"), Selectivity = c("Logistic", "Logistic"),
                   Selectivity_index = c(1L, 1L), Fleet_code = c(1L, 2L),
                   stringsAsFactors = FALSE)
  # Fleet 2 follows fleet 1 (same Selectivity_index, and fleet 1 leads the
  # group) -> the shared block would be penalized twice.
  lt2 <- data.frame(process = "sel", param = "inf_asc", fleet = 2L,
                    prior_family = "normal", stringsAsFactors = FALSE)
  testthat::expect_error(Rceattle:::.check_sel_linkage_support(lt2, fc),
                         "share a Selectivity_index")
  # The lead fleet, i.e. the group's first fleet, carries the block prior.
  lt1 <- lt2; lt1$fleet <- 1L
  testthat::expect_silent(Rceattle:::.check_sel_linkage_support(lt1, fc))
})

testthat::test_that("a covariate (non-prior) sel linkage on a mirror fleet is allowed", {
  fc <- data.frame(Fleet_name = c("F1", "F2"), Selectivity = c("Logistic", "Logistic"),
                   Selectivity_index = c(1L, 1L), Fleet_code = c(1L, 2L),
                   stringsAsFactors = FALSE)
  # prior_family "none" => an additive covariate offset, which is fine on a mirror.
  lt <- data.frame(process = "sel", param = "inf_asc", fleet = 2L,
                   prior_family = "none", stringsAsFactors = FALSE)
  testthat::expect_silent(Rceattle:::.check_sel_linkage_support(lt, fc))
})


testthat::test_that("a prior on a limb the fleet's curve lacks is rejected", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # A Logistic fleet has no descending limb, so sel_inf[2] / log_sel_slp[2] never
  # reach selectivity-at-age: they sit at their build defaults. A prior on them is
  # still added to the objective, so it shifts the reported likelihood by a
  # constant that moves with an unrelated default while doing nothing to the fit.
  # That is how a reconciliation against another model picks up an unexplained
  # offset, so refuse it rather than accept it silently.
  d  <- Rceattle::BS2017SS
  fc <- Rceattle::switch_check(Rceattle::clean_data(d))$fleet_control
  logistic <- fc$Fleet_name[fc$Selectivity == "Logistic"][1]
  testthat::skip_if(is.na(logistic))

  run <- function(selfun) suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, file = NULL, inits = NULL, estimateMode = 3,
    random_rec = FALSE, msmMode = 0, selFun = selfun,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))

  testthat::expect_error(run(Rceattle::build_selectivity(linkages = list(
    inf_desc = Rceattle::linkage_spec(~ 1, fleet = logistic,
                                      priors = list(intercept = normal(10, 3)))))),
    "does not use those")

  testthat::expect_error(run(Rceattle::build_selectivity(linkages = list(
    slp_desc = Rceattle::linkage_spec(~ 1, fleet = logistic,
                                      priors = list(intercept = lognormal(-1, 1.5)))))),
    "does not use those")

  # The limb the curve does have is unaffected.
  testthat::expect_s3_class(run(Rceattle::build_selectivity(linkages = list(
    inf_asc = Rceattle::linkage_spec(~ 1, fleet = logistic,
                                     priors = list(intercept = normal(0, 3)))))),
    "Rceattle")
})
