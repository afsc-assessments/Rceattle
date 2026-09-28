# One-step-ahead (OSA) residuals -- Phase 1 (aggregate catch + index).
#
# A gold-standard "joint nll is numerically unchanged by the obsvec/keep
# refactor" check (comparing against a DLL built from the pre-change source) is
# kept locally -- it is too slow for routine testing because it recompiles a
# second DLL. Here we test the parts that guarantee
# correctness and are fast/CI-appropriate: the obsvec/obs_ctl construction and
# the end-to-end osa_residuals() machinery.

testthat::test_that("build_osa_data lays observations into obsvec correctly", {
  testthat::skip_if_not_installed("Rceattle")

  dat <- make_test_data(nyrs = 8, nages = 5, seed = 123)
  dl  <- Rceattle::rearrange_data(dat, build_osa = TRUE)

  # obs_ctl must be a DATA FRAME (regression guard: rearrange_data's
  # data.frame->matrix coercion must not sweep it up).
  testthat::expect_s3_class(dl$obs_ctl, "data.frame")
  testthat::expect_true(all(c("obs_pos", "source", "data_row", "fleet_code",
                              "year", "is_last_bin") %in% names(dl$obs_ctl)))

  # Defaults and types.
  testthat::expect_identical(dl$osa_mode, 0L)
  testthat::expect_type(dl$obsvec, "double")
  testthat::expect_type(dl$index_obsvec_idx, "integer")
  testthat::expect_type(dl$catch_obsvec_idx, "integer")

  endyr <- dl$endyr
  ft    <- dl$flt_type

  # Index: included rows are exactly those passing the TMB guard, and obsvec
  # stores log(observation) at the mapped position.
  inc_i <- which(dl$index_ctl[, 3] > 0 & dl$index_ctl[, 3] <= endyr &
                   ft[dl$index_ctl[, 1]] > 0 & dl$index_obs[, 1] > 0)
  testthat::expect_true(all(dl$index_obsvec_idx[inc_i] >= 0))
  testthat::expect_equal(dl$obsvec[dl$index_obsvec_idx[inc_i] + 1L],
                         log(dl$index_obs[inc_i, 1]))

  # Catch: same, with the fishery-only guard (flt_type == 1).
  inc_c <- which(dl$catch_ctl[, 3] > 0 & dl$catch_ctl[, 3] <= endyr &
                   ft[dl$catch_ctl[, 1]] == 1 & dl$catch_obs[, 1] > 0)
  testthat::expect_equal(dl$obsvec[dl$catch_obsvec_idx[inc_c] + 1L],
                         log(dl$catch_obs[inc_c, 1]))

  # One obs_ctl row per included aggregate observation.
  testthat::expect_equal(sum(dl$obs_ctl$source %in% c("index", "catch")),
                         length(inc_i) + length(inc_c))

  # Composition: bin counts = (proportion + 1e-5) * Sample_size, with one
  # decomposition group per row and the final bin flagged is_last_bin.
  if (nrow(dl$comp_obs) > 0) {
    r1 <- which(dl$comp_obsvec_idx >= 0)[1]
    sp <- dl$comp_ctl[r1, 2]
    n_comp <- dl$nages[sp]   # make_test_data: single-sex age comps
    start  <- dl$comp_obsvec_idx[r1]
    expected <- (as.numeric(dl$comp_obs[r1, seq_len(n_comp)]) + 0.00001) *
      dl$comp_n[r1, 2]
    testthat::expect_equal(dl$obsvec[start + seq_len(n_comp)], expected)

    comp_ctl_rows <- dl$obs_ctl[dl$obs_ctl$source == "comp", ]
    testthat::expect_true("length" %in% names(dl$obs_ctl))
    # exactly one last-bin (dropped) per composition group
    testthat::expect_equal(sum(comp_ctl_rows$is_last_bin),
                           length(unique(comp_ctl_rows$group_id)))
  }
})


testthat::test_that("default fit skips composition OSA data (fast path)", {
  testthat::skip_if_not_installed("Rceattle")

  dat  <- make_test_data(nyrs = 8, nages = 5, seed = 123)
  fast <- Rceattle::rearrange_data(dat)                    # build_osa = FALSE (default)
  full <- Rceattle::rearrange_data(dat, build_osa = TRUE)

  # Fast path: aggregate index/catch entries are still built (the template
  # always reads them), but no composition / caal / diet metadata.
  testthat::expect_true(all(fast$comp_obsvec_idx < 0))
  testthat::expect_true(all(fast$caal_obsvec_idx < 0))
  testthat::expect_false(any(fast$obs_ctl$source %in% c("comp", "caal", "diet")))

  # Aggregate obsvec entries are identical between the fast and full builds, so
  # the fitted objective is unaffected by the toggle.
  ii <- fast$index_obsvec_idx[fast$index_obsvec_idx >= 0] + 1L
  testthat::expect_equal(fast$obsvec[ii], full$obsvec[ii])

  # The full build does produce composition entries when comp data are present.
  if (nrow(full$comp_obs) > 0) {
    testthat::expect_true(any(full$comp_obsvec_idx >= 0))
    testthat::expect_true(any(full$obs_ctl$source == "comp"))
  }
})


testthat::test_that("build_osa_data lays CAAL observations into obsvec correctly", {
  testthat::skip_if_not_installed("Rceattle")

  # Parametric growth populates conditional age-at-length (CAAL) data.
  dat <- make_test_data(nyrs = 8, nages = 5, seed = 123,
                        growth = "vonBertalanffy")
  dl  <- Rceattle::rearrange_data(dat, build_osa = TRUE)
  testthat::skip_if(nrow(dl$caal_obs) == 0, "no CAAL data in fixture")

  testthat::expect_true("length" %in% names(dl$obs_ctl))

  r1     <- which(dl$caal_obsvec_idx >= 0)[1]
  sp     <- dl$caal_ctl[r1, 2]
  n_caal <- dl$nages[sp]                 # CAAL bins span ages (no joint-sex)
  start  <- dl$caal_obsvec_idx[r1]
  expected <- (as.numeric(dl$caal_obs[r1, seq_len(n_caal)]) + 0.00001) *
    dl$caal_n[r1, 1]
  testthat::expect_equal(dl$obsvec[start + seq_len(n_caal)], expected)

  caal_rows <- dl$obs_ctl[dl$obs_ctl$source == "caal", ]
  testthat::expect_equal(sum(caal_rows$is_last_bin),
                         length(unique(caal_rows$group_id)))
  testthat::expect_true(all(!is.na(caal_rows$length)))           # conditioning length
  testthat::expect_true(all(!is.na(caal_rows$age_length_bin)))    # age bin
})


testthat::test_that("build_osa_data lays diet observations into obsvec correctly", {
  testthat::skip_if_not_installed("Rceattle")

  # Minimal data_list with two stomachs of predator species 1 (suitMode > 0):
  # stomach 0 has 2 prey items, stomach 1 has 1. Aggregate/comp/caal are empty.
  empty_i <- function(nc) matrix(integer(0), 0, nc)
  empty_n <- function(nc) matrix(numeric(0), 0, nc)
  dl <- list(
    endyr = 5, flt_type = c(1L, 2L), nages = c(3L, 3L), nlengths = c(3L, 3L),
    index_ctl = empty_i(3), index_obs = empty_n(2),
    catch_ctl = empty_i(3), catch_obs = empty_n(2),
    comp_ctl = empty_i(5), comp_obs = empty_n(3), comp_n = empty_n(2),
    caal_ctl = empty_i(5), caal_obs = empty_n(3), caal_n = empty_n(1),
    n_stomach_obs = 2L,
    stomach_id = c(0L, 0L, 1L),
    diet_ctl = matrix(c(1, 1, 0, 0, 1, 1, 3,
                        1, 2, 0, 0, 1, 1, 3,
                        1, 1, 0, 0, 2, 1, 4), nrow = 3, byrow = TRUE),
    diet_obs = matrix(c(50, 0.6, 50, 0.3, 40, 0.7), nrow = 3, byrow = TRUE),
    suitMode = c(2L, 0L))

  dl <- Rceattle:::build_osa_data(dl, build_osa = TRUE)

  testthat::expect_length(dl$diet_obsvec_idx, 2L)
  testthat::expect_true(all(dl$diet_obsvec_idx >= 0))

  diet_rows <- dl$obs_ctl[dl$obs_ctl$source == "diet", ]
  # stomach 0: 2 prey + other = 3 bins; stomach 1: 1 prey + other = 2 bins.
  testthat::expect_equal(nrow(diet_rows), 5L)
  testthat::expect_equal(sum(diet_rows$is_last_bin), 2L)   # one "other prey" per stomach
  testthat::expect_setequal(unique(diet_rows$stomach_id), c(0L, 1L))

  # Counts = (proportions + other + 1e-5) normalized, scaled to the sample size.
  v <- c(0.6, 0.3, 1 - 0.9) + 0.00001
  expected0 <- v / sum(v) * 50
  start0 <- dl$diet_obsvec_idx[1]
  testthat::expect_equal(dl$obsvec[start0 + seq_len(3L)], expected0)
})


testthat::test_that("obsvec/keep refactor leaves the fitted objective finite", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  dat <- make_test_data(nyrs = 8, nages = 5, seed = 123)
  res <- Rceattle::fit_mod(data_list = dat, inits = NULL, estimateMode = 3,
                           fit_control = fit_control(phase = FALSE, verbose = 0))

  testthat::expect_false(is.null(res$obj))
  testthat::expect_true(is.finite(res$obj$fn()))
  jc <- res$obj$report(res$obj$env$last.par.best)$jnll_comp
  testthat::expect_true(is.finite(sum(jc[1, ])))   # index slot
  testthat::expect_true(is.finite(sum(jc[2, ])))   # catch slot

  # obs_ctl is carried onto the fitted object as a data frame.
  testthat::expect_s3_class(res$obs_ctl, "data.frame")
})


testthat::test_that("fitted objective is unchanged (golden jnll on BS2017SS)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Total negative log-likelihood at the initial parameters of a fixed bundled
  # model. This guards that the fitted objective is unchanged -- in particular
  # that the OSA obsvec/keep machinery leaves normal fitting bit-for-bit
  # identical. estimateMode = 3 builds the model and reports jnll_comp without
  # optimizing, so the value is deterministic. Update the golden value ONLY when
  # the likelihood is intentionally changed. Uses the default comp_offset = 1e-5
  # (the historical composition proportion offset); set comp_offset = 0 for a
  # WHAM-style multinomial.
  data("BS2017SS", package = "Rceattle")
  fit <- Rceattle::fit_mod(BS2017SS, estimateMode = 3, msmMode = 0,
                           fit_control = fit_control(phase = FALSE, verbose = 0))
  jc <- fit$obj$report(fit$obj$par)$jnll_comp
  testthat::expect_equal(sum(jc[is.finite(jc)]), 1537036.2876293703,
                         tolerance = 1e-6)
})


testthat::test_that("process_residuals() runs on a converging model", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOApollock", package = "Rceattle")
  fit <- Rceattle::fit_mod(
    data_list = GOApollock, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(verbose = 0, phase = TRUE, getsd = TRUE))

  pr <- process_residuals(fit, process = "recruitment")
  testthat::expect_s3_class(pr, "rceattle_osa")
  testthat::expect_equal(unique(pr$source), "recruitment")
  testthat::expect_true(all(is.finite(pr$residual)))
  # One recruitment residual per hindcast year.
  testthat::expect_equal(nrow(pr), length(GOApollock$styr:GOApollock$endyr))
  testthat::expect_setequal(range(pr$year), c(GOApollock$styr, GOApollock$endyr))

  # Deterministic with a fixed seed.
  pr2 <- process_residuals(fit, process = "recruitment")
  testthat::expect_equal(pr$residual, pr2$residual)

  # Does not clobber the caller's global RNG: the posterior draw is seeded on a
  # local stream and .Random.seed is restored on exit.
  set.seed(101L)
  before <- get(".Random.seed", envir = .GlobalEnv)
  invisible(process_residuals(fit, process = "recruitment", seed = 7))
  testthat::expect_identical(before, get(".Random.seed", envir = .GlobalEnv))

  # process = "all" returns every supported process present, all finite. On
  # GOApollock the catchability deviates use a random-walk prior, so the
  # marginal-SD standardization is approximate and emits a warning.
  testthat::expect_warning(pr_all <- process_residuals(fit, process = "all"),
                           "approximate")
  testthat::expect_true("recruitment" %in% pr_all$source)
  testthat::expect_true(all(is.finite(pr_all$residual)))

  # Recruitment-only residuals use an exact iid prior -> no approximation warning.
  testthat::expect_no_warning(process_residuals(fit, process = "recruitment"))

  # residuals(type = "process") reshapes into the common residual schema.
  r <- residuals(fit, type = "process")
  testthat::expect_true(all(c("Source", "Species", "Year", "Residual") %in% names(r)))

  # Composition OSA data is built on demand, so comp residuals work from any fit
  # -- no fit_control(osa = TRUE) needed.
  osa_comp <- osa_residuals(fit, source = "comp")
  testthat::expect_s3_class(osa_comp, "rceattle_osa")
  testthat::expect_true(all(is.finite(osa_comp$residual)))
})


testthat::test_that("osa_residuals() runs end-to-end on a converging model", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOApollock", package = "Rceattle")
  fit <- Rceattle::fit_mod(
    data_list = GOApollock, inits = NULL, file = NULL,
    estimateMode = 1, random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(verbose = 0, phase = TRUE))

  osa <- osa_residuals(fit, source = c("index", "catch", "comp"))

  testthat::expect_s3_class(osa, "rceattle_osa")
  testthat::expect_true(all(c("source", "fleet", "year", "residual") %in% names(osa)))
  testthat::expect_true(all(is.finite(osa$residual)))
  testthat::expect_setequal(unique(osa$source), c("index", "catch", "comp"))

  # Full observation map, regenerated on demand (as osa_residuals() does
  # internally); the stored fit$obs_ctl carries only the aggregate series.
  ctl <- Rceattle:::build_osa_data(fit$obj$env$data, build_osa = TRUE)$obs_ctl

  # Aggregate: one residual per included catch/index observation.
  agg <- osa[osa$source %in% c("index", "catch"), ]
  testthat::expect_equal(nrow(agg),
                         sum(ctl$source %in% c("index", "catch")))

  # Composition: each composition drops its sum-to-N last bin, so there is one
  # fewer residual than bins per composition (here per fleet x year group).
  n_comp_bins   <- sum(ctl$source == "comp")
  n_comp_groups <- length(unique(ctl$group_id[ctl$source == "comp"]))
  testthat::expect_equal(sum(osa$source == "comp"), n_comp_bins - n_comp_groups)
  testthat::expect_true(all(is.finite(osa$residual[osa$source == "comp"])))

  # Deterministic with a fixed seed (cheap aggregate path).
  a1 <- osa_residuals(fit, source = "catch")
  a2 <- osa_residuals(fit, source = "catch")
  testthat::expect_equal(a1$residual, a2$residual)

  # Diagnostics return one row per source plus an overall row, with the null
  # intervals populated.
  diag <- osa_diagnostics(osa)
  testthat::expect_true(nrow(diag) >= 2)
  testthat::expect_true(all(c("sdnr", "sdnr_lo", "sdnr_hi") %in% names(diag)))
  testthat::expect_true(is.finite(diag$sdnr[diag$group == "all"]))

  # residuals(type = "osa") reshapes into the common residual schema; `source`
  # selects the data source(s) (aggregate subset to keep the test cheap).
  r <- residuals(fit, type = "osa", source = c("index", "catch"))
  testthat::expect_true(all(c("Source", "Fleet_code", "Year", "Residual") %in% names(r)))
  testthat::expect_equal(nrow(r), sum(ctl$source %in% c("index", "catch")))

  # glm-style residual kinds: response / pearson cover all sources by default.
  rr <- residuals(fit, type = "response")
  testthat::expect_true(all(c("index", "catch", "comp") %in% rr$Source))
  testthat::expect_true(all(is.finite(rr$Residual)))
  rp <- residuals(fit, type = "pearson", source = "index")
  testthat::expect_setequal(unique(rp$Source), "index")
  testthat::expect_true(all(is.finite(rp$Residual)))

  # `type` selects the residual kind only; passing a data-source name (use
  # `source` for that) is an error.
  testthat::expect_error(residuals(fit, type = "index"))

  # species filter (GOApollock is a single-species model -> species 1).
  testthat::expect_setequal(unique(residuals(fit, species = 1)$Species), 1)
  testthat::expect_equal(nrow(residuals(fit, species = 999)), 0L)

  # plot() builds a separate aggregate (Q-Q only) and composition (Q-Q + OSA and
  # Pearson bubbles) figure. The osa object carries fleet_name and the Pearson
  # residuals needed for the composition bubbles.
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    testthat::expect_true("fleet_name" %in% names(osa))
    testthat::expect_false(is.null(attr(osa, "pearson")))
    pf <- tempfile(fileext = ".pdf"); grDevices::pdf(pf)
    pl    <- plot(osa)
    p_idx <- plot(osa, source = "index")                # source filter
    p_sep <- plot(osa, source = "comp", combine = FALSE) # split age/length
    p_sp  <- plot(osa, species = 1)                     # species filter
    pc    <- plot_comp(fit)                             # ggplot2 comp figures
    grDevices::dev.off()
    testthat::expect_true(all(c("aggregate", "composition") %in% names(pl)))
    testthat::expect_equal(names(p_idx), "aggregate")
    testthat::expect_true("composition_age" %in% names(p_sep))
    testthat::expect_false("composition" %in% names(p_sep))
    testthat::expect_true(length(p_sp) >= 1L)
    testthat::expect_true("pearson" %in% names(pc))
    testthat::expect_true(any(grepl("^annual", names(pc))))
    testthat::expect_true(any(grepl("^aggregated", names(pc))))
  }

  # OSA must refuse a debug (estimateMode >= 3) fit.
  fit_dbg <- Rceattle::fit_mod(data_list = GOApollock, inits = NULL,
                               estimateMode = 3,
                               fit_control = fit_control(phase = FALSE, verbose = 0))
  testthat::expect_error(osa_residuals(fit_dbg), "estimateMode")
})


# Reference diet Pearson residual, written out on the scale the stomach-content
# likelihood works in (ceattle.cpp:4815-4860): each stomach carries an "other
# prey" balance that enters the density but not the residual frame, observed and
# fitted are each normalized by their own offset-inflated total over prey + other,
# and the counts are formed AFTER that normalization, so the total is N_s.
# Before v5.29.0 the residual divided by the raw multinomial variance at N_s,
# ignoring all of this and the Diet_comp_weights effective sample size.
#
# This checks residuals()' wiring -- grouping, normalization, family lookup --
# not the variance itself, which is checked against simulation in
# test-likelihood-pearson-effective-n.R.
diet_pearson_ref <- function(dd, hat, off = 1e-5, family = 0L, weight = 1) {
  sid    <- paste(dd$Pred, dd$Pred_sex, dd$Pred_age, dd$Year)
  obs    <- dd$Stomach_proportion_by_weight
  nbin   <- as.numeric(stats::ave(rep(1, length(obs)), sid, FUN = sum)) + 1
  sum_o  <- as.numeric(stats::ave(obs, sid, FUN = sum))
  sum_h  <- as.numeric(stats::ave(hat, sid, FUN = sum))
  S_o    <- sum_o + (1 - pmin(sum_o, 1)) + off * nbin
  S_h    <- sum_h + pmax(1 - sum_h, 1e-5) + off * nbin
  p_o    <- (obs + off) / S_o
  p_h    <- (hat + off) / S_h
  n      <- as.numeric(dd$Sample_size)
  fam    <- rep_len(family, length(obs)); w <- rep_len(weight, length(obs))
  ifelse(fam == 1L,
         (p_o - p_h) / sqrt(p_h * (1 - p_h) * (n + n * exp(w)) / (n * (1 + n * exp(w)))),
         (p_o - p_h) / sqrt(p_h * (1 - p_h) / (w * n)))
}


testthat::test_that("residuals(source = 'diet') computes diet Pearson residuals", {
  testthat::skip_if_not_installed("Rceattle")

  # Minimal fitted-model stand-in carrying just the diet pieces residuals() uses.
  dd <- data.frame(Pred = c(1L, 1L, 2L), Pred_sex = 0L, Prey = c(1L, 2L, 1L),
                   Prey_sex = 0L, Pred_age = c(3L, 3L, 4L), Prey_age = c(1L, 2L, 1L),
                   Year = c(2000L, 2000L, 2001L),
                   Stomach_proportion_by_weight = c(0.6, 0.3, 0.5),
                   Sample_size = c(50, 50, 40))
  fit <- structure(list(
    data_list  = list(diet_data = dd, spnames = c("A", "B")),
    quantities = list(diet_hat = cbind(NA_real_, c(0.55, 0.35, 0.45)))),
    class = "Rceattle")

  r <- residuals(fit, type = "pearson", source = "diet")
  testthat::expect_setequal(unique(r$Source), "diet")
  testthat::expect_true(all(c("Pred_sex", "Prey", "Prey_sex", "Pred_age", "Prey_age",
                              "Observed", "Fitted", "Residual") %in% names(r)))
  hat <- fit$quantities$diet_hat[, 2]
  # No Diet_distribution and no fitted weights on this stand-in, so the family
  # falls back to the schema default, the multinomial, at weight 1.
  testthat::expect_equal(r$Residual, diet_pearson_ref(dd, hat))
  testthat::expect_equal(residuals(fit, type = "response", source = "diet")$Residual,
                         dd$Stomach_proportion_by_weight - hat)
  # diet uses a predator/prey schema, so it must be requested on its own.
  testthat::expect_error(residuals(fit, source = c("diet", "comp")), "on its own")
  # species filter acts on the predator species.
  testthat::expect_setequal(unique(residuals(fit, source = "diet", species = 1)$Species),
                            1L)
})


testthat::test_that("a Dirichlet-multinomial predator gets the DM variance", {
  # Regression: the family was read off `bioenergetics_control`, a workbook
  # sheet that is never an element of a data_list, so every predator scored as a
  # multinomial and `Diet_comp_weights` -- a LOG under a DM -- was used as a
  # natural-scale multiplier on the sample size. Both defects are silent: the
  # residuals come back finite and plausible, just standardized by the wrong
  # variance.
  dd <- data.frame(Pred = c(1L, 1L, 2L), Pred_sex = 0L, Prey = c(1L, 2L, 1L),
                   Prey_sex = 0L, Pred_age = c(3L, 3L, 4L), Prey_age = c(1L, 2L, 1L),
                   Year = c(2000L, 2000L, 2001L),
                   Stomach_proportion_by_weight = c(0.6, 0.3, 0.5),
                   Sample_size = c(50, 50, 40))
  hat <- c(0.55, 0.35, 0.45)
  # Predator 1 is Dirichlet-multinomial, predator 2 multinomial, so one call
  # covers the branch and its neighbour.
  mk <- function(dist, wt) structure(list(
    data_list = list(diet_data = dd, spnames = c("A", "B"), comp_offset = 1e-5,
                     Diet_distribution = dist, Diet_comp_weights = wt),
    quantities = list(diet_hat = cbind(NA_real_, hat))), class = "Rceattle")

  # A positive log weight, so the multinomial reading of it is still a valid
  # variance -- the contrast below then comes from the family, not from a
  # negative effective sample size turning the comparison into NaN against NaN.
  # theta = 20 rather than 2: at 2 the DM's effective N (34.0) and the
  # multinomial misreading (34.7) very nearly collide, so the contrast would
  # rest on all.equal()'s default tolerance instead of on the families differing.
  lw  <- log(20)
  fit <- mk(c(1L, 0L), c(lw, 1))
  r   <- residuals(fit, type = "pearson", source = "diet")$Residual
  testthat::expect_equal(
    r, diet_pearson_ref(dd, hat, family = c(1L, 1L, 0L),
                        weight = c(lw, lw, 1)))
  testthat::expect_true(all(is.finite(r)))

  # The DM rows must differ from what the multinomial fallback would have given;
  # equality there is the bug returning, not a coincidence.
  mn <- residuals(mk(c(0L, 0L), c(lw, 1)), type = "pearson",
                  source = "diet")$Residual
  testthat::expect_true(all(is.finite(mn)))
  # A real separation, not one resting on all.equal()'s tolerance. At theta = 20
  # the DM assumes an effective N of 47.7 (= n(1+A)/(n+A), A = n*theta = 1000)
  # where the misreading assumes log(20) * 50 = 149.8. A residual scales with
  # the square root of that, so the two differ by sqrt(149.8/47.7) = 1.77x.
  testthat::expect_gt(min(abs(mn[1:2] / r[1:2])), 1.7)
  testthat::expect_equal(r[3], mn[3])
})


testthat::test_that("diet residuals and plot_diet_comp run on a fitted diet model", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Small multispecies model with diet (same fixture as test-diet-likelihood.R);
  # estimateMode = 3 builds the report (diet_hat) without optimizing.
  nyrs <- 10L; nspp <- 2L
  Fmort  <- c(seq(0.02, 0.3, length.out = nyrs / 2), seq(0.3, 0.05, length.out = nyrs / 2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  sim <- make_msm_test_data(
    years = seq_len(nyrs),
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),
    gam_a = c(1, 0.1), gam_b = rep(0.15, nspp),
    log_phi = matrix(c(-5, 0.5, -10, -2), nspp, nspp, byrow = TRUE))
  fit <- Rceattle::fit_mod(sim$data_list, estimateMode = 3, msmMode = 1,
                           suitMode = 4, niter = 5, initMode = "NonEquilibrium",
                           fit_control = fit_control(phase = FALSE, verbose = 0))

  testthat::expect_false(is.null(fit$quantities$diet_hat))
  r <- residuals(fit, type = "pearson", source = "diet")
  testthat::expect_setequal(unique(r$Source), "diet")
  testthat::expect_equal(nrow(r), nrow(fit$data_list$diet_data))
  # The diet Pearson matches the likelihood's own construction (finite rows).
  dd  <- fit$data_list$diet_data
  hat <- fit$quantities$diet_hat[, 2]
  # Read off the data_list itself. `bioenergetics_control` is a workbook sheet
  # name, not a data_list element, so reading the switch off one resolves to
  # NULL and scores every predator as a multinomial. Keeping this reference on
  # the real field is what stops it agreeing with that defect in residuals().
  fam <- fit$data_list$Diet_distribution
  testthat::expect_false(is.null(fam))
  fam <- as.integer(fam)[dd$Pred]
  w   <- fit$estimated_params$diet_comp_weights
  w   <- if (is.null(w)) 1 else as.numeric(w)[dd$Pred]
  fin <- is.finite(r$Residual)
  testthat::expect_equal(
    r$Residual[fin],
    diet_pearson_ref(dd, hat, off = fit$data_list$comp_offset,
                     family = fam, weight = w)[fin])

  # plot_diet_comp() now sources its residuals from residuals(source = "diet").
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    pf <- tempfile(fileext = ".pdf"); grDevices::pdf(pf)
    testthat::expect_error(plot_diet_comp(fit), NA)
    grDevices::dev.off()
  }
})


# A small multispecies diet fixture reused by the two OSA-diet tests below.
# Diet is multinomial (Diet_distribution = 0) with unit weights, so the OSA
# (conditional-binomial) decomposition of the diet likelihood must reproduce the
# ordinary multinomial diet likelihood slot exactly (see test (a)).
.make_diet_osa_fixture <- function() {
  nyrs <- 10L; nspp <- 2L
  Fmort  <- c(seq(0.02, 0.3, length.out = nyrs / 2),
              seq(0.3, 0.05, length.out = nyrs / 2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  make_msm_test_data(
    years   = seq_len(nyrs),
    Fmort   = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),
    gam_a   = c(1, 0.1), gam_b = rep(0.15, nspp),
    log_phi = matrix(c(-5, 0.5, -10, -2), nspp, nspp, byrow = TRUE))
}


testthat::test_that("OSA mode reproduces the diet fitting likelihood (decomposition invariant)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Diet has no WHAM reference (WHAM has no diet data), so we validate the diet
  # OSA path against the model itself. With keep == 1 the conditional-binomial
  # decomposition (dmultinom_osa) sums back to the joint multinomial density, so
  # building the model in OSA mode (osa_mode = 1, reading the diet composition
  # from obsvec) must reproduce the diet likelihood slot of ordinary fitting
  # (osa_mode = 0, reading diet_obs). This is a deterministic check that the diet
  # obsvec plumbing -- diet_obsvec_idx, the n_prey + 1 "other prey" segment, and
  # the obs<->pred alignment -- is wired correctly. The fixture is multinomial
  # diet with unit weights, so the weighted slot equals the unweighted joint.
  sim <- .make_diet_osa_fixture()

  # estimateMode = 3 builds and reports jnll_comp without optimizing; the
  # decomposition identity holds at any parameter values, so no fit is needed.
  # The diet obsvec segment (read when osa_mode = 1) is rebuilt on demand.
  fit <- Rceattle::fit_mod(sim$data_list, estimateMode = 3, msmMode = 1,
                           suitMode = 4, niter = 5, initMode = "NonEquilibrium",
                           fit_control = fit_control(phase = FALSE, verbose = 0))
  ctl <- Rceattle:::build_osa_data(fit$obj$env$data, build_osa = TRUE)$obs_ctl
  testthat::skip_if(!any(ctl$source == "diet"), "no diet observations")

  diet_row <- 19L                       # jnll_comp C++ slot 18 (diet) -> R row 19
  jc0  <- fit$obj$report(fit$obj$par)$jnll_comp                 # osa_mode = 0
  obj1 <- Rceattle:::.osa_build_obj(fit)                        # rebuild osa_mode = 1
  jc1  <- obj1$report(obj1$par)$jnll_comp

  testthat::expect_gte(nrow(jc0), diet_row)
  # The diet likelihood is active (non-trivial) ...
  testthat::expect_true(sum(jc0[diet_row, ]) != 0)
  # ... and the OSA decomposition reproduces it to ~machine precision.
  testthat::expect_equal(sum(jc1[diet_row, ]), sum(jc0[diet_row, ]),
                         tolerance = 1e-8)
})


testthat::test_that("osa_residuals(source = 'diet') runs end-to-end on a fitted diet model", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # End-to-end exercise of oneStepPredict() over the diet obsvec positions (never
  # run elsewhere in the suite). osa_residuals requires an optimized
  # (estimateMode < 3) fit; the simulated data are generated from the model with
  # a large effective sample size, so estimateMode = 1 converges. The diet
  # composition obsvec is rebuilt on demand by osa_residuals(). getsd = FALSE
  # skips sdreport -- OSA residuals use the fitted object, not the standard
  # errors, and this tiny fixture has a non-positive-definite Hessian that would
  # otherwise warn.
  sim <- .make_diet_osa_fixture()
  fit <- Rceattle::fit_mod(sim$data_list, estimateMode = 1, msmMode = 1,
                           suitMode = 4, niter = 5, initMode = "NonEquilibrium",
                           fit_control = fit_control(phase = FALSE, verbose = 0,
                                                     getsd = FALSE))
  ctl <- Rceattle:::build_osa_data(fit$obj$env$data, build_osa = TRUE)$obs_ctl
  testthat::skip_if(!any(ctl$source == "diet"), "no diet observations")

  osa <- osa_residuals(fit, source = "diet")

  testthat::expect_s3_class(osa, "rceattle_osa")
  testthat::expect_setequal(unique(osa$source), "diet")

  # Each stomach drops its sum-to-one "other prey" last bin, so the residual
  # count equals the kept (non-last) diet bins in obs_ctl.
  n_keep <- sum(ctl$source == "diet" & !ctl$is_last_bin)
  testthat::expect_equal(nrow(osa), n_keep)
  testthat::expect_true(all(is.finite(osa$residual)))

  # Deterministic for a fixed seed (the default continuous method has no RNG).
  osa2 <- osa_residuals(fit, source = "diet")
  testthat::expect_equal(osa$residual, osa2$residual)

  # Diagnostics run on the diet source and return a finite SDNR. The fixture's
  # diet proportions are noiseless (expected values, not multinomial draws), so
  # an SDNR ~ 1 calibration is not expected here -- that needs sampled data and
  # belongs in an offline validation; here we only guard that the machinery runs
  # and produces a finite statistic.
  dg <- osa_diagnostics(osa)
  testthat::expect_true(is.finite(dg$sdnr[dg$group == "all"]))
})

testthat::test_that("OSA supports an MVN index fleet without warning (per-family coverage in test-likelihood-osa-index-families.R)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # OSA residuals now cover every index family: the MVN covariance block is
  # whitened by the Cholesky of its covariance so oneStepPredict() can residualize
  # it (the correctness oracle lives in test-likelihood-osa-index-families.R). Here
  # we only pin that build_osa_data() no longer excludes/ warns for an MVN fleet
  # and that its index observations are laid into the OSA obsvec.
  nyrs <- 8; nages <- 5
  dat  <- make_test_data(nyrs = nyrs, nages = nages, seed = 42)
  sds  <- rep(20, nyrs); Rho <- matrix(0.3, nyrs, nyrs); diag(Rho) <- 1
  Sigma <- diag(sds) %*% Rho %*% diag(sds)
  srv  <- dat$fleet_control$Fleet_name == "Survey"       # Survey == Fleet_code 1
  dat$fleet_control$Index_distribution[srv] <- "MVN"
  dat$fleet_control$Catchability[srv]       <- "AnalyticalArith"
  dat$index_cov <- list(Survey = Sigma)

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0,
    fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

  # The OSA build no longer warns, and the MVN Survey's fitted observations get
  # whitened obsvec entries (one per fitted year).
  testthat::expect_no_warning(
    osa_dat <- Rceattle:::build_osa_data(fit$obj$env$data, build_osa = TRUE))
  testthat::expect_equal(sum(osa_dat$obs_ctl$source == "index"), nyrs)
  testthat::expect_true(all(osa_dat$index_obsvec_idx >= 0L))

  osa <- suppressWarnings(osa_residuals(fit, source = c("index", "comp")))
  testthat::expect_s3_class(osa, "rceattle_osa")
  testthat::expect_true(all(is.finite(osa$residual)))
  testthat::expect_true(any(osa$source == "index"))
})


testthat::test_that("a negative expected composition count warns", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Issue #108: "Sometimes the expected counts are negative. This only appears to
  # occur when the sample size is very small." Reproduced by sweeping
  # Sample_size: with the default method, `predicted` goes negative at N <= 2 and
  # is clean by N = 5 (with "oneStepGeneric" it persists to N = 5). The bin
  # counts are (proportion + comp_offset) * N and so already fractional; at a
  # tiny N they sit near zero and oneStepPredict()'s numerical conditional mean
  # can land just below it. An expected count cannot be negative, so say so
  # rather than let it be read as a fitted value.
  d <- make_test_data(nyrs = 20, nages = 10, seed = 123)
  d$fleet_control$Comp_distribution <- "Multinomial"
  d$comp_data$Sample_size <- 2
  fit <- suppressMessages(suppressWarnings(fit_mod(
    d, file = NULL, estimateMode = 1, random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

  testthat::expect_warning(
    osa <- suppressMessages(Rceattle::osa_residuals(fit, source = "comp")),
    "expected count cannot be negative")
  testthat::expect_true(any(osa$predicted < 0))

  # An adequate sample size does not warn.
  d$comp_data$Sample_size <- 100
  fit2 <- suppressMessages(suppressWarnings(fit_mod(
    d, file = NULL, estimateMode = 1, random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
  osa2 <- suppressMessages(Rceattle::osa_residuals(fit2, source = "comp"))
  testthat::expect_false(any(osa2$predicted < 0))
})


testthat::test_that("a failed parallel OSA loop rebuilds the object before retrying", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # oneStepPredict(parallel = TRUE) forks. A worker that aborts leaves the TMB
  # object's external pointers unusable in the parent, so the serial retry used
  # to re-enter the same object and take the whole R session with it ("An
  # irrecoverable exception occurred"). Reported against GOA pollock 2025, whose
  # ar1/rw catchability linkages carry 109 of its 164 random effects.
  #
  # The crash cannot be reproduced in a unit test -- it kills the process -- so
  # pin the property that prevents it: the retry builds a fresh object rather
  # than reusing the one the failed fork touched.
  d <- make_test_data(nyrs = 12, nages = 6, seed = 7)
  fit <- suppressMessages(suppressWarnings(fit_mod(
    d, file = NULL, estimateMode = 1, random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

  # Record the objects handed back, not just the number of calls: the builder
  # short-circuits to fit$obj for a model already in OSA mode, so counting calls
  # would pass even if the retry got the same object straight back.
  built <- list()
  real_build <- Rceattle:::.osa_build_obj
  testthat::local_mocked_bindings(
    .osa_build_obj = function(...) {
      o <- real_build(...)
      built[[length(built) + 1L]] <<- o
      o
    },
    .package = "Rceattle")

  # Fail the parallel attempt exactly as a dead worker does, and let the serial
  # retry through. `openmp()` has to be mocked as well: osa_residuals() probes
  # it and falls straight to serial when it errors, which is what happens once
  # several TMB models are loaded -- as they are when the whole suite runs.
  # Without this the parallel branch is never taken and the retry never fires,
  # so the test passes alone and fails in the suite.
  real_osp <- TMB::oneStepPredict
  testthat::local_mocked_bindings(
    openmp = function(...) 1L,
    oneStepPredict = function(..., parallel = FALSE) {
      if (isTRUE(parallel)) stop("non-numeric argument to mathematical function")
      real_osp(..., parallel = FALSE)
    },
    .package = "TMB")

  osa <- suppressWarnings(suppressMessages(
    Rceattle::osa_residuals(fit, source = "index", parallel = TRUE)))

  # One build up front, one more on the retry ...
  testthat::expect_gte(length(built), 2L)
  # ... and the retry's object is genuinely a different one, not the same
  # object handed back. This is what makes the retry safe.
  testthat::expect_false(identical(built[[1]], built[[length(built)]]))
  testthat::expect_true(nrow(osa) > 0)
  testthat::expect_true(all(is.finite(osa$residual)))
})

testthat::test_that(".osa_build_obj(force = TRUE) rebuilds an OSA-mode model", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Without `force`, a model already in OSA mode gets its own object back. That
  # is right for the first build and wrong for the retry, where reusing the
  # object the failed parallel attempt touched is the whole problem.
  d <- make_test_data(nyrs = 12, nages = 6, seed = 7)
  fit <- suppressMessages(suppressWarnings(fit_mod(
    d, file = NULL, estimateMode = 1, random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

  # Put the fit into OSA mode so the short-circuit applies.
  osa_dat <- Rceattle:::build_osa_data(fit$obj$env$data, build_osa = TRUE)
  fit$obj <- Rceattle:::.osa_build_obj(fit, osa_dat)
  testthat::expect_equal(fit$obj$env$data$osa_mode, 1)

  same  <- Rceattle:::.osa_build_obj(fit, osa_dat)
  fresh <- Rceattle:::.osa_build_obj(fit, osa_dat, force = TRUE)
  testthat::expect_true(identical(same, fit$obj))     # short-circuit
  testthat::expect_false(identical(fresh, fit$obj))   # genuinely rebuilt
})

# CAAL is the one composition source whose bins are ages within a length bin, so
# it is the case where the OSA and Pearson frames could disagree on index_label.
testthat::test_that("CAAL residuals label both frames as age bins", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data(nyrs = 8, nages = 5, seed = 1, growth = "vonBertalanffy")
  testthat::expect_gt(nrow(d$caal_data), 0)

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    growthFun = Rceattle::build_growth(fun = "vonBertalanffy"),
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = FALSE))))
  osa  <- suppressWarnings(Rceattle::osa_residuals(fit, source = "caal"))
  pear <- attr(osa, "pearson")

  testthat::expect_true(all(osa$source == "caal"))
  testthat::expect_identical(unique(osa$index_label), "age")
  testthat::expect_false(is.null(pear))
  testthat::expect_identical(unique(pear$index_label), "age")

  # The two frames differ only where the residual definitions differ: OSA has a
  # conditional `sd`, Pearson the sample size and the `assumed_sd` its own
  # likelihood assumes. The two sd columns are different quantities, so they
  # keep separate names rather than one being renamed onto the other.
  testthat::expect_identical(setdiff(names(osa), names(pear)), "sd")
  testthat::expect_setequal(setdiff(names(pear), names(osa)),
                            c("sample_size", "assumed_sd"))

  # One age per length group is fixed by sum-to-N, so OSA is that much shorter.
  n_len <- length(unique(pear$length))
  testthat::expect_equal(nrow(osa), nrow(pear) - n_len)
})


# A Gaussian method is continuous-only, so `discrete = TRUE` sends the
# composition rows to oneStepGeneric. The fallback happened but was neither
# announced nor written to the `method` attribute, so the attribute -- and
# print() -- named a method no composition row had used.
testthat::test_that("discrete = TRUE announces and records its method override", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  dat <- make_test_data(nyrs = 6, nages = 4, seed = 11)
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0,
    fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

  testthat::expect_message(
    osa <- suppressWarnings(osa_residuals(
      fit, source = "comp", discrete = TRUE,
      method = "oneStepGaussianOffMode")),
    "cannot be scored by a Gaussian")

  # Recorded as its own entry, leaving `default` as the method that was asked for.
  m <- attr(osa, "method")
  testthat::expect_equal(m[["default"]], "oneStepGaussianOffMode")
  testthat::expect_equal(m[["DiscreteComposition"]], "oneStepGeneric")
})
