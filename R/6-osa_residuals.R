# The one-step-ahead methods TMB::oneStepPredict() offers, split by whether they
# approximate the conditional distribution as Gaussian (so they need a
# continuous observation and produce a conditional mean and sd) or read the
# conditional CDF the template supplies.
.OSA_GAUSSIAN_METHODS <- c("oneStepGaussianOffMode", "oneStepGaussian",
                           "fullGaussian")
.OSA_METHODS <- c(.OSA_GAUSSIAN_METHODS, "oneStepGeneric", "cdf")

# Where a Dirichlet-multinomial composition goes when method = "cdf" is asked
# for: its conditional beta-binomial has no closed-form CDF, so the template
# supplies no CDF term and those bins have to be residualized some other way.
.OSA_CDF_FALLBACK <- "oneStepGaussianOffMode"

# Where a `method = "cdf"` residual saturates. The template squeezes each CDF to
# [2*DBL_EPSILON, 1 - 2*DBL_EPSILON] (osa_squeeze_cdf() in comp_osa.hpp), and
# oneStepPredict recovers F in double precision, so no residual can exceed this.
# 8.04; the absolute ceiling, from a CDF read straight off, is 8.21.
.OSA_CDF_CEILING <- stats::qnorm(1 - 2 * .Machine$double.eps)

# How many times to rebuild the model and redo the tail of a one-step-ahead
# sequence after a non-finite residual. Each retry costs the rows from the
# failure onwards, so this bounds the worst case at a few times the original call
# rather than leaving it unbounded. In practice the progress guard stops after
# one attempt whenever the retry does not help, which is every case measured.
.OSA_CDF_MAX_RETRY <- 5L


#' One-step-ahead (OSA) residuals for an Rceattle model
#'
#' @description
#' One-step-ahead residuals, also called forecast or quantile residuals
#' (Thygesen et al. 2017), computed via [TMB::oneStepPredict()]. Unlike Pearson
#' residuals they are iid standard normal under a correctly specified model even
#' when observations are correlated through composition bins or the model holds
#' random effects, so they support objective goodness-of-fit testing (Trijoulet
#' et al. 2023; Stewart and Monnahan 2025). The residualization runs inside the
#' assessment, so it also accounts for correlation induced by the model's own
#' random effects.
#'
#' @details
#' Computed *post hoc* and expensive -- TMB re-optimizes the random effects as
#' each observation is added -- so [fit_mod()] does not produce them, and the fit
#' must have been optimized at `estimateMode < 3`. Composition data are
#' decomposed into univariate conditional residuals (binomial / beta-binomial;
#' Trijoulet et al. 2023), and the last bin of each composition is fixed by the
#' sum-to-N constraint, so it has no residual (`NA`).
#'
#' Choosing between the methods, what each costs, how they behave under random
#' effects, the ceiling on `|residual|` under `"cdf"`, and where compositions at
#' scale defeat it are in `vignette("model-diagnostics")`, with the measurements
#' behind each.
#'
#' @param object A fitted `Rceattle` object, from [fit_mod()].
#' @param fit Deprecated name for `object`, still accepted; supplying both is an
#'   error.
#' @param source Observation sources to residualize: any of `"ecov"`, `"index"`,
#'   `"catch"`, `"comp"`, `"caal"`, `"diet"` or `"all"`, defaulting to the five
#'   non-diet sources. Mirrors the `source` argument of [residuals.Rceattle()];
#'   a source with no observations is skipped silently, `"diet"` is opt-in
#'   because it needs a multispecies model and is expensive, and `"ecov"` is the
#'   state-space covariate, residualized first against its own series as in
#'   WHAM's `make_osa_residuals()`.
#' @param method How the conditional distribution is read: one of
#'   `"oneStepGaussianOffMode"` (default, the WHAM/SAM choice),
#'   `"oneStepGaussian"`, `"fullGaussian"`, `"oneStepGeneric"` or `"cdf"`,
#'   passed to [TMB::oneStepPredict()].
#' @param discrete Whether to treat composition observations as the discrete
#'   counts they are. `NULL` (default) picks `TRUE` under `method = "cdf"`, where
#'   the residual is not standard normal without it, and `FALSE` otherwise.
#'   `TRUE` makes composition residuals randomized quantile residuals (Dunn and
#'   Smyth 1996) and so stochastic -- set `seed`. A Gaussian `method` cannot
#'   score a discrete observation, so that group falls back to
#'   `"oneStepGeneric"`; `"cdf"` with `discrete = FALSE` is allowed and says in a
#'   message that those residuals are biased up.
#' @param parallel Compute the per-observation loop with
#'   \code{\link[parallel]{mclapply}}, default `TRUE`, which is the main speedup
#'   for a model with random effects; set `options(mc.cores = )` to choose how
#'   many. Where a forked worker aborts instead of returning, the loop
#'   recomputes serially and prints the worker's own "irrecoverable exception"
#'   message, which comes from C and cannot be suppressed -- it does not mean the
#'   call failed.
#' @param seed Seed passed to [TMB::oneStepPredict()] for reproducible
#'   randomized-quantile residuals. Default `123`.
#' @param trace Print [TMB::oneStepPredict()] progress. Default `FALSE`.
#' @param ... Further arguments passed to [TMB::oneStepPredict()].
#'
#' @return A data frame of class `rceattle_osa`, one row per residualized
#'   observation, with columns `source`, `fleet`, `fleet_name`, `species`,
#'   `sex`, `year`, `age_length_bin`, `accumulated` (the bin folds neighbouring
#'   ages, so it stands for a range rather than the age named), `length` (the
#'   conditioning length bin for caal, `NA` otherwise), `index_label`,
#'   `observed`, `predicted`, `sd` and `residual`.
#'
#'   `observed` and `predicted` are on the residualization scale: log for
#'   lognormal catch and index, natural scale for `"Normal"` and
#'   `"TruncatedNormal"`, whitened by the lower Cholesky of the survey
#'   covariance for `"MVN"`/`"MVNORM"`, and bin counts for compositions. Three
#'   cases read differently, and none is an error:
#'
#'   * Under `method = "cdf"`, `predicted` is `NA` on every row because the
#'     method forms no conditional mode. `sd` is `NA` under the default method
#'     too; only `"oneStepGaussian"` returns one.
#'   * On a `"TruncatedNormal"` fleet under a Gaussian method, `predicted` is the
#'     truncated conditional mean `E[x | x > 0]`, which sits above the fitted
#'     index (163.8 against a fitted 100.0 on the package fixture), and `sd` is
#'     `NA`. The residual is unaffected.
#'   * A composition `predicted` can go slightly negative where a bin holds
#'     almost no fish, and the function warns, naming the count and the years.
#'     The bin's own count drives it, not the composition's sample size -- median
#'     observed count 0.05 on those rows against 4.9 elsewhere -- so the warning
#'     is not by itself evidence of thin data. Those residuals are biased
#'     positive.
#'
#'   Adding a `"TruncatedNormal"` fleet can move the *other* fleets' residuals on
#'   a random-effects model, because its rows interleave with the other index
#'   fleets by year instead of forming a contiguous block: a contiguous group
#'   reproduces a single call to 1.5e-14, while an interleaved split moved a
#'   residual by 5.8e-2. Both orderings are valid probability-integral-transform
#'   sequences, so neither is wrong. Fixed-effect models are unaffected.
#'
#'   Attributes carry `seed` and `method` -- the string that was passed, or a
#'   named vector where a group of rows ran under its own, such as
#'   `TruncatedNormal = "oneStepGeneric"`. Where composition types are present a
#'   `"pearson"` attribute holds the matching Pearson residuals, on proportions
#'   with the sample size in `sample_size` rather than the bin counts above, so
#'   the two are not comparable column by column.
#'
#' @references
#' Thygesen, U.H., et al. 2017. Validation of ecological state space models using
#'   the Laplace approximation. Environ. Ecol. Stat. 24:317-339.
#'
#' Trijoulet, V., et al. 2023. Model validation for compositional data in stock
#'   assessment models. Fish. Res. 257:106487.
#'
#' Stewart, I.J., and Monnahan, C.C. 2025. Diagnosing common sources of lack of
#'   fit to composition data using one-step-ahead residuals. Can. J. Fish. Aquat.
#'   Sci. 82:1-13.
#'
#' @seealso [osa_diagnostics()], [plot.rceattle_osa()], [process_residuals()],
#'   `vignette("model-diagnostics")`
#' @examples
#' \dontrun{
#' data(BS2017SS)
#' fit <- fit_mod(BS2017SS, estimateMode = "Hindcast")
#' osa <- osa_residuals(fit, source = c("index", "comp"))
#' plot(osa)
#' }
#' @export
osa_residuals <- function(object = NULL,
                          source   = c("ecov", "index", "catch", "comp", "caal"),
                          method   = "oneStepGaussianOffMode",
                          discrete = NULL,
                          parallel = TRUE,
                          seed     = 123,
                          trace    = FALSE,
                          ..., fit = NULL) {
  # `fit` was the old name for `object`; see R/0-deprecate.R.
  if (!missing(fit))
    object <- .rce_deprecated_arg(fit, !missing(object), "fit", "object", "osa_residuals")


  # ---- Validate the fit ----
  if (!inherits(object, "Rceattle")) {
    stop("`object` must be a fitted Rceattle model (from fit_mod()).")
  }
  if (is.null(object$obj)) {
    stop("`object` has no TMB object ($obj); OSA residuals require the fitted ",
         "model object.")
  }
  em <- object$data_list$estimateMode
  if (!is.null(em) && em >= 3) {
    stop("OSA residuals require a model optimized with estimateMode < 3 ",
         "(the returned objective for estimateMode >= 3 is a debug placeholder ",
         "that oneStepPredict cannot differentiate).")
  }
  if (is.null(object$obj$env$last.par.best)) {
    stop("`object` does not appear to have been optimized; OSA residuals require ",
         "a converged fit.")
  }

  # "diet" is supported but opt-in: it applies only to multispecies models with
  # estimated suitability and can be expensive, so it is not in the default set.
  # "all" is a synonym for every source including diet.
  valid_sources <- c("ecov", "index", "catch", "comp", "caal", "diet", "all")
  source <- match.arg(source, choices = valid_sources, several.ok = TRUE)
  if ("all" %in% source) source <- c("ecov", "index", "catch", "comp", "caal", "diet")

  # Check the method here rather than letting TMB reject it one observation group
  # at a time: the group split below reads it, so a typo would otherwise pick the
  # wrong split before failing.
  .method_defaulted <- missing(method)
  method <- match.arg(method, choices = .OSA_METHODS)

  # On composition data the package default is the method its own scoring table
  # rejects: residualized at the parameters that simulated the data it fails the
  # KS test on every replicate, where method = "cdf" passes (the method
  # comparison is in vignette("model-diagnostics")). The default stays put
  # because "cdf" returns non-finite
  # residuals in bulk on a deeply nested random-effects model, but a caller who
  # never chose a method should be told which one they got.
  if (.method_defaulted && any(c("comp", "caal", "diet") %in% source)) {
    message("osa_residuals(): composition residuals are being computed with the default ",
            "method = \"", method, "\", which is biased on composition data. ",
            "method = \"cdf\" is the only one that passes a self-test there; it can return ",
            "non-finite residuals on a model with many random effects. ",
            "See \"Choosing the one-step-ahead method\" in ",
            "vignette(\"model-diagnostics\").")
  }

  # TRUE under method = "cdf", FALSE otherwise, which leaves every other method
  # behaving as it always has. A composition bin holds a count, so its
  # conditional CDF is a step function and qnorm(F(x)) inherits the step,
  # E[F(X)] = 1/2 + sum(p^2)/2, biased up; randomizing over it -- qnorm(F(x) -
  # U f(x)), Dunn and Smyth (1996), as Trijoulet et al. (2023) prescribe --
  # removes it. The three options are scored in vignette("model-diagnostics")
  # (tools/verify/verify-osa-cdf.R).
  discrete_default <- is.null(discrete)
  if (discrete_default) discrete <- identical(method, "cdf")
  if (!is.logical(discrete) || length(discrete) != 1L || is.na(discrete)) {
    stop("`discrete` must be TRUE, FALSE, or NULL to choose per method.",
         call. = FALSE)
  }

  # Build the full OSA observation data (comp / caal / diet segments) on demand,
  # so nothing has to be set at fit time; any fit optimized at estimateMode < 3
  # will do, which the guards above enforce.
  # build_osa_data() reads only the *_ctl / *_obs arrays the model already
  # carries, and `obs_ctl` maps each obsvec position back to its source. The
  # same regenerated data is reused by .osa_build_obj().
  osa_dat <- build_osa_data(object$obj$env$data, build_osa = TRUE)
  obs_ctl <- osa_dat$obs_ctl

  # ---- Select the obsvec positions to residualize ----
  sel <- obs_ctl[obs_ctl$source %in% source & !obs_ctl$is_last_bin, , drop = FALSE]
  if (nrow(sel) == 0) {
    stop("No observations of source(s) ", paste(source, collapse = ", "),
         " are available for OSA residuals in this model.")
  }

  # ---- CONDITIONING ORDER (this is the one-step-ahead sequence) ----------------
  # `subset` order = the order oneStepPredict() conditions each observation on
  # the previously-residualized ones. It leaves the joint likelihood alone, but
  # it moves the per-bin composition residuals, whose within-multinomial
  # conditional binomials are sequenced in this order. A fixed-effects fit is
  # invariant to it, its observations being independent given the parameters.
  #
  # `fleet_code` is a tie-break, so the sequence is fully determined rather than
  # falling to the incidental row order of obs_ctl where several fleets report in
  # one (source, year). Under random effects that order shifts individual
  # residuals, though not the N(0,1) validity conclusion (Trijoulet et al. 2023).
  # Ascending fleet_code also matches WHAM's within-stage ordering (index fleets
  # before the fishery), so the WHAM cross-check is undisturbed.

  # Order by source, then year, fleet, bin. This reproduces WHAM's
  # make_osa_residuals() conditioning: the covariate (ecov) is residualized first
  # and standalone, then index/catch -> comp -> CAAL each conditional on the
  # earlier types.
  sel <- sel[order(match(sel$source, source), sel$year, sel$fleet_code,
                    sel$bin_index, na.last = TRUE), , drop = FALSE]

  # Rebuild at the fitted parameters with osa_mode >= 1, so the composition
  # (comp/caal) likelihoods read their counts from obsvec and use the unweighted
  # proper density oneStepPredict needs. The aggregate catch/index likelihood is
  # unchanged, so one rebuilt object serves every observation type -- including
  # mode 2, which only ADDS the conditional CDF terms method = "cdf" reads, and
  # so still serves the groups routed to a Gaussian method for want of a CDF.
  osa_mode <- if (identical(method, "cdf")) 2L else 1L
  obj_osa <- .osa_build_obj(object, osa_dat, osa_mode = osa_mode)

  # ---- Compute the OSA residuals ----
  # oneStepPredict() takes one `method` and one `discrete` per call, and neither
  # suits every observation type at once: aggregate index/catch are continuous,
  # while composition (comp/caal/diet) bins hold counts. So the rows are split
  # into groups sharing both and each group gets its own call. `subset` is
  # 1-based into obsvec (obs_pos is 0-based); `.row` maps a result back to its
  # `sel` row.
  get_col <- function(df, nm) {
    if (!is.null(df[[nm]])) df[[nm]] else rep(NA_real_, nrow(df))
  }
  is_comp      <- sel$source %in% c("comp", "caal", "diet")
  sel_discrete <- ifelse(is_comp, isTRUE(discrete), FALSE)

  # ---- The method each observation is actually residualized with --------------
  # Everything runs under `method`. Three groups cannot, because the method the
  # caller passed does not describe them, and each is split into its own
  # oneStepPredict() call. All three are announced and recorded on the returned
  # object's `method` attribute.
  fc_osa   <- object$data_list$fleet_control
  dat_osa  <- object$obj$env$data
  ill_osa  <- dat_osa$index_ll_type
  sel_method <- rep(method, nrow(sel))

  # (a) Dirichlet-multinomial compositions under method = "cdf". Their
  # conditional is a beta-binomial, whose CDF is not elementary and cannot be
  # summed at a fractional count, so there is no CDF term (comp_osa.hpp). The
  # miss is silent, not loud: both tails come back equal, so Fx = 0.5 and every
  # bin gets a residual of exactly 0. Hence these fleets take the Gaussian
  # default. Each composition source has its own likelihood-family vector,
  # indexed differently: comp and caal by fleet, diet by predator species. `1`
  # is Dirichlet-multinomial in all three.
  .is_dm <- function(rows, llt, key) {
    if (is.null(llt) || !length(llt)) return(rep(FALSE, length(rows)))
    # `key >= 1` is not decoration: a 0 index drops the element rather than
    # returning NA, which would shift every later lookup by one silently.
    hit <- rows & !is.na(key) & key >= 1L & key <= length(llt)
    hit & c(as.integer(llt), NA_integer_)[ifelse(hit, key, length(llt) + 1L)] == 1L
  }
  sel_dm <- rep(FALSE, nrow(sel))
  if (identical(method, "cdf")) {
    sel_dm <-
      .is_dm(sel$source == "comp", dat_osa$comp_ll_type, sel$fleet_code) |
      .is_dm(sel$source == "caal", dat_osa$caal_ll_type, sel$fleet_code) |
      .is_dm(sel$source == "diet", dat_osa$diet_ll_type, sel$species)
    sel_method[sel_dm] <- .OSA_CDF_FALLBACK
    # ... on the continuous approximation, which is exactly the treatment those
    # fleets get when `method` is left alone. Carrying `discrete` across would
    # send them to oneStepGeneric's numerical integration instead -- a different
    # method again, several times slower, and not one this change has measured.
    sel_discrete[sel_dm] <- FALSE
    if (any(sel_dm)) {
      dm_fleets <- unique(as.character(
        fc_osa$Fleet_name[match(sel$fleet_code[sel_dm & sel$source != "diet"],
                                fc_osa$Fleet_code)]))
      dm_fleets <- c(dm_fleets, if (any(sel_dm & sel$source == "diet"))
        paste("diet species", sort(unique(sel$species[sel_dm & sel$source == "diet"]))))
      message("osa_residuals(): ", paste(stats::na.omit(dm_fleets), collapse = ", "),
              " use a Dirichlet-multinomial, whose conditional beta-binomial has ",
              "no closed-form CDF. Those ", sum(sel_dm), " observation(s) are ",
              "residualized with method = \"", .OSA_CDF_FALLBACK,
              "\" instead of \"cdf\". See ?osa_residuals.")
    }
  }

  # Overriding the default to FALSE under "cdf" leaves the step in qnorm(F(x)),
  # which is the worst-calibrated of the three composition options measured
  # (mean +0.610, sd 1.262 where the answer is known to be standard normal). It
  # is a legitimate thing to ask for, so it runs; every other resolution here
  # that moves numbers announces itself, and so does this one. Counted after the
  # Dirichlet-multinomial split above, whose rows are continuous either way.
  if (identical(method, "cdf") && !discrete_default && !isTRUE(discrete) &&
      any(is_comp & !sel_dm)) {
    message("osa_residuals(): discrete = FALSE with method = \"cdf\" leaves the ",
            "step in the composition transform. Those ", sum(is_comp & !sel_dm),
            " residual(s) are biased up and over-dispersed (mean +0.61, sd 1.26 ",
            "on a self-test). See ?osa_residuals.")
  }

  # (b) Discrete compositions. The Gaussian methods are continuous-only, so
  # `discrete = TRUE` needs a CDF-based method: "cdf" already is one, and any
  # Gaussian choice falls back to the generic (numerically integrated) one.
  sel_gauss_disc <- sel_discrete & sel_method %in% .OSA_GAUSSIAN_METHODS
  if (any(sel_gauss_disc)) {
    message("osa_residuals(): discrete = TRUE cannot be scored by a Gaussian ",
            "method, which is continuous-only. Those ", sum(sel_gauss_disc),
            " composition observation(s) are residualized with method = ",
            "\"oneStepGeneric\" instead of \"", method, "\". See ?osa_residuals.")
  }
  sel_method[sel_gauss_disc] <- "oneStepGeneric"

  # (c) `Index_distribution = "TruncatedNormal"` (index_ll_type 4) under a
  # Gaussian method is residualized on its own with oneStepGeneric over a
  # (0, Inf) range, which normalizes the density by its integral and so gives the
  # truncated CDF F(x) = [Phi((x-mu)/sd) - Phi(-mu/sd)] / Phi(mu/sd), whose
  # qnorm is standard normal. The Gaussian methods cannot see the truncation --
  # it enters the density only through log Phi(mu/sd), a function of the
  # prediction and not the observation -- and return the untruncated residual
  # wherever truncation carries real mass. method = "cdf" needs none of this: the
  # template returns that same F in closed form, so the family runs exactly in
  # the main call.
  #
  # The range belongs to the FAMILY, so only these rows may carry it: "Normal" is
  # genuinely untruncated, and a lognormal fleet's obsvec entry is log(obs),
  # negative for a small index. Keyed off index_ll_type, the same vector
  # build_osa_data() read to decide obs vs log(obs), so the two cannot disagree.
  sel_trunc <- rep(FALSE, nrow(sel))
  if (!is.null(ill_osa) && length(ill_osa)) {
    sel_trunc <- sel$source == "index" & ill_osa[sel$fleet_code] %in% 4L
    sel_trunc[is.na(sel_trunc)] <- FALSE
    sel_trunc <- sel_trunc & sel_method != "cdf"
  }
  sel_method[sel_trunc] <- "oneStepGeneric"
  sel_group <- paste0(sel_method, "|", sel_discrete, "|", sel_trunc)

  # Each group gets its own oneStepPredict() call, handed everything before its
  # FIRST row as `conditional` (see .run_osp()), which reproduces a single call
  # only while the group is a CONTIGUOUS block of `sel`. Nothing above makes it
  # one: `sel` follows the caller's `source` order, so source = c("comp",
  # "index", "caal") leaves the index rows in a hole inside the discrete group,
  # whose first row is then row 1 and whose `conditional` comes back empty --
  # silently dropping the index data terms from the compositions' conditioning.
  # Interleaved Dirichlet-multinomial and multinomial fleets do the same by year.
  #
  # So sort the groups contiguous, keeping each group's internal order and
  # ordering the groups by where they first appear. That moves the one-step-ahead
  # SEQUENCE where a group would have been interleaved; every ordering is a valid
  # probability-integral-transform sequence (Trijoulet et al. 2023), and only a
  # contiguous one has a split that can be made invisible. `sel_pos` records each
  # row's place before the reorder, so the frame a caller gets keeps the order it
  # was requested in; only the CONDITIONING changes.
  sel_pos <- seq_len(nrow(sel))
  if (length(unique(sel_group)) > 1L) {
    ord <- order(match(sel_group, unique(sel_group)))   # stable: keeps within-group order
    if (!identical(ord, seq_along(ord))) {
      sel          <- sel[ord, , drop = FALSE]
      sel_group    <- sel_group[ord]
      sel_method   <- sel_method[ord]
      sel_discrete <- sel_discrete[ord]
      sel_trunc    <- sel_trunc[ord]
      sel_dm       <- sel_dm[ord]
      sel_pos      <- sel_pos[ord]
    }
  }
  stopifnot(!any(diff(match(sel_group, unique(sel_group))) < 0))

  # Say so. `method` is overridden for these rows whatever the caller passed, and
  # a silent override is exactly the kind of thing that makes a Q-Q plot hard to
  # account for later. Also flag the cost: the exact integration is several times
  # slower per observation than a Gaussian approximation.
  if (any(sel_trunc)) {
    trunc_fleets <- unique(as.character(
      fc_osa$Fleet_name[match(sel$fleet_code[sel_trunc], fc_osa$Fleet_code)]))
    message("osa_residuals(): fleet(s) ", paste(trunc_fleets, collapse = ", "),
            " use Index_distribution = \"TruncatedNormal\", whose truncation a ",
            "Gaussian method cannot see. Those ", sum(sel_trunc), " observation(s) ",
            "are residualized with method = \"oneStepGeneric\" over the ",
            "truncated support instead of \"", method,
            "\" -- exact, and slower. See ?osa_residuals.")
  }

  # oneStepPredict(parallel = TRUE) calls TMB::openmp(), which can only resolve
  # the active model when a single TMB DLL is loaded; with several loaded (e.g.
  # Rceattle alongside WHAM, or the full test suite) it errors with "Multiple TMB
  # models loaded". Probe openmp() once and silently fall back to serial when it
  # is unavailable, so parallel = TRUE stays a safe default everywhere.
  parallel_ok <- isTRUE(parallel) &&
    !inherits(try(TMB::openmp(), silent = TRUE), "try-error")
  if (isTRUE(parallel) && !parallel_ok) {
    message("osa_residuals(): parallel unavailable (multiple TMB models loaded?); ",
            "computing one-step-ahead residuals serially.")
  }

  # Upper limit for the truncated group's exact integration. c(0, Inf) is the
  # support, but integrate()'s infinite-limit substitution probes far out in the
  # tail, where a random-effects model's Laplace inner problem does not converge:
  # the integrand comes back non-finite, integrate() errors, and TMB's try()
  # writes NA. On a 15-year fixture with random recruitment that lost 5 of 15
  # residuals, silently, under a warning blaming the fit's convergence. So bound
  # it: ten standard deviations above every fitted index in the group and above
  # every observation being residualized, which discards under 1e-23 of the mass.
  # oneStepGeneric normalizes over the range it is given, so the range need only
  # hold the observation and effectively all the density, not reach infinity.
  .trunc_range <- function(rows) {
    dr <- sel$data_row[rows]
    mu <- suppressWarnings(as.numeric(object$quantities$index_hat[dr]))
    sg <- suppressWarnings(as.numeric(.observation_sd(object$quantities, "index")[dr]))
    ob <- suppressWarnings(as.numeric(osa_dat$obsvec[sel$obs_pos[rows] + 1L]))
    hi <- suppressWarnings(max(c(mu + 10 * sg, ob * 1.5), na.rm = TRUE))
    if (!is.finite(hi) || hi <= 0) hi <- Inf   # nothing usable; fall back
    c(0, hi)
  }

  .run_osp <- function(rows, dsc, meth, trunc = FALSE, spline = FALSE) {
    # `meth` is the group's own method, already resolved above -- a Gaussian
    # method never reaches a discrete group, and the families that cannot use the
    # caller's choice have been split off. Which groups run in parallel is
    # decided at `want_par` below.
    osp <- function(par) {
      args <- list(
        obj                 = obj_osa,
        observation.name    = "obsvec",
        data.term.indicator = "keep",
        method              = meth,
        # Only the truncated family restricts the support; every other group
        # keeps TMB's default (-Inf, Inf).
        range               = if (trunc) .trunc_range(rows) else c(-Inf, Inf),
        subset              = sel$obs_pos[rows] + 1L,
        # Everything earlier in the sequence than this group's first observation
        # is CONDITIONED ON, not discarded: oneStepPredict otherwise marks every
        # row it is not residualizing as unconditional and zeroes its data term,
        # so the group split would change the latent states' conditional
        # distribution and move the residuals. On a 21-random-effect fixture at
        # discrete = FALSE, 1.5e-14 from a single call with this, 5.8e-2 without.
        #
        # The contiguous sort above is what makes `min(rows) - 1` the whole of
        # what precedes the group. It matches a single call on the CONDITIONING
        # only: under discrete = TRUE oneStepPredict re-seeds and draws its
        # randomizing uniforms per call, so a group's residuals also depend on
        # its row count. Fixed-effect models are unaffected either way, their
        # observations being independent given the parameters.
        conditional         = if (min(rows) > 1L) sel$obs_pos[seq_len(min(rows) - 1L)] + 1L
                              else numeric(0),
        discrete            = dsc,
        parallel            = par,
        seed                = seed,
        trace               = trace)
      dots <- list(...)
      if (length(dots) && (is.null(names(dots)) || !all(nzchar(names(dots))))) {
        stop("osa_residuals(): all arguments passed through `...` to ",
             "oneStepPredict() must be named; an unnamed one would be matched ",
             "positionally against an argument this function already sets.",
             call. = FALSE)
      }
      clash <- intersect(names(dots), names(args))
      if (length(clash)) {
        warning("osa_residuals(): ignoring ", paste(clash, collapse = ", "),
                " passed through `...` -- ", if (length(clash) > 1) "these are"
                else "this is", " set per observation group by this function ",
                "(see ?osa_residuals). Use the `method` argument to choose a ",
                "method.", call. = FALSE)
      }
      args <- c(args, dots[setdiff(names(dots), names(args))])
      # `range` binds the integration limits only where oneStepGeneric integrates
      # the density itself. Under TMB's default splineApprox = TRUE it splines a
      # tmbprofile() slice and integrates over whatever range that slice covered,
      # which does not respect the (0, Inf) support: on a fixture at mu = x =
      # 100, sd = 150 it returns -0.790 where the exact transform is -0.437.
      # The exact path is what makes this group worth splitting out; `spline`
      # exists so the group loop can retry with the approximation when it fails.
      if (trunc) args$splineApprox <- isTRUE(spline)
      do.call(TMB::oneStepPredict, args)
    }

    # oneStepPredict(parallel = TRUE) forks via mclapply, and some
    # model/observation combinations abort the child rather than return, leaving
    # an error object where a gradient should be -- which surfaces from deep in
    # TMB as "non-numeric argument to mathematical function". Retry serially:
    # slower, but it returns residuals instead of a message pointing nowhere near
    # the cause. The retry must build a FRESH object; re-entering the one the
    # failed attempt used ends the R session ("An irrecoverable exception
    # occurred") rather than recovering. On GOA pollock 2025 (164 random effects,
    # 109 of them from an ar1 and a rw catchability linkage) forked workers abort
    # at any width above one and mclapply normally absorbs it -- a direct call
    # returned all 90 index observations at 1, 2, 4, 8 and 11 cores -- so the
    # case that surfaces here looks load-dependent rather than a property of the
    # model. Serial never fails.
    #
    # method = "cdf" is reproducible in parallel even when discrete: the forked
    # workers only evaluate the objective, and oneStepPredict draws the
    # randomizing uniforms serially under `seed` once they have all returned. The
    # generic method is left serial, where that has not been measured.
    want_par <- parallel_ok && (!dsc || identical(meth, "cdf"))
    res <- if (!want_par) osp(FALSE) else
      tryCatch(osp(TRUE), error = function(e) {
        message("osa_residuals(): the parallel one-step-ahead loop failed (",
                conditionMessage(e), "); recomputing serially.")
        obj_osa <<- .osa_build_obj(object, osa_dat, osa_mode = osa_mode, force = TRUE)
        osp(FALSE)
      })
    # oneStepGeneric and "cdf" return no `observation` column, so those groups
    # would otherwise report NA for a value that is simply the obsvec entry that
    # was residualized.
    obs_col <- get_col(res, "observation")
    if (all(is.na(obs_col))) {
      obs_col <- as.numeric(osa_dat$obsvec[sel$obs_pos[rows] + 1L])
    }
    # `predicted` is a conditional MODE, which only the Gaussian methods form, so
    # under method = "cdf" it is NA on every row -- including the
    # Dirichlet-multinomial fleets sent to a Gaussian method, which would
    # otherwise be the only rows carrying one, and are also the rows where an
    # expected count goes negative. One column must not mean a conditional mode
    # on half an object and nothing on the rest.
    blank <- identical(method, "cdf")
    data.frame(.row = rows, observed = obs_col,
               predicted = if (blank) NA_real_ else get_col(res, "mean"),
               sd        = if (blank) NA_real_ else get_col(res, "sd"),
               residual  = res$residual)
  }
  osa <- do.call(rbind, lapply(unique(sel_group), function(g) {
    rows  <- which(sel_group == g)
    trunc <- sel_trunc[rows[1]]
    meth  <- sel_method[rows[1]]
    res   <- .run_osp(rows, dsc = sel_discrete[rows[1]], meth = meth, trunc = trunc)

    # The exact integration evaluates the Laplace marginal at arbitrary values of
    # the observation, where a random-effects model's inner Newton problem does
    # not always converge: integrate() errors, TMB writes NA, and the rows vanish
    # -- worse than an approximate residual, since osa_diagnostics() then passes
    # verdict on whatever survived without saying how many it lost. Retry under
    # TMB's spline approximation, robust but blind to the (0, Inf) support, so
    # approximate in the same direction the Gaussian methods are. The whole group
    # is retried, not the failed rows, so a fleet keeps one method.
    if (trunc && any(!is.finite(res$residual))) {
      n_bad <- sum(!is.finite(res$residual))
      warning("osa_residuals(): exact integration of the TruncatedNormal ",
              "density failed for ", n_bad, " of ", nrow(res), " observation(s) ",
              "-- the Laplace inner problem does not converge across the whole ",
              "support, which happens on models with random effects. The fleet ",
              "was recomputed with TMB's spline approximation, so ITS RESIDUALS ",
              "ARE APPROXIMATE and do not carry the truncation exactly. Treat ",
              "them as indicative; the other fleets are unaffected.",
              call. = FALSE)
      res <- .run_osp(rows, dsc = sel_discrete[rows[1]], meth = meth,
                      trunc = TRUE, spline = TRUE)
      attr(res, "approx") <- TRUE
    }

    # Retry the tail from the first non-finite residual. It does NOT fix the
    # bulk failure (.osa_retry_tail() records the measurement), which is driven
    # by the depth of the conditioning rather than a poisoned warm start; it is
    # kept for a genuinely transient failure and costs nothing unless one occurs.
    # "cdf" only: `oneStepGaussianOffMode` walks the sequence REVERSED
    # (`reverse = (method == "oneStepGaussianOffMode")` is oneStepPredict()'s
    # formal default), so a tail retry would redo the wrong end, and
    # `oneStepGaussian` and `fullGaussian` never override the warm start at all.
    if (identical(meth, "cdf") && length(object$obj$env$random) &&
        any(!is.finite(res$residual))) {
      res <- .osa_retry_tail(res, function(from)
        .run_osp(rows[from:length(rows)], dsc = sel_discrete[rows[1]],
                 meth = meth, trunc = FALSE))
    }
    res
  }))
  osa <- osa[order(osa$.row), , drop = FALSE]   # restore chronological 'sel' order
  index_label <- c("age", "length")[sel$comp_type + 1L]   # NA for aggregates
  fc <- object$data_list$fleet_control                        # fleet code -> name
  fleet_name <- if (!is.null(fc)) {
    fc$Fleet_name[match(sel$fleet_code, fc$Fleet_code)]
  } else NA_character_

  out <- data.frame(
    source         = sel$source,
    fleet          = sel$fleet_code,
    fleet_name     = fleet_name,
    species        = sel$species,
    sex            = sel$sex,
    year           = sel$year,
    age_length_bin = sel$age_length_bin,
    accumulated    = if (is.null(sel$accumulated)) FALSE else sel$accumulated,
    length         = sel$length,
    index_label    = index_label,
    observed      = osa$observed,
    predicted     = osa$predicted,
    sd            = osa$sd,
    residual      = osa$residual,
    stringsAsFactors = FALSE)

  # Undo the grouping reorder, so the rows come back in the source/year/fleet/bin
  # order they were requested in. The reorder above buys a correct conditioning
  # sequence; it should not also rearrange the answer.
  out <- out[order(sel_pos), , drop = FALSE]

  rownames(out) <- NULL
  class(out) <- c("rceattle_osa", "data.frame")
  # Record what was actually used, not what was asked for: three groups are
  # residualized with their own method regardless of `method`, so a single
  # string would misdescribe those rows. Stays a plain string when nothing was
  # overridden.
  attr(out, "method") <- if (any(sel_trunc) || any(sel_dm) || any(sel_gauss_disc)) {
    c(default = method,
      if (any(sel_trunc))     c(TruncatedNormal = "oneStepGeneric"),
      if (any(sel_dm))        c(DirichletMultinomial = .OSA_CDF_FALLBACK),
      if (any(sel_gauss_disc)) c(DiscreteComposition = "oneStepGeneric"))
  } else method
  attr(out, "seed")   <- seed
  # Randomized composition residuals carry a draw, so whether they were
  # randomized is part of reading them, and `discrete` resolves per method rather
  # than being whatever the caller passed. Named the way `method` is when a
  # family was treated differently, so a mixed model cannot report a single flag
  # that is wrong for half its fleets. Gated on what was actually randomized: a
  # cdf call over `source = "index"` randomizes nothing however `discrete`
  # resolved, and nor does a model whose every composition fleet is
  # Dirichlet-multinomial, those rows going to a Gaussian method.
  attr(out, "discrete") <- if (!any(sel_discrete)) {
    FALSE
  } else if (any(sel_dm)) {
    c(default = TRUE, DirichletMultinomial = FALSE)
  } else TRUE
  # Per-species bin counts, so plot() can split joint-sex (Sex == 3) composition
  # bins onto a single age/length axis (males are stored as bins nbin+1 .. 2*nbin).
  attr(out, "nages")    <- object$data_list$nages
  attr(out, "nlengths") <- object$data_list$nlengths

  # Attach the matching Pearson residuals for composition sources so the
  # plot() method can show OSA and Pearson bubbles side by side.
  comp_types <- intersect(unique(out$source), c("comp", "caal"))
  if (length(comp_types) > 0) {
    pear <- tryCatch(
      stats::residuals(object, type = "pearson", source = comp_types),
      error = function(e) NULL)
    # residuals() is a general-purpose method and names its columns in the
    # data-sheet style (Fleet_code, Year, Observed, ...); this data frame names
    # them in the style of the object it is attached to. Carrying both
    # conventions on one object means a reader has to know which half they are
    # holding, so rename to match here. `Sd` becomes `assumed_sd`, not `sd`:
    # the OSA frame's `sd` is the conditional sd of the OSA residual.
    if (!is.null(pear)) {
      nm <- c(Source = "source", Fleet_code = "fleet", Fleet_name = "fleet_name",
              Species = "species", Sex = "sex", Year = "year",
              Age0_Length1 = "index_label", Bin = "age_length_bin",
              Length = "length", Sample_size = "sample_size",
              Accumulated = "accumulated",
              Observed = "observed", Fitted = "predicted",
              Residual = "residual", Sd = "assumed_sd")
      hit <- names(pear) %in% names(nm)
      names(pear)[hit] <- unname(nm[names(pear)[hit]])

      # Same form as the OSA frame, not a second encoding under a second name.
      if (!is.null(pear$index_label)) {
        pear$index_label <- ifelse(pear$source == "caal", "age",
                                   ifelse(!is.na(pear$index_label) &
                                            pear$index_label == 1,
                                          "length", "age"))
      }
    }
    attr(out, "pearson") <- pear
  }

  # A `method = "cdf"` residual is qnorm of a CDF read in double precision, so it
  # is censored at 8.04 rather than merely large, and osa_diagnostics()
  # summarises the censored values. The warning names oneStepGaussian rather than
  # "a Gaussian method" because it alone reports such a row uncensored -- the
  # package default returns NaN there and oneStepGeneric compresses it further;
  # measured in vignette("model-diagnostics") and
  # tools/verify/verify-osa-cdf-accuracy.R.
  n_sat <- if (identical(method, "cdf")) {
    sum(abs(out$residual) >= .OSA_CDF_CEILING - 1e-6, na.rm = TRUE)
  } else 0L
  if (n_sat > 0) {
    warning(n_sat, " of ", nrow(out), " one-step-ahead residual(s) are at the ",
            "method = \"cdf\" ceiling of ", formatC(.OSA_CDF_CEILING, digits = 4,
                                                   format = "f"),
            ", where the conditional CDF reaches the last double below 0 or 1. ",
            "Those observations are further from the model than this method can ",
            "report, so treat them as censored -- osa_diagnostics() summarises ",
            "the censored values. method = \"oneStepGaussian\" reports them ",
            "uncensored; the package default returns NaN there and ",
            "\"oneStepGeneric\" compresses them further. See ?osa_residuals.",
            call. = FALSE)
  }

  n_bad <- sum(!is.finite(out$residual))
  if (n_bad > 0) {
    # The "cdf" clause is shown only where it can apply. Under a Gaussian method
    # the same fit returns every residual finite, so leading with it there would
    # explain away a result that does point at the data or the fit.
    warning(n_bad, " of ", nrow(out), " OSA residual(s) are non-finite. Check ",
            "model convergence and the sparsest compositions before ",
            "interpreting the residuals: non-finite values usually point to ",
            "very sparse or degenerate compositions (common for conditional ",
            "age-at-length), or to a poorly converged fit.",
            if (identical(method, "cdf") && length(object$obj$env$random))
              paste0(" On a random-effects fit with a large composition data ",
                     "set this can instead be the documented limitation of ",
                     "method = \"cdf\" -- the Laplace inner problem fails on ",
                     "the depth of conditioning -- which says nothing about ",
                     "the fit.")
            else "",
            " See ?osa_residuals.")
  }

  # A composition `predicted` is an expected bin count, so it cannot be negative.
  # It goes slightly negative where a bin holds almost no fish: observations
  # enter as counts, (proportion + comp_offset) * N, and oneStepPredict()'s
  # conditional mean is a numerical step from the observation that overshoots
  # below zero when the count is near it. The BIN's count drives it, not the
  # composition's sample size -- on EBS pollock the negative rows' median
  # observed count is 0.05 against 4.9 elsewhere, while their sample sizes span
  # the same 1 to 821 as the rest and 69 of 353 sit above 100. Reported rather
  # than clamped: the negative value is the signal that the bin is too sparse to
  # decompose, and clamping hides it.
  is_comp <- out$source %in% c("comp", "caal", "diet")
  bad <- is_comp & is.finite(out$predicted) & out$predicted < 0
  if (any(bad)) {
    yrs <- sort(unique(out$year[bad]))
    warning(sum(bad), " composition `predicted` value(s) are negative, in ",
            "year(s) ", paste(utils::head(yrs, 5), collapse = ", "),
            if (length(yrs) > 5) ", ..." else "",
            ". An expected count cannot be negative: those bins hold too few ",
            "fish for the one-step-ahead decomposition to describe, so treat ",
            "`predicted` there as uninformative and their residuals as biased ",
            "positive. See ?osa_residuals.", call. = FALSE)
  }
  out
}


#' Redo the tail of a one-step-ahead sequence after a failed Laplace solve
#'
#' @description
#' `method = "cdf"` on a random-effects model can return non-finite residuals in
#' bulk: measured on BS2017SS with random recruitment, 1879 of 4538 composition
#' bins, against 0 for the Gaussian method. The failures are a contiguous tail,
#' and the same 1880 rows residualized ON THEIR OWN give 1 failure, so they are
#' not intrinsically bad observations.
#'
#' The hypothesis this helper was written for was a warm-start cascade: TMB's cdf
#' loop is `nll <- fn(observation(k)); lp <- env$last.par; ...; env$last.par <- lp`
#' (TMB 1.9.21), capturing the warm start AFTER the evaluation, so a NaN solve
#' would be restored as the start for everything after it. **Measurement refuted
#' that.** Redoing the tail on a fresh call, which is what clears any poisoned
#' warm start, recovers nothing: 1879 before, 1879 after. What differs between
#' the failing call and the successful isolated one is not the warm start but the
#' DEPTH OF CONDITIONING (2658 prior observations against none), so the Laplace
#' inner problem is failing on the conditioning itself.
#'
#' The retry is kept because it is cheap in the common case (it runs only when
#' something is already non-finite) and it does recover a genuinely transient
#' failure, but **no such case has been observed**: see `?osa_residuals` for
#' the limitation this leaves, and `inst/dev/SESSION_HANDOFF.md` for where to
#' take it next.
#'
#' @param res Result frame from one group, with a `residual` column.
#' @param rerun Function of a starting row index, returning a frame of the same
#'   shape for rows `from:nrow(res)`.
#' @param max_try Attempts before giving up. A retry can fail again further
#'   along, so each round must strictly reduce the number of bad rows.
#' @return `res`, with the recomputed tail spliced in. Messages only if it
#'   recovered something, or if it could not.
#' @keywords internal
#' @noRd
.osa_retry_tail <- function(res, rerun, max_try = .OSA_CDF_MAX_RETRY) {
  n_bad0 <- sum(!is.finite(res$residual))
  if (n_bad0 == 0L) return(res)
  for (attempt in seq_len(max_try)) {
    bad <- which(!is.finite(res$residual))
    if (!length(bad)) break
    from <- min(bad)
    again <- rerun(from)
    # Require progress: a retry that fails at the same place or earlier is not
    # going to converge by being repeated, and would loop to the bound.
    if (sum(!is.finite(again$residual)) >= length(bad)) break
    res[from:nrow(res), ] <- again
  }
  n_left <- sum(!is.finite(res$residual))
  if (n_bad0 > n_left) {
    message("osa_residuals(): a Laplace solve failed partway through the ",
            "one-step-ahead sequence, which under method = \"cdf\" costs every ",
            "observation after it. The block from the first failure was ",
            "recomputed, recovering ", n_bad0 - n_left, " of ", n_bad0,
            " observation(s)",
            if (n_left) paste0("; ", n_left, " remain non-finite") else "", ".")
  } else {
    # Recovering none rules out a poisoned warm start: the failure is in the
    # solve itself, not in the sequence. On a large composition data set that is
    # the known "cdf" limitation; on a small one it still points at the fit.
    message("osa_residuals(): ", n_bad0, " observation(s) are non-finite and ",
            "recomputing the block from the first one recovered none, so the ",
            "Laplace inner problem fails on the conditioning itself rather than ",
            "on a poisoned warm start. On a random-effects fit with a large ",
            "composition data set that is the measured limitation of ",
            "method = \"cdf\" rather than a sign of a bad fit; on a small one, ",
            "check convergence and the sparsest compositions. See ?osa_residuals.")
  }
  res
}


#' Rebuild a fitted Rceattle TMB object in OSA mode at the fitted parameters
#'
#' @description
#' Returns a TMB ADFun object equivalent to `fit$obj` but with `osa_mode >= 1`,
#' built at the fitted parameter values and the same map / random-effect
#' structure, ready for [TMB::oneStepPredict()]. In OSA mode the composition
#' likelihoods read their counts from `obsvec` and use unweighted proper
#' densities; the aggregate catch/index likelihood is identical in both modes,
#' so this single object serves every observation type.
#'
#' @param fit A fitted `Rceattle` object.
#' @param osa_dat Optional pre-built OSA observation data (the list returned by
#'   [build_osa_data()] with `build_osa = TRUE`) to reuse instead of rebuilding
#'   it. `NULL` (the default) rebuilds it from `fit$obj$env$data`.
#' @param osa_mode 1 (default) for the keep-gated densities the Gaussian and
#'   generic methods need; 2 to add the conditional CDF terms
#'   `method = "cdf"` reads. Mode 2 evaluates a `pbeta` per composition bin, so
#'   it is asked for rather than assumed.
#' @param force Build a new object even when `fit` is already in the requested
#'   OSA mode, where this otherwise returns `fit$obj` itself. The retry after a
#'   failed parallel one-step-ahead loop needs a genuinely new one.
#' @return A TMB ADFun object with the requested `osa_mode`.
#' @keywords internal
#' @noRd
.osa_build_obj <- function(fit, osa_dat = NULL, osa_mode = 1L, force = FALSE) {
  obj <- fit$obj
  osa_mode <- as.integer(osa_mode)[1]
  # A model already fitted in the requested OSA mode can use its own object --
  # except when the caller needs a genuinely new one. The retry after a failed
  # parallel loop does: handing back `fit$obj` there would reuse the object the
  # failure touched, which is the thing that ends the R session.
  if (!force && !is.null(obj$env$data$osa_mode) &&
      obj$env$data$osa_mode == osa_mode) {
    return(obj)
  }
  # Regenerate the full OSA observation vector (comp / CAAL / diet segments) on
  # demand, so residuals need nothing set at fit time. build_osa_data() reads
  # only the *_ctl / *_obs arrays the model already carries in obj$env$data, so
  # the regenerated obsvec is what the fit would have carried had it been
  # pre-built. `osa_dat` lets the caller pass a copy to avoid recomputing it.
  data2 <- if (is.null(osa_dat)) build_osa_data(obj$env$data, build_osa = TRUE) else osa_dat
  data2$obs_ctl <- NULL   # R-side metadata table, not a TMB DATA input
  # obj$env$data is already sanitized (stored as double with a 'check.passed'
  # attribute). Overwrite osa_mode as a double (DATA_INTEGER reads it via
  # CppAD::Integer) and drop 'check.passed' so MakeADFun re-sanitizes cleanly.
  data2$osa_mode <- as.numeric(osa_mode)
  attr(data2, "check.passed") <- NULL
  random_names <- if (length(obj$env$random)) {
    unique(names(obj$env$par)[obj$env$random])
  } else {
    NULL
  }
  # parList(last.par.best) returns the best (MLE) parameters as a list structured
  # to match the map -- the reliable way to rebuild at the fitted values. Use
  # last.par.best (not the parList() default last.par) so OSA residuals are
  # computed at the optimum even when getsd = FALSE (no sdreport re-eval).
  obj2 <- TMB::MakeADFun(
    data       = data2,
    parameters = obj$env$parList(par = obj$env$last.par.best),
    map        = obj$env$map,
    random     = random_names,
    DLL        = fit$TMBfilename,
    silent     = TRUE)
  obj2$fn(obj2$par)   # evaluate once so last.par is populated
  obj2
}


#' Statistical diagnostics for OSA residuals
#'
#' @description
#' Computes the Stewart and Monnahan (2025) statistical diagnostics for a set of
#' OSA residuals: the standard deviation of the normalized residuals (SDNR) and
#' the lower/upper tail statistics, each with the 95% interval expected under
#' the standard normal null hypothesis (so departures can be judged objectively
#' rather than by eye). Computed per data source (type x fleet) and overall.
#'
#' Under a correctly specified model OSA residuals are already iid standard
#' normal, so the SDNR is simply their sample standard deviation. Its null
#' interval follows the chi-square result for the sample standard deviation of
#' `n` standard normals (Francis 2014). Each tail statistic is the `r`-th
#' order statistic of the residuals, at index `round(p * (n + 1))`, and its
#' null interval is exact: the `r`-th order statistic of `n` uniforms is
#' `Beta(r, n - r + 1)`. Nothing is simulated.
#'
#' @param osa An `rceattle_osa` object from [osa_residuals()], or a data frame
#'   with `residual` and (optionally) `type`/`fleet` columns.
#' @param nsim Ignored from 5.29.0, when the tail null intervals became exact;
#'   retained so existing calls keep working, and warns when supplied.
#' @param probs Lower/upper tail probabilities. Default `c(0.025, 0.975)`.
#' @param seed Ignored from 5.29.0, as `nsim` is; the residuals themselves are
#'   randomized-quantile, so seed [osa_residuals()] instead.
#'
#' @return A data frame (class `"rceattle_osa_diagnostics"`, so it prints as a
#'   compact severity-tagged summary; every column is still there and `$` works
#'   as before) with one row per data source plus an `"all"` row, with columns:
#'   `group` (the `"<source> fleet <n>"` label), `source`, `fleet`, `n`, `sdnr`,
#'   `sdnr_lo`, `sdnr_hi`, `lower`, `lower_lo`, `lower_hi`, `upper`, `upper_lo`,
#'   `upper_hi`, the order statistic each tail was read at (`lower_r`, `upper_r`)
#'   and its exact nominal probability `r/(n + 1)` (`lower_p`, `upper_p`,
#'   which differ from `probs` at small `n`), and the logical flags `sdnr_ok`,
#'   `lower_ok`, `upper_ok` (TRUE when the statistic is inside its null
#'   interval). On the `"all"` row `source` and `fleet` are `NA`.
#'
#' @references
#' Francis, R.I.C.C. 2014. Replacing the multinomial in stock assessment models:
#'   a first step. Fish. Res. 151:70-84.
#'
#' Stewart, I.J., and Monnahan, C.C. 2025. Can. J. Fish. Aquat. Sci. 82:1-13.
#'
#' @seealso [osa_residuals()]
#' @export
osa_diagnostics <- function(osa, nsim = 10000, probs = c(0.025, 0.975),
                            seed = 123) {

  if (is.null(osa$residual)) {
    stop("'osa' must have a 'residual' column (e.g. output of osa_residuals()).")
  }
  # The tail null intervals are exact from v5.29.0, so nothing here simulates.
  # Both arguments are kept so existing calls keep working.
  if (!missing(nsim) || !missing(seed)) {
    warning("'nsim' and 'seed' are ignored: the tail null intervals are now ",
            "exact (the r-th order statistic is Beta(r, n - r + 1)), so ",
            "osa_diagnostics() no longer simulates. Note the OSA residuals ",
            "themselves are still randomized-quantile residuals -- pass a seed ",
            "to osa_residuals() to control those.", call. = FALSE)
  }

  has_groups <- !is.null(osa$source) && !is.null(osa$fleet)
  groups <- if (has_groups) {
    split(osa, list(osa$source, osa$fleet), drop = TRUE)
  } else {
    list(all = osa)
  }

  rows <- lapply(names(groups), function(g) {
    grp <- groups[[g]]
    stat <- .osa_sdnr_tails(grp$residual, probs = probs)
    data.frame(
      group  = if (has_groups) paste(grp$source[1], "fleet", grp$fleet[1]) else "all",
      source = if (has_groups) grp$source[1] else NA_character_,
      fleet  = if (has_groups) grp$fleet[1] else NA_integer_,
      stat,
      stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)

  # Overall row across all residuals.
  overall <- .osa_sdnr_tails(osa$residual, probs = probs)
  out <- rbind(out, data.frame(group = "all", source = NA_character_,
                               fleet = NA_integer_, overall,
                               stringsAsFactors = FALSE))
  rownames(out) <- NULL
  # Still a data frame -- every column and `$` access is unchanged. The class
  # only adds a print method, so the twenty columns stop wrapping across three
  # screen-widths with no verdict; see print.rceattle_osa_diagnostics().
  class(out) <- c("rceattle_osa_diagnostics", "data.frame")
  out
}


#' @export
print.rceattle_osa_diagnostics <- function(x, ...) {
  df <- as.data.frame(x)
  per <- df[df$group != "all" | !is.na(df$source), , drop = FALSE]
  all_row <- df[df$group == "all" & is.na(df$source), , drop = FALSE]

  # SDNR is the headline statistic, so it carries WARN; a tail outside its null
  # interval with an acceptable SDNR is a NOTE. Both intervals are exact under
  # the standard-normal null -- chi-square for SDNR, Beta for the tail order
  # statistics -- so "outside" already means "further than chance", and neither
  # is a FAIL: these are diagnostics on fit, not a broken model.
  sev <- rep("OK", nrow(per))
  sev[!is.na(per$lower_ok) & !per$lower_ok] <- "NOTE"
  sev[!is.na(per$upper_ok) & !per$upper_ok] <- "NOTE"
  sev[!is.na(per$sdnr_ok)  & !per$sdnr_ok]  <- "WARN"

  n_bad <- sum(sev != "OK")
  .rce_diag_header(
    "OSA diagnostics", .rce_worst(sev),
    paste0(.rce_n_of(n_bad, nrow(per)),
           " source(s) outside a null interval",
           if (nrow(all_row)) paste0("; overall SDNR ",
                                     formatC(all_row$sdnr[1], format = "f", digits = 2),
                                     " (", formatC(all_row$sdnr_lo[1], format = "f", digits = 2),
                                     "-", formatC(all_row$sdnr_hi[1], format = "f", digits = 2),
                                     ")") else ""))

  # Worst first, so a long fleet list opens with what needs attention.
  ord <- order(match(sev, .CONV_SEVERITY), decreasing = TRUE)
  per <- per[ord, , drop = FALSE]; sev <- sev[ord]
  per$.tag <- .rce_sev_tag(sev)
  per$.sdnr_int <- sprintf("%.2f-%.2f", per$sdnr_lo, per$sdnr_hi)

  .rce_diag_table(per, c(" " = ".tag", "source" = "source", "fleet" = "fleet",
                         "n" = "n", "sdnr" = "sdnr", "null" = ".sdnr_int"))
  cat("  tails and their intervals are in the data frame",
      "(lower/upper, *_lo, *_hi)\n")
  invisible(x)
}


#' Exact null interval for an order statistic of standard normals
#'
#' The `r`-th order statistic of `n` uniforms is `Beta(r, n - r + 1)`. What is
#' compared against this must be `sort(resid)[r]`: pairing it with `quantile()`'s
#' type-7 interpolation drops coverage to about 0.86 near n = 50. `r` is clamped
#' to `[1, n]`, unclamped, `qbeta(p, n + 1, 0)` is 1 and the upper tail check
#' passes for every series with n <= 19.
#'
#' @param q Nominal tail probability.
#' @param n Number of residuals.
#' @param probs Lower/upper probabilities for the interval itself.
#' @return A list: the order statistic index `r`, its exact nominal probability
#'   `r/(n+1)`, and the two-element `null` interval.
#' @noRd
.osa_tail_null <- function(q, n, probs = c(0.025, 0.975)) {
  r <- as.integer(min(n, max(1L, round(q * (n + 1)))))
  list(r = r, nominal = r / (n + 1),
       null = stats::qnorm(stats::qbeta(probs, r, n - r + 1)))
}


#' SDNR and tail statistics with standard-normal null intervals
#'
#' @param resid Numeric vector of residuals (assumed standard normal under H0).
#' @param probs See [osa_diagnostics()].
#' @return A one-row data frame of statistics and their null intervals.
#' @noRd
.osa_sdnr_tails <- function(resid, probs = c(0.025, 0.975)) {
  resid <- resid[is.finite(resid)]
  n <- length(resid)
  if (n < 2) {
    return(data.frame(n = n, sdnr = NA_real_, sdnr_lo = NA_real_,
                      sdnr_hi = NA_real_, lower = NA_real_, lower_lo = NA_real_,
                      lower_hi = NA_real_, upper = NA_real_, upper_lo = NA_real_,
                      upper_hi = NA_real_, lower_r = NA_integer_,
                      upper_r = NA_integer_, lower_p = NA_real_,
                      upper_p = NA_real_, sdnr_ok = NA, lower_ok = NA,
                      upper_ok = NA))
  }

  # SDNR and its chi-square null interval (Francis 2014).
  sdnr    <- stats::sd(resid)
  df      <- n - 1
  sdnr_lo <- sqrt(stats::qchisq(probs[1], df) / df)
  sdnr_hi <- sqrt(stats::qchisq(probs[2], df) / df)

  # Tail statistics as order statistics, against their exact nulls, so the
  # observed value and its interval come from the same estimator.
  lo <- .osa_tail_null(probs[1], n, probs)
  hi <- .osa_tail_null(probs[2], n, probs)
  s     <- sort(resid)
  lower <- s[lo$r]
  upper <- s[hi$r]

  data.frame(
    n        = n,
    sdnr     = sdnr,    sdnr_lo  = sdnr_lo,    sdnr_hi  = sdnr_hi,
    lower    = lower,   lower_lo = lo$null[1], lower_hi = lo$null[2],
    upper    = upper,   upper_lo = hi$null[1], upper_hi = hi$null[2],
    lower_r  = lo$r,    upper_r  = hi$r,
    lower_p  = lo$nominal, upper_p = hi$nominal,
    sdnr_ok  = sdnr  >= sdnr_lo     & sdnr  <= sdnr_hi,
    lower_ok = lower >= lo$null[1] & lower <= lo$null[2],
    upper_ok = upper >= hi$null[1] & upper <= hi$null[2])
}


#' @export
print.rceattle_osa <- function(x, ...) {
  cat("Rceattle one-step-ahead (OSA) residuals\n")
  # `method` is a named vector when a family overrode the caller's choice.
  m <- attr(x, "method")
  m_txt <- if (length(m) > 1L && !is.null(names(m))) {
    paste(paste0(names(m), " = ", m), collapse = ", ")
  } else as.character(m)
  # Only when it is on: a randomized composition residual carries a draw, so the
  # seed is part of the answer rather than a formality. The flag is a named
  # vector when only some fleets were randomized.
  d <- attr(x, "discrete")
  cat("  method:", m_txt, " seed:", attr(x, "seed"),
      if (isTRUE(d)) "  randomized (discrete compositions)"
      else if (length(d) > 1L) "  randomized (discrete compositions, except the D-M fleets)"
      else "",
      "\n")
  cat("  ", nrow(x), " residuals across ",
      length(unique(paste(x$source, x$fleet))), " data source(s)\n", sep = "")
  print(utils::head(as.data.frame(x), ...))
  invisible(x)
}


#' @export
summary.rceattle_osa <- function(object, ...) {
  osa_diagnostics(object, ...)
}
