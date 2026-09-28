#' Process residuals for an Rceattle model's random-effect processes
#'
#' @description
#' One-sample process residuals for the model's random-effect (or penalized)
#' process deviations, in the style of SAM's `procres()` (Nielsen and Berg
#' 2014). They validate the *process* model, whether the deviations behave like
#' their assumed iid normal process, a complement to the observation-based
#' [osa_residuals()].
#'
#' Each set of deviations holds a Gaussian process prior in the model. The
#' posterior *mode* of the deviations is shrunk toward that prior, so, following
#' SAM, a single draw is taken from the joint posterior of the deviations (from
#' the joint precision when they are random effects, or the fixed-effect
#' covariance when they are penalized fixed effects) and standardized by the
#' process standard deviation. Under a correctly specified process these are
#' approximately iid N(0, 1). Because a random draw is used the residuals are
#' stochastic; set `seed` for reproducibility.
#'
#' Supported processes and their deviations:
#' \describe{
#'   \item{`"recruitment"`}{`rec_dev`, prior `N(-bias_adjust_proc * R_sd^2/2, R_sd^2)` per species.}
#'   \item{`"initial"`}{`init_dev`, prior `N(-bias_adjust_proc * R_sd^2/2, R_sd^2)` per species.}
#'   \item{`"catchability"`}{`index_q_dev`, prior `N(0, q_dev_sd^2)` per index.}
#' }
#' `"all"` returns every supported process present in the fit. Selectivity and
#' natural-mortality deviations (which use random-walk / 2D-AR1 priors) are not
#' yet supported. Recruitment and initial-abundance residuals are exact (their
#' priors are iid). Catchability residuals are exact only for the iid deviate
#' prior (`Time_varying_q = 1` or `2`); for a random-walk or AR1 catchability
#' prior the marginal-SD standardization ignores the prior correlation, so those
#' residuals are approximate and a warning is emitted.
#'
#' @param object A fitted `Rceattle` model. The targeted deviations must be estimated, as random effects (e.g. `random_rec = TRUE`) or as penalized fixed effects, with a usable covariance.
#' @param fit deprecated name for `object`, still accepted so existing
#'   scripts keep working. Supplying both is an error.
#' @param process One of `"recruitment"`, `"initial"`, `"catchability"`, or
#'   `"all"`.
#' @param seed Seed for the posterior draw. Default 123.
#'
#' @return A data frame of class `rceattle_osa` (so it can be passed to
#'   [osa_diagnostics()] and [plot.rceattle_osa()]) with one row per deviation:
#'   columns `source` (the process), `fleet`, `species`, `year`,
#'   `age_length_bin`, and `residual`.
#'
#' @references
#' Nielsen, A., and Berg, C.W. 2014. Estimation of time-varying selectivity in
#'   stock assessments using state-space models. Fish. Res. 158:96-101.
#'
#' @seealso [osa_residuals()], [osa_diagnostics()]
#' @export
process_residuals <- function(object = NULL,
                              process = c("recruitment", "initial",
                                          "catchability", "all"),
                              seed = 123, fit = NULL) {
  # `fit` was the old name for `object`; see R/0-deprecate.R.
  if (!missing(fit))
    object <- .rce_deprecated_arg(fit, !missing(object), "fit", "object", "process_residuals")

  if (!inherits(object, "Rceattle")) {
    stop("`object` must be a fitted Rceattle model (from fit_mod()).")
  }
  if (is.null(object$obj)) stop("`object` has no TMB object ($obj).")
  process <- match.arg(process)

  specs <- list(recruitment  = "rec_dev",
                initial      = "init_dev",
                catchability = "index_q_dev")
  todo <- if (process == "all") names(specs) else process

  # Seed once here -- saving and restoring the caller's RNG -- and let the
  # processes draw from a single advancing stream. Seeding inside the per-process
  # loop instead would make equal-length deviation blocks (e.g. recruitment vs
  # initial) draw identical standard normals (spurious cross-process correlation)
  # and would clobber the caller's global .Random.seed.
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)

  out <- lapply(todo, function(p) {
    res <- tryCatch(.process_residual_one(object, p, specs[[p]]),
                    error = function(e) {
                      if (process == "all") NULL else stop(e)
                    })
    res
  })
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  if (is.null(out) || nrow(out) == 0) {
    stop("No supported process deviations were found in this fit.")
  }
  rownames(out) <- NULL
  class(out) <- c("rceattle_osa", "data.frame")
  attr(out, "method") <- "procres"
  attr(out, "seed")   <- seed
  out
}


#' Process residuals for a single deviation parameter
#' @keywords internal
.process_residual_one <- function(fit, process, par_name) {
  obj <- fit$obj

  # ---- Posterior mode + covariance of the estimated deviations ----
  # Random effects: use the joint precision (the deviation block of its inverse
  # is the marginal posterior covariance). Penalized fixed effects: use the
  # fixed-effect covariance (inverse Hessian) directly.
  jp <- fit$sdrep$jointPrecision
  if (!is.null(jp) && par_name %in% rownames(jp)) {
    idx   <- which(rownames(jp) == par_name)
    lp    <- obj$env$last.par.best
    mode  <- lp[names(lp) == par_name]
    Sigma <- as.matrix(Matrix::solve(jp))[idx, idx, drop = FALSE]
  } else {
    sdr <- fit$sdrep
    if (is.null(sdr) || is.null(sdr$cov.fixed)) sdr <- TMB::sdreport(obj)
    cf <- sdr$cov.fixed
    if (is.null(cf) || !(par_name %in% colnames(cf))) {
      stop("'", par_name, "' is not estimated in this fit (no usable ",
           "covariance). Refit with the matching option (e.g. random_rec for ",
           "recruitment) so the deviations are estimated.")
    }
    cols  <- which(colnames(cf) == par_name)
    mode  <- sdr$par.fixed[names(sdr$par.fixed) == par_name]
    Sigma <- cf[cols, cols, drop = FALSE]
  }
  n <- length(mode)
  if (n == 0) stop("No estimated '", par_name, "' deviations.")

  # ---- Prior mean / SD and labels for the estimated elements ----
  spec <- .process_prior_spec(fit, process, par_name, n)
  if (length(spec$sd) != n) {
    stop("Could not align the prior for '", par_name, "' (", length(spec$sd),
         " priors vs ", n, " estimated deviations).")
  }

  # ---- One posterior draw, standardized by the prior ----
  # The RNG is seeded once by the caller (process_residuals()); drawing from that
  # advancing stream keeps each process's normals distinct and reproducible.
  L <- t(chol(Sigma))
  draw  <- as.numeric(mode + L %*% stats::rnorm(n))
  resid <- (draw - spec$mean) / spec$sd

  data.frame(
    source         = process,
    fleet          = spec$fleet,
    species        = spec$species,
    sex            = NA_integer_,
    year           = spec$year,
    age_length_bin = spec$age,
    length         = NA_integer_,
    index_label    = NA_character_,
    observed      = NA_real_,
    predicted     = NA_real_,
    sd            = NA_real_,
    residual      = resid,
    stringsAsFactors = FALSE)
}


#' Prior mean/SD and labels for the estimated elements of a deviation parameter
#'
#' Maps the estimated (non-fixed) elements of `par_name`, in the column-major
#' order used by TMB and the covariance matrices, to their species/fleet,
#' year/age labels, and the process prior mean and SD.
#' @keywords internal
.process_prior_spec <- function(fit, process, par_name, n) {
  d    <- fit$data_list
  obj  <- fit$obj
  full <- obj$env$parList()[[par_name]]          # full parameter array
  mapf <- fit$map$mapFactor[[par_name]]
  # Estimated (non-fixed) elements, in column-major order -- the same order as
  # the covariance matrix rows. If this does not match the covariance dimension
  # the element-to-(species/year/age) mapping would be guesswork, so error
  # rather than mislabel.
  est <- if (is.null(mapf)) seq_along(full) else which(!is.na(as.integer(mapf)))
  if (length(est) != n) {
    stop("Could not map the ", n, " estimated '", par_name, "' deviations to ",
         "model elements (found ", length(est), " estimated positions). ",
         "Process residuals are not supported for this parameter's ",
         "configuration.")
  }
  ai <- arrayInd(est, dim(full))                 # one row per estimated element
  styr <- d$styr

  if (process %in% c("recruitment", "initial")) {
    R_sd <- exp(as.numeric(obj$env$parList()$R_log_sd))   # per species
    # rec_dev / init_dev prior mean is the lognormal bias correction
    # -bias_adjust_proc * sigma^2/2 (matches ceattle.cpp slots 9 & 10).
    ba   <- obj$env$data$bias_adjust_proc
    if (is.null(ba)) ba <- 1
    sp   <- ai[, 1]
    list(species = sp, fleet = sp,
         year = if (process == "recruitment") styr + ai[, 2] - 1L else NA_integer_,
         age  = if (process == "initial") ai[, 2] + 1L else NA_integer_,
         mean = -as.numeric(ba) * R_sd[sp]^2 / 2, sd = R_sd[sp])
  } else if (process == "catchability") {
    q_sd <- exp(as.numeric(obj$env$parList()$index_q_dev_log_sd))  # per fleet
    flt  <- ai[, 1]
    # The (mean 0, marginal SD) standardization below is exact only for the iid
    # catchability deviate prior (Time_varying_q = 1 or 2). Warn when an involved
    # fleet uses a random walk (index_varying_q == 4), because ignoring the prior
    # correlation makes those residuals approximate.
    ivq  <- obj$env$data$index_varying_q
    uflt <- unique(flt)
    corr <- rep(FALSE, length(uflt))
    if (!is.null(ivq)) corr <- corr | (ivq[uflt] %in% 4L)
    if (any(corr, na.rm = TRUE)) {
      warning("process = 'catchability': index fleet(s) ",
              paste(uflt[which(corr)], collapse = ", "),
              " use a random-walk catchability deviate prior; ",
              "their residuals standardize by the marginal SD and are therefore ",
              "approximate. They are exact only for the iid prior ",
              "(Time_varying_q = 1 or 2).", call. = FALSE)
    }
    list(species = NA_integer_, fleet = flt,
         year = styr + ai[, 2] - 1L, age = NA_integer_,
         mean = rep(0, n), sd = q_sd[flt])
  } else {
    stop("Process '", process, "' is not supported.")
  }
}
