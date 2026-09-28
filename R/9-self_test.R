#' Self-test simulation analysis
#'
#' @description Simulates data from an Rceattle model and refits the model to the simulated data, to check that the fitting procedure recovers the operating-model parameters. TODO add process variation (i.e. random recruitment deviations) to the simulation.
#'
#' @inheritParams rceattle-refit-args
#' @param seed random number seed. Each simulation \code{i} uses \code{seed + i}
#'   so results are reproducible under both sequential and parallel execution.
#' @param nsim number of simulations
#' @param simulate passed to \code{\link{sim_mod}}. If \code{TRUE} (default),
#'   data are simulated with observation error; if \code{FALSE}, expected
#'   values from the model are used.
#' @param getsd whether each refit runs \code{TMB::sdreport}. Self-test compares
#'   the refit point estimates to the operating model, so \code{FALSE} is faster
#'   with no effect on that comparison. Default \code{NULL} inherits the input
#'   model's setting (\code{TRUE} only if it holds an \code{sdrep}).
#' @param phase as in \code{\link{fit_mod}}. Under the default
#'   \code{start = "initial"} each refit covers the same ground the original fit
#'   did, so a model that needed phasing to fit its real data needs it again for
#'   every simulated one, without it such a model's refits can end many orders
#'   of magnitude from a zero gradient and be dropped as non-converged. Default
#'   \code{NULL} reads the setting \code{fit_mod()} recorded on the source fit
#'   (\code{fit$run_config$fit_control$phase}), so a model fitted under the
#'   package default of \code{phase = FALSE} is refitted unphased; pass
#'   \code{TRUE} for a model that needs phasing but was not fitted with it.
#' @param start which of the input model's parameter sets each refit starts
#'   from. \code{"initial"} (default) uses \code{initial_params}, the values the
#'   original fit itself started from, so the estimator has to travel the same
#'   distance to the optimum on simulated data that it did on the real data.
#'   \code{"estimated"} starts from \code{estimated_params} instead: much faster
#'   and far more likely to converge, but the fixed effects, and, with
#'   \code{random_rec = TRUE}, the inner Laplace problem too, begin at the
#'   generating values, so on a multimodal or weakly identified surface the
#'   optimizer never leaves the basin containing them and recovery is close to
#'   guaranteed by construction. Read it as optimistic about recovery, not
#'   merely less powerful. (Nor is it a complete warm start: \code{fit_mod()}
#'   resets \code{log_Ftarget}, \code{proj_F_prop}, and the stock-recruit
#'   \eqn{\alpha}/\eqn{\beta} from the model's own specification under either
#'   setting.) Non-identifiability shows up in the curvature and so is visible
#'   either way, via \code{$convergence}'s Hessian conditioning and
#'   estimability checks, it is \emph{reachability} that a warm start stops
#'   testing.
#' @param debug return every simulation rather than the converged ones. The
#'   dropped runs are the interesting ones when a self-test comes back short, and
#'   each holds its own \code{$convergence} diagnostics. See \strong{Value}.
#' @param timeout elapsed-second limit per simulation, \code{Inf} (default) for
#'   none. The optimizer runs with no iteration cap, so a replicate that wanders
#'   somewhere pathological can stall the whole run, a hang that no convergence
#'   check can catch, because the fit never returns. One that exceeds the limit
#'   is stopped, counted as non-converged and reported separately. Approximate:
#'   the limit is checked when control returns to R, so it fires between the
#'   optimizer's function evaluations rather than inside one.
#' @param process passed to \code{\link{sim_mod}}. \code{FALSE} (default) keeps
#'   the fitted process deviations, so the test measures whether the estimator
#'   recovers its own parameters from new observations. Naming a process,
#'   \code{"recruitment"}, \code{"M"}, \code{"growth"}, \code{"dynamics"},
#'   \code{TRUE}, ..., redraws it too, so the test instead measures whether the
#'   estimator recovers a process it has not been shown. The deviations behind
#'   each replicate come back in \code{attr(result, "process_sim")}; see
#'   \code{Value}.
#'
#' @return A list of Rceattle models named \code{Sim_1}, \code{Sim_2}, ....
#'   By default only the converged simulations, renumbered contiguously; a
#'   message reports how many were dropped.
#'
#'   The list holds class \code{"Rceattle_selftest"} and the number of
#'   simulations attempted in \code{attr(, "nsim")}, so printing it reports the
#'   convergence rate; see \code{\link{print.Rceattle_selftest}}. It is
#'   otherwise the list it always was, \code{sims[["Sim_1"]]},
#'   \code{length(sims)}, \code{lapply()} and
#'   \code{plot_biomass(c(sims, list(fit)))} are all unchanged, and \code{c()}
#'   and \code{[} return a plain list of fits.
#'
#'   With \code{debug = TRUE}, every simulation, with \code{Sim_i} being
#'   simulation \code{i} (so it pairs with the seed \code{seed + i}), and a
#'   logical vector of the convergence verdicts in \code{attr(, "converged")}.
#'   Inspect a failure with \code{sims[[j]]$convergence}. A simulation whose
#'   refit errored outright is returned as the condition object rather than a
#'   model, so it cannot abort the run.
#'
#'   When \code{process} redrew something, \code{attr(, "process_sim")} holds the
#'   deviations that generated each replicate's data, a list keyed by the same
#'   \code{Sim_i} names, so \code{attr(x, "process_sim")[["Sim_1"]]} belongs to
#'   \code{x[["Sim_1"]]}, subset and renumbered alongside the models. Each entry
#'   is a named list of whichever of \code{rec_dev}, \code{init_dev},
#'   \code{log_M1_dev} and \code{beta_linkage_re} were drawn, each with a
#'   same-shaped \code{_drawn} logical marking the cells the draw touched,
#'   restrict any recovery statistic to those, since the rest are fitted values
#'   (see \code{\link{sim_mod}}). Compare estimates against these, not against
#'   the operating model: its fitted deviations are no longer what generated the
#'   data.
#'
#' @section Interpreting the spread:
#' \code{\link{sim_mod}} redraws the observations only, indices, catch,
#' compositions, CAAL and stomach contents. Some rows are deliberately left
#' alone, and \code{\link{sim_mod}} warns about each: a predator whose
#' suitability is empirical rather than estimated has no predicted diet to draw
#' from, and a covariance (MVN/MVNORM) survey fleet has no covariance outside
#' the years it is fitted to. Those data are held fixed across every replicate,
#' so recovery of whatever they inform is optimistic.
#'
#' By default it does not redraw recruitment, so with
#' \code{random_rec = TRUE} every replicate shares the operating model's single
#' recruitment realization, and that realization is its shrunk empirical-Bayes
#' modes rather than a draw from N(0, sigmaR). Two consequences: the spread#' across replicates reflects observation error only and is a lower bound on
#' estimation uncertainty in SSB and recruitment (do not read it against the
#' model's own uncertainty bands, which include process error); and sigmaR is
#' re-estimated from deviations that were shrunk toward zero the same way in
#' every replicate, a downward bias that averaging over simulations does not
#' remove. Pass \code{process = "recruitment"} (or \code{"dynamics"}, or
#' \code{TRUE}) to redraw it and remove both, at the cost of asking a different
#' question, see \code{process} above.
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(data_list = BS2017SS,
#'     inits = NULL, file = NULL,
#'     estimateMode = 0, random_rec = FALSE,
#'     msmMode = 0, avgnMode = 0,
#'     phase = FALSE, verbose = 0)
#' sims <- self_test(ss_run, nsim = 10)
#' }
#' @export
self_test <- function(object = NULL, nsim = 50, simulate = TRUE, seed = 123, cores = NULL, getsd = NULL, phase = NULL, start = c("initial", "estimated"), debug = FALSE, timeout = Inf, process = FALSE, fit_control = NULL, Rceattle = NULL) {
  # `Rceattle` was the old name for `object`; see R/0-deprecate.R.
  if (!missing(Rceattle))
    object <- .rce_deprecated_arg(Rceattle, !missing(object), "Rceattle", "object", "self_test")
  # Cleared once consumed: .parallel_lapply() exports every binding in this frame
  # to the PSOCK workers, and two names for one fitted model send it twice.
  Rceattle <- NULL

  # `process` stays after the arguments that predate it so positional calls keep
  # their meaning. The deprecated `Rceattle` formal sits last; see R/0-deprecate.R.
  if (!inherits(object, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }

  start <- match.arg(start)
  if (is.null(object[[paste0(start, "_params")]])) {
    stop("`start = \"", start, "\"` needs the model's ", start,
         "_params, which this fit does not carry.", call. = FALSE)
  }

  # rho/self-test read point estimates, so getsd = FALSE is faster and neutral;
  # default inherits the input model's setting (matches retrospective/jitter).
  if (is.null(getsd)) getsd <- !is.null(object$sdrep)

  # Phasing likewise inherits the input model's setting. Under the default
  # `start = "initial"` a refit covers the same ground the original fit did, so
  # if that fit needed phasing, so does every refit.
  #
  # Read the value fit_mod() recorded, so a custom phase list carries over as
  # the list rather than collapsing to TRUE and being rebuilt from set_phases()
  # defaults -- that would phase the refit on a different schedule than the fit
  # it is testing. Fall back to whether `phase_params` was attached (fit_mod()
  # does that only when it phased) for models predating `run_config`.
  if (is.null(phase)) {
    phase <- object$run_config$fit_control$phase
    if (is.null(phase)) phase <- !is.null(object$phase_params)
  }

  # fit_control() overrides `phase` and `getsd`, but only where the caller
  # changed them from fit_control()'s own defaults; see .rce_refit_control().
  ctl <- .rce_refit_control(fit_control, "self_test")
  if (!is.null(ctl$phase)) phase <- ctl$phase
  if (!is.null(ctl$getsd)) getsd <- ctl$getsd

  # Cross-platform parallel via parallel::parLapply on a PSOCK cluster
  # (same approach as run_mse). Respect the CRAN core limit
  # ('_R_CHECK_LIMIT_CORES_' is set during R CMD check;
  # parallel::makeCluster errors if we exceed 2 cores then).
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  if (is.null(cores)) {
    cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  } else {
    cores <- max(1L, as.integer(cores))
    if (cran_cap) cores <- min(cores, 2L)
  }
  use_parallel <- nsim > 1L && cores > 1L

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Per-simulation closure ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  run_one_sim <- function(i) {

    set.seed(seed + i) # unique seed per sim for reproducibility under parallel

    # * Simulate new data
    sim_data <- Rceattle::sim_mod(object, simulate = simulate, process = process)
    data_list <- sim_data
    # The deviations that generated this replicate, when process error was
    # redrawn. Carried through so the caller compares estimates against what
    # actually generated the data; the source model's fitted deviations are no
    # longer the truth once a process has been redrawn.
    truth <- attr(sim_data, "process_sim")

    # * Adjust initial values ----
    inits <- switch(start,
                    initial   = object$initial_params,
                    estimated = object$estimated_params)


    # * Refit ----
    # Bounded and error-trapped: a replicate that errors or hangs would otherwise
    # abort self_test() -- and, under a cluster, every other replicate with it --
    # which is the opposite of what `debug` is for.
    newmod <- .refit_with_timeout(
      suppressMessages(
        suppressWarnings(
          # Refit the simulated data set from `start`, under the source model's
          # phasing (see `phase` above).
          .refit_like(
            data_list        = data_list,
            inits            = inits,
            estimateMode     = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode),
            phase            = phase,
            getsd            = getsd,
            srr_mse_switchyr = min(data_list$srr_mse_switchyr, data_list$endyr),
            suit_endyr       = pmin(data_list$suit_endyr, data_list$endyr))
        )
      ),
      timeout = timeout)

    # Report the verdict beside the model rather than dropping it here: the
    # dispatcher filters below, so `debug = TRUE` can hand back the runs that
    # failed with their $convergence diagnostics intact.
    if (inherits(newmod, "condition")) {
      return(list(model = newmod, converged = FALSE, error = conditionMessage(newmod),
                  truth = truth))
    }
    list(model = newmod, converged = .refit_converged(newmod), truth = truth)
  } # End run_one_sim closure


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Dispatch sims (parallel via PSOCK or sequential) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if (use_parallel) {
    mod_list <- .parallel_lapply(1:nsim, run_one_sim, min(cores, nsim), environment())
  } else {
    mod_list <- lapply(1:nsim, run_one_sim)
  }

  # Split the verdict from the models ----
  # Sim_i is simulation i here, before any filtering, so a debug caller can line
  # a failed run up with the seed (`seed + i`) that produced it.
  converged <- vapply(mod_list, function(x) isTRUE(x$converged), logical(1))
  errs      <- vapply(mod_list, function(x) x$error %||% NA_character_, character(1))
  truths    <- lapply(mod_list, function(x) x$truth)
  mod_list  <- lapply(mod_list, function(x) x$model)
  # Unconditional, so `names(which(!attr(sims, "converged")))` works at every
  # length -- including the empty case, where a guarded assignment would have
  # left the attribute unnamed.
  names(mod_list) <- names(converged) <- paste0("Sim_", seq_along(mod_list))
  names(truths) <- names(mod_list)
  .report_dropped(sum(!converged), length(converged), "simulation")
  .report_errors(errs, "simulation")


  # Return ----
  # debug: every simulation, `converged` naming which are which, so a run that
  # failed can be read through its own $convergence rather than inferred from a
  # short list. Otherwise the converged runs only, renumbered Sim_1..Sim_k --
  # a self-test is read as a distribution over runs, so the gaps carry no
  # meaning, and plot_biomass(model_names = names(sims)) stays contiguous.
  # `process_sim` is subset and renumbered alongside the models so that
  # attr(sims, "process_sim")[[k]] is always the truth behind sims[[k]].
  if (isTRUE(debug)) {
    attr(mod_list, "converged") <- converged
    if (any(!vapply(truths, is.null, logical(1)))) {
      attr(mod_list, "process_sim") <- truths
    }
    return(.as_selftest(mod_list, nsim))
  }

  keep_truth <- truths[converged]
  mod_list <- mod_list[converged]
  if (length(mod_list) > 0) {
    names(mod_list) <- names(keep_truth) <- paste0("Sim_", seq_along(mod_list))
  }
  if (any(!vapply(keep_truth, is.null, logical(1)))) {
    attr(mod_list, "process_sim") <- keep_truth
  }
  return(.as_selftest(mod_list, nsim))
}


# Still the same list of models: `sims[["Sim_1"]]`, `length(sims)`,
# `plot_biomass(c(sims, list(fit)))` and compare_sim() all read it exactly as
# before -- c() and [ drop the class, which is right, because what they return
# is a plain list of fits. Only class() gains an entry and print() gets an
# opinion.
#
# `nsim` is carried because the returned runs cannot say how many were
# attempted: the non-converged ones are dropped, and the fraction that reached
# an optimum is the first thing a self-test has to report. Same reason
# jitter() carries `njitter`.
.as_selftest <- function(mod_list, nsim) {
  attr(mod_list, "nsim") <- nsim
  class(mod_list) <- "Rceattle_selftest"
  mod_list
}


#' Print method for a self-test
#'
#' @description Reports what a self-test is run to find out: how many
#' simulations the estimator brought back to an optimum, and under what the
#' replicates were generated. The returned list cannot say the first on its own
#', non-converged runs are dropped before it is returned, so the number of
#' fits returned is not the number attempted.
#'
#' @details
#' The status is the convergence RATE, since that is the quantity the run
#' produces: `FAIL` if nothing converged, `WARN` below `rate`, `NOTE` if any
#' replicate was dropped, `OK` otherwise. The table beneath tallies each
#' returned fit's own `$convergence$status`, which is a separate question,
#' a replicate can reach a zero gradient and still hold a `NOTE` or `WARN`,
#' and those are not folded into the header for that reason.
#'
#' The last line says what generated the replicates, because it decides what the
#' spread means. With processes held fixed (the default) `sim_mod()` redraws the
#' observations alone, so the spread reflects observation error only and is a
#' lower bound on estimation uncertainty, do not read it against the model's
#' own uncertainty bands, which include process error. See
#' **Interpreting the spread** in [self_test()].
#'
#' @param x A `"Rceattle_selftest"` object from [self_test()].
#' @param rate Fraction of simulations that must converge before the run counts
#'   as clean. Default `0.9`, matching [jitter()]'s.
#' @param ... Currently unused.
#' @return `x`, invisibly.
#' @export
print.Rceattle_selftest <- function(x, rate = 0.9, ...) {
  tried <- attr(x, "nsim") %||% NA_integer_
  # debug = TRUE returns every simulation, so the count of entries is not the
  # count that converged; the attribute it carries is.
  conv <- attr(x, "converged")
  kept <- if (is.null(conv)) length(x) else sum(conv)

  sev <- if (kept == 0L) "FAIL"
         else if (!is.na(tried) && kept < rate * tried) "WARN"
         else if (!is.na(tried) && kept < tried) "NOTE"
         else "OK"

  .rce_diag_header(
    "self-test", sev,
    paste0(kept, " of ", if (is.na(tried)) length(x) else tried,
           " simulation(s) converged",
           if (!is.null(conv)) " (debug: every simulation returned)" else ""))

  # Each returned fit's own four-level verdict. Reported rather than folded into
  # the status above: a gradient at zero and a clean convergence report are
  # different claims, and a self-test that converges everywhere while every
  # replicate carries a WARN is worth seeing.
  st <- vapply(x, function(m) {
    s <- m$convergence$status
    if (is.null(s) || !length(s)) NA_character_ else as.character(s)[1]
  }, character(1))
  st <- st[!is.na(st)]
  if (length(st)) {
    tally <- as.data.frame(table(status = factor(st, levels = .CONV_SEVERITY)),
                           stringsAsFactors = FALSE)
    tally <- tally[tally$Freq > 0, , drop = FALSE]
    # table() returns the grouping as a factor whatever stringsAsFactors says,
    # and sprintf() on one prints its level code.
    tally$status <- as.character(tally$status)
    tally$.tag <- .rce_sev_tag(tally$status)
    .rce_diag_table(tally, stats::setNames(
      c(".tag", "status", "Freq"), c(" ", "convergence", "n")))
  }

  proc <- attr(x, "process_sim")
  # Each deviation comes back beside a same-shaped `_drawn` mask (see sim_mod()),
  # which is bookkeeping rather than a process; naming it here would report
  # "rec_dev, rec_dev_drawn" as two things redrawn.
  drawn <- if (is.null(proc)) character(0) else
    grep("_drawn$", unique(unlist(lapply(proc, names))), value = TRUE,
         invert = TRUE)
  if (length(drawn)) {
    cat("  processes redrawn: ", paste(sort(drawn), collapse = ", "),
        " -- compare against attr(x, \"process_sim\"), not the operating model\n",
        sep = "")
  } else {
    cat("  observations redrawn, processes held at their fitted values:",
        "the spread is a lower bound\n")
  }
  invisible(x)
}
