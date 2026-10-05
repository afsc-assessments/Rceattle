# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("Test retrospective", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  library(Rceattle)


  # Data
  data(BS2017SS) # ?BS2017SS for more information on the data
  BS2017SS$projyr <- 2020
  BS2017SS$fleet_control$Proj_F_proportion <-rep(1,7)


  # Single-species with fixed M
  ss_run <- fit_mod(data_list = BS2017SS,
                    inits = NULL, # Initial parameters = 0
                    file = NULL, # Don't save
                    estimateMode = 1, # Estimate hindcast only
                    random_rec = FALSE, # No random recruitment
                    msmMode = 0, # Single species mode
                    fit_control = fit_control(getsd = FALSE, 
                      phase = TRUE,
                      verbose = 1))

  # Retro
  ret <-retrospective(ss_run, peels = 5)

  # By hand
  retro_list <- list()
  for(y in 1:5){
    dat_tmp <- BS2017SS
    dat_tmp$endyr <- dat_tmp$endyr - y
    retro_list[[y]] <- Rceattle::fit_mod(data_list = dat_tmp,
                                         inits = NULL, # Initial parameters = 0
                                         file = NULL, # Don't save
                                         estimateMode = 1, # Estimate hindcast only
                                         random_rec = FALSE, # No random recruitment
                                         msmMode = 0, # Single species mode
                                         fit_control = fit_control(getsd = FALSE, 
                                           phase = TRUE,
                                           loopnum = 5,
                                           verbose = 1))
  }
  retro_list <- rev(retro_list)


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Tests ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  yrs <- ss_run$data_list$styr:ss_run$data_list$endyr
  nyrs <- length(yrs)

  testthat::expect_equal(6, length(ret$Rceattle_list))
  testthat::expect_equal(2012:2017, as.numeric(sapply(ret$Rceattle_list, function(x) x$data_list$endyr_peel)))

  for(i in 1:5){
    # Matches
    yrs <- retro_list[[i]]$data_list$styr:retro_list[[i]]$data_list$endyr
    nyrs <- length(yrs)
    testthat::expect_equal(retro_list[[i]]$quantities$biomass[,1:nyrs], ret$Rceattle_list[[i]]$quantities$biomass[,1:nyrs], tolerance = 0.01)

    # Endyr of peel
    testthat::expect_equal(ret$Rceattle_list[[i]]$data_list$endyr_peel, 2017-rev(1:5)[i])

    # Data removed
    # * Index
    index_tmp <- ret$Rceattle_list[[i]]$data_list$index_data %>% filter(Year > 2017-rev(1:5)[i])
    testthat::expect_equal(nrow(index_tmp), 0)

    # * Comp
    comp_tmp <- ret$Rceattle_list[[i]]$data_list$comp_data %>% filter(Year > 2017-rev(1:5)[i])
    testthat::expect_equal(nrow(comp_tmp), 0)

    # * CAAL
    caal_tmp <- ret$Rceattle_list[[i]]$data_list$caal_data %>% filter(Year > 2017-rev(1:5)[i])
    testthat::expect_equal(nrow(caal_tmp), 0)
  }
})

testthat::test_that("self_test resolves getsd and completes (regression)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Regression: self_test()'s per-sim closure run_one_sim references `getsd`, but
  # (unlike retrospective/jitter/profile) self_test() never defined it -- no
  # argument, no default -- so every refit died with "object 'getsd' not found",
  # in both the sequential and the parallel-cluster dispatch. Before the fix the
  # self_test() call below errors (never returns), so reaching the assertions at
  # all proves the closure now resolves getsd. A small synthetic fixture keeps
  # this fast; the refits' convergence count is not the point.
  d <- make_test_data()
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, file = NULL, estimateMode = 1,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  # Returns a (possibly empty, if none converged) list of refit Rceattle models
  # named Sim_1, Sim_2, ... -- NOT a $Rceattle_list/$nll wrapper (that is jitter's).
  st <- suppressMessages(suppressWarnings(self_test(fit, nsim = 2, cores = 1)))
  testthat::expect_type(st, "list")
  if (length(st) > 0) {
    testthat::expect_true(all(vapply(st, inherits, logical(1), "Rceattle")))
    testthat::expect_true(all(grepl("^Sim_", names(st))))
  }

  # getsd is an accepted argument (the crux of the fix).
  testthat::expect_no_error(
    suppressMessages(suppressWarnings(self_test(fit, nsim = 1, cores = 1, getsd = FALSE))))
})


testthat::test_that(".parallel_lapply falls back to sequential when a worker dies", {
  # Each worker fits whole models, so several large ones at once can exhaust
  # memory and the OS kills one. All the parent sees is
  # "Error in unserialize(node$con) : error reading from connection", and every
  # peel or simulation finished up to that point is thrown away. That is a
  # non-Windows-only failure: retrospective(), jitter() and run_mse() take a FORK
  # cluster there and PSOCK on Windows, so the platforms do not fail alike.
  #
  # Stub the cluster constructor rather than trying to starve a real worker --
  # the point is what .parallel_lapply() does with the error, not how the OS
  # produced it.
  boom <- function(...) stop("error reading from connection")

  testthat::local_mocked_bindings(
    makeCluster = boom, .package = "parallel")

  out <- NULL
  testthat::expect_message(
    out <- Rceattle:::.parallel_lapply(1:3, function(i) i^2, 2, environment()),
    "computing the 3 tasks sequentially")

  # Falling back has to give the same answer as the cluster would have.
  testthat::expect_equal(out, list(1, 4, 9))
})


testthat::test_that(".parallel_lapply agrees with lapply on a real cluster", {
  testthat::skip_on_cran()

  # Nothing else in the suite passes cores > 1, so without this the cluster
  # itself -- makeCluster, the FORK/PSOCK branch, parLapply -- is never run, and
  # only the fallback around it would be covered. Kept to a pure function so it
  # stays fast and cannot flake on a model fit.
  items <- 1:4
  f <- function(i) list(sq = i^2, chr = letters[i])
  testthat::expect_equal(
    Rceattle:::.parallel_lapply(items, f, 2, environment()),
    lapply(items, f))
})


testthat::test_that(".parallel_lapply still surfaces a genuine error in fun()", {
  # The fallback must not turn a real bug into a silent serial re-run that
  # reports the cluster as the problem. A failure inside the worker closure
  # propagates out of the sequential retry as itself.
  testthat::local_mocked_bindings(
    makeCluster = function(...) stop("error reading from connection"),
    .package = "parallel")

  testthat::expect_error(
    suppressMessages(
      Rceattle:::.parallel_lapply(1:2, function(i) stop("bug in the peel"), 2,
                                  environment())),
    "bug in the peel")
})


testthat::test_that("retrospective survives dropped peels (regression)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # retrospective() named its output `Year_(endyr-peels):endyr` -- always
  # peels + 1 labels -- and looped `1:(length(mod_list) - 1)`. Both assume every
  # peel survives. Dropping was effectively unreachable while the keep/drop gate
  # compared a string that is never assigned (see .refit_converged), so the
  # assumption held by accident; making the gate work exposes it. With one peel
  # dropped the names assignment errors ("'names' attribute [3] must be the same
  # length as the vector [2]"), and with all of them dropped `1:0` counts DOWN
  # and indexes mod_list[[2]] ("subscript out of bounds").
  #
  # Drive it through the real function by stubbing the gate, so this pins the
  # surviving code paths rather than a hand-built list.
  d <- make_test_data()
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, file = NULL, estimateMode = 0,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  endyr <- fit$data_list$endyr

  # Key the stub off the model's OWN terminal year, not a call counter: the gate
  # is called twice per peel but `hindcast_converged && .refit_converged(newmod)`
  # short-circuits, so a peel that fails the first check never makes the second
  # and any count-based scheme silently drifts out of step.
  keep_only <- function(keep_endyrs) {
    function(newmod, max_grad = 1) {
      isTRUE(newmod$data_list$endyr_peel %in% keep_endyrs)
    }
  }

  # 1 of 2 peels dropped: labels must follow the surviving peel's terminal year,
  # not a fixed endyr-peels:endyr sequence.
  local({
    testthat::local_mocked_bindings(.refit_converged = keep_only(endyr - 1),
                                    .package = "Rceattle")
    r <- suppressMessages(suppressWarnings(
      retrospective(fit, peels = 2, cores = 1, getsd = FALSE)))
    testthat::expect_length(r$Rceattle_list, 2L)      # survivor + full model
    endyrs <- vapply(r$Rceattle_list,
                     function(x) as.numeric(x$data_list$endyr_peel), numeric(1))
    testthat::expect_equal(names(r$Rceattle_list), paste0("Year_", endyrs))
    testthat::expect_equal(unname(endyrs), c(endyr - 1, endyr))
    testthat::expect_false(anyNA(names(r$Rceattle_list)))
  })

  # Every peel dropped: only the original model remains, and Mohn's rho has
  # nothing to average -- it must warn and return, not error.
  local({
    testthat::local_mocked_bindings(.refit_converged = keep_only(numeric(0)),
                                    .package = "Rceattle")
    # Two warnings now: the drop itself, then rho having nothing to average.
    testthat::expect_warning(
      testthat::expect_warning(
        r <- suppressMessages(retrospective(fit, peels = 2, cores = 1, getsd = FALSE)),
        "2 of 2 peels dropped"),
      "Mohn's rho is undefined")
    testthat::expect_length(r$Rceattle_list, 1L)
    testthat::expect_equal(names(r$Rceattle_list), paste0("Year_", endyr))
  })

  # And the drop WARNS rather than messages: rho is averaged over the peels
  # that survive, so a drop changes a reported number, and `suppressMessages()`
  # is common enough in the assessment scripts to hide a message.
  local({
    testthat::local_mocked_bindings(.refit_converged = keep_only(endyr - 1),
                                    .package = "Rceattle")
    testthat::expect_warning(
      r <- suppressMessages(retrospective(fit, peels = 2, cores = 1, getsd = FALSE)),
      "1 of 2 peel dropped")
    # The count asked for rides on the object, so print() can still say a peel
    # was lost once the warning has scrolled away.
    testthat::expect_equal(r$peels_requested, 2)
    testthat::expect_true(any(grepl(
      "1 peel\\(s\\) dropped", utils::capture.output(print(r)))))
  })
})


testthat::test_that("each peel reports its own terminal year, and rho is unmoved", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Every plot builds its year axis per model as `styr:endyr`, and nothing
  # outside retrospective() reads `endyr_peel` -- so while each peel reported the
  # FULL model's `endyr`, all the peels were drawn to the same terminal year and
  # were indistinguishable, which is the opposite of what a retrospective plot
  # is for.
  #
  # The assignment must happen AFTER both refits. `endyr` sizes the model, and
  # the forecast refit turns F back on over `(nyrs_peel+1):nyrs` against the FULL
  # nyrs, so peeling `endyr` any earlier indexes off the end of log_F. The
  # projection-width check below is what catches a regression that moves it.
  d <- make_test_data()
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, file = NULL, estimateMode = 0,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  r <- suppressMessages(suppressWarnings(
    retrospective(fit, peels = 2, cores = 1, getsd = FALSE)))

  full_endyr <- fit$data_list$endyr
  endyrs <- vapply(r$Rceattle_list,
                   function(x) as.numeric(x$data_list$endyr), numeric(1))
  peels  <- vapply(r$Rceattle_list,
                   function(x) as.numeric(x$data_list$endyr_peel), numeric(1))

  # The point of the change.
  testthat::expect_equal(unname(endyrs), unname(peels))

  # The peels really do differ, or plotting gains nothing; and the unpeeled
  # model (last, after the rev()) keeps the true terminal year.
  testthat::expect_true(any(endyrs < full_endyr))
  testthat::expect_equal(unname(endyrs[length(endyrs)]), full_endyr)

  # `endyr` is output metadata only: each peel still carries its retrospective
  # forecast, so its quantities span the full projection. This is the assertion
  # that fails if the assignment is moved ahead of the refits.
  for (m in r$Rceattle_list) {
    testthat::expect_equal(ncol(m$quantities$ssb),
                           m$data_list$projyr - m$data_list$styr + 1L)
  }

  # Mohn's rho reads `endyr_peel` and the full model's `endyr` from the calling
  # frame, never a peel's `data_list$endyr`, so it is still computed over every
  # surviving peel rather than collapsing to NaN.
  scored <- r$mohns[r$mohns$N > 0, , drop = FALSE]
  testthat::expect_true(nrow(scored) > 0)
  testthat::expect_false(any(is.na(unlist(scored[, -(1:3)]))))

  # Setting `endyr` to the peel makes it duplicate `endyr_peel`, which would
  # otherwise lose the unpeeled terminal year -- and with it the boundary
  # between the retrospective forecast and the true projection. `endyr_full`
  # keeps that recoverable from a peel on its own.
  fulls <- vapply(r$Rceattle_list,
                  function(x) as.numeric(x$data_list$endyr_full), numeric(1))
  testthat::expect_equal(unname(fulls), rep(full_endyr, length(fulls)))
  testthat::expect_true(all(fulls >= peels))
  testthat::expect_true(all(fulls <= vapply(r$Rceattle_list,
    function(x) as.numeric(x$data_list$projyr), numeric(1))))
})


testthat::test_that("self_test inherits the source model's phasing", {
  # self_test() refits from `initial_params`, so it has to cover the same ground
  # the original fit did -- but it was the only refitting diagnostic that never
  # phased (retrospective fixes phase = TRUE, jitter exposes it), taking
  # .refit_like()'s phase = FALSE default with no way to change it. On a model
  # that needs phasing the refits then land many orders of magnitude from the
  # optimum and are all dropped. (Since 5.11.0 retrospective() exposes `phase`
  # as an argument too, defaulting to the TRUE it used to hardcode.)
  #
  # Phasing is a property of the refit, not of the fixture, so assert on the
  # resolved default rather than on a fit: `phase_params` is attached by
  # fit_mod() only when it phased, exactly as `sdrep` is attached only when it
  # ran sdreport. Both defaults read the source model the same way.
  fml <- formals(self_test)
  testthat::expect_true("phase" %in% names(fml))
  testthat::expect_null(eval(fml$phase))          # NULL = inherit

  # The new arguments are appended, so documented positional calls keep meaning.
  # The fitted model is `object` since 5.11.0; `Rceattle` still works and is
  # declared last, which is what keeps positional calls binding to `object`.
  testthat::expect_equal(names(fml)[1:6],
                         c("object", "nsim", "simulate", "seed", "cores", "getsd"))
  testthat::expect_equal(names(fml)[length(fml)], "Rceattle")

  phased   <- structure(list(phase_params = list(dummy = 0)), class = "Rceattle")
  unphased <- structure(list(), class = "Rceattle")
  testthat::expect_true(!is.null(phased$phase_params))
  testthat::expect_false(!is.null(unphased$phase_params))
})


testthat::test_that("self_test passes phase and start through to the refit", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Asserting that the arguments EXIST is not the same as asserting they are
  # used: the switch could be inverted, or `phase` never forwarded, and the
  # signature tests above would still pass. Capture what .refit_like() actually
  # receives, and stop before any real fitting.
  d <- make_test_data()
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, file = NULL, estimateMode = 3,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  fit$initial_params   <- list(tag = "INITIAL")
  fit$estimated_params <- list(tag = "ESTIMATED")

  seen <- NULL
  capture <- function(...) {
    args <- list(...)
    seen <<- c(seen, list(list(phase = args$phase, tag = args$inits$tag)))
    NULL   # dropped by the gate; we only care about the call
  }

  grab <- function(...) {
    seen <<- NULL
    # .package explicitly: without it local_mocked_bindings() falls back to
    # pkgload's dev package, which does not exist when the tests run against an
    # INSTALLED package -- i.e. under R CMD check, which is what CI runs.
    testthat::local_mocked_bindings(.refit_like = capture, .package = "Rceattle")
    suppressMessages(suppressWarnings(self_test(fit, nsim = 1, cores = 1, ...)))
    seen[[1]]
  }

  # phase: an explicit value always wins.
  testthat::expect_true(grab(phase = TRUE)$phase)
  testthat::expect_false(grab(phase = FALSE)$phase)

  # NULL inherits, primarily from the value fit_mod() recorded on the fit -- so a
  # custom phase LIST carries over as the list instead of collapsing to TRUE.
  testthat::expect_false(grab()$phase)                 # fixture recorded FALSE
  fit$run_config$fit_control$phase <- TRUE
  testthat::expect_true(grab()$phase)
  fit$run_config$fit_control$phase <- list(dummy = 1)
  testthat::expect_equal(grab()$phase, list(dummy = 1))

  # Fallback for fits predating run_config: whether phasing left phase_params.
  fit$run_config <- NULL
  testthat::expect_false(grab()$phase)
  fit$phase_params <- list(dummy = 0)
  testthat::expect_true(grab()$phase)

  # start: selects which parameter set is handed over, and defaults to initial.
  testthat::expect_equal(grab(start = "initial")$tag,   "INITIAL")
  testthat::expect_equal(grab(start = "estimated")$tag, "ESTIMATED")
  testthat::expect_equal(grab()$tag,                    "INITIAL")
})


testthat::test_that("self_test start= is validated before any fitting", {
  # match.arg + the carries-it check both run before the first simulation, so a
  # typo or a stripped model fails immediately instead of after nsim refits.
  fml <- formals(self_test)
  testthat::expect_equal(eval(fml$start), c("initial", "estimated"))
  testthat::expect_false(eval(fml$debug))

  full <- structure(list(initial_params = list(a = 1), estimated_params = list(a = 2)),
                    class = "Rceattle")
  testthat::expect_error(self_test(full, start = "mle"), "should be one of")

  # model_average() strips both parameter sets; asking to start from one of them
  # should say so rather than pass NULL inits to fit_mod().
  stripped <- structure(list(initial_params = NULL, estimated_params = NULL),
                        class = "Rceattle")
  testthat::expect_error(self_test(stripped, start = "initial"), "initial_params")
  testthat::expect_error(self_test(stripped, start = "estimated"), "estimated_params")

  # Not an Rceattle is still caught first.
  testthat::expect_error(self_test(list()), "not of class")
})


testthat::test_that("self_test debug= returns every simulation with its verdict", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data()
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, file = NULL, estimateMode = 1,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  dbg <- suppressMessages(suppressWarnings(
    self_test(fit, nsim = 2, cores = 1, debug = TRUE)))

  # Every simulation comes back, indexed by simulation number (Sim_i <-> seed+i),
  # never filtered -- that is the point of debug.
  testthat::expect_length(dbg, 2L)
  testthat::expect_equal(names(dbg), c("Sim_1", "Sim_2"))
  testthat::expect_true(all(vapply(dbg, inherits, logical(1), "Rceattle")))

  verdict <- attr(dbg, "converged")
  testthat::expect_type(verdict, "logical")
  testthat::expect_equal(names(verdict), names(dbg))

  # The default return is exactly the converged subset, renumbered contiguously.
  plain <- suppressMessages(suppressWarnings(self_test(fit, nsim = 2, cores = 1)))
  testthat::expect_null(attr(plain, "converged"))
  testthat::expect_length(plain, sum(verdict))
  if (length(plain) > 0) {
    testthat::expect_equal(names(plain), paste0("Sim_", seq_along(plain)))
  }

  # A dropped run still carries the diagnostics that explain why it was dropped.
  testthat::expect_true(all(vapply(
    dbg, function(x) is.null(x$convergence) || inherits(x$convergence, "Rceattle_convergence"),
    logical(1))))
})


testthat::test_that(".refit_converged judges the gradient, not the string", {
  # The gate the four refitting diagnostics share. It used to compare
  # opt$Convergence_check against TMBhelper's non-invertible-Hessian string,
  # which is unreachable when getsd = TRUE (fit_tmb returns before assigning it)
  # and never assigned at all when getsd = FALSE -- so the check either dropped
  # runs incidentally, via the enclosing is.null() guard, or dropped nothing.
  mk <- function(max_gradient, sd_requested = FALSE, sd_present = FALSE) {
    structure(list(.conv_hindcast = list(max_gradient = max_gradient,
                                         sd_requested = sd_requested,
                                         sd_present   = sd_present)),
              class = "Rceattle")
  }

  testthat::expect_true(Rceattle:::.refit_converged(mk(1e-4)))
  testthat::expect_true(Rceattle:::.refit_converged(mk(0.5)))   # under the FAIL tier

  # The case the old gate kept whenever getsd = FALSE.
  testthat::expect_false(Rceattle:::.refit_converged(mk(1e13)))
  testthat::expect_false(Rceattle:::.refit_converged(mk(1.0001)))
  # At the boundary the gate must agree with .check_optimizer(), which FAILs on
  # `mg > 1` -- so exactly 1 is a WARN the gate keeps, not a FAIL it drops.
  testthat::expect_true(Rceattle:::.refit_converged(mk(1)))
  testthat::expect_false(Rceattle:::.refit_converged(mk(NA_real_)))
  testthat::expect_false(Rceattle:::.refit_converged(mk(numeric(0))))
  testthat::expect_false(Rceattle:::.refit_converged(NULL))

  # A requested sdreport that did not come back still drops the run: that is the
  # condition the old string test was reaching for, and it is kept.
  testthat::expect_false(
    Rceattle:::.refit_converged(mk(1e-4, sd_requested = TRUE, sd_present = FALSE)))
  testthat::expect_true(
    Rceattle:::.refit_converged(mk(1e-4, sd_requested = TRUE, sd_present = TRUE)))
  # ... but only when one was asked for.
  testthat::expect_true(
    Rceattle:::.refit_converged(mk(1e-4, sd_requested = FALSE, sd_present = FALSE)))

  # No hindcast snapshot (a refit at an estimateMode that does not optimize):
  # fall back to the optimizer's verdict, which is what these call sites did
  # before -- and which is absent for those modes, so they still drop.
  no_snap <- function(cc) structure(list(opt = if (is.null(cc)) NULL else
                                           list(Convergence_check = cc)),
                                    class = "Rceattle")
  testthat::expect_false(Rceattle:::.refit_converged(no_snap(NULL)))
  testthat::expect_false(
    Rceattle:::.refit_converged(no_snap("The model is definitely not converged")))
  testthat::expect_true(
    Rceattle:::.refit_converged(no_snap("There is no evidence that the model is not converged")))
})


testthat::test_that(".refit_with_timeout bounds a run and returns rather than throws", {
  # The failure this exists for is a hang, so the contract is: stop it, and hand
  # the condition back instead of throwing (one replicate must not abort a run).
  fast <- Rceattle:::.refit_with_timeout(41 + 1, timeout = 30)
  testthat::expect_equal(fast, 42)

  # Inf / missing means no limit, and must not leave one set behind either.
  testthat::expect_equal(Rceattle:::.refit_with_timeout(42), 42)
  testthat::expect_equal(Rceattle:::.refit_with_timeout(42, timeout = Inf), 42)

  # An ordinary error comes back as the condition, not thrown.
  boom <- Rceattle:::.refit_with_timeout(stop("kaboom"), timeout = Inf)
  testthat::expect_s3_class(boom, "error")
  testthat::expect_match(conditionMessage(boom), "kaboom")
  testthat::expect_false(Rceattle:::.is_timeout(boom))

  # A run that overruns is stopped and identified as a timeout, not as a model
  # error -- the two get different advice.
  slow <- Rceattle:::.refit_with_timeout(
    { repeat invisible(Sys.time()) }, timeout = 1)
  testthat::expect_s3_class(slow, "condition")
  testthat::expect_true(Rceattle:::.is_timeout(slow))

  # The limit must not leak into whatever runs next: a per-replicate limit set
  # with transient = FALSE stays in force until it is reset, so a missed reset
  # would kill every later replicate.
  testthat::expect_equal(Rceattle:::.refit_with_timeout(42, timeout = Inf), 42)
  later <- Rceattle:::.refit_with_timeout({ Sys.sleep(1.5); "done" }, timeout = Inf)
  testthat::expect_equal(later, "done")
})


testthat::test_that(".report_dropped and .report_errors say the right thing", {
  testthat::expect_silent(Rceattle:::.report_dropped(0, 5, "peel"))
  testthat::expect_message(Rceattle:::.report_dropped(1, 5, "peel"),
                           "1 of 5 peel dropped")
  testthat::expect_message(Rceattle:::.report_dropped(2, 5, "peel"),
                           "2 of 5 peels dropped.*3 returned")

  testthat::expect_silent(Rceattle:::.report_errors(c(NA_character_, NA_character_), "jitter"))
  # A timeout is reported as a timeout (raise the limit), not as an error.
  testthat::expect_message(
    Rceattle:::.report_errors(c(NA, "reached elapsed time limit"), "jitter"),
    "exceeded `timeout`")
  testthat::expect_message(
    Rceattle:::.report_errors(c(NA, "singular matrix"), "jitter"),
    "errored; first: singular matrix")
})


testthat::test_that(".normalize_fit_tmb unwraps fit_tmb's non-PD early return", {
  # With getsd = TRUE and a Hessian that fails chol(), TMBhelper::fit_tmb()
  # returns list(opt = , h = ) instead of the estimates -- a different shape,
  # with no $objective / $max_gradient / $Convergence_check at the top level.
  est <- list(par = c(a = 1), objective = 42, max_gradient = 3,
              Convergence_check = "The model is likely not converged")
  h <- matrix(c(1, 2, 2, 1), 2)

  got <- Rceattle:::.normalize_fit_tmb(list(opt = est, h = h))
  testthat::expect_equal(got$objective, 42)
  testthat::expect_equal(got$max_gradient, 3)
  testthat::expect_equal(got$hessian, h)
  testthat::expect_equal(got$Convergence_check, "The model is definitely not converged")
  # No $SD, so fit_mod() still records sdrep = NULL and sdreport_failed fires.
  testthat::expect_null(got$SD)

  # A no-op on the documented shape, which the guard tells apart by which names
  # are present -- so it has to match names exactly (`[[`, not `$`).
  testthat::expect_identical(Rceattle:::.normalize_fit_tmb(est), est)
  testthat::expect_identical(Rceattle:::.normalize_fit_tmb(c(est, list(SD = "sd"))),
                             c(est, list(SD = "sd")))
})


testthat::test_that("a peel reports standard errors for its own hindcast", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # A peel is fitted twice, and the second refit estimates only the peeled
  # years' F. Reported from that, every hindcast quantity gets a standard error
  # of exactly zero, so a peel drew with no confidence band while the unpeeled
  # model drew one. The peel is now reported from the same parameters with the
  # hindcast free in the map, built without optimizing.
  # BS2017SS rather than make_test_data(): `getsd = TRUE` also turns on the
  # positive-definite Hessian check in .refit_converged(), and the small fixture
  # loses every peel to it, leaving nothing to assert on.
  data(BS2017SS, envir = environment())
  d <- BS2017SS
  d$projyr <- 2020
  d$fleet_control$Proj_F_proportion <- rep(1, 7)

  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, getsd = TRUE, verbose = 0))))

  r <- suppressMessages(suppressWarnings(
    retrospective(fit, peels = 1, cores = 1, getsd = TRUE)))

  # The last entry is the input model passed through; the refits are the rest.
  full_endyr <- fit$data_list$endyr
  is_refit <- vapply(r$Rceattle_list,
                     function(x) x$data_list$endyr_peel < full_endyr, logical(1))
  testthat::expect_true(any(is_refit))

  ssb_se <- function(m) {
    s <- m$sdrep$sd[which(names(m$sdrep$value) == "ssb")]
    nh <- m$data_list$endyr_peel - m$data_list$styr + 1
    matrix(s, nrow = m$data_list$nspp)[1, seq_len(nh)]
  }

  for (m in r$Rceattle_list[is_refit]) {
    se <- ssb_se(m)
    testthat::expect_true(all(is.finite(se)))
    testthat::expect_true(all(se > 0))   # the defect: every one was exactly 0
  }

  # The report model is built at estimateMode 3, which skips the projection
  # stage. If that changed the ADREPORT span, every plotter indexing by year
  # would silently misalign.
  n_ssb <- function(m) sum(names(m$sdrep$value) == "ssb")
  for (m in r$Rceattle_list[is_refit]) {
    testthat::expect_equal(n_ssb(m), n_ssb(fit))
  }

  # Nothing may be re-estimated: a peel must have the same SSB and the same
  # objective whether or not standard errors were asked for. This is what stops
  # the report pass being "improved" into a refit.
  #
  # Compared peel by peel and not run against run: `getsd = TRUE` also turns on
  # the positive-definite Hessian check in .refit_converged(), which can drop a
  # peel that `getsd = FALSE` keeps, so the two lists need not hold the same
  # peels.
  r_no <- suppressMessages(suppressWarnings(
    retrospective(fit, peels = 1, cores = 1, getsd = FALSE)))

  shared <- intersect(names(r$Rceattle_list), names(r_no$Rceattle_list))
  testthat::expect_gt(length(shared), 0)
  for (nm in shared) {
    testthat::expect_equal(r$Rceattle_list[[nm]]$quantities$ssb,
                           r_no$Rceattle_list[[nm]]$quantities$ssb)
    testthat::expect_equal(r$Rceattle_list[[nm]]$opt$objective,
                           r_no$Rceattle_list[[nm]]$opt$objective)
  }

  # And it is not paid for when it was not asked for.
  no_refit <- vapply(r_no$Rceattle_list,
                     function(x) x$data_list$endyr_peel < full_endyr, logical(1))
  testthat::expect_true(all(vapply(r_no$Rceattle_list[no_refit],
                                   function(m) is.null(m$sdrep), logical(1))))
})


# `peels` as a vector of peel depths.
#
# A retrospective used to be a count: `peels = n` meant 1:n. A caller needs to
# be able to ask for 2:10 and skip the first -- each peel costs a model fit, and
# the shallow ones contribute only to the shortest horizons that
# hindcast_skill() scores. Each element is a number of YEARS
# REMOVED from endyr, not a year.
testthat::test_that("retrospective takes a vector of peel depths", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data()
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, file = NULL, estimateMode = 0,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  endyr <- fit$data_list$endyr

  depths <- function(r) {
    sort(endyr - vapply(r$Rceattle_list,
                        function(x) as.numeric(x$data_list$endyr_peel),
                        numeric(1)))
  }

  # The named depths, and only those: no 1-year peel.
  r <- suppressMessages(suppressWarnings(
    retrospective(fit, peels = 2:3, cores = 1, getsd = FALSE)))
  testthat::expect_equal(unname(depths(r)), c(0, 2, 3))  # 0 is the unpeeled model
  testthat::expect_equal(r$peel_depths, 2:3)
  testthat::expect_equal(r$peels_requested, 2L)

  # Out of order and duplicated resolves to the same thing, so the returned list
  # runs deepest-first however the depths were written.
  r_msg <- suppressMessages(suppressWarnings(
    retrospective(fit, peels = c(3, 2, 3), cores = 1, getsd = FALSE)))
  testthat::expect_equal(r_msg$peel_depths, 2:3)
  testthat::expect_equal(names(r_msg$Rceattle_list), names(r$Rceattle_list))

  # A scalar still means 1:n.
  r1 <- suppressMessages(suppressWarnings(
    retrospective(fit, peels = 2, cores = 1, getsd = FALSE)))
  testthat::expect_equal(r1$peel_depths, 1:2)
  testthat::expect_equal(unname(depths(r1)), c(0, 1, 2))

  # A peel reads only the input model, so WHICH other depths were asked for
  # cannot move it. This is what makes a subset request legitimate: the 2-year
  # peel from `peels = 2:3` has to be the same fit as the one from `peels = 2`.
  nm <- paste0("Year_", endyr - 2)
  testthat::expect_true(nm %in% names(r$Rceattle_list))
  testthat::expect_equal(r$Rceattle_list[[nm]]$quantities$ssb,
                         r1$Rceattle_list[[nm]]$quantities$ssb)
  testthat::expect_equal(r$Rceattle_list[[nm]]$opt$objective,
                         r1$Rceattle_list[[nm]]$opt$objective)
})


# A bad `peels` must be refused before any peel is dispatched. Left to
# run_one_peel(), a depth of 0 makes `(endyr_peel + 1):endyr` count DOWN and the
# retained years get fitted as a forecast, and an over-deep one builds a model
# that ends before styr -- both silently, after paying for the fits.
testthat::test_that("retrospective refuses an unusable `peels`", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data()
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, file = NULL, estimateMode = 3,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  nyr <- fit$data_list$endyr - fit$data_list$styr

  testthat::expect_error(retrospective(fit, peels = 0, cores = 1),
                         "must be at least 1")
  testthat::expect_error(retrospective(fit, peels = c(2, 0), cores = 1),
                         "must be at least 1")
  testthat::expect_error(retrospective(fit, peels = -1, cores = 1),
                         "must be at least 1")
  testthat::expect_error(retrospective(fit, peels = 2.5, cores = 1),
                         "whole numbers")
  testthat::expect_error(retrospective(fit, peels = c(1, NA), cores = 1),
                         "containing NA")
  testthat::expect_error(retrospective(fit, peels = integer(0), cores = 1),
                         "positive whole number")
  testthat::expect_error(retrospective(fit, peels = "3", cores = 1),
                         "positive whole number")
  # A list reaches the message builder, where `is.finite()` has no list method
  # and errors outright unless the clause is guarded on is.numeric(). And a
  # character must not be diagnosed "containing Inf" -- is.finite("3") is FALSE,
  # so an unguarded clause fires on every non-numeric atomic. Assert the
  # diagnosis, not just that something was refused.
  testthat::expect_error(retrospective(fit, peels = list(1), cores = 1),
                         "got a list of length 1")
  testthat::expect_error(retrospective(fit, peels = list(1, 2), cores = 1),
                         "got a list of length 2")
  testthat::expect_error(retrospective(fit, peels = "3", cores = 1),
                         "got a character of length 1\\.")
  testthat::expect_error(retrospective(fit, peels = Inf, cores = 1),
                         "containing Inf")
  # Deeper than the model can carry: named with the deepest depth that works,
  # so the message is actionable rather than a convergence failure.
  testthat::expect_error(retrospective(fit, peels = nyr, cores = 1),
                         "deepest peel that leaves")
  testthat::expect_error(retrospective(fit, peels = c(2, nyr + 100), cores = 1),
                         "deepest peel that leaves")
})


# A length-1 `peels` cannot name a single depth -- `peels = 5` and `peels = c(5)`
# are the same object in R, and back-compatibility fixes that to mean 1:5. The
# de-duplicating path makes `c(5, 5)` the way to ask for one peel, which is worth
# pinning: it is the documented workaround, not an accident to be optimized away.
testthat::test_that("a repeated depth asks for that one peel", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data()
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, file = NULL, estimateMode = 0,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  r <- suppressMessages(suppressWarnings(
    retrospective(fit, peels = c(3, 3), cores = 1, getsd = FALSE)))
  testthat::expect_equal(r$peel_depths, 3L)
  testthat::expect_equal(r$peels_requested, 1L)
  testthat::expect_equal(
    sort(fit$data_list$endyr - vapply(r$Rceattle_list,
      function(x) as.numeric(x$data_list$endyr_peel), numeric(1))),
    c(0, 3), ignore_attr = TRUE)
})


# An unrecorded `random_vars` is reported once, by retrospective(), before any
# peel is dispatched.
#
# It used to be raised inside .rce_peel_map(), which runs inside the per-peel
# closure, so it was discarded by the FORK/PSOCK worker that ran the peel -- and
# the default `cores` is detectCores() - 6, so the warning had never reached
# anyone running a multi-peel retrospective. It announces a BIASED Mohn's rho
# (-6.6% on sigma at 5 peels, monotone in peel depth), which is the last thing
# that should be silent.
#
# It also fired on fits that were perfectly fine: `random_vars` started as `c()`
# and `mod_objects$random_vars <- NULL` DELETES the element, so an ordinary
# random_rec = FALSE fit was indistinguishable from one made before 5.10.0
# recorded it. fit_mod() now stores character(0) for "recorded, and there were
# none", and only an absent record warns.
testthat::test_that("retrospective reports an unrecorded random_vars, once, and not an empty one", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data()
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, file = NULL, estimateMode = 1, random_rec = FALSE,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  # A fit with no random effects RECORDS that fact rather than omitting it.
  testthat::expect_true("random_vars" %in% names(fit))
  testthat::expect_identical(fit$random_vars, character(0))

  warns_of <- function(f) {
    w <- character(0)
    withCallingHandlers(
      suppressMessages(try(retrospective(f, peels = c(2, 2), cores = 1,
                                         getsd = FALSE), silent = TRUE)),
      warning = function(cnd) {
        w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
      })
    w
  }

  # Recorded-and-empty: pinning every deviation is what the fit itself did, so
  # there is nothing to report.
  testthat::expect_length(
    grep("does not record which blocks", warns_of(fit)), 0L)

  # Unrecorded: reported, and exactly once however many peels are asked for.
  old <- fit
  old$random_vars <- NULL
  testthat::expect_length(
    grep("does not record which blocks", warns_of(old)), 1L)

  # Once per CALL, not once per peel -- and a subset that never includes peel 1
  # still gets it, which an `i == 1L` gate inside the closure would not.
  w3 <- character(0)
  withCallingHandlers(
    suppressMessages(try(retrospective(old, peels = 2:3, cores = 1,
                                       getsd = FALSE), silent = TRUE)),
    warning = function(cnd) {
      w3 <<- c(w3, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  testthat::expect_length(grep("does not record which blocks", w3), 1L)
})


# forecast_rec = "model" on a fit WITHOUT a DSEM.
#
# Regression: the predicate behind "do the latent states supply the forecast?"
# called a helper `.pin()` that was never defined anywhere in the package. `&&`
# and `||` short-circuit, so it was reached only when forecast_rec = "model" met
# proj_mean_rec = FALSE on a fit with no DSEM -- and every DSEM test took the
# `.has_dsem()` arm before it. That combination errored with "could not find
# function '.pin'", which made hindcast_skill() -- whose forecast_rec DEFAULTS to
# "model" -- unusable on every non-DSEM model.
#
# With rec_dev a fixed effect the forecast comes off the stock-recruit curve (a
# zero deviation), which is the third arm of the precedence order documented on
# `forecast_rec`; with it a random effect the states supply it. Exercise both, so
# the predicate is evaluated either way rather than short-circuited.
testthat::test_that("forecast_rec = 'model' runs on a fit with no DSEM", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data()
  run <- function(rr) {
    fit <- suppressMessages(suppressWarnings(fit_mod(
      data_list = d, file = NULL, estimateMode = 1, random_rec = rr,
      recFun = build_srr(srr_fun = 0, proj_mean_rec = FALSE),
      fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
    # The guard above is only meaningful if the fit really is DSEM-free and
    # really did keep proj_mean_rec = FALSE -- otherwise the predicate is
    # short-circuited and this test proves nothing.
    testthat::expect_null(fit[["dsem"]])
    testthat::expect_false(isTRUE(as.logical(fit[["data_list"]][["proj_mean_rec"]])))
    list(fit = fit, retro = suppressMessages(suppressWarnings(
      retrospective(fit, peels = c(2, 2), cores = 1, getsd = FALSE,
                    forecast_rec = "model"))))
  }

  out <- list()
  fits <- list()
  for (rr in c(FALSE, TRUE)) {
    z <- testthat::expect_no_error(run(rr))
    testthat::expect_s3_class(z[["retro"]], "Rceattle_retro")
    testthat::expect_equal(z[["retro"]][["peel_depths"]], 2L)
    out[[as.character(rr)]]  <- z[["retro"]]
    fits[[as.character(rr)]] <- z[["fit"]]
  }

  # The predicate behind the branch, asserted directly: it must disagree between
  # the two arms, or `forecast_rec = "model"` cannot tell them apart at all. This
  # is what an inverted or constant `pin_block()` would break.
  testthat::expect_identical(fits[["FALSE"]][["random_vars"]], character(0))
  testthat::expect_true("rec_dev" %in% fits[["TRUE"]][["random_vars"]])

  # With rec_dev a FIXED effect the peeled years take the stock-recruit curve's
  # zero deviation -- the third arm of the precedence order on `forecast_rec`.
  endyr <- fits[["FALSE"]][["data_list"]][["endyr"]]
  styr  <- fits[["FALSE"]][["data_list"]][["styr"]]
  fore  <- (endyr - 2 - styr + 2L):(endyr - styr + 1L)
  dev_of <- function(r) {
    m <- r[["Rceattle_list"]][[paste0("Year_", endyr - 2)]]
    as.numeric(m[["estimated_params"]][["rec_dev"]][1, fore])
  }
  testthat::expect_true(all(dev_of(out[["FALSE"]]) == 0))

  # The random-effect arm is NOT asserted by value, deliberately. On this
  # fixture that fit drives sigma_R to zero -- measured: every hindcast rec_dev
  # comes back around 2e-86 -- so the data-free tail the Laplace approximation
  # integrates out sits at the prior mean of 0 as well, and the two arms are
  # numerically indistinguishable here for a legitimate reason. The branch where
  # latent states really do supply a different forecast is covered on a LAGGED
  # DSEM, where it moves recruitment 46.8% (test-dsem-retrospective.R); an
  # IID process with no data has nothing to propagate.
  testthat::expect_true(all(abs(dev_of(out[["TRUE"]])) < 1e-10))
})


# A peel whose sdreport() throws loses its band, and SAYS so.
#
# Two regressions in one handler. It returned `newmod$sdrep` on failure, which
# at that point is the FORECAST refit's sdreport -- built with the whole
# hindcast pinned, so every hindcast standard error in it is exactly zero. That
# is the thing the report pass exists to replace, so a peel that could not be
# reported handed back a zero-width band as though it were real. And it said so
# with a warning() raised inside the per-peel closure, which a FORK/PSOCK worker
# discards, so at the default `cores` nobody was told.
#
# Now: NULL (the state getsd = FALSE already produces, giving NA bands) and the
# message is collected and re-raised by the caller after dispatch.
#
# Driven by stubbing TMB::sdreport to throw, because reaching this for real
# needs a peel whose Hessian fails exactly here -- which is why the zero-band
# fallback survived this long untested.
testthat::test_that("a peel whose sdreport fails returns no band and reports why", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # BS2017SS, not make_test_data(): that fixture cannot produce a
  # positive-definite Hessian, so under getsd = TRUE every peel is dropped by
  # .refit_converged() before this path is reached and the test would pass
  # vacuously on an empty peel list. See TRAPS.md.
  data(BS2017SS, package = "Rceattle")
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = BS2017SS, file = NULL, estimateMode = 1,
    fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0))))

  # Fail ONLY the report pass. .fit_tmb() also calls sdreport() while fitting
  # each peel, with par.fixed / hessian.fixed and the rest; the report pass
  # calls it bare, with just the object. A blanket stub kills the peel's own fit
  # instead, and the error escapes the tryCatch under test.
  real_sdreport <- TMB::sdreport
  testthat::local_mocked_bindings(
    sdreport = function(obj, ...) {
      if (...length() == 0L) stop("boom")
      real_sdreport(obj, ...)
    }, .package = "TMB")

  w <- character(0)
  r <- withCallingHandlers(
    suppressMessages(retrospective(fit, peels = c(2, 2), cores = 1,
                                   getsd = TRUE)),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })

  # The caller is told, once, naming the peel.
  hit <- grep("could not report hindcast standard errors", w, value = TRUE)
  testthat::expect_length(hit, 1L)
  testthat::expect_match(hit, "Peel 2", fixed = TRUE)
  testthat::expect_match(hit, "no uncertainty band", fixed = TRUE)

  # The band is ABSENT, not zero. A zero-valued sdrep would pass an
  # is.null() check only by accident, so assert both.
  peels <- Filter(function(m) m$data_list$endyr_peel < fit$data_list$endyr,
                  r$Rceattle_list)
  testthat::skip_if(length(peels) == 0L,
                    "the peel did not converge, so the report path is unreached")
  testthat::expect_length(peels, 1L)
  testthat::expect_null(peels[[1]][["sdrep"]])

  # And the plumbing does not ride out on the returned object.
  testthat::expect_null(attr(peels[[1]], "peel_warnings"))

  # Point estimates are unaffected: the fit is still there.
  testthat::expect_false(is.null(peels[[1]]$quantities$ssb))
})
