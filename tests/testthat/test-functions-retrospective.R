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
  # is called twice per peel but a peel that fails the first check returns before
  # the second, so any count-based scheme silently drifts out of step.
  keep_only <- function(keep_endyrs) {
    function(newmod, max_grad = 1, hessian = TRUE) {
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


# Mohn's rho has TWO live conventions in this field, differing by exactly the
# number of peels, and r4ss returns both from one function: `SSmohnsrho()` gives
# `$SSB` as sum(rel_error) and `$AFSC_Hurtado_SSB` as sum(rel_error)/n. SAM's
# `mohn.samset()` uses colMeans, i.e. the mean. Rceattle accumulates and divides
# by N, so it is the MEAN form, and that is what makes `print(retro, band = 0.2)`
# legitimate: the Hurtado-Ferro et al. (2015) ranges the band comes from are on
# the mean. Verified 2026-10-12 against both implementations on the same three
# peels of BS2017SS -- agreement to 5.6e-17 or better, and the sum form exactly
# n x larger on every quantity and species.
#
# Pinned because a switch between the two is invisible in a fit: rho would move
# by a factor of the peel count and still look like a rho. At the default
# peels = 5 a clean 0.04 becomes 0.20, the band edge.
testthat::test_that("Mohn's rho is the mean over peels, not the sum", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  data(BS2017SS, envir = environment())
  d <- BS2017SS
  d$projyr <- 2020
  d$fleet_control$Proj_F_proportion <- rep(1, 7)

  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, file = NULL,
    estimateMode = 1,            # hindcast only; rho reads the terminal year
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  n_ask <- 3
  r <- suppressMessages(suppressWarnings(
    retrospective(fit, peels = n_ask, cores = 1, getsd = FALSE)))

  styr <- fit$data_list$styr
  endyr <- fit$data_list$endyr
  nspp <- fit$data_list$nspp
  full <- r$Rceattle_list[[length(r$Rceattle_list)]]
  peels <- r$Rceattle_list[vapply(
    r$Rceattle_list, function(x) x$data_list$endyr_peel < endyr, logical(1))]
  testthat::expect_equal(length(peels), n_ask)

  m <- as.data.frame(r$mohns)
  spp <- setdiff(names(m), c("Object", "Forecast year", "N"))

  for (q in c("ssb", "biomass", "R")) {
    # Each peel's relative error at its OWN terminal year -- the quantity all
    # three implementations agree on before they aggregate.
    re <- t(vapply(peels, function(p) {
      k <- p$data_list$endyr_peel - styr + 1
      (p$quantities[[q]][, k] - full$quantities[[q]][, k]) /
        full$quantities[[q]][, k]
    }, numeric(nspp)))

    row <- m[m$Object == q & m[["Forecast year"]] == 0, , drop = FALSE]
    rce <- suppressWarnings(as.numeric(unlist(row[, spp, drop = FALSE])))

    testthat::expect_equal(row[["N"]], n_ask)
    testthat::expect_equal(rce, unname(colMeans(re)), tolerance = 1e-10)

    # And it is NOT the sum form. With n peels the two differ by exactly n, so
    # this fails loudly if the convention is ever switched.
    testthat::expect_equal(unname(colSums(re)) / rce,
                           rep(n_ask, length(rce)), tolerance = 1e-8)
  }
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
    # The positive case for the pdHess arm of the usability check below: a
    # healthy rebuild must pass it, or a real band would be declined instead.
    testthat::expect_true(isTRUE(m$sdrep$pdHess))
    # Predicted catch is ADREPORTed only inside the correction, not in what it returns.
    testthat::expect_false("catch_hat" %in% names(m$sdrep$value))
  }
  testthat::expect_false("catch_hat" %in% names(fit$sdrep$value))

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


# A peel's report rebuilt outside retrospective(): the full fit's map with the
# peeled years' deviations off, at the peel's estimates, built without
# optimizing. With `hold_F` the peeled F is fixed at the peel's values and the
# peeled catch zeroed, as retrospective()'s report does before its correction;
# without, F is free and the peeled catch stays in the likelihood, as the report
# before 5.66.3 did. `catch_sd` then sets that catch to its own prediction with
# this sd, so as it shrinks the band converges to the catch-conditioned one.
peel_reference_sdreport <- function(fit, peel, ep, hold_F, catch_sd = NULL) {
  styr <- fit$data_list$styr
  map <- fit$map
  for (nm in c("rec_dev", "log_M1_dev", "index_q_dev", "log_sel_slp_dev",
               "sel_inf_dev", "sel_coff_dev")) {
    x <- map$mapList[[nm]]
    x[slice.index(x, length(dim(x))) > ep - styr + 1] <- NA
    map$mapList[[nm]] <- x
    map$mapFactor[[nm]] <- factor(x)
  }
  dl <- peel$data_list
  dl$endyr <- fit$data_list$endyr
  if (hold_F) {
    peeled <- slice.index(map$mapList$log_F, 2) > ep - styr + 1
    map$mapList$log_F[peeled] <- NA
    map$mapFactor$log_F <- factor(map$mapList$log_F)
  }
  ref <- suppressMessages(suppressWarnings(Rceattle:::.refit_like(
    data_list = dl, inits = peel$estimated_params, map = map,
    estimateMode = 3, getsd = FALSE,
    srr_mse_switchyr = min(dl$srr_mse_switchyr, ep),
    srr_hat_endyr = min(dl$srr_hat_endyr, ep),
    suit_endyr = pmin(dl$suit_endyr, ep))))
  obj <- ref$obj
  peeled <- obj$env$data$catch_ctl[, 3] > ep
  if (hold_F) {
    # Zeroed in the TMB data, not data_list, which would also set that F to 0.
    obj$env$data$catch_obs[peeled, 1] <- 0
    obj$retape()
  } else if (!is.null(catch_sd)) {
    # The likelihood reads the observed log catch from obsvec, centred
    # bias_adjust_obs * sd^2 / 2 below the prediction.
    dat <- obj$env$data
    pos <- dat$catch_obsvec_idx[peeled] + 1
    hat <- obj$report()$catch_hat[peeled]
    obj$env$data$obsvec[pos[pos > 0]] <- log(hat[pos > 0]) -
      dat$bias_adjust_obs * catch_sd^2 / 2
    obj$env$data$catch_obs[peeled, 2] <- catch_sd
    obj$retape()
  }
  TMB::sdreport(obj)
}

# The rebuild used to keep the peeled years' catch in its likelihood with their
# F free. Fitted exactly, that catch still entered the Laplace approximation
# through the hindcast's random effects, so the report sat off the optimum (max
# gradient 0.15-4.19 on parameters other than sigma_R across 10 peels of GOA
# pollock 2026; 0.277 on this fixture's 2-year peel). Holding the peeled F and dropping that catch leaves one term: the
# peel fixes its peeled recruitment deviations yet still scores them, at 0 in
# the peeled hindcast and at the projected deviation in this report, so R_log_sd
# alone keeps a gradient of -sum((dev / sigma_R)^2).
testthat::test_that("a peel's sd rebuild is stationary but for R_log_sd under random recruitment", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  data(BS2017SS, envir = environment())
  d <- BS2017SS
  d$projyr <- 2020
  d$fleet_control$Proj_F_proportion <- rep(1, 7)

  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = TRUE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, getsd = TRUE, verbose = 0))))
  r <- suppressMessages(suppressWarnings(
    retrospective(fit, peels = 2, cores = 1, getsd = TRUE)))

  peel <- r$Rceattle_list[[paste0("Year_", fit$data_list$endyr - 2)]]
  testthat::expect_false(is.null(peel$sdrep))
  grad <- peel$sdrep$gradient.fixed
  is_sigma_r <- names(peel$sdrep$par.fixed) == "R_log_sd"
  testthat::expect_lt(max(abs(grad[!is_sigma_r])), 0.01)
  testthat::expect_true(isTRUE(peel$sdrep$pdHess))

  # The band is centred on the peel's own point estimates, peeled years too.
  nspp <- fit$data_list$nspp
  ny <- fit$data_list$endyr - fit$data_list$styr + 1
  ssb_value <- matrix(peel$sdrep$value[names(peel$sdrep$value) == "ssb"],
                      nrow = nspp)
  testthat::expect_equal(unname(ssb_value[, seq_len(ny)]),
                         unname(peel$quantities$ssb[, seq_len(ny)]))

  # Against the held-F report the correction starts from: fitted-year bands do
  # not depend on the peeled F, so they must be equal, which pins the joint
  # covariance of fixed and random effects; a forecast year's band, which
  # depends on that F, must be widened by the correction.
  ep <- fit$data_list$endyr - 2
  ref_sd <- peel_reference_sdreport(fit, peel, ep, hold_F = TRUE)
  fitted_yrs <- seq_len(ep - fit$data_list$styr + 1)
  forecast_yrs <- (ep - fit$data_list$styr + 2):ny
  ssb_sd <- function(s) matrix(s$sd[names(s$value) == "ssb"], nrow = nspp)
  testthat::expect_equal(ssb_sd(peel$sdrep)[, fitted_yrs],
                         ssb_sd(ref_sd)[, fitted_yrs], tolerance = 1e-6)
  testthat::expect_gt(max(ssb_sd(peel$sdrep)[, forecast_yrs] /
                            ssb_sd(ref_sd)[, forecast_yrs]), 1.01)

  # And by the right amount: the old joint report, with the peeled catch at its
  # own prediction and a small sd, converges to the catch-conditioned band.
  limit_sd <- peel_reference_sdreport(fit, peel, ep, hold_F = FALSE,
                                      catch_sd = 0.01)
  testthat::expect_lt(max(abs(ssb_sd(peel$sdrep)[, forecast_yrs] /
                               ssb_sd(limit_sd)[, forecast_yrs] - 1)), 0.01)
})

# Without random effects the old report -- hindcast and peeled F free, peeled
# catch in the likelihood -- sat at its optimum, its band the catch-conditioned
# one plus the catch's observation error. The correction must reproduce it in
# the fitted years, and in the forecast years fall short only by that error,
# which it leaves out (catch known). Two peels, so a forecast year's SSB
# depends on the peeled F: BS2017SS spawns on 1 January, before that year's
# fishing.
testthat::test_that("a peel's catch-conditioned band matches the old joint report without random effects", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  data(BS2017SS, envir = environment())
  d <- BS2017SS
  d$projyr <- 2020
  d$fleet_control$Proj_F_proportion <- rep(1, 7)
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, getsd = TRUE, verbose = 0))))
  r <- suppressMessages(suppressWarnings(
    retrospective(fit, peels = 2, cores = 1, getsd = TRUE)))

  styr <- fit$data_list$styr; endyr <- fit$data_list$endyr; ep <- endyr - 2
  peel <- r$Rceattle_list[[paste0("Year_", ep)]]
  testthat::expect_false(is.null(peel$sdrep))

  ref_sd <- peel_reference_sdreport(fit, peel, ep, hold_F = FALSE)

  ssb_sd <- function(s) matrix(s$sd[names(s$value) == "ssb"],
                              nrow = fit$data_list$nspp)
  fitted_yrs <- seq_len(ep - styr + 1)
  forecast_yrs <- (ep - styr + 2):(endyr - styr + 1)
  testthat::expect_lt(max(abs(ref_sd$gradient.fixed)), 0.01)
  testthat::expect_equal(ssb_sd(peel$sdrep)[, fitted_yrs],
                         ssb_sd(ref_sd)[, fitted_yrs], tolerance = 1e-6)
  ratio <- ssb_sd(peel$sdrep)[, forecast_yrs] / ssb_sd(ref_sd)[, forecast_yrs]
  testthat::expect_true(all(ratio <= 1 + 1e-6))
  testthat::expect_gt(min(ratio), 0.97)
})

# A peel whose standard-error rebuild fails used to fall back to the forecast
# refit's own sdrep -- the object the rebuild exists to replace, whose hindcast
# standard errors are all exactly 0 (the documented defect behind "a peel
# reports standard errors for its own hindcast", which asserts the complement:
# the REBUILT errors are > 0). The band was
# then drawn as a certain quantity rather than an unmeasured one. The peel now
# keeps what TMB::sdreport() returned, or no sdrep if it threw, and records why
# on attr(peel, "sd_report_error"), which is what stops plots drawing it. It
# raises no warning, because the peel converged and Mohn's rho reads point
# estimates only. A rebuild fails four ways -- it throws, cannot factor its
# Hessian, reports nothing finite, or its correction throws -- each covered below.
testthat::test_that("a peel whose sd rebuild fails keeps its sdreport, and records why", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  data(BS2017SS, envir = environment())
  d <- BS2017SS
  d$projyr <- 2020
  d$fleet_control$Proj_F_proportion <- rep(1, 7)

  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, getsd = TRUE, verbose = 0))))

  # Fail only the report build. It is the one .refit_like() call made at
  # estimateMode 3, so the two fitting refits are untouched and the peel still
  # converges -- which is the point: a missing band is not a failed peel.
  real_refit <- Rceattle:::.refit_like
  # Runs retrospective() with the report build forced to throw, collecting
  # every warning rather than letting one fail the test.
  run_throwing_rebuild <- function(peels, cores) {
    res <- NULL
    warns <- character()
    local({
      testthat::local_mocked_bindings(
        .refit_like = function(...) {
          a <- list(...)
          if (identical(a$estimateMode, 3)) stop("forced report-build failure")
          do.call(real_refit, a)
        }, .package = "Rceattle")
      withCallingHandlers(
        res <<- suppressMessages(
          retrospective(fit, peels = peels, cores = cores, getsd = TRUE)),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
    })
    list(fit = res, warns = warns)
  }

  got <- run_throwing_rebuild(peels = 1, cores = 1)
  r <- got$fit
  testthat::expect_false(any(grepl("standard errors", got$warns)))

  full_endyr <- fit$data_list$endyr
  is_refit <- vapply(r$Rceattle_list,
                     function(x) x$data_list$endyr_peel < full_endyr, logical(1))
  testthat::expect_true(any(is_refit))

  # Absent, not zero, and the reason is on the peel.
  for (m in r$Rceattle_list[is_refit]) {
    testthat::expect_null(m$sdrep)
    testthat::expect_match(attr(m, "sd_report_error", exact = TRUE),
                           "forced report-build failure")
  }

  # The unpeeled model is passed through untouched and keeps its own sdrep, so
  # the change cannot be mistaken for "getsd stopped working".
  testthat::expect_false(is.null(r$Rceattle_list[[length(r$Rceattle_list)]]$sdrep))

  # Nothing about the fit moved: a missing band costs no point estimate.
  r_no <- suppressMessages(suppressWarnings(
    retrospective(fit, peels = 1, cores = 1, getsd = FALSE)))
  shared <- intersect(names(r$Rceattle_list), names(r_no$Rceattle_list))
  testthat::expect_gt(length(shared), 0)
  for (nm in shared) {
    testthat::expect_equal(r$Rceattle_list[[nm]]$quantities$ssb,
                           r_no$Rceattle_list[[nm]]$quantities$ssb)
    # The objective too: this is what stops the report pass being "improved"
    # into a refit.
    testthat::expect_equal(r$Rceattle_list[[nm]]$opt$objective,
                           r_no$Rceattle_list[[nm]]$opt$objective)
  }

  # The run above is sequential (`peels = 1`). Two peels on two cores run the
  # rebuild inside .parallel_lapply() workers, and the reason must still come
  # back on each peel.
  #
  # FORK only. A PSOCK worker runs `library(Rceattle)` into a pristine namespace,
  # so the mock does not reach it and the rebuild succeeds.
  testthat::skip_on_os("windows")
  got_p <- run_throwing_rebuild(peels = 2, cores = 2)
  rp <- got_p$fit
  testthat::expect_false(any(grepl("standard errors", got_p$warns)))

  refits_p <- vapply(rp$Rceattle_list,
                     function(x) x$data_list$endyr_peel < full_endyr, logical(1))
  testthat::expect_true(any(refits_p))
  for (m in rp$Rceattle_list[refits_p]) {
    testthat::expect_null(m$sdrep)
    testthat::expect_match(attr(m, "sd_report_error", exact = TRUE),
                           "forced report-build failure")
  }

  # The two non-throwing failures. One stub per arm, each asserting its OWN
  # message: a single stub setting both conditions would stay green if either
  # half of the predicate were deleted. The stub is reached only through the
  # mocked rebuild, so the fitting refits get the real sdreport and the peel
  # still converges.
  real_sdreport <- TMB::sdreport
  # The rebuild zeroes the peeled catch in obj$env$data and retapes before
  # calling sdreport, so the stub carries an empty catch table and a no-op retape.
  stub <- structure(list(
    env = list2env(list(data = list(catch_obs = matrix(numeric(0), 0, 2),
                                    catch_ctl = matrix(integer(0), 0, 3)))),
    retape = function(...) invisible(NULL)), class = "rce_stub_obj")

  run_bad_rebuild <- function(fake) {
    res <- NULL
    warns <- character()
    local({
      testthat::local_mocked_bindings(
        .refit_like = function(...) {
          a <- list(...)
          if (identical(a[["estimateMode"]], 3)) return(list(obj = stub))
          do.call(real_refit, a)
        }, .package = "Rceattle")
      testthat::local_mocked_bindings(
        sdreport = function(obj, ...) {
          if (inherits(obj, "rce_stub_obj")) return(fake)
          real_sdreport(obj, ...)
        }, .package = "TMB")
      withCallingHandlers(
        res <<- suppressMessages(
          retrospective(fit, peels = 1, cores = 1, getsd = TRUE)),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
    })
    list(fit = res, warns = warns)
  }

  check_flagged <- function(got, fake, pattern) {
    testthat::expect_false(any(grepl("standard errors", got$warns)))
    refits_b <- vapply(got$fit$Rceattle_list,
                       function(x) x$data_list$endyr_peel < full_endyr,
                       logical(1))
    # Kept, not dropped: recording the failure on .conv_hindcast$sd_present
    # instead would have made .refit_converged() drop the peel, changing the
    # number of peels Mohn's rho averages over.
    testthat::expect_equal(sum(refits_b), 1L)
    for (m in got$fit$Rceattle_list[refits_b]) {
      testthat::expect_identical(m$sdrep, fake)   # returned as TMB gave it
      testthat::expect_match(attr(m, "sd_report_error", exact = TRUE), pattern)
    }
  }

  # Arm 1: an unfactorable Hessian, but the report itself is finite.
  fake_indef <- structure(list(pdHess = FALSE, sd = 1, value = c(ssb = 0)),
                          class = "sdreport")
  check_flagged(run_bad_rebuild(fake_indef), fake_indef, "not positive definite")

  # Arm 2: a factorable Hessian that still reports nothing finite.
  fake_nan <- structure(list(pdHess = TRUE, sd = NA_real_, value = c(ssb = 0)),
                        class = "sdreport")
  check_flagged(run_bad_rebuild(fake_nan), fake_nan, "no finite standard error")

  # Arm 3: a real sdreport whose correction throws, so forecast-year bands are
  # left at the held-F report and the reason says so.
  local({
    testthat::local_mocked_bindings(
      .rce_catch_conditioned_sd = function(...) stop("forced correction failure"),
      .package = "Rceattle")
    r_c <- suppressMessages(suppressWarnings(
      retrospective(fit, peels = 1, cores = 1, getsd = TRUE)))
    m <- r_c$Rceattle_list[[paste0("Year_", full_endyr - 1)]]
    testthat::expect_true(isTRUE(m$sdrep$pdHess))
    testthat::expect_match(attr(m, "sd_report_error", exact = TRUE),
                           "catch-conditioned correction failed: forced correction failure")
  })
})

# Only a hindcast's non-convergence warns. A forecast refit fits the peeled F to
# catch with the hindcast held fixed, so its failure leaves the terminal-year rho
# standing but its forecast years are not a forecast of that catch: the peel stays,
# silently, and leaves the forecast-year rows of $mohns, which would otherwise
# average in an F that never fitted the catch. Its Hessian covers F alone, so that
# alone neither drops the peel nor takes it out of those rows.
testthat::test_that("only a peel's hindcast non-convergence warns or drops it", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  data(BS2017SS, envir = environment())
  d <- BS2017SS
  d$projyr <- 2020
  d$fleet_control$Proj_F_proportion <- rep(1, 7)
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  endyr <- fit$data_list$endyr
  ep <- endyr - 2

  # Mark one refit of the 2-year peel as unconverged. Its catch tells them apart:
  # the hindcast refit has peeled catch at 0, the forecast refit has it back. The
  # 1-year peel's forecast refit gets only an indefinite Hessian.
  real_refit <- Rceattle:::.refit_like
  run <- function(fail, getsd) {
    warns <- character()
    res <- NULL
    local({
      testthat::local_mocked_bindings(
        .refit_like = function(...) {
          a <- list(...)
          m <- do.call(real_refit, a)
          if (identical(a[["estimateMode"]], 3)) return(m)
          cd <- a$data_list$catch_data
          is_forecast <- any(cd$Catch[cd$Year > a$data_list$endyr_peel] > 0)
          if (identical(a$data_list$endyr_peel, ep) &&
              identical(fail, if (is_forecast) "forecast" else "hindcast")) {
            m$.conv_hindcast$max_gradient <- 5
          } else if (is_forecast && identical(fail, "forecast")) {
            m$.conv_hindcast$pdHess <- FALSE
          }
          m
        }, .package = "Rceattle")
      withCallingHandlers(
        res <<- suppressMessages(
          retrospective(fit, peels = 2, cores = 1, getsd = getsd)),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
    })
    list(fit = res, warns = warns)
  }
  rows <- function(r, yr) r$mohns[r$mohns$Object == "ssb" &
                                    r$mohns[["Forecast year"]] == yr, ]

  ok <- run("none", getsd = FALSE)$fit

  # Forecast refit failed: kept, no warning, out of the forecast-year rows only.
  got <- run("forecast", getsd = TRUE)
  testthat::expect_false(any(grepl("dropped|converge", got$warns)))
  testthat::expect_true(paste0("Year_", ep) %in% names(got$fit$Rceattle_list))
  r <- got$fit
  testthat::expect_match(
    attr(r$Rceattle_list[[paste0("Year_", ep)]], "forecast_error", exact = TRUE),
    "did not fit F to the peeled catch")
  testthat::expect_null(
    attr(r$Rceattle_list[[paste0("Year_", endyr - 1)]], "forecast_error", exact = TRUE))

  testthat::expect_equal(rows(r, 0)$N, 2)
  sp <- fit$data_list$spnames
  testthat::expect_equal(rows(r, 0)[, sp], rows(ok, 0)[, sp], tolerance = 1e-6)
  # Year 1 is the 1-year peel's alone; year 2 only the failed peel reached.
  testthat::expect_equal(rows(r, 1)$N, 1)
  one <- ok$Rceattle_list[[paste0("Year_", endyr - 1)]]$quantities$ssb
  col <- endyr - fit$data_list$styr + 1
  testthat::expect_equal(
    unlist(rows(r, 1)[, sp], use.names = FALSE),
    as.numeric(one[, col] / fit$quantities$ssb[, col] - 1), tolerance = 1e-6)
  testthat::expect_equal(rows(r, 2)$N, 0)
  testthat::expect_true(all(is.nan(unlist(rows(r, 2)[, sp]))))
  testthat::expect_true(any(grepl(
    "1 peel\\(s\\) left out of the forecast-year rows",
    utils::capture.output(print(r)))))
  testthat::expect_false(any(grepl(
    "left out of the forecast-year rows", utils::capture.output(print(ok)))))

  # Hindcast failed: dropped, with the warning.
  got_h <- run("hindcast", getsd = FALSE)
  testthat::expect_true(any(grepl("1 of 2 peel dropped", got_h$warns)))
  testthat::expect_null(
    attr(got_h$fit$Rceattle_list[[paste0("Year_", endyr - 1)]], "forecast_error", exact = TRUE))
  testthat::expect_false(paste0("Year_", ep) %in% names(got_h$fit$Rceattle_list))
})


# A penalty-form peel that ends before srr_hat_styr has no curve fitted to
# recruitment. retrospective() warns once, before the peels run (warnings from peels
# run in parallel are lost), and each such peel averages over its own years.
testthat::test_that("retrospective() warns once about penalty-form peels with no penalty years", {
  testthat::skip_on_cran()
  d <- make_test_data()
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, file = NULL, estimateMode = 0,
    recFun = build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt",
                       srr_hat_styr = d$endyr),
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  w <- character(0)
  r <- withCallingHandlers(
    suppressMessages(retrospective(fit, peels = 2, cores = 1, getsd = FALSE)),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  hit <- grep("contain no stock-recruit penalty years", w, value = TRUE)
  testthat::expect_length(hit, 1L)
  testthat::expect_match(hit, paste0("ending ", d$endyr - 1, ", ", d$endyr - 2), fixed = TRUE)
  # Both peels survive, each projecting the mean ratio to the curve over its own
  # years (styr skipped); peels are picked by terminal year, as the list is by year.
  testthat::expect_length(r$Rceattle_list, 3L)
  peel_mods <- Filter(function(m) m$data_list$endyr_peel < d$endyr, r$Rceattle_list)
  testthat::expect_length(peel_mods, 2L)
  for (m in peel_mods) {
    n  <- m$data_list$endyr_peel - m$data_list$styr + 1
    rd <- m$estimated_params$rec_dev[1, n + 1]
    testthat::expect_true(is.finite(rd))
    testthat::expect_equal(rd, log(mean((m$quantities$R / m$quantities$R_hat)[1, 2:n])),
                           tolerance = 1e-6)
  }
})


# A multispecies hindcast curve has no unfished R0, so each peel projects
# log(mean(R / R_hat)) over its own years. Alpha and beta are held by linkages:
# estimated on this fixture they run toward overflow (TODO-srr-multispecies.md).
testthat::test_that("retrospective() scales a multispecies hindcast curve by its mean ratio to R", {
  testthat::skip_on_cran()
  set.seed(123)
  d  <- make_msm_test_data()$data_list
  fc <- fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  msm_fit <- function(recFun, mode) suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = mode, msmMode = 1,
    suitMode = 0, initMode = "NonEquilibrium", random_rec = FALSE,
    recFun = recFun, fit_control = fc)))

  # A Beverton-Holt curve through each species' mean (SSB, R): alpha = 2 Rbar / Sbar, beta = 1 / Sbar.
  mr   <- msm_fit(build_srr(), 1)
  nh   <- d$endyr - d$styr + 1
  Rbar <- rowMeans(mr$quantities$R[, 1:nh, drop = FALSE])
  Sbar <- rowMeans(mr$quantities$ssb[, 1:nh, drop = FALSE])
  fixed <- function(v) lapply(seq_along(v), function(sp) linkage_spec(
    ~ 1, species = sp, est_phase = 0, init = list(`(Intercept)` = unname(v[sp]))))
  fit <- msm_fit(build_srr(srr_fun = "BevertonHolt", linkages = list(
    alpha = fixed(2 * Rbar / Sbar), beta = fixed(1 / Sbar))), 0)

  r <- suppressWarnings(suppressMessages(retrospective(fit, peels = 1, cores = 1, getsd = FALSE)))
  peel <- Filter(function(m) m$data_list$endyr_peel < d$endyr, r$Rceattle_list)
  testthat::expect_length(peel, 1L)
  m  <- peel[[1]]
  n  <- m$data_list$endyr_peel - m$data_list$styr + 1
  rd <- m$estimated_params$rec_dev[, n + 1]
  testthat::expect_true(all(is.finite(rd)))
  testthat::expect_equal(unname(rd),
                         unname(log(rowMeans((m$quantities$R / m$quantities$R_hat)[, 1:n, drop = FALSE]))),
                         tolerance = 1e-6)
})
