# The 2D-AR1 SELECTIVITY density, pinned against an independent multivariate
# normal -- the companion to test-mortality-2dar1-orientation.R.
#
# SEPARABLE(f, g) applies f to the outermost array dimension. The selectivity
# field is (bin, year), so the YEAR correlation is passed first and the BIN
# correlation second. Swapped, the two exchange silently: the fit stays
# self-consistent and nothing errors, but Sel_curve_pen1 is reported as a
# correlation across bins while it acted across years. The field here is 8 bins
# by 42 years, so the two orientations are not even the same shape.
#
# calc_nll_ar1_2d() (helpers.R) builds kronecker(R_year, R_bin) and is the
# outside check; its `rho_a` argument is the FAST dimension, which is the bin.

testthat::test_that("the 2D-AR1 selectivity density correlates bins by Sel_curve_pen1 and years by Sel_curve_pen2", {
  # No skip_on_cran(): estimateMode = 3 builds without optimizing.
  testthat::skip_if_not_installed("mvtnorm")

  FLT <- 2L
  N_BINS <- 8L
  SIG <- 0.4; RHO_BIN <- 0.15; RHO_YEAR <- 0.81   # deliberately unequal

  data("GOA2018SS")
  d <- GOA2018SS
  # One fleet on 2DAR1, with its own selectivity block, so jnll_comp's column
  # for that fleet holds this density and nothing else.
  d$fleet_control$Selectivity_index <- seq_len(nrow(d$fleet_control))
  d$fleet_control$Selectivity[FLT] <- "2DAR1"
  d$fleet_control$Time_varying_sel_sd[FLT] <- 1
  d$fleet_control$Age_first_selected[FLT] <- 1
  d$fleet_control$N_sel_bins[FLT] <- N_BINS
  d$fleet_control$Sel_curve_pen1[FLT] <- 0
  d$fleet_control$Sel_curve_pen2[FLT] <- 0
  d$fleet_control$Sel_norm_bin[FLT] <- NA          # do not normalize
  d$catch_data$Catch <- 1e6                        # zero catch turns sel devs off

  fit <- suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3,
    random_rec = FALSE, random_sel = TRUE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)))

  nyrs <- d$endyr - d$styr + 1L

  obj <- fit$obj
  par <- obj$env$last.par.best
  if (is.null(par)) { obj$fn(obj$par); par <- obj$env$last.par.best }

  # Write a whole block back through its map. `par` holds one entry per free
  # parameter, ordered by map LEVEL -- and the codes in mapList are not
  # contiguous (build_map() numbers them across all fleets, so a block that only
  # some fleets estimate has gaps). Take the sorted unique codes, not 1:n.
  set_block <- function(par, name, arr) {
    pos <- which(names(par) == name)
    mp  <- as.integer(fit$map$mapList[[name]])
    lv  <- sort(unique(mp[!is.na(mp)]))
    testthat::expect_length(pos, length(lv))
    par[pos] <- vapply(seq_along(pos), function(k) arr[which(mp == lv[k])[1]], numeric(1))
    par
  }

  # sel_curve_pen is [n_flt, 3]: slot 1 across bins, slot 2 across years.
  # rho_trans() is tanh(), so the parameter is atanh(rho).
  pen <- fit$estimated_params$sel_curve_pen
  pen[FLT, 1] <- atanh(RHO_BIN)
  pen[FLT, 2] <- atanh(RHO_YEAR)
  par <- set_block(par, "sel_curve_pen", pen)

  sds <- fit$estimated_params$sel_dev_log_sd
  sds[FLT] <- log(SIG)
  par <- set_block(par, "sel_dev_log_sd", sds)

  # A known deviation field, written straight into the parameter vector and read
  # back with obj$report(), which evaluates at exactly these values -- obj$fn()
  # would run the Laplace inner problem and move the random effects.
  set.seed(23)
  field <- matrix(stats::rnorm(N_BINS * nyrs), N_BINS, nyrs)

  arr <- array(0, dim = dim(fit$estimated_params$sel_coff_dev))
  arr[FLT, 1, seq_len(N_BINS), seq_len(nyrs)] <- field
  par <- set_block(par, "sel_coff_dev", arr)

  JNLL_SEL_DEV <- 5L
  got <- obj$report(par)$jnll_comp[JNLL_SEL_DEV + 1L, FLT]

  want  <- calc_nll_ar1_2d(field, SIG, rho_a = RHO_BIN, rho_y = RHO_YEAR)
  wrong <- calc_nll_ar1_2d(field, SIG, rho_a = RHO_YEAR, rho_y = RHO_BIN)

  testthat::expect_equal(got, want, tolerance = 1e-8)
  # The check only means something if the two orientations actually differ.
  testthat::expect_gt(abs(want - wrong), 1)
})
