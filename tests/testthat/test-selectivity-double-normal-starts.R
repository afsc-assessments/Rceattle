# DoubleNormal (type 8) reads the two shared `sel_inf` slots as a PEAK and a
# right-tail floor, where the logistic family reads them as two inflections.
# The shared defaults are wrong for that reading: the peak started at 0, below
# the first age, and the floor at 10 on the LOGIT scale -- plogis(10) = 0.99996
# -- so the starting curve was flat at ~1 for every age. A flat curve leaves the
# ascending width with no gradient, so the optimizer had no descent direction
# and stayed on that ridge.
#
# Measured on the GOApollock fishery (fleet 8), static selectivity, phased,
# `inits = NULL`:
#
#   starts                      sel-at-age                  objective    AIC
#   peak 0, floor logit 10      constant 0.999996           3085.9797   6611.96
#   peak 5.5, floor logit 0     0.0024 .. 0.9815 .. 0.4407   914.1043   2268.21
#   (DoubleLogistic, same defaults, for scale)               918.1536   2276.31
#
# 2171.88 nats, and the fixed form now edges out the double logistic by 4.05.
# 914.1043 is the same optimum a hand-tuned `inits` reached, so the derived
# starts are not merely better -- they find what someone had to set by hand.
#
# No bundled dataset sets `Selectivity = "DoubleNormal"` -- all 12 use codes
# {0,1,2,3,4,6} -- so none of the four golden references reaches this, and
# `build_params()` output is bit-identical on all four. Four scripts under
# `Rceattle-models/GOA cod/Bridging/` DO set it, but inject their own `sel_inf`
# immediately after, and the live bridges have since moved to
# `DoubleNormalSS3`; so they are inert, for that reason rather than for absence.
# `test-selectivity-double-normal.R` never saw this because it supplies its own
# starting values.
#
# The peak is derived, not chosen: an age-based curve is evaluated at `bin + 1`
# (`selectivity.hpp`), so the x-axis is the 1-based BIN ORDINAL 1..nages and the
# mid-range is (nages + 1) / 2 -- the same expression the `DoubleNormalSS3`
# block uses for the same slot, and the bin-scale analogue of the
# `length_midpoint()` the length branch applies.
#
# The floor change is NOT what escapes the ridge: with the peak at mid-range the
# fit reaches 914.1043 from a floor start of either 0 or 10, and with the peak at
# 0 a floor of 0 still fails (2973.69, peak running to -859). The peak carries
# the fix; logit 0 is kept because it is the logistic's maximum-gradient point.
#
# The fitted floor is NOT identified: it lands near -17.8 with a standard error
# in the thousands, and anywhere in [-21.6, -10.4] across multi-starts with the
# objective unchanged to 4 dp. That is disclosed in
# `vignettes/articles/adding-a-selectivity-form.Rmd` and is a property of the
# form, not of these starts -- but it is the default outcome now, so it is
# named here too.

.dn_data <- function(dimension = NULL) {
  data("GOApollock", package = "Rceattle", envir = environment())
  d <- GOApollock
  d$fleet_control$Selectivity[8]      <- "DoubleNormal"
  d$fleet_control$Time_varying_sel[8] <- "Off"
  if (!is.null(dimension)) d$fleet_control$Selectivity_dimension[8] <- dimension
  d
}

.dn_pars <- function(d) {
  dd <- suppressMessages(Rceattle::switch_check(d))
  list(pars = suppressWarnings(suppressMessages(Rceattle::build_params(dd))),
       data = dd)
}


testthat::test_that("an age-based DoubleNormal starts mid-range with a low floor", {
  got <- .dn_pars(.dn_data())
  p <- got$pars
  d <- got$data
  # A LITERAL, not a recomputation of the implementation's expression. An
  # earlier version of this block recomputed `minage + (nages - 1) / 2`, which
  # asserted self-consistency rather than correctness and so could not see that
  # the formula was on the wrong scale.
  #
  # An age-based curve is evaluated at `bin + 1` (selectivity.hpp), so the
  # x-axis is the 1-based BIN ORDINAL 1..nages and the mid-range is
  # (nages + 1) / 2. GOApollock has nages = 10, hence 5.5.
  testthat::expect_equal(unname(d$nages[1]), 10L)
  testthat::expect_equal(unname(p$sel_inf[1, 8, 1]), 5.5)
  # The right-tail floor is a LOGIT: 0 is a floor of 0.5, where the old default
  # of 10 was a floor of 0.99996 and flattened the whole curve.
  testthat::expect_equal(unname(p$sel_inf[2, 8, 1]), 0)
  testthat::expect_lt(stats::plogis(p$sel_inf[2, 8, 1]), 0.9)
})


testthat::test_that("the peak is a bin ordinal, not an absolute age", {
  # The two conventions coincide only at minage = 1, which every bundled dataset
  # has -- so a formula written on the age scale measures correctly on all of
  # them and silently wrong on a minage = 3 stock. Drive minage directly.
  #
  # `Bin_first_selected` is a 1-based bin ordinal while `Sel_norm_bin` is an
  # absolute age (CLAUDE.md rule 10); this slot follows the former, because
  # selectivity.hpp evaluates the curve at `bin + 1`.
  for (ma in c(1L, 3L)) {
    d <- .dn_data()
    d$minage <- rep(ma, d$nspp)
    p <- suppressWarnings(suppressMessages(Rceattle::build_params(
      suppressMessages(Rceattle::switch_check(d)))))
    # Unchanged by minage: 10 bins always put the mid-range at 5.5.
    testthat::expect_equal(unname(p$sel_inf[1, 8, 1]), 5.5,
                           info = paste("minage", ma))
  }
})


testthat::test_that("the other fleets' starts are untouched", {
  # The block keys on Selectivity == DoubleNormal, so every other form must keep
  # the shared defaults. This is the assertion that would catch a `which()`
  # matching too widely, which is the way a start change leaks into a fit that
  # never asked for it.
  got <- .dn_pars(.dn_data())
  p <- got$pars
  others <- setdiff(seq_len(nrow(got$data$fleet_control)), 8L)
  testthat::expect_true(all(p$sel_inf[2, others, 1] == 10))
  # Fleet 1 is DescendingLogistic and age-based: 0 is an inflection there, and
  # the curve still has a gradient, so it is deliberately left alone.
  testthat::expect_equal(unname(p$sel_inf[1, 1, 1]), 0)
})


testthat::test_that("a length-based DoubleNormal keeps its length midpoint", {
  # The length branch runs first and sets sel_inf[1] to a length. Clobbering it
  # with an AGE midpoint would be the same defect in the other direction, so the
  # DoubleNormal block sets the peak only on age-based fleets.
  age_p <- .dn_pars(.dn_data())$pars
  len_p <- .dn_pars(.dn_data(dimension = "Length"))$pars
  testthat::expect_false(isTRUE(all.equal(unname(len_p$sel_inf[1, 8, 1]),
                                          unname(age_p$sel_inf[1, 8, 1]))))
  # Still a sensible length, and still with the low floor.
  testthat::expect_gt(len_p$sel_inf[1, 8, 1], 0)
  testthat::expect_equal(unname(len_p$sel_inf[2, 8, 1]), 0)
})


testthat::test_that("the fitted curve is a dome rather than flat at one", {
  testthat::skip_on_cran()   # a real phased fit_mod() optimization
  m <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = .dn_data(), inits = NULL, estimateMode = "Hindcast",
    msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = TRUE, verbose = 0,
                                        getsd = FALSE))))
  s <- m$quantities$sel_at_age
  # Assert the SHAPE before reading the value: `sel` is not a reported quantity
  # name, and all() over a zero-length vector is TRUE, so an empty slice would
  # otherwise satisfy every check below.
  testthat::expect_equal(dim(s)[1], 8L)
  v <- s[8, 1, , dim(s)[4]]
  testthat::expect_length(v, 10L)

  # The defect: a constant curve. 0.999996 at every age, min == max.
  testthat::expect_gt(max(v) - min(v), 0.5)
  testthat::expect_lt(min(v), 0.1)
  # And it reaches the optimum a hand-tuned inits reached, not the flat ridge.
  testthat::expect_lt(m$opt$objective, 1000)
  testthat::expect_equal(m$opt$objective, 914.1043, tolerance = 1e-3)
})
