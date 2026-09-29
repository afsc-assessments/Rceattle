# DoubleNormalSS3 (type 15): Stock Synthesis size-selectivity pattern 24, with
# its six parameters in their own array (sel_dn6) and time variation only
# through selectivity linkages. The R implementation below follows SS3 3.30's
# SS_selex.tpl (pattern 24, Apical_Selex = 1): peak P1, logit top width P2, log
# ascending / descending widths P3 / P4, logit initial / final selectivity
# P5 / P6, with -999 on an end switching its scaling off. The curve is not
# normalized.

testthat::skip_on_cran()

ss3_pattern24 <- function(x, P, init_on = TRUE, final_on = TRUE, w = diff(x)[1]) {
  n <- length(x)
  peak2  <- P[1] + w + (0.99 * x[n] - P[1] - w) / (1 + exp(-P[2]))
  up     <- exp(P[3]); dn <- exp(P[4])
  point1 <- 1 / (1 + exp(-P[5])); point2 <- 1 / (1 + exp(-P[6]))
  t1min  <- exp(-(x[1] - P[1])^2 / up)
  t2min  <- exp(-(x[n] - peak2)^2 / dn)
  t1 <- x - P[1]; t2 <- x - peak2
  j1 <- 1 / (1 + exp(-(20 * t1 / (1 + abs(t1)))))
  j2 <- 1 / (1 + exp(-(20 * t2 / (1 + abs(t2)))))
  asc <- if (init_on) point1 + (1 - point1) * (exp(-t1^2 / up) - t1min) / (1 - t1min) else exp(-t1^2 / up)
  dsc <- if (final_on) 1 + (point2 - 1) * (exp(-t2^2 / dn) - 1) / (t2min - 1) else exp(-t2^2 / dn)
  asc * (1 - j1) + j1 * ((1 - j2) + dsc * j2)
}

dn6_build <- function(d, inits = NULL, selFun = build_selectivity()) {
  suppressWarnings(suppressMessages(fit_mod(
    data_list = d, inits = inits, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    growthFun = build_growth(fun = "vonBertalanffy"), selFun = selFun,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
}

set.seed(21)
d <- make_msm_test_data()$data_list
fl <- which(d$fleet_control$Species == 1 & d$fleet_control$Fleet_type %in% c(1, "Fishery"))[1]
d$fleet_control$Selectivity[fl]           <- "DoubleNormalSS3"
d$fleet_control$Selectivity_dimension[fl] <- "Length"
d$fleet_control$Time_varying_sel[fl]      <- "Off"
edges <- sort(unique(d$caal_data$Length[d$caal_data$Species == 1]))
mids  <- edges + diff(edges)[1] / 2

testthat::test_that("the curve is SS3 pattern 24, unnormalized, with both ends scaled", {
  m0 <- dn6_build(d)
  P  <- c(55, -3, 4.5, 5.5, -2, 0.5)
  ip <- m0$estimated_params
  ip$sel_dn6[, fl, ] <- P
  m  <- dn6_build(d, inits = ip)
  got <- m$quantities$sel_at_length[fl, 1, seq_along(edges), 1]
  testthat::expect_equal(unname(got), ss3_pattern24(mids, P), tolerance = 1e-10)
  # Age selectivity is the length curve through the fleet's age-length key
  gm <- m$quantities$growth_matrix[2 * d$nspp + fl, 1, , seq_along(edges), 1]
  testthat::expect_equal(unname(m$quantities$sel_at_age[fl, 1, , 1]),
                         unname(as.numeric(gm %*% got)), tolerance = 1e-10)
})

testthat::test_that("an end at SS3's -999 is unscaled and its parameter fixed", {
  m0 <- dn6_build(d)
  # The default start leaves both ends at -999
  testthat::expect_equal(unname(m0$estimated_params$sel_dn6[5:6, fl, 1]), c(-999, -999))
  testthat::expect_true(all(is.na(m0$map$mapList$sel_dn6[5:6, fl, ])))
  testthat::expect_false(any(is.na(m0$map$mapList$sel_dn6[1:4, fl, ])))
  P  <- c(60, -4, 5, 5, -999, -999)
  ip <- m0$estimated_params
  ip$sel_dn6[, fl, ] <- P
  m  <- dn6_build(d, inits = ip)
  testthat::expect_equal(unname(m$quantities$sel_at_length[fl, 1, seq_along(edges), 1]),
                         ss3_pattern24(mids, P, init_on = FALSE, final_on = FALSE),
                         tolerance = 1e-10)
})

# The configuration the AI cod bridge actually uses: control.ss gives its fishery
# Size_DblN_start_logit = -999 (unscaled) and Size_DblN_end_logit = 4 (scaled).
# Only both-on and both-off were covered, so the mixed case the bridge relies on
# went unexercised.
testthat::test_that("one end scaled and one at -999 is SS3 pattern 24", {
  m0 <- dn6_build(d)
  P  <- c(55, -3, 4.5, 5.5, -999, 0.5)
  ip <- m0$estimated_params
  ip$sel_dn6[, fl, ] <- P
  m  <- dn6_build(d, inits = ip)

  # The unscaled end's parameter is fixed, the scaled one stays estimated.
  testthat::expect_true(is.na(m$map$mapList$sel_dn6[5, fl, 1]))
  testthat::expect_false(is.na(m$map$mapList$sel_dn6[6, fl, 1]))
  testthat::expect_equal(unname(m$quantities$sel_at_length[fl, 1, seq_along(edges), 1]),
                         ss3_pattern24(mids, P, init_on = FALSE, final_on = TRUE),
                         tolerance = 1e-10)
})

# One end flag serves both sexes, so a two-sex species whose sexes disagree would
# have had one sex fitted with the other's curve shape -- its own end parameter
# estimated, entering nothing, at zero gradient.
testthat::test_that("the sexes must agree on whether an end is scaled", {
  a <- array(0, dim = c(6, 2, 2))
  a[5:6, 1, ] <- -999                      # both ends off, both sexes
  a[5, 2, ] <- -999; a[6, 2, ] <- -3       # start off, end on: allowed
  nm <- c("A", "B")

  # Mixed ENDS is SS3's own configuration and stays allowed.
  testthat::expect_equal(
    Rceattle:::.rce_dn6_ends(a, c(2, 2), c(1, 1), nm),
    matrix(c(0L, 0L, 0L, 1L), ncol = 2, byrow = TRUE))

  # Mixed SEXES on the same end is refused.
  b <- a; b[5, 2, 2] <- -4
  testthat::expect_error(Rceattle:::.rce_dn6_ends(b, c(2, 2), c(1, 1), nm),
                         "sexes disagree")
  # A one-sex species has no second sex to disagree, so it is unaffected.
  testthat::expect_equal(
    Rceattle:::.rce_dn6_ends(b, c(1, 1), c(1, 1), nm),
    matrix(c(0L, 0L, 0L, 1L), ncol = 2, byrow = TRUE))
})

# Fleets sharing a Selectivity_index estimate one sel_dn6 block, but the end flag
# is per fleet and read off each fleet's own start. adjust_map_shared_params()
# shares the map and the -999 pass then NAs the follower's own cell, so the two
# would share a parameter and still get different curves.
testthat::test_that("a shared Selectivity_index must agree on its ends", {
  a <- array(0, dim = c(6, 2, 1))
  a[5:6, 1, ] <- -999                     # fleet A: both ends off
  a[5, 2, ] <- -999; a[6, 2, ] <- -3      # fleet B: end on
  nm <- c("A", "B")
  both <- c(TRUE, TRUE)

  testthat::expect_error(
    Rceattle:::.rce_dn6_ends(a, 1, c(1, 1), nm, is_dn6 = both,
                             sel_index = c(1L, 1L)),
    "share Selectivity_index")

  # Separate blocks estimate separate parameters, so they may differ.
  testthat::expect_equal(
    Rceattle:::.rce_dn6_ends(a, 1, c(1, 1), nm, is_dn6 = both,
                             sel_index = c(1L, 2L)),
    matrix(c(0L, 0L, 0L, 1L), ncol = 2, byrow = TRUE))

  # A follower on another form never reads sel_dn6, so it is not compared.
  testthat::expect_equal(
    Rceattle:::.rce_dn6_ends(a, 1, c(1, 1), nm, is_dn6 = c(TRUE, FALSE),
                             sel_index = c(1L, 1L)),
    matrix(c(0L, 0L, 0L, 1L), ncol = 2, byrow = TRUE))
})

# sel_dn6 is the only array case 15 reads, and the other forms never read it, so a
# linkage named the wrong way round is estimated and changes nothing. `peak`
# (sel_inf) against `dn_peak` (sel_dn6) is the pair to watch: both the Doxygen and
# parameter_dictionary() call form 15's P1 "peak".
testthat::test_that("a linkage parameter must match the fleet's form", {
  fc <- data.frame(Fleet_name = c("Fsh15", "Srv1"),
                   Selectivity = c("DoubleNormalSS3", "Logistic"),
                   stringsAsFactors = FALSE)
  row <- function(param, fleet) {
    data.frame(process = "sel", param = param, fleet = fleet,
               stringsAsFactors = FALSE)
  }
  chk <- Rceattle:::.check_sel_linkage_support

  # Each name on the form that reads it.
  testthat::expect_silent(chk(row("dn_peak", 1L), fc))
  testthat::expect_silent(chk(row("peak", 2L), fc))
  # And each on the form that does not.
  testthat::expect_error(chk(row("peak", 1L), fc), "reads only its own six")
  testthat::expect_error(chk(row("dn_peak", 2L), fc), "only that form reads it")
  # An unstratified row targets every fleet, as the form check above treats it.
  testthat::expect_error(chk(row("dn_peak", NA_integer_), fc),
                         "only that form reads it")
})

testthat::test_that("block linkages replace a parameter, SS3 Blk_Fxn 2", {
  yrs <- d$styr:d$endyr
  brk <- c(-Inf, yrs[4] - 0.5, Inf)                   # two blocks
  selFun <- build_selectivity(linkages = list(
    dn_peak = linkage_spec(~ cut(Year, breaks = brk), fleet = fl, link = "identity")))
  m0 <- dn6_build(d, selFun = selFun)
  tbl <- m0$data_list$linkage_table
  j <- which(tbl$process == "sel" & tbl$design_col != "(Intercept)")
  testthat::expect_length(j, 1L)
  # Second block's peak is 8 cm above the base (first block) peak
  ip <- m0$estimated_params
  ip$beta_linkage[j] <- 8
  m <- dn6_build(d, inits = ip, selFun = selFun)
  P  <- m$estimated_params$sel_dn6[, fl, 1]
  Pb <- P; Pb[1] <- P[1] + 8
  testthat::expect_equal(unname(m$quantities$sel_at_length[fl, 1, seq_along(edges), 1]),
                         ss3_pattern24(mids, P, FALSE, FALSE), tolerance = 1e-10)
  testthat::expect_equal(unname(m$quantities$sel_at_length[fl, 1, seq_along(edges), length(yrs)]),
                         ss3_pattern24(mids, Pb, FALSE, FALSE), tolerance = 1e-10)
})

testthat::test_that("the new form is registered where the drift guards look", {
  testthat::expect_equal(unname(sel_map["DoubleNormalSS3"]), 15)
  testthat::expect_true(all(.SEL_DN6_PARAMS %in% SEL_LINKAGE_PARAMS))
  testthat::expect_setequal(unique(unname(LINKAGE_PARAM_CODES$sel[.SEL_DN6_PARAMS])), 6:11)
  d2 <- d; d2$fleet_control$Time_varying_sel[fl] <- "IID"
  testthat::expect_error(data_check(switch_check(d2)), "must be 'Off'")
})
