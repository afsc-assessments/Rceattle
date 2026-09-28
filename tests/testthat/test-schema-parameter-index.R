# parameter_index() locates every estimated parameter in the model's own
# coordinates, so a convergence diagnostic can name the fleet, bin and year it
# means rather than repeating the parameter block.
#
# Provenance: the estimability check reported "(sel_coff_dev, log_sel_slp,
# sel_inf)" for 49 flagged parameters, which does not say where to look. The
# recovery goes through TMB's own parList(), so it tracks build_map() exactly;
# these tests pin that against obj$par rather than against a re-derivation.

testthat::test_that("the dictionary's declared rank matches the built array", {
  testthat::skip_on_cran()
  # The dictionary names what each dimension MEANS; parameter_index() reads the
  # extent off the array. That only works if the two agree on how many
  # dimensions a block has. M1_dev_log_sd was declared [nspp, nsex, 2] against a
  # [nspp, nsex] array, which is the drift this catches.
  fit <- Rceattle::fit_mod(data_list = Rceattle::BS2017SS, estimateMode = 3,
                           msmMode = 0,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0))
  pl   <- fit$obj$env$parList(x = seq_along(fit$obj$par))
  dict <- Rceattle::parameter_dictionary()

  strip <- function(s) gsub("]", "", gsub("[", "", s, fixed = TRUE), fixed = TRUE)
  for (block in names(pl)) {
    dd <- dict$dims[match(block, dict$internal)]
    testthat::expect_false(is.na(dd),
      info = paste(block, "is not in parameter_dictionary()"))
    toks <- trimws(strsplit(strip(dd), ",")[[1]])
    toks <- toks[nzchar(toks)]
    actual <- dim(pl[[block]])
    if (is.null(actual)) actual <- length(pl[[block]])
    testthat::expect_equal(length(toks), length(actual),
      info = paste0(block, ": dims '", dd, "' declares ", length(toks),
                    " dimension(s), array has ", length(actual)))
  }
})


testthat::test_that("every shortened age axis states how many ages it omits", {
  # A token like `nages-1` does not say which end it drops. A new one added to
  # the dictionary but not to .PAR_AXIS_OFFSET would label its ages one bin
  # adrift.
  dict <- Rceattle::parameter_dictionary()
  strip <- function(s) gsub("]", "", gsub("[", "", s, fixed = TRUE), fixed = TRUE)
  toks <- unlist(lapply(stats::na.omit(dict$dims),
                        function(x) trimws(strsplit(strip(x), ",")[[1]])))
  short <- unique(toks[grepl("-", toks, fixed = TRUE)])
  testthat::expect_true(all(short %in% names(Rceattle:::.PAR_AXIS_OFFSET)),
    info = paste("no .PAR_AXIS_OFFSET entry for:",
                 paste(setdiff(short, names(Rceattle:::.PAR_AXIS_OFFSET)),
                       collapse = ", ")))
})


testthat::test_that("the dictionary's declared extent matches the built array", {
  testthat::skip_on_cran()
  # Rank alone does not catch a token naming the wrong axis: log_F was declared
  # [n_fsh, nyrs] against an [n_flt, nyrs_hind] array, and log_pop_scalar
  # [nspp, nyrs] against an age axis -- which would have labelled ages as
  # calendar years.
  d <- Rceattle::BS2017SS
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0))
  pl   <- fit$obj$env$parList(x = seq_along(fit$obj$par))
  dl   <- fit$data_list
  dict <- Rceattle::parameter_dictionary()

  # Arrays are allocated to the widest species, so nsex/nages take the maximum.
  ext <- c(nspp = dl$nspp, nsex = max(dl$nsex), n_flt = nrow(dl$fleet_control),
           n_sel = nrow(dl$fleet_control),
           nyrs = dl$projyr - dl$styr + 1L, nyrs_hind = dl$endyr - dl$styr + 1L,
           nages = max(dl$nages),
           n_sel_bins = max(dl$fleet_control$N_sel_bins, na.rm = TRUE))
  # init_dev estimates nages-1 deviations inside a [nspp, nages] array, so the
  # declared extent is the estimated portion. The only such block.
  skip_extent <- list(init_dev = 2L)

  strip <- function(s) gsub("]", "", gsub("[", "", s, fixed = TRUE), fixed = TRUE)
  for (block in names(pl)) {
    dd <- dict$dims[match(block, dict$internal)]
    if (is.na(dd)) next
    toks <- trimws(strsplit(strip(dd), ",")[[1]])
    toks <- toks[nzchar(toks)]
    actual <- dim(pl[[block]])
    if (is.null(actual)) actual <- length(pl[[block]])
    if (length(toks) != length(actual)) next          # the rank test owns this
    for (k in seq_along(toks)) {
      if (identical(skip_extent[[block]], k)) next
      e <- ext[toks[k]]
      if (is.na(e)) next                              # extent not knowable here
      testthat::expect_equal(unname(e), actual[k],
        info = paste0(block, " dim ", k, ": token '", toks[k], "' is ", e,
                      ", array is ", actual[k]))
    }
  }
})


testthat::test_that("every estimated parameter is located exactly once", {
  testthat::skip_on_cran()
  fit <- Rceattle::fit_mod(data_list = Rceattle::BS2017SS, estimateMode = 3,
                           msmMode = 0,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0))
  idx <- Rceattle::parameter_index(fit)

  # One row per element of obj$par, and the blocks agree with TMB's own names.
  testthat::expect_equal(nrow(idx), length(fit$obj$par))
  testthat::expect_equal(sort(idx$par_index), seq_along(fit$obj$par))
  testthat::expect_equal(idx$block[order(idx$par_index)],
                         unname(names(fit$obj$par)))
})


testthat::test_that("coordinates are the model's own labels, not indices", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0))
  idx <- Rceattle::parameter_index(fit)

  # Fleets by name, years as calendar years -- an index would be no better than
  # the block name this replaces.
  f <- idx[idx$block == "log_F", ]
  testthat::expect_true(all(f$fleet %in% d$fleet_control$Fleet_name))
  yr <- as.numeric(stats::na.omit(f$year))
  testthat::expect_true(all(yr >= d$styr & yr <= d$projyr))

  # rec_dev is indexed by species, and BS2017SS names three.
  r <- idx[idx$block == "rec_dev", ]
  testthat::expect_true(all(r$species %in% d$spnames))
})


# ceattle.cpp reads init_dev(sp, age - 1) over age = 1 .. nages-1, so its ages
# run minage + 1 .. minage + nages - 1. minage is 1 in every bundled dataset, so
# this test is the only thing holding an age apart from an array position here.
testthat::test_that("init_dev ages start one above minage", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0))
  idx <- Rceattle::parameter_index(fit)
  i <- idx[idx$block == "init_dev", ]
  testthat::expect_gt(nrow(i), 0)

  for (sp in seq_along(d$spnames)) {
    ages <- sort(as.numeric(i$age[i$species == d$spnames[sp]]))
    testthat::expect_equal(
      ages, seq(d$minage[sp] + 1L, d$minage[sp] + d$nages[sp] - 1L),
      info = d$spnames[sp])
  }

  # The oldest initial age is named; reading the axis from minage dropped it.
  testthat::expect_true(any(as.numeric(i$age) == max(d$minage + d$nages - 1L)))
})


testthat::test_that("a mirrored parameter is one row naming every fleet it drives", {
  testthat::skip_on_cran()
  # Fleets sharing a Selectivity_index share ONE parameter block, so reporting
  # them as several parameters would overstate how many things are unidentified.
  d <- Rceattle::BS2017SS
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0))
  idx <- Rceattle::parameter_index(fit)

  testthat::expect_true(all(idx$n_cells >= 1))
  # However many cells a parameter drives, it is still one parameter.
  testthat::expect_equal(sum(duplicated(idx$par_index)), 0L)
  mirrored <- idx[idx$n_cells > 1, ]
  if (nrow(mirrored) > 0) {
    testthat::expect_true(all(grepl(" \\+ ", mirrored$fleet) |
                              grepl(" \\+ ", mirrored$species)))
  }
})


testthat::test_that("an axis the model does not distinguish is not printed", {
  testthat::skip_on_cran()
  # GOApollock is one species and one sex, so naming either on every line says
  # nothing. The structured column still carries the species.
  d <- Rceattle::GOApollock
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0))
  idx <- Rceattle::parameter_index(fit)

  testthat::expect_equal(d$nsex, 1)
  testthat::expect_true(all(is.na(idx$sex)))          # single-sex: no sex exists
  testthat::expect_false(any(grepl("combined", idx$label)))

  # The species is constant, so it is not a component of any label. Test the
  # components, not the string: every GOApollock fleet is named "Pollock_..."
  # and a substring match would find the species inside the fleet.
  AX <- c("species", "fleet", "sex", "age", "bin", "year", "slot")
  testthat::expect_false("species" %in% Rceattle:::.rce_varying_axes(idx, AX))
  parts <- unlist(strsplit(idx$label, ", ", fixed = TRUE))
  testthat::expect_false(d$spnames[1] %in% parts)
  # ... but the structured column still carries it.
  testthat::expect_true(any(idx$species == d$spnames[1], na.rm = TRUE))
})


testthat::test_that("the summary collapses ordinal axes and counts the group", {
  testthat::skip_on_cran()
  fit <- Rceattle::fit_mod(data_list = Rceattle::GOApollock, estimateMode = 3,
                           msmMode = 0,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0))
  idx <- Rceattle::parameter_index(fit)
  lf  <- idx$par_index[idx$block == "log_F"]
  out <- Rceattle:::.rce_par_summary(lf, idx)

  testthat::expect_length(out, 1L)
  # Years are ordinal, so they read as a range rather than 49 values.
  testthat::expect_match(out, "log_F (F)", fixed = TRUE)
  testthat::expect_match(out, "\\d{4}-\\d{4}")
  testthat::expect_match(out, sprintf("\\(%d\\)$", length(lf)))
})


testthat::test_that("a selectivity slot is named from the fleet's Selectivity", {
  testthat::skip_on_cran()
  # Slot 2 of sel_inf is a descending inflection only for the double-logistic
  # family. DoubleNormal (8) reuses it for the right-tail floor and LogisticPM
  # (11) for the free age-1 selectivity, so a fixed label names the wrong
  # quantity on those fleets.
  d <- Rceattle::BS2017SS
  dbl <- Rceattle:::.rce_sel_slot_labels("sel_inf", d, 1L)

  d8 <- d; d8$fleet_control$Selectivity[1] <- 8
  d11 <- d; d11$fleet_control$Selectivity[1] <- 11
  testthat::expect_equal(Rceattle:::.rce_sel_slot_labels("sel_inf", d8, 1L)[2],
                         "logit right-tail floor")
  testthat::expect_equal(Rceattle:::.rce_sel_slot_labels("sel_inf", d11, 1L)[2],
                         "age-1 log-selectivity")
  testthat::expect_equal(dbl[2], "descending")
  # The deviation blocks carry the slots of the parameter they deviate from.
  testthat::expect_equal(Rceattle:::.rce_sel_slot_labels("sel_inf_dev", d8, 1L),
                         Rceattle:::.rce_sel_slot_labels("sel_inf", d8, 1L))
  # An unrecognised or out-of-range fleet falls back rather than erroring.
  testthat::expect_equal(Rceattle:::.rce_sel_slot_labels("sel_inf", d, NA), dbl)
})


testthat::test_that("sel_curve_pen slots are the AR1 correlations, not penalties", {
  testthat::skip_on_cran()
  # build_map_selectivity() maps sel_curve_pen out entirely and re-enables it
  # only under 2DAR1 (slots 1-2) and 3DAR1 (slots 1-3), where the slots hold
  # logit-scale AR1 correlations across bins, years and cohorts. An estimated
  # element is therefore never the fleet_control penalty weight.
  testthat::expect_equal(Rceattle:::.PAR_SLOT_LABELS$sel_curve_pen,
                         c("bin correlation", "year correlation",
                           "cohort correlation"))
})


testthat::test_that("rec_pars slot 2 is the SRR alpha", {
  testthat::skip_on_cran()
  # ceattle.cpp forms alpha(sp, yr) = exp(rec_pars(sp, 1) + linkage offset).
  # Steepness is derived from alpha and SPR0; it is not the parameter.
  testthat::expect_equal(Rceattle:::.PAR_SLOT_LABELS$rec_pars[2], "alpha")
  testthat::expect_false(any(grepl("steepness",
                                   Rceattle:::.PAR_SLOT_LABELS$rec_pars)))
})


testthat::test_that("the index covers a fit carrying random effects", {
  testthat::skip_on_cran()
  # parameter_index() tags the FIXED effects only, through parList()'s
  # fixed/random split. With random recruitment deviations the random block must
  # be absent from the index and obj$par alike.
  fit <- Rceattle::fit_mod(data_list = Rceattle::BS2017SS, estimateMode = 3,
                           msmMode = 0, random_rec = TRUE,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0))
  idx <- Rceattle::parameter_index(fit)
  testthat::expect_equal(sort(idx$par_index), seq_along(fit$obj$par))
  testthat::expect_equal(idx$block[order(idx$par_index)],
                         unname(names(fit$obj$par)))
  testthat::expect_false("rec_dev" %in% idx$block)   # random, so not a fixed par
})


testthat::test_that("the bounds table's rows line up with the index", {
  testthat::skip_on_cran()
  # parameters_on_bounds indexes .capture_opt_convergence()'s parameter vector,
  # which is assembled by walking the map rather than by reading obj$par. The
  # coordinate labels are only right if the two orders agree.
  fit <- Rceattle::fit_mod(data_list = Rceattle::BS2017SS, estimateMode = 1,
                           msmMode = 0,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0,
                                                               phase = FALSE))
  idx <- Rceattle::parameter_index(fit)
  ch  <- fit$.conv_hindcast
  testthat::skip_if(is.null(ch$par), "no bounds vector captured")
  testthat::expect_equal(length(ch$par), nrow(idx))
  testthat::expect_equal(unname(names(ch$par)), idx$block[order(idx$par_index)])
})


# Under estimateMode = "Estimate" with an estimating HCR -- the default, and what
# every Tier 3 assessment runs -- fit_mod() rebuilds the model for the
# projection, so fit$obj holds log_Ftarget / log_Flimit alone while the
# estimability and bounds checks still index the 584-parameter hindcast vector.
# Labelling one from the other named log_Flimit where rec_pars was flagged.
testthat::test_that("an HCR projection does not relabel the hindcast parameters", {
  testthat::skip_on_cran()

  d <- Rceattle::Atka2022
  fsh <- which(as.character(d$fleet_control$Fleet_type) %in% c("Fishery", "1"))
  d$fleet_control$Proj_F_proportion[fsh] <- 1 / length(fsh)  # a Tier 3 HCR needs a split

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, msmMode = 0, estimateMode = "Estimate",
    HCR = Rceattle::build_hcr(HCR = 5, Ftarget = 0.4, Flimit = 0.35,
                              Plimit = 0.2, Alpha = 0.05),
    fit_control = Rceattle::fit_control(verbose = 0))))

  hind <- names(fit$.conv_hindcast$par)
  testthat::expect_true(length(fit$obj$par) < length(hind))   # the two differ
  testthat::expect_setequal(unique(names(fit$obj$par)),
                            c("log_Flimit", "log_Ftarget"))

  # The hindcast index was captured before the projection replaced obj, and it
  # is the one the hindcast-indexed checks resolve to.
  idx <- Rceattle:::.conv_par_index(fit)
  testthat::expect_equal(nrow(idx$hindcast), length(hind))
  testthat::expect_identical(Rceattle:::.conv_index_for(idx, hind), idx$hindcast)

  # The projection index describes a different vector and is refused for it, so
  # a check falls back to block names rather than printing a wrong coordinate.
  testthat::expect_null(
    Rceattle:::.conv_index_for(list(final = idx$final), hind))
  testthat::expect_false(Rceattle:::.conv_index_ok(idx$final, hind))

  # convergence_diagnostics() runs on such a fit without inventing coordinates.
  cd <- Rceattle::convergence_diagnostics(fit)
  testthat::expect_true(cd$status %in% c("OK", "WARN", "FAIL"))
})
