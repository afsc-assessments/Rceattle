# Fleets sharing a Catchability_index share ONE catchability (index_log_q). TMB
# collapses a shared parameter to the mean of its members' starting values
# (updateMap(): tapply(par, map, mean)), and index_log_q is
# log(Catchability_init), so the group starts at the GEOMETRIC MEAN of the
# members' inits -- no fleet keeps the value in its own row.
#
# This is worth its own test because a shared q that starts at the mean scales a
# survey's whole predicted index by a constant factor, and no residual pattern
# distinguishes that from a real change in abundance. On GOA Pacific cod a fleet
# left at the default dragged the survey's q to sqrt(1.4964) = 1.2233, costing
# 18% of the index and 6.23 nats with nothing else out of place.
#
# data_check() already reports a shared group whose Catchability forms or
# Time_varying_q differ, and a Fixed lead whose inits differ. The gap this covers
# is the ordinary case: forms agreeing, inits not, a q estimated for the group.
#
# Kept apart from test-selectivity-shared-sigma.R, which covers the same helper
# for the two deviation SDs: those fire only under random_sel / random_q, where
# the consequence is a starting value, while Catchability_init is also a prior
# centre and can be left non-positive, so the cases and the fixtures differ.

make_shared_q_data <- function(init1 = 0.5, init2 = 2.0, form = "Estimated") {
  d <- make_test_data()
  fc <- d$fleet_control
  # Make the second fleet a survey as well, so both carry an estimated q, and put
  # the two in one catchability group. Everything else about them matches, so the
  # inits are the only thing they disagree on.
  fc$Fleet_name          <- c("Survey", "Survey2")
  fc$Fleet_type          <- c("Survey", "Survey")
  fc$Catchability        <- rep(form, 2)
  fc$Catchability_index  <- c(1, 1)
  fc$Catchability_init   <- c(init1, init2)
  # Every other index setting is copied from the first fleet, so the two are
  # identical apart from the init. data_check() requires each fitted index fleet
  # to carry its own.
  for (col in c("Catchability_prior_sd", "Estimate_index_sd", "Index_sd_prior",
                "Index_distribution", "Time_varying_q", "Time_varying_q_sd",
                "Estimate_q"))
    if (col %in% names(fc)) fc[[col]] <- rep(fc[[col]][1], 2)
  d$fleet_control <- fc

  # The second fleet needs index rows of its own, or its index_log_q is mapped
  # out and it contributes nothing to the mean.
  idx <- d$index_data
  second <- idx[idx$Fleet_code == 1, , drop = FALSE]
  second$Fleet_code <- 2
  if ("Fleet_name" %in% names(second)) second$Fleet_name <- "Survey2"
  d$index_data <- rbind(idx, second)
  d
}

build_shared_q <- function(d) {
  suppressMessages(Rceattle::fit_mod(d, estimateMode = 3, msmMode = 0))
}


test_that("a shared estimated catchability warns when the inits differ", {
  w <- collect_warnings(fit <- build_shared_q(make_shared_q_data(0.5, 2.0)))
  hit <- grep("Catchability_init", w, value = TRUE)
  expect_length(hit, 1)
  expect_match(hit, "Fleets sharing Catchability_index 1")
  expect_match(hit, "estimates one catchability")
  # The value reported is the geometric mean, not either fleet's own.
  expect_match(hit, "\\(1\\)", fixed = FALSE)   # sqrt(0.5 * 2) = 1
  expect_match(hit, "no fleet keeps the value in its own row")
})

test_that("the shared catchability really starts at the geometric mean", {
  init1 <- 0.5; init2 <- 2.0
  fit <- suppressWarnings(build_shared_q(make_shared_q_data(init1, init2)))
  q <- fit$obj$par[names(fit$obj$par) == "index_log_q"]
  # One shared parameter for the group, and it is the mean on the log scale --
  # which is neither fleet's Catchability_init.
  expect_length(q, 1)
  expect_equal(as.numeric(q), mean(log(c(init1, init2))), tolerance = 1e-10)
  expect_equal(exp(as.numeric(q)), sqrt(init1 * init2), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(exp(as.numeric(q)), init1)))
  expect_false(isTRUE(all.equal(exp(as.numeric(q)), init2)))
})

test_that("matching inits raise nothing, so the warning is not a false alarm", {
  w <- collect_warnings(build_shared_q(make_shared_q_data(0.8, 0.8)))
  expect_length(grep("Catchability_init", w), 0)
})

test_that("a Fixed group is left to data_check, which reports it differently", {
  # Nothing is estimated, so there is no shared parameter to collapse and each
  # fleet sits at its own init. build_map() must stay quiet; data_check() carries
  # this case and says the inits are used separately.
  w <- collect_warnings(build_shared_q(make_shared_q_data(0.5, 2.0, form = "Fixed")))
  expect_length(grep("estimates one catchability", w), 0)
  expect_gt(length(grep("each fleet uses its own", w)), 0)
})

# The rest exercise .warn_shared_block_start() directly. A shared block is seeded
# from log(Catchability_init), and two configurations can leave that column
# non-positive or blank for a fleet that still joins the block: Analytical and
# AnalyticalArith solve q from the data, so data_check() exempts them from
# requiring a positive value, and a fleet with no fitted index rows is never asked
# for one. Either way the group's start becomes -Inf or NA rather than a mean, so
# the warning must say that instead of naming a geometric mean the fit will never
# start from.
fake_q_group <- function(init, name = c("A", "B"), form = "Estimated") {
  list(fleet_control = data.frame(
    Fleet_name = name, Fleet_code = seq_along(name),
    Catchability_index = rep(1, length(name)),
    Catchability = rep(form, length(name)),
    Catchability_init = init, stringsAsFactors = FALSE))
}

# `note_when` mirrors the live call site: the caveat applies only where a prior
# is actually scored, which ceattle.cpp gates on est_index_q == 2.
warn_q <- function(init, est = c(1, 1), note = NULL, form = "Estimated") {
  collect_warnings(
    Rceattle:::.warn_shared_block_start(
      list(index_log_q = est), fake_q_group(init, form = form),
      "Catchability_index", "Catchability_init", "index_log_q",
      what = "catchability", note = note,
      note_when = function(fc, est)
        any(Rceattle:::.canon_switch(fc$Catchability[est], Rceattle:::q_map) ==
              "Estimated-with-prior")))
}

test_that("a zero init reports a non-finite start, not a geometric mean of 0", {
  w <- warn_q(c(0.5, 0))
  expect_length(w, 1)
  expect_match(w, "whose log is not finite")
  expect_match(w, "starts at -Inf")
  expect_match(w, "cannot fit")
  # The old text would have claimed the group starts at 0, which it does not.
  expect_false(grepl("geometric mean", w))
})

test_that("a blank init reports a non-finite start", {
  w <- warn_q(c(0.5, NA))
  expect_length(w, 1)
  expect_match(w, "<blank>")
  expect_match(w, "starts at NA")
  expect_false(grepl("geometric mean", w))
})

test_that("two usable inits still report the geometric mean", {
  w <- warn_q(c(0.5, 2.0))
  expect_match(w, "geometric mean of those values \\(1\\)")
  expect_false(grepl("not finite", w))
})

test_that("the catchability call carries the prior-centre caveat", {
  # Catchability_init is read twice: as the shared start, and under
  # Estimated-with-prior as the LEAD fleet's prior centre. Differing values there
  # move the objective, so "no fleet keeps its own value" is not the whole story
  # and the message has to say so -- but only where a prior is scored.
  w <- warn_q(c(0.5, 2.0), note = "PRIOR NOTE HERE",
              form = "Estimated-with-prior")
  expect_match(w, "PRIOR NOTE HERE")
})

test_that("a fleet whose block slot is NA contributes nothing", {
  # Only one estimated member left, so there is no shared parameter to collapse.
  expect_length(warn_q(c(0.5, 2.0), est = c(1, NA)), 0)
})

# Round 2 found these: a factor column defeated the non-finite guard, an equally
# unusable pair was silent, and the caveat fired where no prior is scored.
test_that("a factor column does not defeat the non-finite guard", {
  # as.numeric() on a factor returns LEVEL CODES, all positive, so a zero would
  # have been reported as a real geometric mean (1.414) of values whose true
  # shared start is -Inf.
  w <- warn_q(factor(c("0.5", "0")))
  expect_length(w, 1)
  expect_match(w, "whose log is not finite")
  expect_false(grepl("geometric mean", w))
  expect_false(grepl("1.414", w, fixed = TRUE))
})

test_that("members equally unusable still warn", {
  # The geometric-mean message only matters when the inits DIFFER, but "seeds the
  # group at -Inf and cannot fit" does not depend on that.
  w <- warn_q(c(0, 0))
  expect_length(w, 1)
  expect_match(w, "A, B carry 0, 0")
  expect_match(w, "cannot fit")
})

test_that("an empty string reads as blank, and plurals agree", {
  expect_match(warn_q(c("0.5", "")), "B carries <blank>")
  expect_match(warn_q(c(0, NA)), "A, B carry 0, <blank>")
})

test_that("a negative init does not leak log()'s own NaN warning", {
  w <- warn_q(c(0.5, -1))
  expect_length(w, 1)
  expect_false(any(grepl("NaNs produced", w)))
})

test_that("the prior caveat appears only where a prior is scored", {
  # Catchability code 6 ("AR1") is REMOVED and data_check() errors on it, so the
  # caveat must not mention it, and plain Estimated scores no prior at all.
  expect_false(any(grepl("NOTE", warn_q(c(0.5, 2), note = "NOTE"))))
  expect_true(any(grepl("NOTE",
    warn_q(c(0.5, 2), note = "NOTE", form = "Estimated-with-prior"))))
})
