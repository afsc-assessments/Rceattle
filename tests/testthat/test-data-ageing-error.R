# A fleet picks its ageing-error matrix through fleet_control's
# Ageing_error_index, so one species can carry several: GOA cod reads its
# pre-2007 survey otoliths about a year older than its post-2007 ones. The array
# is filled per index, and these are the checks that must be keyed the same way.
#
# Two of them were not. Coverage was checked per SPECIES, so a complete matrix
# vouched for an incomplete sibling whose missing true ages stayed 0 and were
# renormalized away by the template -- a silently wrong age composition. And the
# index was checked only for existence, so with no index column on age_error
# (where the index IS the species) a fleet could be fitted with another species'
# ageing error, over another species' age range.

e <- new.env(parent = asNamespace("Rceattle"))
for (h in c("helpers.R", "helpers-make-msm-data.R")) {
  sys.source(testthat::test_path(h), envir = e)
}

ae_data <- local({
  d <- e$make_msm_test_data()$data_list
  list(d = d, nages = d$nages[1])
})

ae_check <- function(d) {
  suppressWarnings(Rceattle:::data_check(Rceattle::switch_check(d)))
}

testthat::test_that("age_error coverage is per matrix, not per species", {
  d <- ae_data$d
  testthat::skip_if(d$nspp < 2)
  ae <- as.data.frame(d$age_error)

  # Species 1 keeps its complete matrix (index 1) and gains a second one that
  # stops short of the oldest ages -- the GOA cod shape, typed by hand.
  short <- ae[ae$Species == 1 & ae$True_age <= 3, , drop = FALSE]
  ae$Ageing_error_index <- ae$Species
  short$Ageing_error_index <- d$nspp + 1L
  d$age_error <- rbind(ae, short)

  # Named by matrix, so the incomplete sibling shows behind the complete one.
  testthat::expect_message(ae_check(d),
                           paste("age_error` matrix", d$nspp + 1L))
})

testthat::test_that("one Ageing_error_index cannot span two species", {
  d <- ae_data$d
  testthat::skip_if(d$nspp < 2)
  ae <- as.data.frame(d$age_error)
  # A matrix's rows are sized and offset by its species' nages / minage.
  ae$Ageing_error_index <- 1L
  d$age_error <- ae

  testthat::expect_error(suppressMessages(ae_check(d)),
                         "more than one Species")
})

testthat::test_that("a fleet cannot borrow another species' ageing error", {
  d <- ae_data$d
  testthat::skip_if(d$nspp < 2)
  # No Ageing_error_index on age_error, so the index IS the species: asking for
  # 2 on a species-1 fleet reads species 2's matrix, which existed and passed.
  flt <- which(as.integer(d$fleet_control$Species) == 1L)[1]
  testthat::skip_if(is.na(flt))
  d$fleet_control$Ageing_error_index <- NA_integer_
  d$fleet_control$Ageing_error_index[flt] <- 2L

  testthat::expect_error(suppressMessages(ae_check(d)),
                         "belonging to another species")

  # The fleet's own matrix is fine.
  d$fleet_control$Ageing_error_index[flt] <- 1L
  testthat::expect_false(inherits(
    tryCatch(suppressMessages(ae_check(d)), error = function(x) x), "error"))
})

testthat::test_that("a blank Ageing_error_index falls back to its species", {
  # The only case here that builds a TMB object; the data_check cases above are
  # fast and stay unguarded.
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  d <- ae_data$d
  testthat::skip_if(d$nspp < 2)
  ae <- as.data.frame(d$age_error)
  # Adding a matrix by typing indices on the new rows leaves the original rows
  # blank; those still mean "the species' own matrix".
  extra <- ae[ae$Species == 1, , drop = FALSE]
  ae$Ageing_error_index <- NA_integer_
  extra$Ageing_error_index <- d$nspp + 1L
  d$age_error <- rbind(ae, extra)

  m <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 0,
    random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
  testthat::expect_true(is.finite(m$quantities$jnll))
})
