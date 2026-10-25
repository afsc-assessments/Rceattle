# The TMB diet likelihood reads each stomach's prey as the contiguous run of
# diet_ctl where stomach_id equals the stomach number (ceattle.cpp section 13.6),
# while build_osa_data() matches with which(stomach_id == i). The two agree only
# while the ids are numbered 0, 1, ... with each stomach's rows consecutive.
# clean_data() guarantees that; these tests pin that data_check() enforces it,
# because breaking it is silent -- every stomach after an out-of-order id drops
# out of the likelihood, with no warning and a lower jnll.

testthat::test_that("clean_data() numbers stomach_id contiguously from 0", {
  testthat::skip_on_cran()

  data("BS2017MS", package = "Rceattle", envir = environment())
  # The bundled diet data carries no stomach_id at all -- clean_data() builds it.
  testthat::expect_false("stomach_id" %in% colnames(BS2017MS$diet_data))

  cleaned <- suppressMessages(suppressWarnings(Rceattle::clean_data(BS2017MS)))
  sid <- as.integer(cleaned$diet_data$stomach_id)

  testthat::expect_gt(length(sid), 0)
  testthat::expect_identical(sort(unique(sid)), 0:max(sid))
  # One run per stomach == every stomach's rows are consecutive.
  testthat::expect_identical(length(rle(sid)$lengths), length(unique(sid)))
})


testthat::test_that("data_check() rejects non-contiguous or gapped stomach_id", {
  testthat::skip_on_cran()

  data("BS2017MS", package = "Rceattle", envir = environment())
  # fit_mod() runs clean_data() -> build_hcr() -> switch_check() -> data_check();
  # data_check() needs the defaults those fill in, so prepare the same way.
  prep <- function(x) {
    x$HCR <- 0                       # build_hcr() would set this
    suppressMessages(suppressWarnings(
      Rceattle::switch_check(Rceattle::clean_data(x))))
  }
  cleaned <- prep(BS2017MS)
  quiet_check <- function(x) suppressMessages(suppressWarnings(data_check(x)))

  # The cleaned data must pass -- otherwise the check is rejecting valid models.
  testthat::expect_no_error(quiet_check(cleaned))

  # Rows of one stomach split apart: move the last row of stomach 0 to the end,
  # so its id reappears after every other stomach.
  scrambled <- cleaned
  first0 <- which(scrambled$diet_data$stomach_id == 0)
  testthat::skip_if(length(first0) < 2, "need a stomach with >1 prey row")
  moved <- utils::tail(first0, 1)
  scrambled$diet_data <- rbind(scrambled$diet_data[-moved, ],
                               scrambled$diet_data[moved, ])
  testthat::expect_error(quiet_check(scrambled), "sorted by 'stomach_id'")

  # Blocks intact but out of ascending order -- what re-sorting diet_data by
  # predator age or year does. Each stomach's rows are still together, so a check
  # that only tested grouping would pass this, while the C++ cursor starts at
  # stomach 0 and runs past nearly every block: on BS2017MS, reversing the rows
  # leaves 1 of 45 stomachs in the likelihood.
  reversed <- cleaned
  reversed$diet_data <- reversed$diet_data[rev(seq_len(nrow(reversed$diet_data))), ]
  testthat::expect_true(
    all(rle(as.integer(reversed$diet_data$stomach_id))$lengths ==
          rle(as.integer(rev(cleaned$diet_data$stomach_id)))$lengths))  # blocks intact
  testthat::expect_error(quiet_check(reversed), "sorted by 'stomach_id'")

  # And a single block moved to the front, which keeps every block intact too.
  moved_block <- cleaned
  blk <- which(moved_block$diet_data$stomach_id == max(moved_block$diet_data$stomach_id))
  moved_block$diet_data <- rbind(moved_block$diet_data[blk, ],
                                 moved_block$diet_data[-blk, ])
  testthat::expect_error(quiet_check(moved_block), "sorted by 'stomach_id'")

  # A gap in the numbering: renumber the last stomach past the end.
  gapped <- cleaned
  mx <- max(gapped$diet_data$stomach_id)
  gapped$diet_data$stomach_id[gapped$diet_data$stomach_id == mx] <- mx + 5L
  testthat::expect_error(quiet_check(gapped), "no gaps")

  # NA ids.
  na_ids <- cleaned
  na_ids$diet_data$stomach_id[1] <- NA_integer_
  testthat::expect_error(quiet_check(na_ids), "contains NA")
})


testthat::test_that("clean_data() repairs a reordered diet table", {
  testthat::skip_on_cran()

  data("BS2017MS", package = "Rceattle", envir = environment())
  cleaned <- suppressMessages(suppressWarnings(Rceattle::clean_data(BS2017MS)))

  # Shuffle the rows, drop the stale ids, and re-clean: the invariant returns.
  shuffled <- cleaned
  set.seed(1)
  shuffled$diet_data <- shuffled$diet_data[sample(nrow(shuffled$diet_data)), ]
  shuffled$diet_data$stomach_id <- NULL
  shuffled$diet_data$stratum_id <- NULL

  shuffled$HCR <- 0
  recleaned <- suppressMessages(suppressWarnings(
    Rceattle::switch_check(Rceattle::clean_data(shuffled))))
  sid <- as.integer(recleaned$diet_data$stomach_id)
  testthat::expect_identical(sort(unique(sid)), 0:max(sid))
  testthat::expect_identical(length(rle(sid)$lengths), length(unique(sid)))
  testthat::expect_no_error(
    suppressMessages(suppressWarnings(data_check(recleaned))))
})
