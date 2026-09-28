# plot() filtered to the finite residuals and warned only when EVERY residual
# was non-finite, so at the documented 1879-of-4538 loss under method = "cdf" it
# drew a clean Q-Q panel, with an SDNR annotation, on a time-biased subset. The
# warning counts OSA rows only: a Pearson panel beside it is a different
# population, with its own exclusions and its own finite filter.

osa_frame <- function(residual, method = "oneStepGaussianOffMode") {
  out <- data.frame(source = "index", fleet = 1L, fleet_name = "Survey",
                    year = seq_along(residual), age_length_bin = 1L,
                    residual = residual, stringsAsFactors = FALSE)
  class(out) <- c("rceattle_osa", "data.frame")
  attr(out, "method") <- method
  out
}

testthat::test_that("plot() warns when it drops non-finite residuals", {
  testthat::skip_if_not_installed("ggplot2")
  set.seed(1)
  x <- osa_frame(c(stats::rnorm(20), NA, NaN, Inf))
  w <- testthat::capture_warnings(plot(x))
  # The dropped count against the total, then what the panel actually describes.
  testthat::expect_true(any(grepl("3 of 23 OSA residual", w, fixed = TRUE)))
  testthat::expect_true(any(grepl("describe only the 20 shown", w,
                                  fixed = TRUE)))
})

testthat::test_that("the contiguous-tail note is only for method = 'cdf'", {
  testthat::skip_if_not_installed("ggplot2")
  set.seed(2)
  r <- c(stats::rnorm(10), NA)

  # A Gaussian method's failures are not a tail, so the note would misdescribe
  # them.
  gauss <- testthat::capture_warnings(plot(osa_frame(r)))
  testthat::expect_false(any(grepl("contiguous tail", gauss, fixed = TRUE)))

  cdf <- testthat::capture_warnings(plot(osa_frame(r, method = "cdf")))
  testthat::expect_true(any(grepl("contiguous tail", cdf, fixed = TRUE)))

  # The attribute is a named vector whenever a family was residualized with its
  # own method, so the note must be found through it, not by a scalar compare.
  nm <- c(default = "cdf", DirichletMultinomial = "oneStepGaussianOffMode")
  named <- testthat::capture_warnings(plot(osa_frame(r, method = nm)))
  testthat::expect_true(any(grepl("contiguous tail", named, fixed = TRUE)))
})

testthat::test_that("a wholly non-finite panel warns and draws nothing", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::expect_warning(plot(osa_frame(c(NA, NaN, Inf))),
                           "No finite residuals to plot")
  out <- suppressWarnings(plot(osa_frame(c(NA, NaN, Inf))))
  testthat::expect_null(out)
})
