# Under the model OSA residuals are i.i.d. N(0, 1), so a fixed |resid| > 3 flag
# expects 0.0027 n outliers per panel: 0.3 at n = 100, 13.5 at n = 5,000, so a
# long length panel always looked worse than a short age panel. The OSA panel
# flags above qnorm(1 - 0.05 / (2 n)) per panel (Bonferroni; 0.05 expected flags
# whatever n). The Pearson panel keeps 3: those residuals are not N(0, 1).

bubble_frame <- function(n, source, r) {
  data.frame(source = source, year = seq_len(n), age_length_bin = 1L,
             residual = c(r, rep(0.1, n - 1)))
}

testthat::test_that("the OSA flag is per panel and does not grow with the panel's size", {
  osa <- rbind(bubble_frame(30, "small", 3.5), bubble_frame(3000, "large", 3.5))
  p <- Rceattle:::.osa_bubble_plot(osa)
  shape <- p$data$shape
  testthat::expect_equal(shape[p$data$source == "small"][1], "outlier")   # 3.5 > 3.14
  testthat::expect_equal(shape[p$data$source == "large"][1], "normal")    # 3.5 < 4.31
  # The cut is two-sided: 3.0 at n = 30 sits below 3.144 (a one-sided cut, 2.935, flags it).
  p30 <- Rceattle:::.osa_bubble_plot(bubble_frame(30, "small", 3.0))
  testthat::expect_equal(p30$data$shape[1], "normal")
  # Exactly at the Bonferroni cut for n = 100 and 5,000: just above flags, just below does not.
  for (n in c(100, 5000)) {
    cut <- stats::qnorm(1 - 0.05 / (2 * n))
    above <- Rceattle:::.osa_bubble_plot(bubble_frame(n, "p", cut + 1e-6))
    below <- Rceattle:::.osa_bubble_plot(bubble_frame(n, "p", cut - 1e-6))
    testthat::expect_equal(above$data$shape[1], "outlier", info = n)
    testthat::expect_equal(below$data$shape[1], "normal", info = n)
  }
})

testthat::test_that("the Pearson panel keeps the fixed 3, and NA residuals do not count", {
  osa <- bubble_frame(3000, "large", 3.5)
  osa$residual[2:1500] <- NA
  p <- Rceattle:::.osa_bubble_plot(osa, outlier = "fixed")
  testthat::expect_equal(p$data$shape[1], "outlier")
  p2 <- Rceattle:::.osa_bubble_plot(osa)                      # n = 1501 finite
  testthat::expect_equal(p2$data$shape[1], "normal")           # cut 4.15
  # A residual beyond the bubble scale is flagged on its untruncated value.
  osa$residual[1] <- 8
  p3 <- suppressWarnings(Rceattle:::.osa_bubble_plot(osa))   # warns that 8 is truncated
  testthat::expect_equal(p3$data$shape[1], "outlier")
})
