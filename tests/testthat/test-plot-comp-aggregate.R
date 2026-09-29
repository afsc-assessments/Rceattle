# The aggregated composition pools counts across years rather than averaging
# proportions, so the aggregate is a sum of independent draws and its variance
# is the exact sum of the per-year variances. afscOSA's qbinom band at the
# pooled sample size was deliberately not ported: it assumes a distribution the
# old mean-of-proportions statistic did not have. Checked by simulation.

set.seed(90210)

# One panel's worth of long-format rows, as .comp_resid_long() produces.
make_long <- function(obs, p, N, sd = NULL) {
  Y <- nrow(p); B <- ncol(p)
  data.frame(
    Fleet = 1L, Fleet_name = "f", Species = 1L, Sex = 0L, comp_type = 0L,
    Year = rep(seq_len(Y), times = B), N = rep(N, times = B),
    bin = rep(seq_len(B), each = Y),
    obs = as.numeric(obs), hat = as.numeric(p),
    Sd = if (is.null(sd)) as.numeric(sqrt(p * (1 - p) / N)) else as.numeric(sd),
    sex_grp = "combined", bin_lab = "Age", panel = "f", type_lab = "age",
    stringsAsFactors = FALSE)
}


test_that("the aggregate pools counts rather than averaging proportions", {
  # Two years with very different sample sizes and different compositions. The
  # pooled statistic must follow the big year; a mean of proportions would sit
  # halfway between them.
  p   <- rbind(c(0.8, 0.1, 0.1), c(0.1, 0.1, 0.8))
  obs <- p
  N   <- c(10, 990)

  agg <- .comp_aggregate(make_long(obs, p, N))
  expect_equal(agg$obs[agg$bin == 1], (0.8 * 10 + 0.1 * 990) / 1000,
               tolerance = 1e-10)
  expect_equal(agg$ISS[1], 1000)
  # A mean of proportions would have given 0.45 in bin 1.
  expect_false(isTRUE(all.equal(agg$obs[agg$bin == 1], 0.45, tolerance = 1e-6)))
})


test_that("the aggregate band covers at its nominal rate", {
  skip_on_cran()
  Y <- 25L; B <- 8L; nrep <- 1500
  p <- t(replicate(Y, { v <- exp(-0.5 * ((1:B - 4) / 1.5)^2) + 0.02; v / sum(v) }))
  N <- round(exp(stats::rnorm(Y, log(150), 0.6)))

  inside <- matrix(NA, nrep, B)
  for (i in seq_len(nrep)) {
    obs <- t(vapply(seq_len(Y),
                    function(y) stats::rmultinom(1, N[y], p[y, ])[, 1] / N[y],
                    numeric(B)))
    a <- .comp_aggregate(make_long(obs, p, N))
    a <- a[order(a$bin), ]
    inside[i, ] <- a$obs >= a$lwr & a$obs <= a$upr
  }
  cover <- colMeans(inside)
  # Normal approximation on an exact variance: close to nominal, and never
  # wildly conservative the way a pooled-binomial band would be.
  expect_gt(min(cover), 0.90)
  expect_lt(max(cover), 0.99)
})


test_that("projection rows are excluded from the aggregate and its ISS", {
  p   <- rbind(c(0.5, 0.5), c(0.5, 0.5), c(0.5, 0.5))
  obs <- p
  N   <- c(100, 100, 100)
  long <- make_long(obs, p, N)
  long$Year <- rep(c(1, 2, 99), times = 2)   # year 99 is past endyr

  expect_equal(.comp_aggregate(long, endyr = 2)$ISS[1], 200)
  expect_equal(.comp_aggregate(long)$ISS[1], 300)   # ungated, for comparison
})


test_that("a bin with no assumed sd keeps its composition, losing only the band", {
  # aggregate()'s formula interface drops whole rows on any NA, so this bin
  # would otherwise vanish from the aggregated composition entirely.
  p   <- rbind(c(0.5, 0.3, 0.2), c(0.5, 0.3, 0.2))
  obs <- p
  N   <- c(100, 100)
  long <- make_long(obs, p, N)
  long$Sd[long$bin == 2] <- NA_real_

  agg <- .comp_aggregate(long)
  expect_equal(sort(agg$bin), 1:3)
  expect_equal(agg$obs[agg$bin == 2], 0.3, tolerance = 1e-10)
  expect_true(is.na(agg$lwr[agg$bin == 2]))
  expect_true(all(is.finite(agg$lwr[agg$bin != 2])))
})


test_that("a wider assumed sd widens the band proportionally", {
  # The band reads its variance from the same Sd the residuals were divided by,
  # so an overdispersed fleet gets a wider band without any change here.
  p   <- rbind(c(0.4, 0.6), c(0.4, 0.6))
  obs <- p
  N   <- c(100, 100)
  narrow <- .comp_aggregate(make_long(obs, p, N))
  wide   <- .comp_aggregate(make_long(obs, p, N,
                                      sd = sqrt(p * (1 - p) / N) * 3))
  expect_equal((wide$upr - wide$lwr) / (narrow$upr - narrow$lwr), rep(3, 2),
               tolerance = 1e-8)
})


# The annotated ESS is recovered from the assumed sd rather than from the weight
# column: sd^2 = p(1-p)/N_eff holds for every family, so one expression covers a
# weighted multinomial and a Dirichlet-multinomial alike. Handing the function an
# sd built by hand and asserting it inverts that sd is algebra returned to
# sender, so the sds below come from .rce_comp_pearson() itself -- the same code
# the residuals and the band use -- and the target is measured by simulation.
test_that("the annotated ESS is the effective N the likelihood assumed", {
  skip_on_cran()
  Y <- 4L; B <- 8L; nrep <- 6000
  p <- t(replicate(Y, { v <- exp(-0.5 * ((1:B - 4) / 1.6)^2); v / sum(v) }))
  N <- c(30, 150, 400, 1200)
  off <- 1e-5

  # sd exactly as the fleet's own likelihood assumes it, via the real helper.
  real_sd <- function(family, weight) {
    t(vapply(seq_len(Y), function(y)
      attr(.rce_comp_pearson(p[y, ], p[y, ], N[y], rep(1L, B),
                             rep(family, B), rep(weight, B), off), "sd"),
      numeric(B)))
  }
  # Effective N measured from the spread of draws, bin by bin, pooled over years
  # exactly as .comp_aggregate() pools counts.
  measured_ess <- function(draw) {
    o <- replicate(nrep, colSums(t(vapply(seq_len(Y),
      function(y) draw(y) * N[y], numeric(B)))))
    tot  <- sum(N)
    pbar <- rowMeans(o) / tot
    mean(pbar * (1 - pbar) / (apply(o / tot, 1L, stats::var)))
  }

  # A weighted multinomial: N_eff is w * N, year by year.
  w   <- 2.5
  agg <- .comp_aggregate(make_long(p, p, N, sd = real_sd(0L, w)))
  expect_equal(agg$ESS[1], measured_ess(function(y)
    stats::rmultinom(1, round(w * N[y]), p[y, ])[, 1] / round(w * N[y])),
    tolerance = 0.05)

  # A Dirichlet-multinomial: N_eff is below N and not proportional to it, so a
  # formula that only rescaled the input N would pass the multinomial case and
  # fail here.
  theta <- 0.4
  agg_dm <- .comp_aggregate(make_long(p, p, N, sd = real_sd(1L, log(theta))))
  expect_equal(agg_dm$ESS[1], measured_ess(function(y) {
    a <- N[y] * p[y, ] * theta
    g <- stats::rgamma(B, a, 1)
    stats::rmultinom(1, N[y], g / sum(g))[, 1] / N[y]
  }), tolerance = 0.05)
  expect_lt(agg_dm$ESS[1], agg_dm$ISS[1])

  # The offset is carried in the sd but not in the raw hat this reads, so the
  # recovery is approximate. Pin how approximate, since the roxygen claims it.
  expect_equal(agg$ESS[1] / (w * sum(N)), 1, tolerance = 1e-3)
})


test_that("each observation carries its own effective N, not each year", {
  # Regression: grouping on Year alone pooled a female-only and a male-only row
  # from the same fleet-year into one observation, halving ESS while ISS
  # correctly counted both -- the figure then read as a fleet downweighted to
  # 0.5 when it was at weight 1. data_check() treats (Fleet, Species, Sex, Year)
  # as the uniqueness key, so that is the key an effective N belongs to.
  pv <- c(0.4, 0.3, 0.2, 0.1)

  # One observation: one (Fleet, Species, Sex, Year), its bins, its sample size.
  one_obs <- function(sex, year, n, grp = "combined") {
    d <- make_long(matrix(pv, 1L), matrix(pv, 1L), n,
                   sd = sqrt(pv * (1 - pv) / n))
    d$Sex <- sex; d$Year <- year; d$sex_grp <- grp
    d
  }

  # A female-only and a male-only row in the SAME year: two observations, which
  # data_check() explicitly allows. ESS must count both.
  d <- rbind(one_obs(1L, 2000L, 200), one_obs(2L, 2000L, 370))
  agg <- .comp_aggregate(d)
  expect_equal(agg$ESS[1], 570, tolerance = 1e-6)
  expect_equal(agg$ISS[1], 570)

  # Joint-sex (Sex == 3) is the opposite case: ONE multinomial spans both sexes,
  # so its two sex_grp halves must NOT be counted as two observations.
  dj <- rbind(one_obs(3L, 2000L, 200, "female"),
              one_obs(3L, 2000L, 200, "male"))
  expect_equal(.comp_aggregate(dj)$ESS[1], 200, tolerance = 1e-6)

  # And two years of one sex still sum, so the new key did not lose the year.
  dy <- rbind(one_obs(1L, 2000L, 200), one_obs(1L, 2001L, 370))
  expect_equal(.comp_aggregate(dy)$ESS[1], 570, tolerance = 1e-6)
})


test_that("a missing assumed sd gives NA rather than a partial effective N", {
  # A silently short ESS would read as a downweighted fleet.
  Y <- 3L; B <- 4L
  p <- t(replicate(Y, c(0.4, 0.3, 0.2, 0.1)))
  N <- c(50, 60, 70)
  sd <- sqrt(p * (1 - p) / N)
  sd[2, ] <- NA_real_

  agg <- .comp_aggregate(make_long(p, p, N, sd = sd))
  expect_true(is.na(agg$ESS[1]))
  expect_equal(agg$ISS[1], sum(N))

  # One unusable bin costs that observation precision, not the panel its
  # annotation: the largest usable bin still carries the effective N.
  sd2 <- sqrt(p * (1 - p) / N); sd2[2, 1] <- NA_real_
  expect_equal(.comp_aggregate(make_long(p, p, N, sd = sd2))$ESS[1], sum(N),
               tolerance = 1e-6)

  # Comp_weights = 0 drops the fleet's composition from the likelihood, giving
  # an infinite assumed sd. That is an effective sample size of zero, which is
  # known, not unknown.
  sd0 <- matrix(Inf, Y, B)
  expect_equal(.comp_aggregate(make_long(p, p, N, sd = sd0))$ESS[1], 0)
})


# A fleet can carry more than one observation structure -- sexes combined in the
# early years, disaggregated later. The likelihood scores one density per
# structure, so pooling them into one panel divides both by a shared total:
# neither series sums to 1 and the two are not on comparable axes. Reported by an
# IPHC user running a two-sex model, 2026-09-28.
test_that("a fleet mixing sex structures gets one panel per structure", {
  # Two rows for one fleet and comp type: a sexes-combined year and a joint-sex
  # year, built as .comp_resid_long() labels them.
  long <- data.frame(
    Fleet = 1L, Fleet_name = "f", Species = 1L, comp_type = 0L,
    Sex = c(0L, 0L, 3L, 3L), Year = c(1L, 1L, 2L, 2L), N = 100,
    bin = c(1L, 2L, 1L, 2L), obs = c(0.5, 0.5, 0.5, 0.5),
    hat = c(0.5, 0.5, 0.5, 0.5), Sd = 0.05,
    sex_grp = c("combined", "combined", "female", "male"),
    norm_grp = c("sexes combined", "sexes combined", "joint-sex", "joint-sex"),
    bin_lab = "Age", type_lab = "age", stringsAsFactors = FALSE)
  long$panel <- paste0(long$Fleet_name, " - ", long$type_lab, " comp (",
                       long$norm_grp, ")")

  # Two panels, and each aggregates to a composition that sums to 1 on its own.
  panels <- unique(long$panel)
  expect_length(panels, 2L)
  for (nm in panels) {
    agg <- Rceattle:::.comp_aggregate(long[long$panel == nm, , drop = FALSE])
    expect_equal(sum(agg$obs), 1)
    expect_equal(sum(agg$hat), 1)
  }

  # Pooled into one panel, as before the split, neither structure sums to 1:
  # that is the defect, and it is what the panel key now prevents.
  both <- Rceattle:::.comp_aggregate(long)
  expect_equal(sum(both$obs), 1)                       # across BOTH structures
  expect_lt(sum(both$obs[both$sex_grp == "combined"]), 1)
})

test_that("a fleet with one sex structure keeps its original panel name", {
  # The overwhelming majority, and every bundled data set: the label must not
  # gain a suffix, or every existing figure's title changes.
  fit <- readRDS(testthat::test_path("fixtures", "fit_baseline.rds"))
  long <- Rceattle:::.comp_resid_long(fit)
  skip_if(is.null(long) || !nrow(long), "no composition rows in the fixture")
  by_fleet <- split(long$norm_grp, paste(long$Fleet_name, long$type_lab))
  skip_if(any(vapply(by_fleet, function(x) length(unique(x)) > 1L, logical(1))),
          "fixture has a mixed-structure fleet")
  expect_false(any(grepl("\\(", long$panel, fixed = FALSE) &
                     grepl("only|combined|joint", long$panel)))
})
