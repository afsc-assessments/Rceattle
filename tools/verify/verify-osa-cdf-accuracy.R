# tools/verify/verify-osa-cdf-accuracy.R
# The two claims about method = "cdf" that are NOT about its null distribution,
# each measured rather than asserted. `verify-osa-cdf.R` covers the
# distributional self-test; this covers the two properties the documentation
# states as trade-offs, so every number in ?osa_residuals has a script behind it.
#
#   1. cdf is exact only WITHOUT random effects. With them, oneStepPredict
#      integrates the CDF over the latent states by Laplace, and that integrand
#      is a Gaussian times a sigmoid rather than a density. Scored against the
#      exact Kalman innovations of a linear-Gaussian state space model, where the
#      truth is computable in closed form.
#
#   2. |residual| is censored at 8.04 under cdf, in both directions and
#      deliberately. What that costs a reported statistic, and which method
#      actually reports the uncensored value -- which is NOT the package default.
#
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript tools/verify/verify-osa-cdf-accuracy.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))

# ---- 1. How exact is "cdf" under random effects? ----------------------------
# A minimal AR1-latent / Gaussian-observation model, compiled standalone. Its
# one-step-ahead residuals are exactly the standardized Kalman innovations, so
# every method can be scored against a known answer. Rceattle's own dynamics are
# nonlinear in the recruitment deviations, which is why the oracle needs its own
# model rather than a fixture.
cat("== 1. accuracy under random effects, against an exact Kalman oracle ==\n\n")

src <- '
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(obsvec);
  DATA_VECTOR_INDICATOR(keep, obsvec);
  DATA_SCALAR(sd_y); DATA_SCALAR(sd_u); DATA_SCALAR(rho);
  PARAMETER_VECTOR(u);
  int n = obsvec.size();
  Type nll = 0;
  nll -= dnorm(u(0), Type(0), sd_u / sqrt(Type(1) - rho * rho), true);
  for (int t = 1; t < n; ++t) nll -= dnorm(u(t), rho * u(t - 1), sd_u, true);
  for (int t = 0; t < n; ++t) {
    nll -= keep(t) * dnorm(obsvec(t), u(t), sd_y, true);
    Type F = pnorm((obsvec(t) - u(t)) / sd_y);
    nll -= keep.cdf_lower(t) * log(F);
    nll -= keep.cdf_upper(t) * log(Type(1) - F);
  }
  return nll;
}
'
# Build in a tempdir, but do NOT use on.exit() to get back: registered at a
# script's top level it has no function frame to attach to, so under source() it
# fires immediately and the compile then writes ssm.cpp/.o/.so into the package
# root. Same trap as tools/verify/verify-safebounds.R. Wrap it in a function.
dir <- tempfile("ssm"); dir.create(dir)
local({
  owd <- setwd(dir); on.exit(setwd(owd), add = TRUE)
  writeLines(src, "ssm.cpp")
  TMB::compile("ssm.cpp")
})
dyn.load(file.path(dir, TMB::dynlib("ssm")))

# Exact one-step-ahead standardized residuals of the same model, by hand.
kalman <- function(y, rho, sd_u, sd_y) {
  n <- length(y); r <- numeric(n)
  m <- 0; P <- sd_u^2 / (1 - rho^2)                 # stationary prior
  for (t in seq_len(n)) {
    if (t > 1) { m <- rho * m; P <- rho^2 * P + sd_u^2 }
    S <- P + sd_y^2
    r[t] <- (y[t] - m) / sqrt(S)
    K <- P / S; m <- m + K * (y[t] - m); P <- (1 - K) * P
  }
  r
}

set.seed(7)
n <- 20; rho <- 0.7; sd_u <- 0.8
u <- as.numeric(arima.sim(list(ar = rho), n, sd = sd_u))
cat(sprintf("%8s %7s | %-12s %-12s %-12s\n", "sd_y", "ratio",
            "fullGaussian", "oneStepGauss", "cdf"))
for (sd_y in c(2.0, 1.0, 0.5)) {
  y <- u + rnorm(n, 0, sd_y)
  truth <- kalman(y, rho, sd_u, sd_y)
  obj <- TMB::MakeADFun(data = list(obsvec = y, sd_y = sd_y, sd_u = sd_u, rho = rho),
                        parameters = list(u = rep(0, n)), random = "u",
                        DLL = "ssm", silent = TRUE)
  obj$fn(obj$par)
  err <- vapply(c("fullGaussian", "oneStepGaussian", "cdf"), function(m)
    max(abs(suppressWarnings(TMB::oneStepPredict(
      obj, "obsvec", "keep", method = m, discrete = FALSE,
      parallel = FALSE, trace = FALSE)$residual) - truth)), 1)
  cat(sprintf("%8.2f %7.2f | %-12.1e %-12.1e %-12.1e\n",
              sd_y, sqrt(sd_u^2 / (1 - rho^2)) / sd_y, err[1], err[2], err[3]))
}
cat("\n  On a LINEAR-Gaussian model the Gaussian methods integrate a density and\n")
cat("  are exact; cdf integrates a Gaussian times a sigmoid and is not.\n\n")

# ---- 1b. ... but does that carry over to Rceattle? ---------------------------
# It does not follow. Those methods are exact when the one-step-ahead PREDICTIVE
# is Gaussian, i.e. when the model is linear in the random effects. Rceattle's
# index and catch are exp() of cumulated log recruitment deviations pushed
# through the population dynamics, which is not. The test: fullGaussian and
# oneStepGaussian are algebraically identical for a Gaussian conditional, so any
# disagreement between THEM measures how far from Gaussian the conditional is.
cat("== 1b. the same question on Rceattle, where the dynamics are nonlinear ==\n\n")
source("tests/testthat/helpers.R", local = TRUE)
dre <- make_test_data(nyrs = 12, nages = 6, seed = 51)
fre <- suppressWarnings(suppressMessages(fit_mod(
  dre, file = NULL, estimateMode = 1, msmMode = 0, random_rec = TRUE,
  fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
cat("  random effects:", length(fre$obj$env$random), "\n")
cat(sprintf("  %-6s %8s | %-28s %-24s\n", "source", "sd",
            "fullGauss vs oneStepGauss", "cdf vs oneStepGauss"))
for (src in c("index", "catch")) {
  gg <- function(m) suppressWarnings(suppressMessages(osa_residuals(
    fre, source = src, method = m, parallel = FALSE)))$residual
  fg <- gg("fullGaussian"); os <- gg("oneStepGaussian"); cd <- gg("cdf")
  cat(sprintf("  %-6s %8.3f | %-28.4f %-24.4f\n",
              src, stats::sd(os), max(abs(fg - os)), max(abs(cd - os))))
}
cat("\n  Two methods that CANNOT differ for a Gaussian conditional differ by more\n")
cat("  than cdf differs from either. So no method is exact for Rceattle index or\n")
cat("  catch under random effects, and the linear-Gaussian result above does not\n")
cat("  license a recommendation for them. It does license one for `ecov`, whose\n")
cat("  conditional -- a Gaussian measurement of an AR1 latent, with every other\n")
cat("  data term unconditional -- genuinely is linear-Gaussian.\n\n")

# ---- 2. What the 8.04 ceiling costs, and which method escapes it ------------
# One survey observation driven far past what any CDF-reading method can report.
# Fixed effects, so the ONLY thing separating the methods is how each handles it.
cat("== 2. the |residual| ceiling, on an observation past it ==\n\n")
source("tests/testthat/helpers.R", local = TRUE)
dat <- make_test_data(nyrs = 12, nages = 6, seed = 61)
srv <- dat$index_data$Fleet_code == 1
dat$index_data$Observation[srv][6] <- dat$index_data$Observation[srv][6] * 200
fit <- suppressWarnings(suppressMessages(fit_mod(
  dat, file = NULL, estimateMode = 1, msmMode = 0,
  fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))

cat(sprintf("%-24s %9s %9s %12s %8s\n", "method", "min", "max", "non-finite", "SDNR"))
for (m in c("oneStepGaussianOffMode", "oneStepGaussian", "oneStepGeneric", "cdf")) {
  o <- try(suppressWarnings(suppressMessages(
    osa_residuals(fit, source = "index", method = m, parallel = FALSE))),
    silent = TRUE)
  if (inherits(o, "try-error")) { cat(sprintf("%-24s  ERROR\n", m)); next }
  r <- o$residual
  cat(sprintf("%-24s %9.3f %9.3f %12d %8s\n", m,
              suppressWarnings(min(r[is.finite(r)])),
              suppressWarnings(max(r[is.finite(r)])), sum(!is.finite(r)),
              formatC(stats::sd(r), format = "f", digits = 3)))
}
cat("\n  Only oneStepGaussian reports the uncensored value. The package default\n")
cat("  returns NaN on that row -- so its SDNR is unusable, not merely large --\n")
cat("  and oneStepGeneric compresses the outlier harder than cdf's ceiling does.\n")
cat("\ndone.\n")
