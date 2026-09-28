# tools/verify/verify-osa-cdf.R
# Does osa_residuals(method = "cdf") return residuals with the distribution it
# claims? "No more negative expected counts" is true by construction and proves
# nothing, so this measures the thing that matters instead.
#
# The exact test. Simulate observations from a fitted model, then residualize
# them AT THE PARAMETERS THAT GENERATED THEM (a build with estimateMode = 3 and
# inits = the fit's own estimates sits exactly there, and returns the real
# objective). Under a correct conditional CDF the probability integral transform
# makes the residuals exactly iid standard normal -- no estimation error, no
# asymptotics. So per-replicate Kolmogorov-Smirnov p-values must look uniform,
# not merely "not tiny", and the within-composition lag-1 autocorrelation must
# sit at zero.
#
# It also reports the Gaussian default over the same replicates, because the
# comparison is the point: the two methods answer the same question and the
# question has a known answer here.
#
# Reports, per method:
#   pooled mean / sd, share of KS p-values below 0.05, and the mean lag-1
#   autocorrelation within a fleet-year composition with its null standard error.
#
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript tools/verify/verify-osa-cdf.R [nrep]

suppressMessages(pkgload::load_all(".", quiet = TRUE))

nrep <- {
  a <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)[1]))
  if (is.na(a)) 20L else a
}

fit_bs <- function() {
  data(BS2017SS, package = "Rceattle", envir = environment())
  suppressWarnings(suppressMessages(fit_mod(
    BS2017SS, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
}

# A fitted object rebuilt on `dat` at `fit`'s own parameter estimates. Mode 3
# builds without optimizing, so the object sits at `inits` exactly; the mode is
# then relabelled because osa_residuals() refuses mode >= 3 -- a guard against an
# unoptimized model, which this deliberately is not.
at_truth <- function(fit, dat) {
  sf <- suppressWarnings(suppressMessages(fit_mod(
    dat, file = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    inits = fit$estimated_params,
    fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
  sf$data_list$estimateMode <- 1
  sf
}

# Lag-1 autocorrelation of the residual sequence within each composition, pooled
# over compositions. Under the null the conditional binomials are independent, so
# this is 0 with standard error 1/sqrt(number of adjacent pairs).
within_acf1_parts <- function(osa) {
  g <- split(osa$residual, list(osa$source, osa$fleet, osa$year), drop = TRUE)
  num <- 0; den <- 0; npair <- 0
  for (r in g) {
    r <- r[is.finite(r)]
    if (length(r) < 3) next
    num <- num + sum(r[-1] * r[-length(r)])
    den <- den + sum(r^2)
    npair <- npair + length(r) - 1
  }
  c(num = num, den = den, npair = npair)
}

# Pool those parts over replicates, so the ratio and its null standard error
# both describe the whole run rather than one replicate of it.
pool_acf1 <- function(parts) {
  if (is.null(parts)) return(c(acf1 = NA_real_, se = NA_real_))
  num <- sum(parts[, "num"]); den <- sum(parts[, "den"])
  np  <- sum(parts[, "npair"])
  c(acf1 = if (den > 0) num / den else NA_real_, se = 1 / sqrt(max(np, 1)))
}

summarize <- function(label, res, ks_p, acf) {
  r <- res[is.finite(res)]
  cat(sprintf(
    "  %-26s n=%6d  mean %+7.4f  sd %6.4f  KS p<0.05 in %2d/%2d reps  acf1 %+7.4f (se %.4f)  non-finite %d\n",
    label, length(r), mean(r), stats::sd(r), sum(ks_p < 0.05), length(ks_p),
    acf["acf1"], acf["se"], sum(!is.finite(res))))
}

cat("== OSA method = \"cdf\": distributional self-test ==\n")
cat("BS2017SS, single species, fixed effects; ", nrep, " simulation replicates\n\n", sep = "")

fit <- fit_bs()
cat(sprintf("fitted objective %.6f\n", fit$opt$objective))

# The relabelling above is only sound if a mode-3 build reproduces the fitted
# hindcast objective on the SAME data. Check it rather than assume it.
chk <- at_truth(fit, fit$data_list)
cat(sprintf("mode-3 rebuild at the estimates, same data: objective %.6f (fitted %.6f, diff %.3e)\n\n",
            chk$obj$fn(chk$obj$par), fit$opt$objective,
            chk$obj$fn(chk$obj$par) - fit$opt$objective))

configs <- list(
  list(label = "oneStepGaussianOffMode", method = "oneStepGaussianOffMode", discrete = FALSE),
  list(label = "cdf (continuous)",       method = "cdf",                    discrete = FALSE),
  list(label = "cdf + discrete = TRUE",  method = "cdf",                    discrete = TRUE))

for (src in list(c("index", "catch"), "comp")) {
  cat("-- source: ", paste(src, collapse = " + "), " --\n", sep = "")
  acc <- lapply(configs, function(x) list(res = numeric(0), ks = numeric(0)))
  for (rep in seq_len(nrep)) {
    set.seed(1000 + rep)
    dsim <- suppressWarnings(suppressMessages(sim_mod(fit, simulate = TRUE)))
    sf   <- at_truth(fit, dsim)
    for (i in seq_along(configs)) {
      cfg <- configs[[i]]
      o <- tryCatch(suppressWarnings(suppressMessages(osa_residuals(
        sf, source = src, method = cfg$method, discrete = cfg$discrete,
        parallel = FALSE, seed = 2000 + rep))), error = function(e) e)
      if (inherits(o, "error")) {
        cat("    replicate ", rep, " ", cfg$label, " ERROR: ",
            conditionMessage(o), "\n", sep = "")
        next
      }
      r <- o$residual[is.finite(o$residual)]
      acc[[i]]$res <- c(acc[[i]]$res, o$residual)
      acc[[i]]$ks  <- c(acc[[i]]$ks, stats::ks.test(r, "pnorm")$p.value)
      # Accumulate the autocorrelation across EVERY replicate. Keeping the last
      # replicate's object and summarising that reported a one-replicate
      # statistic beside twenty-replicate means and sds, which understates the
      # null standard error by a factor of sqrt(nrep).
      acc[[i]]$acf <- rbind(acc[[i]]$acf, within_acf1_parts(o))
    }
  }
  for (i in seq_along(configs)) {
    a <- acc[[i]]
    summarize(configs[[i]]$label, a$res, a$ks, pool_acf1(a$acf))
  }
  cat("\n")
}

# ---- Agreement where the Gaussian approximation is exact --------------------
# The aggregate series are genuinely normal on the scale obsvec stores them in,
# so the Gaussian method's conditional mode IS the mean and the two methods must
# return the same number. Any disagreement beyond the |residual| ceiling is an
# implementation bug, not a better method.
cat("-- cdf vs oneStepGaussianOffMode on the aggregate (exactly Gaussian) series --\n")
g <- suppressWarnings(suppressMessages(osa_residuals(fit, source = c("index", "catch"), parallel = FALSE)))
c_ <- suppressWarnings(suppressMessages(osa_residuals(fit, source = c("index", "catch"), method = "cdf", parallel = FALSE)))
d <- abs(g$residual - c_$residual)
cat(sprintf("  n = %d   max |diff| = %.3e   at |gaussian residual| = %.3f\n",
            length(d), max(d), abs(g$residual)[which.max(d)]))
cat(sprintf("  max |diff| over rows with |residual| < 8: %.3e\n",
            max(c(0, d[abs(g$residual) < 8]))))
cat("\n")

# ---- The same on compositions, at a sample size that is not sparse ----------
# The Gaussian approximation is good where every bin holds plenty of fish, so the
# two methods should agree closely there and diverge only as the bins empty.
cat("-- composition agreement against composition sample size --\n")
for (n in c(2, 20, 200, 2000)) {
  dat <- fit$data_list
  dat$comp_data$Sample_size <- n
  f <- suppressWarnings(suppressMessages(fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
  og <- suppressWarnings(suppressMessages(osa_residuals(f, source = "comp", parallel = FALSE)))
  oc <- suppressWarnings(suppressMessages(osa_residuals(f, source = "comp", method = "cdf",
                                                        discrete = TRUE, parallel = FALSE)))
  ok <- is.finite(og$residual) & is.finite(oc$residual)
  cat(sprintf("  Sample_size %5d : cor %.4f   median |diff| %.4f   negative predicted (gaussian) %4d/%4d\n",
              n, stats::cor(og$residual[ok], oc$residual[ok]),
              stats::median(abs(og$residual[ok] - oc$residual[ok])),
              sum(is.finite(og$predicted) & og$predicted < 0), nrow(og)))
}

# ---- With random effects ----------------------------------------------------
# method = "cdf" evaluates the objective three times per observation instead of
# once plus a gradient, and each evaluation re-solves the Laplace inner problem.
# This is the check that it finishes and returns finite residuals there. The two
# methods are NOT expected to agree exactly: the population dynamics are
# nonlinear in the recruitment deviations, so the marginal conditional is only
# approximately Gaussian and the approximation is what "cdf" avoids.
#
# The aggregate series only, because the composition loop on a random-effects
# BS2017SS is minutes per method. The composition path under random effects is
# covered instead by test-likelihood-osa-cdf.R ("splitting a call into groups
# conditions on what came before it"), which residualizes compositions under cdf
# on a 21-random-effect fixture.
cat("-- random recruitment (159 random effects), aggregate series --\n")
data(BS2017SS, package = "Rceattle", envir = environment())
fit_re <- suppressWarnings(suppressMessages(fit_mod(
  BS2017SS, file = NULL, estimateMode = 1, msmMode = 0, random_rec = TRUE,
  fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
tg <- system.time(gr <- suppressWarnings(suppressMessages(osa_residuals(
  fit_re, source = c("index", "catch"), parallel = TRUE))))[["elapsed"]]
tc <- system.time(cr <- suppressWarnings(suppressMessages(osa_residuals(
  fit_re, source = c("index", "catch"), method = "cdf", parallel = TRUE))))[["elapsed"]]
cat(sprintf("  n = %d   gaussian sd %.4f (%.0fs)   cdf sd %.4f (%.0fs)   cor %.4f   non-finite %d / %d\n",
            nrow(gr), stats::sd(gr$residual), tg, stats::sd(cr$residual), tc,
            stats::cor(gr$residual, cr$residual),
            sum(!is.finite(gr$residual)), sum(!is.finite(cr$residual))))
cat("\n")

# ---- The same distributional test, WITH random effects -----------------------
# The self-test above is on a fixed-effect model, where "cdf" is exact. It is not
# exact once there are random effects: TMB integrates exp(-h(u))*F(x|u) by
# Laplace, and that integrand is a Gaussian times a sigmoid rather than a
# density. Against an exact Kalman oracle on a linear-Gaussian state space model
# the error runs 7e-4 to 3.6e-2 as the latent state gets more informative, while
# fullGaussian and oneStepGaussian are exact to machine precision. So the
# question is whether it is large enough to overturn the recommendation for
# COMPOSITIONS, where the Gaussian methods are wrong for a bigger reason.
#
# The recruitment deviations are REDRAWN. A one-step-ahead residual marginalizes
# over the random effects, so its null distribution assumes they came from their
# prior; holding them at their fitted values would score the methods against data
# the model did not generate.
cat("-- composition self-test WITH random effects --\n")
source("tests/testthat/helpers.R", local = TRUE)
dat_re <- make_test_data(nyrs = 15, nages = 8, seed = 77)
fit_c <- suppressWarnings(suppressMessages(fit_mod(
  dat_re, file = NULL, estimateMode = 1, msmMode = 0, random_rec = TRUE,
  fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
cat("  random effects:", length(fit_c$obj$env$random), "\n")
acc_re <- lapply(configs, function(x) list(res = numeric(0), ks = numeric(0)))
nrep_re <- nrep * 5L    # this fixture yields ~14 residuals a replicate
cat("  replicates:", nrep_re, "\n")
for (rep in seq_len(nrep_re)) {
  set.seed(5000 + rep)
  ds <- try(suppressWarnings(suppressMessages(
    sim_mod(fit_c, simulate = TRUE, process = "recruitment"))), silent = TRUE)
  if (inherits(ds, "try-error")) next
  sf <- try(suppressWarnings(suppressMessages(fit_mod(
    ds, file = NULL, estimateMode = 3, msmMode = 0, random_rec = TRUE,
    inits = fit_c$estimated_params,
    fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE)))), silent = TRUE)
  if (inherits(sf, "try-error")) next
  sf$data_list$estimateMode <- 1
  for (i in seq_along(configs)) {
    o <- try(suppressWarnings(suppressMessages(osa_residuals(
      sf, source = "comp", method = configs[[i]]$method,
      discrete = configs[[i]]$discrete, parallel = FALSE, seed = 9000 + rep))),
      silent = TRUE)
    if (inherits(o, "try-error")) next
    r <- o$residual[is.finite(o$residual)]
    if (length(r) < 3) next
    acc_re[[i]]$res <- c(acc_re[[i]]$res, r)
    acc_re[[i]]$ks  <- c(acc_re[[i]]$ks, stats::ks.test(r, "pnorm")$p.value)
    acc_re[[i]]$acf <- rbind(acc_re[[i]]$acf, within_acf1_parts(o))
  }
}
for (i in seq_along(configs)) {
  a <- acc_re[[i]]
  summarize(configs[[i]]$label, a$res, a$ks, pool_acf1(a$acf))
}

cat("\ndone.\n")
