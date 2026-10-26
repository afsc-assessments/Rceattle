# tools/verify/verify-sim-centering.R
# Are the simulated observations CENTRED where the likelihood expects? An sd
# check sails straight past a biased mean, and a centring error is a bias no
# number of replicates averages away -- it looks exactly like a real estimation
# failure. So check the location separately from the spread.
#
# ceattle.cpp fits log(obs) ~ N(log(hat) - bias_adjust_obs*sigma^2/2, sigma) for
# both the index (section 13.1, lognormal family) and total catch (section 13.2). The
# likelihood's expectation of log(obs) is therefore log(hat) - ba*sigma^2/2, and
# E[obs] = hat*exp((1-ba)*sigma^2/2) -- equal to hat only when ba = 1. The
# quantity to watch is the log offset against exp(-ba*sigma^2/2); it must not
# drift with sigma.
#
# Use an inflated sd. At the fixture default the offset is within Monte Carlo
# noise, so a real error would pass unnoticed.
#
# Usage:
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript tools/verify/verify-sim-centering.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))
e <- new.env(parent = asNamespace("Rceattle"))
for (f in list.files("tests/testthat", "^helper", full.names = TRUE)) sys.source(f, e)

for (cv in c(0.1, 0.2, 0.4)) {
  d <- e$make_test_data(nyrs = 20, nages = 5, seed = 123)
  d$index_data$Log_sd[d$index_data$Fleet_name == "Survey"] <- cv
  om <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = FALSE))))
  srv <- om$data_list$index_data$Fleet_name == "Survey"
  hat <- om$quantities$index_hat[srv]
  sdv <- om$quantities$log_index_sd[srv]

  set.seed(11)
  reps <- replicate(3000, Rceattle::sim_mod(om, simulate = TRUE)$
                      index_data$Observation[srv])

  cat(sprintf("\n--- nominal Log_sd %.2f ---\n", cv))
  cat(sprintf("  model's log_index_sd            : %.4f\n", sdv[1]))
  cat(sprintf("  mean(sim)/hat                   : %.4f  (1.0000 = mean-unbiased)\n",
              mean(rowMeans(reps) / hat)))
  cat(sprintf("  exp(mean(log(sim)) - log(hat))  : %.4f  (exp(-sd^2/2) = %.4f)\n",
              exp(mean(rowMeans(log(reps)) - log(hat))), exp(-sdv[1]^2 / 2)))
  cat(sprintf("  realised log sd                 : %.4f\n",
              mean(apply(log(reps), 1, stats::sd))))

  ## And what the OM itself thinks: refit the expected data, look at the index fit
  cat(sprintf("  ratio mean(index_obs)/mean(hat) in the ORIGINAL data: %.4f\n",
              mean(om$data_list$index_data$Observation[srv]) / mean(hat)))
}


## ---- Total catch (drawn by the template's SIMULATE block) ------------------
## Same test on the other lognormal observation. The template must use the mean
## the dnorm beside it uses, bias term included; if it does not, the log offset
## below walks away from exp(-ba*sd^2/2) as sd grows. Rows with catch_hat == 0
## (the projection rows clean_data() appends) are excluded -- their draw is
## exactly 0 by construction and carries no information about centring.

cat("\n\n================ total catch ================\n")
for (sd_c in c(0.1, 0.2, 0.4)) {
  d <- e$make_test_data(nyrs = 20, nages = 5, seed = 123)
  d$catch_data$Log_sd <- sd_c
  om <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = FALSE))))
  fitted <- om$quantities$catch_hat > 0
  hat <- om$quantities$catch_hat[fitted]
  sdv <- om$quantities$log_catch_sd[fitted]

  ba <- om$data_list$bias_adjust_obs
  if (is.null(ba)) ba <- om$obj$env$data$bias_adjust_obs
  if (is.null(ba)) ba <- 1
  ba <- as.numeric(ba)[1]

  set.seed(11)
  reps <- replicate(3000, suppressWarnings(
    Rceattle::sim_mod(om, simulate = TRUE)$catch_data$Catch[fitted]))

  cat(sprintf("\n--- nominal catch Log_sd %.2f (bias_adjust_obs %.1f) ---\n",
              sd_c, ba))
  cat(sprintf("  model's log_catch_sd            : %.4f\n", sdv[1]))
  cat(sprintf("  mean(sim)/hat                   : %.4f  (1.0000 = mean-unbiased)\n",
              mean(rowMeans(reps) / hat)))
  cat(sprintf("  exp(mean(log(sim)) - log(hat))  : %.4f  (exp(-ba*sd^2/2) = %.4f)\n",
              exp(mean(rowMeans(log(reps)) - log(hat))), exp(-ba * sdv[1]^2 / 2)))
  cat(sprintf("  realised log sd                 : %.4f\n",
              mean(apply(log(reps), 1, stats::sd))))
  cat(sprintf("  rows drawn: %d fitted, %d at catch_hat == 0 (excluded)\n",
              sum(fitted), sum(!fitted)))
}
