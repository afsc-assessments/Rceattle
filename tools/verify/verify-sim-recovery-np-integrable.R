# Simulation recovery of the selectivity-deviation SD under NonParametricIntegrable.
#
# The form exists so that random_sel = TRUE integrates a complete density and
# the reported SD means the SD of the deviations. This harness fits Atka2022's
# fishery as NonParametricIntegrable with the SD held at 0.35, draws iid deviates at a
# known SD for the estimated coefficient bins (sim_mod() has no draw for the
# Time_varying_sel deviates yet; see inst/dev/CLEANUP_BACKLOG.md), simulates the
# observations from that operating model, and refits with random_sel = TRUE
# from an SD of 0.2. Reported: the estimated SD's mean, empirical SD and mean
# reported SE across replicates.
#
# Measured on 5.40.0: the estimate is biased LOW, 0.24 to 0.26 at a true SD
# of 0.35 (every replicate below the truth), 13% at 0.70, and 9% at 0.35 with the
# input sample sizes multiplied by ten. The density is complete and every
# estimated cell is scored; this is the Laplace marginal likelihood's known
# downward bias for a variance component on small multinomial samples, so read
# the reported SD as a lower bound unless the compositions are well sampled.
#
# Usage: Rscript tools/verify/verify-sim-recovery-np-integrable.R [n_reps] [true_sd] [seed_offset]
# Run time: two to three minutes per replicate.

args    <- commandArgs(trailingOnly = TRUE)
n_reps  <- if (length(args) >= 1) as.integer(args[1]) else 12L
true_sd <- if (length(args) >= 2) as.numeric(args[2]) else 0.35
seed0   <- if (length(args) >= 3) as.integer(args[3]) else 5000L

suppressMessages(pkgload::load_all(".", compile = FALSE, quiet = TRUE))
d <- Rceattle::Atka2022
d$fleet_control$Selectivity[2]      <- "NonParametricIntegrable"
d$fleet_control$Time_varying_sel[2] <- "IID"

qfit <- function(dat, inits = NULL, mode = "Hindcast", random_sel = FALSE, sd = FALSE, phase = FALSE) {
  suppressMessages(suppressWarnings(fit_mod(
    data_list = dat, inits = inits, msmMode = 0, estimateMode = mode, random_sel = random_sel,
    fit_control = fit_control(phase = phase, verbose = 0, getsd = sd))))
}

base  <- qfit(d, phase = TRUE)
cells <- which(!is.na(base$map$mapList$sel_coff_dev[2, 1, , ]))
cat(sprintf("base fit: objective %.4f, max|grad| %.1e, %d deviate cells\n",
            base$opt$objective, max(abs(base$obj$gr())), length(cells)))

res <- data.frame()
for (r in seq_len(n_reps)) {
  set.seed(seed0 + r)
  ini <- base$estimated_params
  ini$sel_coff_dev[2, 1, , ][cells] <- stats::rnorm(length(cells), 0, true_sd)
  om  <- qfit(d, inits = ini, mode = "DebugBuild")
  sim <- sim_mod(om, simulate = TRUE)
  ini$sel_dev_log_sd[2] <- log(0.2)
  em  <- qfit(sim, inits = ini, random_sel = TRUE, sd = TRUE)
  sm  <- tryCatch(summary(em$sdrep), error = function(e) NULL)
  se  <- if (is.matrix(sm) && "sel_dev_log_sd" %in% rownames(sm)) sm[rownames(sm) == "sel_dev_log_sd", 2][1] else NA_real_
  sd_hat <- exp(em$estimated_params$sel_dev_log_sd[2])
  res <- rbind(res, data.frame(rep = r, sd_hat = sd_hat, log_sd = log(sd_hat), se_log_sd = as.numeric(se),
                               maxgr = max(abs(em$obj$gr()))))
  cat(sprintf("rep %2d  sd_hat %.3f  se(log sd) %.3f  max|grad| %.1e\n", r, sd_hat, res$se_log_sd[r], res$maxgr[r]))
}
ok <- res[res$maxgr < 1e-2, ]
cat(sprintf("\n%d of %d converged. true sd %.3f  mean %.3f  empirical sd %.3f  |  log scale: mean %.3f (true %.3f), sd %.3f, mean SE %.3f\n",
            nrow(ok), n_reps, true_sd, mean(ok$sd_hat), sd(ok$sd_hat),
            mean(ok$log_sd), log(true_sd), sd(ok$log_sd), mean(ok$se_log_sd, na.rm = TRUE)))
