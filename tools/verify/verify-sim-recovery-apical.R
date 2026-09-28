# Simulation recovery of the per-sex apical selectivity offset (5.38.0).
#
# The offset is informed only by joint compositions (comp_data$Sex = 3), so its
# estimate rests on the sex ratio in the catch. This harness sets a known male
# offset on GOAatf's fishery (Logistic, across-sex normalization), simulates
# every observation with sim_mod(), refits from the truth, and reports the mean
# estimate, the empirical SD and the mean reported SE. A mean more than two
# standard errors of the mean from the truth is the number to investigate.
#
# Usage: Rscript tools/verify/verify-sim-recovery-apical.R [n_reps] [true_log_offset] [seed_offset]
# Run time: about 30 s per replicate.

args <- commandArgs(trailingOnly = TRUE)
n_reps  <- if (length(args) >= 1) as.integer(args[1]) else 40L
true_ap <- if (length(args) >= 2) as.numeric(args[2]) else log(2)
seed0   <- if (length(args) >= 3) as.integer(args[3]) else 1000L

suppressMessages(pkgload::load_all(".", compile = FALSE, quiet = TRUE))
d <- Rceattle::GOAatf
d$fleet_control$Selectivity_index[3] <- 3L
d$fleet_control$Sel_norm_scope[3]    <- "AcrossSexes"
d$fleet_control$Sel_norm_bin[3]      <- "Max"
spec <- build_selectivity(linkages = list(
  apical = linkage_spec(~ 1, by = ~ fleet + sex, fleet = 3L, sex = "male")))

qfit <- function(dat, inits = NULL, mode = "Hindcast", sd = FALSE) {
  suppressMessages(suppressWarnings(fit_mod(
    data_list = dat, inits = inits, msmMode = 0, estimateMode = mode, selFun = spec,
    fit_control = fit_control(phase = FALSE, getsd = sd, verbose = 0))))
}

# Operating model: the no-offset MLE with the offset imposed.
base <- qfit(d)
ini  <- base$estimated_params
ini$log_sel_apical[3, 2] <- true_ap
om   <- qfit(d, inits = ini, mode = "DebugBuild")

res <- data.frame()
for (s in seq_len(n_reps)) {
  set.seed(seed0 + s)
  sim <- sim_mod(om, simulate = TRUE)
  em  <- qfit(sim, inits = ini, sd = TRUE)
  sm  <- tryCatch(summary(em$sdrep), error = function(e) NULL)
  se  <- if (!is.null(sm) && "log_sel_apical" %in% rownames(sm))
    sm[rownames(sm) == "log_sel_apical", 2][1] else NA_real_
  res <- rbind(res, data.frame(rep = s, est = em$estimated_params$log_sel_apical[3, 2],
                               se = as.numeric(se), maxgr = max(abs(em$obj$gr()))))
  cat(sprintf("rep %2d  est %.3f  se %.3f  max|grad| %.1e\n", s, res$est[s], res$se[s], res$maxgr[s]))
}
ok <- res[res$maxgr < 1e-2, ]
cat(sprintf("\n%d of %d converged. true %.4f  mean %.4f  empirical sd %.4f  mean SE %.4f\n",
            nrow(ok), n_reps, true_ap, mean(ok$est), sd(ok$est), mean(ok$se, na.rm = TRUE)))
cat(sprintf("z of the mean: %.2f  (share below truth %.2f)\n",
            (mean(ok$est) - true_ap) / (sd(ok$est) / sqrt(nrow(ok))), mean(ok$est < true_ap)))
