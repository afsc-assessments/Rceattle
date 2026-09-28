# Simulation recovery of a Beverton-Holt curve under predation (msmMode > 0).
#
# Spawning biomass per recruit is undefined under predation, so nothing ties
# alpha and beta to an equilibrium: a curve the data do not inform runs to a
# flat or a linear ridge. This harness asks whether the data CAN inform one.
# The two-species fixture is fished with a strong pulse (SSB spans about 7x
# for species 1, 5x for species 2), a Beverton-Holt curve with its bend inside
# that range is imposed on a mean-recruitment fit, every observation and the
# recruitment process are redrawn with sim_mod(), and the curve is refitted in
# the hindcast. Reported per species: log alpha, log beta and the asymptote
# alpha/beta, as mean, empirical SD and mean reported SE.
#
# The refit starts at the truth by default, which tests local recovery only.
# `start = dispersed` starts log alpha and log beta at truth +/- 1 (alternating
# by replicate), which is what can reach the flat or linear ridge the
# convergence check reports.
#
# Usage: Rscript tools/verify/verify-sim-recovery-srr-msm.R [n_reps] [seed_offset] [truth|dispersed]
# Run time: about a minute per replicate.

args   <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 30L
seed0  <- if (length(args) >= 2) as.integer(args[2]) else 3000L
start  <- if (length(args) >= 3) args[3] else "truth"

suppressMessages(pkgload::load_all(".", compile = FALSE, quiet = TRUE))
e <- new.env(parent = asNamespace("Rceattle"))
sys.source("tests/testthat/helpers-make-msm-data.R", envir = e)

set.seed(123)
Fm <- matrix(c(seq(0.02, 0.9, length.out = 15), seq(0.9, 0.05, length.out = 15),
               seq(0.02, 0.6, length.out = 30)), 2, 30, byrow = TRUE)
d <- e$make_msm_test_data(Fmort = Fm)$data_list

qfit <- function(dat, recFun, mode = "Hindcast", inits = NULL, sd = FALSE) {
  suppressMessages(suppressWarnings(fit_mod(
    data_list = dat, inits = inits, estimateMode = mode, msmMode = 1, suitMode = 0,
    initMode = "NonEquilibrium", random_rec = FALSE, recFun = recFun,
    fit_control = fit_control(phase = FALSE, loopnum = 3, verbose = 0, getsd = sd))))
}

# Operating model: the mean-recruitment fit with a curve imposed whose bend
# (predicted/asymptote = 0.5) sits at the median hindcast SSB, and whose
# recruitment at that SSB equals the fitted mean recruitment.
base <- qfit(d, build_srr(srr_fun = "mean"))
nyh  <- d$endyr - d$styr + 1
ssb  <- base$quantities$ssb[, 1:nyh]
Rfit <- base$quantities$R[, 1:nyh]
s_med <- apply(ssb, 1, stats::median)
beta_true  <- 1 / s_med
alpha_true <- 2 * rowMeans(Rfit) / s_med
bh <- build_srr(srr_fun = "BevertonHolt", srr_alpha_init = alpha_true, srr_beta_init = beta_true)
ini <- base$estimated_params
ini$rec_pars[, 2] <- log(alpha_true)
ini$rec_pars[, 3] <- log(beta_true)
om  <- qfit(d, bh, mode = "DebugBuild", inits = ini)
cat(sprintf("true log alpha %s, log beta %s, asymptote %s\n",
            paste(round(log(alpha_true), 3), collapse = "/"),
            paste(round(log(beta_true), 3), collapse = "/"),
            paste(round(alpha_true / beta_true, 1), collapse = "/")))

res <- data.frame()
for (r in seq_len(n_reps)) {
  set.seed(seed0 + r)
  sim <- sim_mod(om, simulate = TRUE, process = "recruitment")
  # build_srr()'s srr_*_init override supplied inits, so the start is set there.
  shift  <- if (start == "dispersed") (-1)^r else 0
  bh_em  <- build_srr(srr_fun = "BevertonHolt",
                      srr_alpha_init = alpha_true * exp(shift),
                      srr_beta_init  = beta_true * exp(-shift))
  em  <- qfit(sim, bh_em, inits = ini, sd = TRUE)
  sm  <- tryCatch(summary(em$sdrep), error = function(e) NULL)
  rp  <- em$estimated_params$rec_pars
  se  <- if (is.matrix(sm)) sm[rownames(sm) == "rec_pars", 2] else rep(NA_real_, 3 * nrow(em$estimated_params$rec_pars))
  # rec_pars is [nspp, 3] column-major in obj$par: R0 (all species), alpha, beta.
  nspp <- nrow(rp)
  row <- data.frame(rep = r, maxgr = max(abs(em$obj$gr())), obj = em$opt$objective)
  for (sp in seq_len(nspp)) {
    row[[paste0("la", sp)]] <- rp[sp, 2]; row[[paste0("lb", sp)]] <- rp[sp, 3]
    row[[paste0("se_la", sp)]] <- se[nspp + sp]; row[[paste0("se_lb", sp)]] <- se[2 * nspp + sp]
  }
  res <- rbind(res, row)
  cat(sprintf("rep %2d  max|grad| %.1e  log alpha %s  log beta %s\n", r, row$maxgr,
              paste(round(rp[, 2], 3), collapse = "/"), paste(round(rp[, 3], 3), collapse = "/")))
}
report <- function(set, label) {
  cat(sprintf("\n%s (%d replicates)\n", label, nrow(set)))
  for (sp in seq_len(nrow(om$estimated_params$rec_pars))) {
    la <- set[[paste0("la", sp)]]; lb <- set[[paste0("lb", sp)]]
    cat(sprintf("species %d  log alpha: true %.3f mean %.3f sd %.3f mean SE %.3f | log beta: true %.3f mean %.3f sd %.3f mean SE %.3f | asymptote: true %.1f mean %.1f\n",
                sp, log(alpha_true[sp]), mean(la), sd(la), mean(set[[paste0("se_la", sp)]], na.rm = TRUE),
                log(beta_true[sp]), mean(lb), sd(lb), mean(set[[paste0("se_lb", sp)]], na.rm = TRUE),
                alpha_true[sp] / beta_true[sp], mean(exp(la - lb))))
  }
}
# The fixture's refits stop with gradients near 1e-2 without phasing, so both
# the whole set and the better-converged subset are reported.
report(res, "all replicates")
report(res[res$maxgr < 5e-2, ], "max|grad| < 5e-2")
