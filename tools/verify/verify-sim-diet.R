# tools/verify/verify-sim-diet.R
# Are simulated stomach contents drawn from the distribution the diet likelihood
# assumes? Diet was the last observation type still not simulated, so there is no
# earlier implementation to compare against -- the moments are the whole check.
#
# Two things this harness exists to get right, both of which look like bugs when
# measured the wrong way:
#
# 1. The draw is centred on the proportions the DENSITY uses, not on diet_hat.
#    The likelihood appends an "other prey" bin, adds comp_offset and
#    renormalizes before taking the density. Against diet_hat, a bin whose
#    prediction is near the offset shows an apparent 2x bias that is not there.
#
# 2. The Dirichlet-multinomial must be checked at a SMALL concentration. This
#    model fits diet_comp_weights to about 16 (theta ~ 1e7), where the
#    Dirichlet-multinomial is numerically a multinomial and its expected sd
#    inflation is 1.0000 -- a comparison run there passes whatever the simulator
#    does. Every parameter is therefore held at the fitted values and only theta
#    is varied, so the sd inflation is the single thing that moves.
#
# Expected: realised/expected sd ~ 1 in all three cases, and centring sum ratio
# ~ 1. The Dirichlet-multinomial at theta = 0.1 carries an expected inflation of
# about 2.7, which is what gives the check its power.
#
# Usage:
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript tools/verify/verify-sim-diet.R

pkgload::load_all(".", quiet = TRUE)
suppressPackageStartupMessages(library(Rceattle))
data(BS2017MS); data(BS2017SS)
ss <- Rceattle::fit_mod(BS2017SS, inits=NULL, msmMode=0, estimateMode=1, phase=FALSE, getsd=FALSE, verbose=0)

# Fit once with an estimated suitability so diet_hat is concentrated.
mk <- function(fam) { d <- BS2017MS; d$Diet_distribution <- rep(fam, d$nspp); d }
fit1 <- Rceattle::fit_mod(mk("DirichletMultinomial"), inits=ss$estimated_params,
                          msmMode=1, suitMode=4, estimateMode=1, niter=3,
                          phase=FALSE, getsd=FALSE, verbose=0)

# Hold every parameter at that fit and vary only the DM concentration, so
# diet_hat is unchanged and theta is the single thing that moves.
rebuild <- function(fam, w) {
  inits <- fit1$estimated_params
  inits$diet_comp_weights[] <- w
  Rceattle::fit_mod(mk(fam), inits=inits, msmMode=1, suitMode=4, estimateMode=3,
                    niter=3, phase=FALSE, getsd=FALSE, verbose=0)
}

target_props <- function(mod) {
  re <- Rceattle::rearrange_data(mod$data_list); sid <- re$stomach_id
  hat <- mod$quantities$diet_hat[, 2]; out <- rep(NA_real_, length(hat))
  for (i in unique(sid)) {
    idx <- which(sid == i); p <- c(hat[idx], 1 - sum(hat[idx]))
    p <- (p + 1e-5)/sum(p + 1e-5); out[idx] <- p[seq_along(idx)]
  }
  out
}

N <- 20
# The weight is on each family's own scale (ceattle.cpp, section 13.6): the
# multinomial multiplies N by it raw, the Dirichlet-multinomial uses theta = exp(w).
for (spec in list(list(f="Multinomial",          w=1),
                  list(f="DirichletMultinomial", w=log(0.1)),
                  list(f="DirichletMultinomial", w=log(1)))) {
  m <- rebuild(spec$f, spec$w)
  mn <- spec$f == "Multinomial"
  theta <- exp(spec$w); phi <- N*theta
  infl <- if (mn) 1 else sqrt((N + phi)/(1 + phi))
  tgt <- target_props(m)
  big <- which(tgt > 0.05)
  set.seed(7)
  reps <- replicate(400, Rceattle::sim_mod(m, simulate=TRUE)$diet_data$Stomach_proportion_by_weight)
  sd_emp <- apply(reps[big, , drop=FALSE], 1, stats::sd)
  sd_exp <- sqrt(tgt[big]*(1-tgt[big])/N) * infl
  cat(sprintf("\n%-22s %-12s expected sd inflation %.3f (%d bins)\n", spec$f,
              if (mn) paste0("weight=", spec$w) else paste0("theta=", format(theta, digits=3)),
              infl, length(big)))
  cat(sprintf("  realised/expected sd : %.4f\n", mean(sd_emp/sd_exp)))
  cat(sprintf("  centring sum ratio   : %.4f\n",
              sum(rowMeans(reps)[tgt>1e-4])/sum(tgt[tgt>1e-4])))
}
