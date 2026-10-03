# Cold-start convergence for the four golden references.
#
#   export PATH=/usr/bin:$PATH
#   NOT_CRAN=true Rscript tools/verify/verify-golden-cold-start.R
#
# This is the half of the old golden check that test-golden-regression.R no
# longer does. It is a harness rather than a test because its failure mode is
# platform-dependent and the golden job must not be: `goa_ss` has a second local
# minimum 52.9 units above the reference, and a one-ULP change in one `log_F`
# gradient element is enough for `nlminb` to reach it. A cold fit reached it at
# 5.48.0 (12920.1030998153 at max|gradient| 5.7e-11, so the higher minimum
# polishes as well as the lower one) -- but on a HALF-APPLIED length-comp fix,
# and an incoherent likelihood is its own way in. On coherent code the cold path
# holds: the preceding release reproduces the reference exactly and the completed
# fix lands in the same basin. `goa_ms` warm-started from those upper-basin MLEs
# still reached its own reference objective, so it does not always inherit the
# basin. Which minimum a cold fit finds is the flake this was kept out of CI to
# avoid; the regeneration script pins the lower of two starts for the same
# reason.
#
# What it does assert is the direction that cannot flake: a cold fit must not
# land BELOW the pinned reference. A lower objective is not a basin lottery, it
# is a likelihood that changed.
#
# Run it when touching the optimizer, the phases, the starting values or the
# bounds -- the class test-golden-regression.R stopped covering when it stopped
# re-optimizing. Expect ~3-4 minutes locally.

devtools::load_all(".", quiet = TRUE)

# The same literals test-golden-regression.R pins.
ref <- c(ss     = 10241.0304272585,
         ms     = 10267.2478324443,
         goa_ss = 12866.8457276232,
         goa_ms = 12931.8602763619)

fc <- function(...) Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                          newtonsteps = 3, ...)

ss <- Rceattle::fit_mod(
  data_list = Rceattle::BS2017SS, file = NULL, inits = NULL, estimateMode = 0,
  random_rec = FALSE, msmMode = 0, fit_control = fc(phase = TRUE))
goa_ss <- Rceattle::fit_mod(
  data_list = Rceattle::GOA2018SS, file = NULL, inits = NULL, estimateMode = 0,
  random_rec = FALSE, msmMode = 0, fit_control = fc(phase = TRUE))
ms <- Rceattle::fit_mod(
  data_list = Rceattle::BS2017MS, inits = ss$estimated_params, file = NULL,
  estimateMode = 0, niter = 5, random_rec = FALSE, msmMode = 1, suitMode = 0,
  fit_control = fc())
goa_ms <- Rceattle::fit_mod(
  data_list = Rceattle::GOA2018SS, inits = goa_ss$estimated_params, file = NULL,
  estimateMode = 0, niter = 3, random_rec = FALSE, msmMode = 1, suitMode = 0,
  fit_control = fc(phase = TRUE))

fits <- list(ss = ss, ms = ms, goa_ss = goa_ss, goa_ms = goa_ms)
bad  <- character(0)

for (m in names(fits)) {
  f     <- fits[[m]]
  grad  <- max(abs(f$obj$gr(f$opt$par)))
  delta <- f$opt$objective - ref[[m]]

  # A cold fit BELOW the reference means the likelihood changed. Above it may
  # only mean this machine took the other basin, so it is reported, not failed.
  if (delta < -1e-6)
    bad <- c(bad, sprintf("%s reached %.10f, BELOW the pinned %.10f by %.3g",
                          m, f$opt$objective, ref[[m]], -delta))
  if (!is.finite(grad) || grad > 1e-4)
    message(sprintf("  NOTE %s did not polish to a stationary point: max|gradient| = %.3g",
                    m, grad))

  message(sprintf("%-7s objective %.10f  (%+.3g vs reference)  max|gradient| %.3g%s",
                  m, f$opt$objective, delta, grad,
                  if (abs(delta) > 1) "   <- a different basin" else ""))
}

if (length(bad)) {
  stop("cold start reached a lower objective than the reference:\n  ",
       paste(bad, collapse = "\n  "), call. = FALSE)
}
message("\nNo cold fit landed below its reference. A fit ~52.9 above on goa_ss ",
        "is the documented second minimum, not a regression -- see ",
        "inst/dev/TRAPS.md.")
