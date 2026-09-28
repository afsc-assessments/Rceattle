# Regenerate tests/testthat/fixtures/golden-reference.rds
#
# The golden check asks one question -- is the likelihood still the one the
# reference was fit to -- by evaluating it AT the reference parameters rather
# than by re-optimizing. That needs the parameters themselves on disk, which is
# what this script writes.
#
# Run it only when a model change is INTENDED, and say in the commit which
# reference moved and by how much. Re-running it after an unintended change
# launders the very regression the check exists to catch, exactly as editing the
# pinned objectives by hand would.
#
#   Rscript tools/verify/regenerate-golden-reference.R
#
# The four fits are the ones test-golden-regression.R pins: BS2017SS single and
# multispecies, GOA2018SS single and multispecies, each Newton-polished at
# newtonsteps = 3 so the reference is a stationary point and not wherever
# nlminb's relative tolerance happened to bite. Takes a few minutes.

devtools::load_all(".", quiet = TRUE)

fc <- function(...) Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                          newtonsteps = 3, ...)

ss <- Rceattle::fit_mod(
  data_list = Rceattle::BS2017SS, file = NULL, inits = NULL, estimateMode = 0,
  random_rec = FALSE, msmMode = 0, fit_control = fc(phase = TRUE))
ms <- Rceattle::fit_mod(
  data_list = Rceattle::BS2017MS, inits = ss$estimated_params, file = NULL,
  estimateMode = 0, niter = 5, random_rec = FALSE, msmMode = 1, suitMode = 0,
  fit_control = fc())
goa_ss <- Rceattle::fit_mod(
  data_list = Rceattle::GOA2018SS, file = NULL, inits = NULL, estimateMode = 0,
  random_rec = FALSE, msmMode = 0, fit_control = fc(phase = TRUE))
goa_ms <- Rceattle::fit_mod(
  data_list = Rceattle::GOA2018SS, inits = goa_ss$estimated_params, file = NULL,
  estimateMode = 0, niter = 3, random_rec = FALSE, msmMode = 1, suitMode = 0,
  fit_control = fc(phase = TRUE))

fits <- list(ss = ss, ms = ms, goa_ss = goa_ss, goa_ms = goa_ms)

# Refuse to write a reference that is not a stationary point: the whole check
# rests on the gradient at these parameters being ~0, so a fit that stopped
# early would pin a number the next run cannot reproduce for the right reason.
for (m in names(fits)) {
  g <- max(abs(fits[[m]]$obj$gr(fits[[m]]$opt$par)))
  if (!is.finite(g) || g > 1e-4)
    stop(sprintf("%s did not converge: max|gradient| = %.3g", m, g))
  message(sprintf("%-7s objective %.10f  max|gradient| %.3g",
                  m, fits[[m]]$opt$objective, g))
}

ref <- lapply(fits, function(f) list(
  objective = as.numeric(f$opt$objective),
  params    = f$estimated_params))

out <- "tests/testthat/fixtures/golden-reference.rds"
saveRDS(ref, out, version = 2)
message("wrote ", out, " (", round(file.size(out) / 1024), " KB)")
message("Rceattle ", as.character(utils::packageVersion("Rceattle")),
        " on ", R.version.string)
