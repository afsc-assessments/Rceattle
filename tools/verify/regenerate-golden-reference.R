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
# pinned objectives by hand would -- so overwriting an existing fixture needs
# RCEATTLE_GOLDEN_REPIN=true, and the objectives stay as literals in
# test-golden-regression.R where a reviewer reads them in the diff.
#
#   RCEATTLE_GOLDEN_REPIN=true Rscript tools/verify/regenerate-golden-reference.R
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

# The parameters are taken at `opt$par`, the Newton-POLISHED point, and not from
# `fit$estimated_params`. Those are two different points: `newtonsteps` moves
# `opt$par` and `opt$objective` only, while `estimated_params` and `quantities`
# stay at TMB's `last.par.best`. On goa_ss they differ by 2.3e-08 in parameter
# space, which is 9.1e-12 in objective and -- the part that matters here -- a
# gradient of 3.2e-06 rather than 1.3e-11. Pinning the pair from one point makes
# the check reproduce exactly and leaves the gradient gate its full margin.
polished <- function(f) f$obj$env$parList(x = f$opt$par,
                                          par = f$obj$env$last.par.best)

# Refuse to write a reference that is not a stationary point: the whole check
# rests on the gradient AT THESE PARAMETERS being ~0, so a fit that stopped early
# would pin a number the next run cannot reproduce for the right reason.
obj_at <- numeric(0)
for (m in names(fits)) {
  f <- fits[[m]]
  g <- max(abs(f$obj$gr(f$opt$par)))
  if (!is.finite(g) || g > 1e-8)
    stop(sprintf("%s is not a stationary point: max|gradient| = %.3g", m, g))
  obj_at[m] <- as.numeric(f$obj$fn(f$opt$par))
  message(sprintf("%-7s objective %.10f  max|gradient| %.3g", m, obj_at[m], g))
}

# Say what moved, against the literals the test pins, so a re-pin cannot be done
# without seeing the size of it.
pinned <- c(ss = 10241.0304272585, ms = 10267.2478324443,
            goa_ss = 12867.9902664788, goa_ms = 12932.7902167145)
for (m in names(obj_at))
  message(sprintf("%-7s vs the literal in test-golden-regression.R: %+.3e", m,
                  obj_at[m] - pinned[[m]]))

out <- "tests/testthat/fixtures/golden-reference.rds"
if (file.exists(out) && !identical(Sys.getenv("RCEATTLE_GOLDEN_REPIN"), "true"))
  stop(out, " exists. Re-pinning is not a routine step: set ",
       "RCEATTLE_GOLDEN_REPIN=true only when the model change above is intended, ",
       "and update the literals in tests/testthat/test-golden-regression.R to the ",
       "objectives printed here, in the same commit.", call. = FALSE)

# Parameters only. The objectives stay as literals in the test file, where a
# reviewer sees a re-pin as changed text rather than as a changed binary.
ref <- lapply(fits, function(f) list(params = polished(f)))
attr(ref, "provenance") <- list(
  rceattle = as.character(utils::packageVersion("Rceattle")),
  commit   = tryCatch(system("git rev-parse --short HEAD", intern = TRUE),
                      error = function(e) NA_character_),
  r        = R.version.string,
  platform = R.version$platform,
  date     = format(Sys.Date()))

saveRDS(ref, out, version = 2)
message("wrote ", out, " (", round(file.size(out) / 1024), " KB)")
message(paste(utils::capture.output(str(attr(ref, "provenance"))), collapse = "\n"))
