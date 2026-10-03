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
# The references are the ones test-golden-regression.R pins: BS2017SS single and
# multispecies, GOA2018SS single and multispecies, each Newton-polished at
# newtonsteps = 3 so the reference is a stationary point and not wherever
# nlminb's relative tolerance happened to bite. Each is fitted twice (see
# below), so expect roughly ten minutes.

devtools::load_all(".", quiet = TRUE)

out <- "tests/testthat/fixtures/golden-reference.rds"

# Each reference is fitted from TWO starts and the LOWER objective is pinned.
# The first is the recipe that created the references -- a phased fit from the
# build defaults for the single-species pair, and their MLEs for the
# multispecies pair, whose predation likelihood is non-convex so that start is
# part of what the reference means. The second, on a re-pin, is the committed
# reference itself.
#
# Both are needed. `goa_ss` has a second local minimum ~52 nats up and a one-ULP
# change in one `log_F` gradient element is enough to send a phased fit there,
# so pinning the recipe alone can move a reference to a worse optimum and report
# it as a model change (it did, at 5.48.0, on a half-applied length-comp fix).
# But pinning the reference's own basin alone would hide the opposite case: a
# change that moves the OPTIMIZER, on the cold path every assessment script
# takes. Taking the lower of the two cannot do either.
prev <- if (file.exists(out)) readRDS(out) else NULL

fc <- function(...) Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                          newtonsteps = 3, ...)

# Phasing walks a cold start into a basin; from the reference the parameters are
# already in one, so that fit estimates them all at once.
lower <- function(m, recipe, from_ref) {
  if (is.null(from_ref)) return(recipe)
  message(sprintf("%-7s recipe %.10f   from reference %.10f", m,
                  recipe$opt$objective, from_ref$opt$objective))
  if (recipe$opt$objective <= from_ref$opt$objective) recipe else from_ref
}
from_ref <- function(m, f) if (is.null(prev)) NULL else f(prev[[m]][["params"]])

bs_ss <- function(inits, phase) Rceattle::fit_mod(
  data_list = Rceattle::BS2017SS, file = NULL, inits = inits,
  estimateMode = 0, random_rec = FALSE, msmMode = 0,
  fit_control = fc(phase = phase))
bs_ms <- function(inits) Rceattle::fit_mod(
  data_list = Rceattle::BS2017MS, file = NULL, inits = inits,
  estimateMode = 0, niter = 5, random_rec = FALSE, msmMode = 1, suitMode = 0,
  fit_control = fc())
goa <- function(inits, phase, msm, niter) Rceattle::fit_mod(
  data_list = Rceattle::GOA2018SS, file = NULL, inits = inits,
  estimateMode = 0, niter = niter, random_rec = FALSE, msmMode = msm,
  suitMode = 0, fit_control = fc(phase = phase))

ss <- lower("ss", bs_ss(NULL, TRUE),
            from_ref("ss", function(p) bs_ss(p, FALSE)))
ms <- lower("ms", bs_ms(ss$estimated_params), from_ref("ms", bs_ms))
goa_ss <- lower("goa_ss", goa(NULL, TRUE, 0, 3),
                from_ref("goa_ss", function(p) goa(p, FALSE, 0, 3)))
goa_ms <- lower("goa_ms", goa(goa_ss$estimated_params, TRUE, 1, 3),
                from_ref("goa_ms", function(p) goa(p, FALSE, 1, 3)))

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
            goa_ss = 12866.8457276232, goa_ms = 12931.8602763619)
for (m in names(obj_at))
  message(sprintf("%-7s vs the literal in test-golden-regression.R: %+.3e", m,
                  obj_at[m] - pinned[[m]]))

# The fixture is a binary, so the parameters move invisibly in a diff. Say how
# far, and in which block: a move far larger than the objective gain implies is
# a flat ridge or a different basin rather than a re-polished optimum.
if (!is.null(prev)) {
  for (m in names(fits)) {
    was <- prev[[m]][["params"]]
    now <- polished(fits[[m]])
    # A release that adds a parameter block leaves the old reference without it,
    # so compare the blocks the two share and name any that are new.
    gained <- setdiff(names(now), names(was))
    shared <- intersect(names(now), names(was))
    shared <- shared[lengths(now[shared]) == lengths(was[shared])]
    d <- vapply(shared, function(k) max(abs(unlist(now[[k]]) - unlist(was[[k]]))),
                numeric(1))
    message(sprintf("%-7s max|dparam| %.3e in %s%s", m, max(d),
                    names(which.max(d)),
                    if (length(gained))
                      paste0("   new block(s): ", paste(gained, collapse = ", "))
                    else ""))
  }
}

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
