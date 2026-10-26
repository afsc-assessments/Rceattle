# Regenerate tests/testthat/fixtures/golden-reference.rds
#
# The golden check asks one question -- is the likelihood still the one the
# reference was fit to -- by evaluating it AT the reference parameters rather
# than by re-optimizing. That needs the parameters themselves on disk, which is
# what this script writes, with the SSB, biomass, recruitment, unfished and
# dynamic B0, and predation mortality the test also pins at them.
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

# The quantities that set catch advice, evaluated at the polished parameters by
# the same call test-golden-regression.R makes. SSB, biomass, SB0, B0 and dynamic
# B0 never enter the objective under mean recruitment, so its pin cannot see them.
# Names and dimnames are dropped so a relabelled output does not read as a move.
pinned_derived <- c("ssb", "biomass", "R", "SB0", "B0", "DynamicSB0",
                    "DynamicB0", "M2_at_age")
eval_spec <- list(
  ss     = list(data = Rceattle::BS2017SS,  msmMode = 0, niter = 3),
  ms     = list(data = Rceattle::BS2017MS,  msmMode = 1, niter = 5),
  goa_ss = list(data = Rceattle::GOA2018SS, msmMode = 0, niter = 3),
  goa_ms = list(data = Rceattle::GOA2018SS, msmMode = 1, niter = 3))
derived_at <- function(m, params) {
  s <- eval_spec[[m]]
  q <- Rceattle::fit_mod(
    data_list = s$data, inits = params, file = NULL,
    estimateMode = 3, niter = s$niter, random_rec = FALSE,
    msmMode = s$msmMode, suitMode = 0,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))$quantities
  stopifnot(all(pinned_derived %in% names(q)))
  lapply(q[pinned_derived], unname)
}

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
    # Compare only the cells that are FINITE on both sides, and report every
    # other kind of move separately. 7 of the 49 blocks hold a non-finite value
    # in both fixtures: log_M1 past a species' nages or in the male slot of a
    # one-sex species is NA, index_log_sd and catch_log_sd are NA throughout,
    # and sel_dev_log_sd is -Inf in 7 of 7 cells on `ss`. A bare max() therefore
    # returned NA for ALL FOUR models -- NA from a missing cell, and NaN from
    # Inf minus Inf, which an na.rm would not have caught either -- so this
    # report, whose whole job is to make a parameter move visible BEFORE a
    # re-pin, measured nothing. That is how a 1.2730 move in goa_ss log_M1 went
    # unreported: the cells are inert, outside (nsex, nages), but the guard
    # could not say so.
    d <- vapply(shared, function(k) {
      a <- unlist(now[[k]]); b <- unlist(was[[k]])
      ok <- is.finite(a) & is.finite(b)
      if (!any(ok)) 0 else max(abs(a[ok] - b[ok]))
    }, numeric(1))
    # A cell moving between finite, NA and +/-Inf is a change no difference can
    # express, so count those rather than letting them vanish.
    kind <- function(x) ifelse(is.na(x), "NA",
                               ifelse(is.finite(x), "finite",
                                      paste0("Inf", sign(x))))
    moved <- vapply(shared, function(k)
      sum(kind(unlist(now[[k]])) != kind(unlist(was[[k]]))), integer(1))
    message(sprintf("%-7s max|dparam| %.3e in %s%s%s", m, max(d),
                    names(which.max(d)),
                    if (any(moved > 0))
                      sprintf("   finite/NA/Inf pattern moved in %s (%d cell(s))",
                              paste(names(moved)[moved > 0], collapse = ", "),
                              sum(moved))
                    else "",
                    if (length(gained))
                      paste0("   new block(s): ", paste(gained, collapse = ", "))
                    else ""))
  }
}

# The derived references, and how far each moved from the committed one.
ref <- lapply(setNames(nm = names(fits)), function(m) {
  p <- polished(fits[[m]])
  list(params = p, derived = derived_at(m, p))
})
for (m in names(ref)) {
  was <- prev[[m]][["derived"]]
  if (is.null(was)) {
    message(sprintf("%-7s derived quantities are new", m))
    next
  }
  d <- vapply(pinned_derived, function(k) {
    now <- ref[[m]]$derived[[k]]
    if (!identical(dim(now), dim(was[[k]]))) return(Inf)
    nz <- was[[k]] != 0
    max(0, abs(now[nz] / was[[k]][nz] - 1), abs(now[!nz]))
  }, numeric(1))
  message(sprintf("%-7s derived max|relative change| %.3e in %s", m, max(d),
                  names(which.max(d))))
}

if (file.exists(out) && !identical(Sys.getenv("RCEATTLE_GOLDEN_REPIN"), "true"))
  stop(out, " exists. Re-pinning is not a routine step: set ",
       "RCEATTLE_GOLDEN_REPIN=true only when the model change above is intended, ",
       "and update the literals in tests/testthat/test-golden-regression.R to the ",
       "objectives printed here, in the same commit.", call. = FALSE)

# Parameters and derived quantities. The objectives stay as literals in the test
# file, where a reviewer sees a re-pin as changed text rather than as a changed
# binary.
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
