# RTMB port gate: the RTMB objective (R/6-rtmb_*.R) against the C++ template on BS2017SS,
# single-species, which is the scope of stage 1.
#
#   Rscript tools/verify/verify-rtmb-oracle.R [<no-FMA DLL path>] [installed=<lib>]
#
# The reference is the C++ template. Build it with -ffp-contract=off (pass its .so path)
# so fused multiply-adds in the report path do not count as disagreement; without one the
# package's own DLL is used. With installed=<lib> the package is loaded from that library
# (byte-compiled, as users run it); otherwise from source, with every .rtmb_* function
# byte-compiled as an install would be.
#
# Precision (every point): objective |df| <= 1e-12 |f|; gradient max|dg| <= 1e-12
# max(|f|, max|g|); every REPORT and every ADREPORT, each on its own scale, elementwise
# |d| <= 1e-12 |c| + 1e-15 max|C| with identical NA/Inf positions; on both sides
# fn(x) = sum(jnll_comp). A non-finite value on one side only fails.
# Derivatives (golden, perturbed, all-free map, jittered points): the Hessian under the
# gradient rule, the ADREPORT Jacobian per quantity on that quantity's own scale;
# standard errors reported, not gated.
# Speed: RTMB retape <= 2x and gr <= 1.2x the C++ build, alternating points.
# Exits non-zero if any gate, or any negative control, fails.

args <- commandArgs(trailingOnly = TRUE)
nofma <- args[!grepl("=", args)][1]
lib <- sub("^installed=", "", grep("^installed=", args, value = TRUE))
if (length(lib)) {
  library(Rceattle, lib.loc = lib)
} else {
  suppressMessages(pkgload::load_all(".", compile = FALSE, quiet = TRUE))
  ns <- asNamespace("Rceattle")
  for (nm in grep("^\\.rtmb_", ls(ns, all.names = TRUE), value = TRUE)) {
    fn <- get(nm, ns)
    if (is.function(fn)) assignInNamespace(nm, compiler::cmpfun(fn), "Rceattle")
  }
}
ref_dll <- "ceattle"
if (!is.na(nofma)) {
  dyn.load(nofma)
  ref_dll <- sub("\\.so$|\\.dll$", "", basename(nofma))
}
cat(sprintf("RTMB ORACLE START %s | reference DLL %s | package from %s\n",
            format(Sys.time()), ref_dll, if (length(lib)) lib else "source"))
make_rtmb <- getFromNamespace(".rtmb_make_adfun", "Rceattle")
data(BS2017SS, package = "Rceattle")
failures <- character(0)
fail <- function(msg) { failures <<- c(failures, msg); cat("  FAIL:", msg, "\n") }

# ---- objects ----------------------------------------------------------------------
quiet_fit <- function(dat, inits = NULL) suppressMessages(suppressWarnings(
  fit_mod(dat, inits = inits, file = NULL, msmMode = 0, estimateMode = 3, niter = 3,
          suitMode = 0, random_rec = FALSE, fit_control = fit_control(getsd = FALSE, verbose = 0))))
pieces <- function(fit) {
  list(data = fit$obj$env$data,
       params = fit$estimated_params[names(fit$obj$env$parameters)],
       map = fit$map$mapFactor)
}
build <- function(pc, ADreport = FALSE, dll = ref_dll) {
  list(cpp = TMB::MakeADFun(pc$data, pc$params, map = pc$map, DLL = dll, silent = TRUE,
                            ADreport = ADreport),
       rtmb = make_rtmb(pc$data, pc$params, map = pc$map, ADreport = ADreport))
}

# ---- comparisons ------------------------------------------------------------------
dims <- function(a) if (is.null(dim(a))) length(a) else dim(a)
# Worst elementwise excess over the report rule (<= 0 passes); Inf on a structural mismatch
report_excess <- function(cc, rr) {
  if (!identical(as.integer(dims(cc)), as.integer(dims(rr)))) return(Inf)
  if (!is.numeric(cc) || !is.numeric(rr)) return(if (identical(cc, rr)) -Inf else Inf)
  cc <- as.vector(cc); rr <- as.vector(rr)
  if (!identical(is.na(cc), is.na(rr))) return(Inf)
  ok <- !is.na(cc)
  if (!identical(is.infinite(cc[ok]), is.infinite(rr[ok])) ||
      any(cc[ok][is.infinite(cc[ok])] != rr[ok][is.infinite(cc[ok])])) return(Inf)
  fin <- ok & is.finite(cc)
  if (!any(fin)) return(-Inf)
  tol <- 1e-12 * abs(cc[fin]) + 1e-15 * max(abs(cc[fin]))
  max(abs(cc[fin] - rr[fin]) - tol)
}
compare_point <- function(o, x, label, oA = NULL) {
  fc <- o$cpp$fn(x); fr <- o$rtmb$fn(x)
  out <- list(df = NA, dg = NA, rep = NA, ad = NA, self = NA, skip = FALSE)
  if (!is.finite(fc)) {
    # A non-finite C++ objective must be the same kind in RTMB; nothing else is compared
    if (!((is.nan(fc) && is.nan(fr)) || identical(fc, fr)))
      fail(sprintf("%s: C++ f = %s, RTMB f = %s", label, fc, fr))
    out$skip <- TRUE
    return(out)
  }
  gc <- o$cpp$gr(x); gr <- o$rtmb$gr(x)
  out$df <- abs(fc - fr) / abs(fc)   # NaN when RTMB alone is non-finite: gate_point fails it
  out$dg <- if (identical(is.finite(gc), is.finite(gr))) {
    max(abs(gc - gr)[is.finite(gc)]) / max(abs(fc), max(abs(gc[is.finite(gc)])))
  } else Inf
  rc <- o$cpp$report(x); rr <- o$rtmb$report(x)
  if (!setequal(names(rc), names(rr))) {
    fail(sprintf("%s: REPORT names differ (%s)", label,
                 paste(c(setdiff(names(rc), names(rr)), setdiff(names(rr), names(rc))), collapse = ", ")))
  }
  ex <- vapply(intersect(names(rc), names(rr)), function(n) report_excess(rc[[n]], rr[[n]]), 0)
  out$rep <- max(ex)
  out$worst <- names(which.max(ex))
  out$self <- max(abs(fc - sum(rc$jnll_comp)) / abs(fc), abs(fr - sum(rr$jnll_comp)) / abs(fr))
  if (!is.null(oA)) out$ad <- adreport_excess(oA, x)
  out$rc <- rc; out$rr <- rr
  out
}
# A NaN or NA comparison fails (!isTRUE(x <= tol)), except at a skipped non-finite point
gate_point <- function(res, label, with_ad = TRUE) {
  if (isTRUE(res$skip)) return(invisible())
  if (!isTRUE(res$df <= 1e-12)) fail(sprintf("%s: objective rel %.3g", label, res$df))
  if (!isTRUE(res$dg <= 1e-12)) fail(sprintf("%s: gradient %.3g", label, res$dg))
  if (!isTRUE(res$rep <= 0)) fail(sprintf("%s: REPORT %s exceeds by %.3g", label, res$worst, res$rep))
  if (with_ad && !isTRUE(res$ad <= 0)) fail(sprintf("%s: ADREPORT %s exceeds by %.3g", label, attr(res$ad, "worst"), res$ad))
  if (!isTRUE(res$self <= 1e-12)) fail(sprintf("%s: fn != sum(jnll_comp), %.3g", label, res$self))
}
# ADREPORTs split by quantity, each compared on its own scale, in the same order both sides
adreport_excess <- function(oA, x) {
  ic <- oA$cpp$env$ADreportIndex(); ir <- oA$rtmb$env$ADreportIndex()
  if (!identical(names(ic), names(ir)) || !identical(lengths(ic), lengths(ir))) return(Inf)
  vc <- oA$cpp$fn(x); vr <- oA$rtmb$fn(x)
  ex <- vapply(names(ic), function(n) report_excess(vc[ic[[n]]], vr[ir[[n]]]), 0)
  structure(max(ex), worst = names(which.max(ex)))
}
# Record every point an nlminb run evaluates (fn and gr calls alike)
path_points <- function(obj, start, max_iter = 1e4) {
  pts <- list()
  rec <- function(x) pts[[length(pts) + 1]] <<- x
  suppressWarnings(nlminb(start, function(x) { rec(x); obj$fn(x) },
                          function(x) { rec(x); obj$gr(x) },
                          control = list(eval.max = 1e4, iter.max = max_iter)))
  unique(pts)
}
summarise <- function(results, label) {
  kept <- Filter(function(r) !isTRUE(r$skip), results)
  get1 <- function(k) max(vapply(kept, function(r) as.numeric(r[[k]]), 0))
  cat(sprintf("  %-28s n=%4d (%d non-finite C++, skipped)  max rel df %.2e  dg %.2e  report excess %.2e  ADREPORT %.2e  self %.2e\n",
              label, length(results), length(results) - length(kept), get1("df"), get1("dg"),
              get1("rep"), get1("ad"), get1("self")))
}
# Coverage, accumulated point by point: REPORTs that never leave 0 or never change,
# jnll_comp cells that are ever non-zero, floored N cells and active decline penalties
cover <- new.env()
cover$n <- 0
update_cover <- function(rr) {
  z <- vapply(rr, function(v) all(v == 0, na.rm = TRUE), TRUE)
  np <- which(cover$np_fleets)
  dec <- sum(vapply(np, function(f) { s <- rr$sel_at_age[f, 1, , 1]; s <- s[s > 0]
    sum(diff(log(s)) < 0) }, 0))
  if (cover$n == 0) {
    cover$first <- rr; cover$zero <- z; cover$const <- setNames(rep(TRUE, length(rr)), names(rr))
    cover$jnz <- rr$jnll_comp != 0; cover$floored <- integer(0); cover$decline <- integer(0)
  } else {
    cover$zero <- cover$zero & z[names(cover$zero)]
    for (n in names(cover$const)[cover$const]) cover$const[n] <- identical(rr[[n]], cover$first[[n]])
    cover$jnz <- cover$jnz | rr$jnll_comp != 0
  }
  N <- rr$N_at_age[, , , 2:H]
  cover$floored <- c(cover$floored, sum(N > 0 & N < 0.001))   # padded ages are exact 0
  cover$decline <- c(cover$decline, dec)
  cover$n <- cover$n + 1
}
run_points <- function(o, oA, pts, label) {
  res <- lapply(seq_along(pts), function(i) {
    r <- compare_point(o, pts[[i]], sprintf("%s #%d", label, i), oA)
    gate_point(r, sprintf("%s #%d", label, i))
    if (!is.null(r$rr)) update_cover(r$rr)
    r[c("df", "dg", "rep", "ad", "self", "skip")]
  })
  summarise(res, label)
  cat(sprintf("     floored N cells per point %d-%d, active decline pairs %d-%d\n",
              min(tail(cover$floored, length(pts))), max(tail(cover$floored, length(pts))),
              min(tail(cover$decline, length(pts))), max(tail(cover$decline, length(pts)))))
  invisible(res)
}

# ---- fixtures ---------------------------------------------------------------------
cat("\n[building fixtures]\n")
fit0 <- quiet_fit(BS2017SS)
golden <- readRDS(file.path("tests", "testthat", "fixtures", "golden-reference.rds"))
fitg <- quiet_fit(BS2017SS, inits = golden$ss$params)
pc0 <- pieces(fit0); pcg <- pieces(fitg)
o <- build(pc0); oA <- build(pc0, ADreport = TRUE)
stopifnot(identical(names(o$cpp$par), names(o$rtmb$par)))
x0 <- o$cpp$par
cover$np_fleets <- pc0$data$flt_sel_type == 2 & pc0$data$flt_sel_lead == 1
H <- pc0$data$endyr - pc0$data$styr + 1; Y <- pc0$data$projyr - pc0$data$styr + 1
xg <- TMB::MakeADFun(pcg$data, pcg$params, map = pcg$map, DLL = ref_dll, silent = TRUE)$par
all_res <- list()

# ---- precision --------------------------------------------------------------------
cat("\n[precision]\n")
p1 <- path_points(o$cpp, x0)
all_res$path <- run_points(o, oA, p1, "1 nlminb path, default start")
set.seed(20261117)
jit <- lapply(1:2, function(k) path_points(o$cpp, x0 + rnorm(length(x0), 0, 0.1), max_iter = 100))
all_res$jitter1 <- run_points(o, oA, jit[[1]], "2 jittered path a (100 iter)")
all_res$jitter2 <- run_points(o, oA, jit[[2]], "2 jittered path b (100 iter)")
logF <- names(xg) == "log_F"
# posfun returns a floored cell in [eps/2, eps), so a straddle point is one with cells
# just under eps AND just over it: maximise the smaller of the two counts
n_near <- function(x) {
  N <- o$cpp$report(x)$N_at_age[, , , 2:H]
  c(below = sum(N > 0 & N < 0.001), under = sum(N >= 0.001 / 1.2 & N < 0.001),
    over = sum(N >= 0.001 & N <= 0.00125))
}
offs <- seq(0, 6, by = 0.1)
nn <- vapply(offs, function(d) { x <- xg; x[logF] <- x[logF] + d; n_near(x) }, numeric(3))
best <- which.max(pmin(nn["under", ], nn["over", ]))
if (best %in% c(1, length(offs))) fail("straddle scan maximised at the end of its grid")
x_deep <- xg; x_deep[logF] <- x_deep[logF] + 3
x_str <- xg; x_str[logF] <- x_str[logF] + offs[best]
all_res$floor <- run_points(o, oA, list(x_deep, x_str),
                            sprintf("3 deep floor, straddle (+%.1f)", offs[best]))
cat(sprintf("     floored N cells: deep %d; straddle %d, with %d in [eps/1.2, eps) and %d in [eps, 1.25 eps]\n",
            n_near(x_deep)[["below"]], nn["below", best], nn["under", best], nn["over", best]))
all_res$golden <- run_points(o, oA, list(xg), "4 golden point")

# Point 5: fixed parameters moved off BS2017SS's identity, zero or constant values
pert <- pcg$params
set.seed(5)
pert$comp_weights[] <- 0.7
pert$proj_F_prop[pc0$data$flt_type == 1] <- 0.5
pert$log_Ftarget[] <- log(0.2); pert$log_Flimit[] <- log(0.3)
pert$log_M1_dev[] <- rnorm(length(pert$log_M1_dev), 0, 0.05)
pert$index_q_dev[] <- rnorm(length(pert$index_q_dev), 0, 0.05)
pert$rec_dev[, (H + 1):Y] <- rnorm(3 * (Y - H), 0, 0.3)
np <- which(pc0$data$flt_sel_type == 2)
pert$sel_coff[np, 1, 6:8] <- pert$sel_coff[np, 1, 6:8] - 0.4 * rep(1:3, each = length(np))
pert$M1_beta[] <- 0.05                                   # environmental effect on M
pert$sel_coff_dev[] <- rnorm(length(pert$sel_coff_dev), 0, 0.05)
pert$log_sel_slp_dev[] <- rnorm(length(pert$log_sel_slp_dev), 0, 0.05)
pert$sel_inf_dev[] <- rnorm(length(pert$sel_inf_dev), 0, 0.05)
pert$log_pop_scalar[] <- 0.1; pert$log_gam_a[] <- 0.2; pert$log_gam_b[] <- -0.3
pert$log_growth_pars[] <- pert$log_growth_pars[] + 0.1
pc5 <- modifyList(pcg, list(params = pert))
o5 <- build(pc5); o5A <- build(pc5, ADreport = TRUE)
all_res$perturbed <- run_points(o5, o5A, list(o5$cpp$par), "5 perturbed fixed parameters")

# Point 6: an all-free map, so derivatives through normally fixed parameters are tested
free <- function(nm, cells = !is.na(pcg$params[[nm]])) {
  f <- rep(NA_integer_, length(pcg$params[[nm]])); f[cells] <- seq_len(sum(cells))
  factor(f)
}
map6 <- pcg$map
for (nm in c("R_log_sd", "comp_weights", "log_M1")) map6[[nm]] <- free(nm)
np_lead <- which(pc0$data$flt_sel_type == 2 & pc0$data$flt_sel_lead == 1)
pen_cells <- array(FALSE, dim(pcg$params$sel_curve_pen)); pen_cells[np_lead, 1:2] <- TRUE
map6$sel_curve_pen <- free("sel_curve_pen", pen_cells)
map6$index_log_q <- free("index_log_q", seq_along(pcg$params$index_log_q) %in% 4:7)
pc6 <- modifyList(pcg, list(map = map6))
o6 <- build(pc6); o6A <- build(pc6, ADreport = TRUE)
# Off the optimum too: R0 is free and confounded with the mean deviation, so at the
# optimum the deviations' centring has zero derivative and a broken one goes unseen
set.seed(6)
x6j <- o6$cpp$par + rnorm(length(o6$cpp$par), 0, 0.05)
# The same jittered point with sigma_R at its starting value, so the deviations' centring
# is unchanged in value and differs only in its derivative
x6r <- x6j; x6r[names(x6r) == "R_log_sd"] <- o6$cpp$par[names(o6$cpp$par) == "R_log_sd"]
all_res$allfree <- run_points(o6, o6A, list(o6$cpp$par, x6j, x6r), "6 all-free map (+ 2 jittered)")

# Point 7: data perturbed off BS2017SS's identity ageing error and January spawning
dat7 <- BS2017SS
for (i in unique(dat7$age_error[, "Species"])) {
  # asymmetric ageing error: 15% read one year young, 5% one year old
  rows <- which(dat7$age_error[, "Species"] == i)
  n <- length(rows); m <- diag(0.8, n)
  m[cbind(2:n, 1:(n - 1))] <- 0.15; m[cbind(1:(n - 1), 2:n)] <- 0.05
  dat7$age_error[rows, 2 + seq_len(n)] <- m / rowSums(m)
}
dat7$spawn_month[] <- 3
pc7 <- pieces(quiet_fit(dat7, inits = golden$ss$params))
o7 <- build(pc7); o7A <- build(pc7, ADreport = TRUE)
stopifnot(identical(names(o7$cpp$par), names(x0)))
all_res$data7 <- run_points(o7, o7A, list(o7$cpp$par, x0), "7 perturbed data")

# Point 8: a survey with compositions but no index and no Catchability form, whose q is
# its starting value; and point 9: a model with no projection years
dat8 <- BS2017SS
dat8$index_data <- dat8$index_data[dat8$index_data$Fleet_code != 4, ]
dat8$fleet_control$Catchability[4] <- NA
pc8 <- pieces(quiet_fit(dat8))
o8 <- build(pc8); o8A <- build(pc8, ADreport = TRUE)
set.seed(8)
all_res$compq <- run_points(o8, o8A, list(o8$cpp$par, o8$cpp$par + rnorm(length(o8$cpp$par), 0, 0.05)),
                            "8 comp-only survey, start q")
dat9 <- BS2017SS
dat9$projyr <- dat9$endyr
pc9 <- pieces(quiet_fit(dat9))
o9 <- build(pc9); o9A <- build(pc9, ADreport = TRUE)
set.seed(9)
all_res$noproj <- run_points(o9, o9A, list(o9$cpp$par, o9$cpp$par + rnorm(length(o9$cpp$par), 0, 0.05)),
                             "9 no projection years")

# ---- scope guard ------------------------------------------------------------------
cat("\n[scope guard: each input the C++ fits differently must be refused]\n")
refusals <- list(
  "estimateMode 4" = function(d) { d$estimateMode <- 4; d },
  "Catchability with prior on a fishery" = function(d) { d$est_index_q[1] <- 2; d },
  "q deviations on a fleet without an index" = function(d) { d$est_index_q[1] <- 1; d$index_varying_q[1] <- 1; d },
  "index row under another species" = function(d) { d$index_ctl[d$index_ctl[, 1] == 4, 2] <- 2; d },
  "catch row under another species" = function(d) { d$catch_ctl[d$catch_ctl[, 1] == 1, 2] <- 2; d },
  "composition row under another species" = function(d) { d$comp_ctl[d$comp_ctl[, 1] == 1, 2] <- 2; d },
  "emp_sel Sex 2" = function(d) { d$emp_sel_ctl[, 3] <- 2; d },
  "Analytical q on a comp-only survey" = function(d) { d$est_index_q[4] <- 3; d })
for (nm in names(refusals)) {
  e <- try(make_rtmb(refusals[[nm]](pc0$data), pc0$params, map = pc0$map), silent = TRUE)
  ok <- inherits(e, "try-error") && grepl("does not yet support", e)
  cat(sprintf("  %-45s %s\n", nm, if (ok) "refused" else "ACCEPTED"))
  if (!ok) fail(paste("scope guard accepted:", nm))
}
e <- try(make_rtmb(pc0$data, pc0$params, map = pc0$map, random = "rec_dev"), silent = TRUE)
cat(sprintf("  %-45s %s\n", "random effects", if (inherits(e, "try-error")) "refused" else "ACCEPTED"))
if (!inherits(e, "try-error")) fail("scope guard accepted random effects")

# ---- derivatives ------------------------------------------------------------------
cat("\n[derivatives: ADREPORT Jacobian, Hessian, standard errors]\n")
deriv <- function(ob, obA, x, label) {
  f <- abs(ob$cpp$fn(x))
  jc <- obA$cpp$gr(x); jr <- obA$rtmb$gr(x)
  hc <- ob$cpp$he(x); hr <- ob$rtmb$he(x)
  # Jacobian rows grouped by ADREPORT quantity, each on its own value and derivative scale
  idx <- obA$cpp$env$ADreportIndex(); val <- obA$cpp$fn(x)
  dj <- max(vapply(idx[lengths(idx) > 0], function(i)
    max(abs(jc[i, , drop = FALSE] - jr[i, , drop = FALSE])) /
      max(abs(val[i]), abs(jc[i, , drop = FALSE]), 1e-300), 0))
  dh <- max(abs(hc - hr)) / max(f, max(abs(hc)))
  if (dj > 1e-12) fail(sprintf("%s: Jacobian %.3g", label, dj))
  if (dh > 1e-12) fail(sprintf("%s: Hessian %.3g", label, dh))
  se <- function(obj, sdr) {
    s <- summary(suppressWarnings(sdr(obj, par.fixed = x)), "report")
    s[rownames(s) %in% c("ssb", "log_ssb", "biomass", "R", "ssb_depletion", "biomass_depletion"), 2]
  }
  sc <- se(ob$cpp, TMB::sdreport); sr <- se(ob$rtmb, RTMB::sdreport)
  kap <- tryCatch(kappa(hc, exact = TRUE), error = function(e) NA)
  same_na <- identical(is.na(sc), is.na(sr))
  ok <- !is.na(sc) & !is.na(sr) & sc != 0
  cat(sprintf("  %-26s Jacobian %.2e  Hessian %.2e  SE max rel %.2e over %d, NaN positions %s  (Hessian condition %.2e)\n",
              label, dj, dh, if (any(ok)) max(abs(sc[ok] - sr[ok]) / abs(sc[ok])) else NA,
              sum(ok), if (same_na) "same" else "DIFFER", kap))
}
deriv(o, oA, xg, "4 golden point")
deriv(o5, o5A, o5$cpp$par, "5 perturbed parameters")
deriv(o6, o6A, o6$cpp$par, "6 all-free map")
deriv(o6, o6A, x6j, "6 all-free map, jittered")
deriv(o, oA, jit[[2]][[length(jit[[2]])]], "2 jittered point b")
deriv(o, oA, jit[[1]][[length(jit[[1]])]], "2 jittered point")

# ---- negative controls ------------------------------------------------------------
cat("\n[negative controls: each deliberate break must fail where named]\n")
with_break <- function(fun, from, to, expr) {
  orig <- getFromNamespace(fun, "Rceattle")
  txt <- deparse(orig, width.cutoff = 500L)
  hit <- grep(from, txt, fixed = TRUE)
  if (length(hit) != 1) stop("negative control: pattern not found once in ", fun, ": ", from)
  txt[hit] <- sub(from, to, txt[hit], fixed = TRUE)
  brk <- eval(parse(text = txt), envir = environment(orig))
  environment(brk) <- environment(orig)
  assignInNamespace(fun, compiler::cmpfun(brk), "Rceattle")
  on.exit(assignInNamespace(fun, orig, "Rceattle"))
  expr
}
detects <- function(ob, x) {
  r <- compare_point(ob, x, "control")
  !isTRUE(r$df <= 1e-12) || !isTRUE(r$dg <= 1e-12) || !isTRUE(r$rep <= 0)
}
sd0 <- exp(pcg$params$R_log_sd[1])
stopifnot(all(exp(pcg$params$R_log_sd) == sd0))
controls <- list(
  list("drop the ZERO_N row (start point)", ".rtmb_objective",
       "jnll_comp[jrow[[\"Zero n-at-age penalty\"]], sp] <- zero_N_pen[sp]",
       "jnll_comp[jrow[[\"Zero n-at-age penalty\"]], sp] <- 0 * zero_N_pen[sp]",
       function() detects(build(pc0), x0)),
  list("transpose age_error (point 7)", ".rtmb_prepare_data",
       "ae_t = t(matrix(ae, na, na))", "ae_t = matrix(ae, na, na)",
       function() detects(build(pc7), o7$cpp$par)),
  list("posfun as a plain comparison (MakeADFun must error)", ".rtmb_posfun",
       "m <- .rtmb_ind_ge(x, eps)", "m <- x >= eps",
       function() {
         e <- try(make_rtmb(pc0$data, pc0$params, map = pc0$map), silent = TRUE)
         inherits(e, "try-error") && grepl("Comparison is generally unsafe", e)
       }),
  list("bias centring from a fixed sigma_R (point 6, derivative only)", ".rtmb_objective",
       "mu <- -pd$bias_adjust_proc * (R_sd[sp] * R_sd[sp])/2",
       sprintf("mu <- -pd$bias_adjust_proc * (%.17g)/2", sd0 * sd0),
       function() detects(build(pc6), x6r)),
  list("off-by-one year in M (point 5)", ".rtmb_objective",
       "p[[\"log_M1_dev\"]][sp, 1, a, ]", "p[[\"log_M1_dev\"]][sp, 1, a, c(2:Y, 1)]",
       function() detects(build(pc5), o5$cpp$par)))
for (k in controls) {
  # Null control first: the same function round-tripped unbroken must pass
  null_ok <- !with_break(k[[2]], k[[3]], k[[3]], k[[5]]())
  caught <- with_break(k[[2]], k[[3]], k[[4]], k[[5]]())
  cat(sprintf("  %-62s %s (unbroken: %s)\n", k[[1]], if (caught) "caught" else "MISSED",
              if (null_ok) "passes" else "FAILS"))
  if (!caught) fail(paste("negative control missed:", k[[1]]))
  if (!null_ok) fail(paste("negative control's unbroken round trip fails:", k[[1]]))
}

# ---- coverage ---------------------------------------------------------------------
cat("\n[coverage]\n")
cat(sprintf("  points compared: %d (path %d, jittered %d, floor 2, golden 1, perturbed 1, all-free 3, data 2, comp-only q 2, no projection 2)\n",
            cover$n, length(p1), sum(lengths(jit))))
empty <- names(Filter(function(v) length(v) == 0, cover$first))
cat("  REPORTs empty (length 0):", paste(sort(empty), collapse = ", "), "\n")
cat("  REPORTs zero at every point:", paste(sort(setdiff(names(which(cover$zero)), empty)), collapse = ", "), "\n")
cat("  REPORTs constant (not zero) at every point:",
    paste(sort(setdiff(names(which(cover$const)), names(which(cover$zero)))), collapse = ", "), "\n")
nz <- which(cover$jnz, arr.ind = TRUE)
rows <- getFromNamespace(".JNLL_ROW_AXIS", "Rceattle")
cat("  jnll_comp cells non-zero at some point:",
    paste(sprintf("%s[%d]", names(rows)[nz[, 1]], nz[, 2]), collapse = "; "), "\n")

# ---- speed ------------------------------------------------------------------------
cat("\n[speed, interleaved; reference = the package's own -O2 DLL]\n")
pcs <- pc0
tape <- matrix(NA, 5, 2, dimnames = list(NULL, c("cpp", "rtmb")))
for (i in 1:5) {
  tape[i, 1] <- system.time(TMB::MakeADFun(pcs$data, pcs$params, map = pcs$map, DLL = "ceattle", silent = TRUE))[["elapsed"]]
  tape[i, 2] <- system.time(make_rtmb(pcs$data, pcs$params, map = pcs$map))[["elapsed"]]
}
os <- build(pcs, dll = "ceattle")
xa <- xg; xb <- xg + 1e-3
time_gr <- function(obj, n = 300) {
  t <- system.time(for (i in seq_len(n)) { x <- if (i %% 2) xa else xb; obj$fn(x); obj$gr(x) })[["elapsed"]]
  t / n
}
gr_only <- function(obj, n = 300) system.time(for (i in seq_len(n)) obj$gr(if (i %% 2) xa else xb))[["elapsed"]] / n
pair <- gro <- matrix(NA, 5, 2)
for (i in 1:5) {   # C++ and RTMB blocks alternate
  pair[i, ] <- c(time_gr(os$cpp), time_gr(os$rtmb))
  gro[i, ] <- c(gr_only(os$cpp), gr_only(os$rtmb))
}
rt <- apply(tape, 2, median); g_c <- median(pair[, 1]); g_r <- median(pair[, 2])
gc1 <- median(gro[, 1]); gr1 <- median(gro[, 2])
cat(sprintf("  retape: C++ %.3f s  RTMB %.3f s  ratio %.2f (gate <= 2)\n", rt[1], rt[2], rt[2] / rt[1]))
cat(sprintf("  fn+gr pair: C++ %.2f ms  RTMB %.2f ms  ratio %.2f (reported)\n", 1e3 * g_c, 1e3 * g_r, g_r / g_c))
cat(sprintf("  gr alone: C++ %.2f ms  RTMB %.2f ms  ratio %.2f (gate <= 1.2)\n", 1e3 * gc1, 1e3 * gr1, gr1 / gc1))
if (rt[2] / rt[1] > 2) fail(sprintf("retape ratio %.2f > 2", rt[2] / rt[1]))
if (gr1 / gc1 > 1.2) fail(sprintf("gr ratio %.2f > 1.2", gr1 / gc1))

# ---- fit --------------------------------------------------------------------------
cat("\n[fit: both optimized from the default start]\n")
fit_one <- function(obj) suppressWarnings(nlminb(x0, obj$fn, obj$gr,
                                                 control = list(eval.max = 1e4, iter.max = 1e4)))
oc <- fit_one(os$cpp); orr <- fit_one(os$rtmb)
rc <- os$cpp$report(oc$par); rr <- os$rtmb$report(orr$par)
rel <- function(a, b) max(abs(a - b) / abs(a))
cat(sprintf("  objective C++ %.10f  RTMB %.10f  (iterations %d, %d)\n", oc$objective, orr$objective,
            oc$iterations, orr$iterations))
cat(sprintf("  max rel difference at the optima: ssb %.2e  biomass %.2e  R %.2e; max|g| %.1e, %.1e\n",
            rel(rc$ssb, rr$ssb), rel(rc$biomass, rr$biomass), rel(rc$R, rr$R),
            max(abs(os$cpp$gr(oc$par))), max(abs(os$rtmb$gr(orr$par)))))
# Golden's convention for a re-optimized objective: 1e-6 relative, and both converged
if (oc$convergence != 0 || orr$convergence != 0) fail("an optimization did not converge")
if (abs(oc$objective - orr$objective) > 1e-6 * abs(oc$objective)) fail("optimized objectives differ by more than 1e-6 relative")

cat(sprintf("\nRTMB ORACLE %s | %d failure(s) | END %s\n",
            if (length(failures)) "FAIL" else "PASS", length(failures), format(Sys.time())))
if (length(failures)) { cat(paste(" -", failures), sep = "\n"); quit(status = 1) }
