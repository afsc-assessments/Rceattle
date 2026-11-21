# RTMB port gate: the RTMB objective (R/6-rtmb_*.R) against the C++ template.
#   fixture=bs2017ss (default): BS2017SS, single-species, the scope of stage 1.
#   fixture=bs2017ms: BS2017MS, multispecies predation (msmMode 1 and 2, suitMode 0 and
#     4), the scope of stage 2; config=m1s0,m2s0,m1s4,m2s4 runs a subset.
#
#   Rscript tools/verify/verify-rtmb-oracle.R [<no-FMA DLL path>] [installed=<lib>]
#           [fixture=bs2017ss|bs2017ms] [config=<list>] [speed=yes|no|only]
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
# BS2017MS only: a C++ `if` on an estimated value is decided once, when the tape is
# recorded, and RTMB's indicators decide at every point; where the two differ (a "flip"),
# the C++ is re-taped at that point and RTMB is gated against the re-taped C++.
# Exits non-zero if any gate, or any negative control, fails.

args <- commandArgs(trailingOnly = TRUE)
nofma <- args[!grepl("=", args)][1]
lib <- sub("^installed=", "", grep("^installed=", args, value = TRUE))
fixture <- sub("^fixture=", "", grep("^fixture=", args, value = TRUE))
if (!length(fixture)) fixture <- "bs2017ss"
stopifnot(fixture %in% c("bs2017ss", "bs2017ms"))
configs <- c("m1s0", "m2s0", "m1s4", "m2s4")
cfg_arg <- sub("^config=", "", grep("^config=", args, value = TRUE))
if (length(cfg_arg)) configs <- intersect(configs, strsplit(cfg_arg, ",")[[1]])
# speed=no skips the speed gate (for parallel runs); speed=only runs nothing else (BS2017MS)
speed <- sub("^speed=", "", grep("^speed=", args, value = TRUE))
if (!length(speed)) speed <- "yes"
stopifnot(speed %in% c("yes", "no", "only"))
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
ver <- function(pkg) as.character(utils::packageVersion(pkg))
cat(sprintf("RTMB ORACLE START %s | fixture %s | reference DLL %s | package from %s | RTMB %s TMB %s R %s\n",
            format(Sys.time()), fixture, ref_dll, if (length(lib)) lib else "source",
            ver("RTMB"), ver("TMB"), paste(R.version$major, R.version$minor, sep = ".")))
make_rtmb <- getFromNamespace(".rtmb_make_adfun", "Rceattle")
data(BS2017SS, package = "Rceattle")
data(BS2017MS, package = "Rceattle")
failures <- character(0)
fail <- function(msg) { failures <<- c(failures, msg); cat("  FAIL:", msg, "\n") }

# ---- objects ----------------------------------------------------------------------
quiet_fit <- function(dat, inits = NULL, msmMode = 0, suitMode = 0, niter = 3, ...) {
  suppressMessages(suppressWarnings(
    fit_mod(dat, inits = inits, file = NULL, msmMode = msmMode, estimateMode = 3, niter = niter,
            suitMode = suitMode, random_rec = FALSE,
            fit_control = fit_control(getsd = FALSE, verbose = 0), ...)))
}
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
  # The C++ reports log(x) through Eigen's vectorized log, which clamps a subnormal x to
  # DBL_MIN; RTMB reports log(x). Those cells must hold exactly those two values.
  out$subnormal <- 0
  for (nm in c("ssb", "biomass", "R", "index_hat")) {
    lg <- paste0("log_", nm); a <- rc[[nm]]
    sub <- !is.na(a) & a > 0 & a < .Machine$double.xmin
    if (!any(sub)) next
    if (!isTRUE(all(rc[[lg]][sub] == log(.Machine$double.xmin))) || !isTRUE(all(rr[[lg]][sub] == log(rr[[nm]][sub]))))
      fail(sprintf("%s: %s at a subnormal %s is neither the C++ clamp nor log(x)", label, lg, nm))
    rc[[lg]][sub] <- rr[[lg]][sub]
    out$subnormal <- out$subnormal + sum(sub)
  }
  if (!setequal(names(rc), names(rr))) {
    fail(sprintf("%s: REPORT names differ (%s)", label,
                 paste(c(setdiff(names(rc), names(rr)), setdiff(names(rr), names(rc))), collapse = ", ")))
  }
  ex <- vapply(intersect(names(rc), names(rr)), function(n) report_excess(rc[[n]], rr[[n]]), 0)
  out$rep <- max(ex)
  out$worst <- names(which.max(ex))
  out$self_cpp <- abs(fc - sum(rc$jnll_comp)) / abs(fc)
  out$self_rtmb <- abs(fr - sum(rr$jnll_comp)) / abs(fr)
  out$self <- max(out$self_cpp, out$self_rtmb)
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

# ---- shared checks ----------------------------------------------------------------
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
finish <- function() {
  cat(sprintf("\nRTMB ORACLE %s | %d failure(s) | END %s\n",
              if (length(failures)) "FAIL" else "PASS", length(failures), format(Sys.time())))
  if (length(failures)) { cat(paste(" -", failures), sep = "\n"); quit(status = 1) }
  quit(status = 0)
}

# ---- BS2017MS (stage 2: predation) ------------------------------------------------
run_bs2017ms <- function(configs) {
  golden <- readRDS(file.path("tests", "testthat", "fixtures", "golden-reference.rds"))
  specs <- list(m1s0 = list(msm = 1, suit = 0, niter = 5, inits = golden$ms$params),
                m2s0 = list(msm = 2, suit = 0, niter = 5, inits = golden$ms$params),
                m1s4 = list(msm = 1, suit = 4, niter = 3, inits = golden$ss$params),
                m2s4 = list(msm = 2, suit = 4, niter = 3, inits = golden$ss$params))
  msfit <- function(sp, inits = sp$inits, dat = BS2017MS, niter = sp$niter, ...) {
    pieces(quiet_fit(dat, inits = inits, msmMode = sp$msm, suitMode = rep(sp$suit, 3),
                     niter = niter, ...))
  }
  H <<- BS2017MS$endyr - BS2017MS$styr + 1
  Y <- BS2017MS$projyr - BS2017MS$styr + 1

  # Every RTMB indicator the C++ decides once at taping, as it comes out at x
  ind_states <- function(rtmb, x) {
    rec <- list()
    orig <- list(le = getFromNamespace(".rtmb_ind_le", "Rceattle"),
                 gt = getFromNamespace(".rtmb_ind_gt", "Rceattle"))
    wrap <- function(f) function(x, a) { v <- f(x, a); rec[[length(rec) + 1L]] <<- v; v }
    assignInNamespace(".rtmb_ind_le", wrap(orig$le), "Rceattle")
    assignInNamespace(".rtmb_ind_gt", wrap(orig$gt), "Rceattle")
    on.exit({
      assignInNamespace(".rtmb_ind_le", orig$le, "Rceattle")
      assignInNamespace(".rtmb_ind_gt", orig$gt, "Rceattle")
    })
    rtmb$report(x)
    rec
  }
  mc <- new.env()
  ms_cover <- function(rr, sid) {
    fl <- 1e-6 / (1 + 1e-6)
    mc$floor <- c(mc$floor, sum(rr$suit_other > 0 & rr$suit_other <= fl * (1 + 1e-9), na.rm = TRUE))
    mc$maxM2 <- c(mc$maxM2, max(rr$M2_at_age, na.rm = TRUE))
    mc$dsbf <- c(mc$dsbf, any(rr$DynamicSBF != rr$DynamicSB0, na.rm = TRUE))
    mc$nan <- c(mc$nan, sum(is.na(rr$M2_at_age)))
    mc$stom_floor <- c(mc$stom_floor, sum(1 - tapply(rr$diet_hat[, 2], sid, sum) < 1e-5))
  }
  # Gate points; where a branch the C++ froze has flipped, re-tape the C++ there
  run_ms <- function(o, oA, pc, pts, label) {
    tape_state <- ind_states(o$rtmb, o$cpp$par)
    sid <- pc$data$stomach_id
    nflip <- 0
    nsub <- 0
    res <- lapply(seq_along(pts), function(i) {
      x <- pts[[i]]
      lab <- sprintf("%s #%d", label, i)
      oo <- o; ooA <- oA
      if (!identical(ind_states(o$rtmb, x), tape_state)) {
        nflip <<- nflip + 1
        pl <- o$cpp$env$parList(x)
        oo$cpp <- TMB::MakeADFun(pc$data, pl, map = pc$map, DLL = ref_dll, silent = TRUE)
        ooA$cpp <- TMB::MakeADFun(pc$data, pl, map = pc$map, DLL = ref_dll, silent = TRUE,
                                  ADreport = TRUE)
      }
      r <- compare_point(oo, x, lab, ooA)
      gate_point(r, lab)
      if (!is.null(r$rr)) { update_cover(r$rr); ms_cover(r$rr, sid) }
      nsub <<- nsub + if (is.null(r$subnormal)) 0 else r$subnormal
      r[c("df", "dg", "rep", "ad", "self", "skip")]
    })
    summarise(res, label)
    cat(sprintf("     branch flips re-taped: %d of %d points%s\n", nflip, length(pts),
                if (nsub > 0) sprintf("; %d subnormal log cells (C++ clamp)", nsub) else ""))
    invisible(nflip)
  }
  free <- function(params, nm, cells = !is.na(params[[nm]])) {
    f <- rep(NA_integer_, length(params[[nm]])); f[cells] <- seq_len(sum(cells))
    factor(f)
  }

  for (cf in configs) {
    sp <- specs[[cf]]
    cat(sprintf("\n==== %s: msmMode %d, suitMode %d, niter %d ====\n", cf, sp$msm, sp$suit, sp$niter))
    cover$n <- 0
    for (k in c("floor", "maxM2", "dsbf", "nan", "stom_floor")) mc[[k]] <- NULL
    pcg <- msfit(sp)
    o <- build(pcg); oA <- build(pcg, ADreport = TRUE)
    stopifnot(identical(names(o$cpp$par), names(o$rtmb$par)))
    cover$np_fleets <- pcg$data$flt_sel_type == 2 & pcg$data$flt_sel_lead == 1
    xg <- o$cpp$par
    pc0 <- msfit(sp, inits = NULL)
    o0 <- build(pc0); o0A <- build(pc0, ADreport = TRUE)
    x0 <- o0$cpp$par
    npts <- 0
    if (speed != "only") {

    cat("\n[precision]\n")
    if (cf == "m1s0") {
      # Reference positive controls: the golden objective, and decision 7's NaN (#324)
      gold <- o$cpp$fn(xg)
      cat(sprintf("  reference: C++ golden objective %.10f (golden-reference 10267.2478324443)\n", gold))
      if (abs(gold - 10267.2478324443) > 1e-10 * gold) fail("reference C++ does not reproduce the golden objective")
      p1 <- path_points(o0$cpp, x0)
      run_ms(o0, o0A, pc0, p1, "1 nlminb path, cold start")
    } else {
      p1 <- path_points(o$cpp, xg)
      run_ms(o, oA, pcg, p1, "1 nlminb path, from inits")
      run_ms(o0, o0A, pc0, list(x0), "1 cold start (its own tape)")
    }
    set.seed(20261117)
    jit <- lapply(1:2, function(k) path_points(o$cpp, xg + rnorm(length(xg), 0, 0.1), max_iter = 100))
    run_ms(o, oA, pcg, jit[[1]], "2 jittered path a (100 iter)")
    run_ms(o, oA, pcg, jit[[2]], "2 jittered path b (100 iter)")
    npts <- npts + length(p1) + sum(lengths(jit))
    logF <- names(xg) == "log_F"
    x_deep <- xg; x_deep[logF] <- x_deep[logF] + 3
    if (cf == "m1s0") {
      n_near <- function(x) {
        N <- o$cpp$report(x)$N_at_age[, , , 2:H]
        c(under = sum(N >= 0.001 / 1.2 & N < 0.001), over = sum(N >= 0.001 & N <= 0.00125))
      }
      offs <- seq(0, 6, by = 0.1)
      nn <- vapply(offs, function(d) { x <- xg; x[logF] <- x[logF] + d; n_near(x) }, numeric(2))
      best <- which.max(pmin(nn["under", ], nn["over", ]))
      if (best %in% c(1, length(offs))) fail("straddle scan maximised at the end of its grid")
      x_str <- xg; x_str[logF] <- x_str[logF] + offs[best]
      run_ms(o, oA, pcg, list(x_deep, x_str), sprintf("3 deep floor, straddle (+%.1f)", offs[best]))
      npts <- npts + 2
    }
    run_ms(o, oA, pcg, list(xg), "4 golden point")

    # Point 5: fixed parameters moved off their identity, zero or symmetric values
    pert <- pcg$params
    set.seed(5)
    pert$comp_weights[] <- 0.7
    pert$proj_F_prop[pcg$data$flt_type == 1] <- 0.5
    pert$log_Ftarget[] <- log(0.2); pert$log_Flimit[] <- log(0.3)
    pert$log_M1_dev[] <- rnorm(length(pert$log_M1_dev), 0, 0.05)
    pert$index_q_dev[] <- rnorm(length(pert$index_q_dev), 0, 0.05)
    pert$rec_dev[, (H + 1):Y] <- rnorm(length(pert$rec_dev[, (H + 1):Y]), 0, 0.3)
    np <- which(pcg$data$flt_sel_type == 2)
    pert$sel_coff[np, 1, 6:8] <- pert$sel_coff[np, 1, 6:8] - 0.4 * rep(1:3, each = length(np))
    pert$M1_beta[] <- 0.05
    pert$diet_comp_weights[] <- c(0.7, 1, 1.3)
    pert$log_phi[] <- c(0.5, -0.3, 1, 0.2, 0.8, -1, 0, 0.4, -0.6)
    pert$log_gam_a[] <- log(c(0.8, 1.2, 1.5)); pert$log_gam_b[] <- log(c(0.6, 0.9, 1.2))
    pc5 <- modifyList(pcg, list(params = pert))
    o5 <- build(pc5); o5A <- build(pc5, ADreport = TRUE)
    run_ms(o5, o5A, pc5, list(o5$cpp$par), "5 perturbed fixed parameters")

    # Point 6: an all-free map; under suitMode 0 log_phi and gam are inert, gradient 0
    map6 <- pcg$map
    for (nm in c("R_log_sd", "comp_weights", "log_M1")) map6[[nm]] <- free(pcg$params, nm)
    np_lead <- which(pcg$data$flt_sel_type == 2 & pcg$data$flt_sel_lead == 1)
    pen_cells <- array(FALSE, dim(pcg$params$sel_curve_pen)); pen_cells[np_lead, 1:2] <- TRUE
    map6$sel_curve_pen <- free(pcg$params, "sel_curve_pen", pen_cells)
    q_cells <- !is.na(pcg$params$index_log_q) & pcg$data$flt_type == 2
    map6$index_log_q <- free(pcg$params, "index_log_q", q_cells)
    if (sp$suit > 0) map6$diet_comp_weights <- free(pcg$params, "diet_comp_weights")
    if (sp$suit == 0) for (nm in c("log_phi", "log_gam_a", "log_gam_b")) map6[[nm]] <- free(pcg$params, nm)
    pc6 <- modifyList(pcg, list(map = map6))
    o6 <- build(pc6); o6A <- build(pc6, ADreport = TRUE)
    set.seed(6)
    x6j <- o6$cpp$par + rnorm(length(o6$cpp$par), 0, 0.05)
    run_ms(o6, o6A, pc6, list(o6$cpp$par, x6j), "6 all-free map (+ 1 jittered)")
    if (sp$suit == 0) {
      inert <- grepl("^log_phi|^log_gam", names(o6$cpp$par))
      g6 <- rbind(o6$cpp$gr(x6j)[inert], o6$rtmb$gr(x6j)[inert])
      if (any(g6 != 0)) fail(sprintf("%s: gradient on inert log_phi/gam is not 0", cf))
    }
    npts <- npts + 6

    # Point 7: data off their identity values: other food by predator, a window per
    # predator, ageing error, spawning month, and multispecies SB0/B0
    if (cf %in% c("m1s0", "m1s4")) {
      dat7 <- BS2017MS
      dat7$other_food <- dat7$other_food * c(0.5, 1, 2)
      for (i in unique(dat7$age_error[, "Species"])) {
        rows <- which(dat7$age_error[, "Species"] == i)
        n <- length(rows); m <- diag(0.8, n)
        m[cbind(2:n, 1:(n - 1))] <- 0.15; m[cbind(1:(n - 1), 2:n)] <- 0.05
        dat7$age_error[rows, 2 + seq_len(n)] <- m / rowSums(m)
      }
      dat7$spawn_month[] <- 3
      pc7 <- msfit(sp, dat = dat7, suit_styr = c(1985, 1979, 1990), suit_endyr = c(2010, 2017, 2005))
      pc7$data$MSSB0[] <- c(1e6, 5e5, 2e5); pc7$data$MSB0[] <- c(2e6, 1e6, 4e5)
      o7 <- build(pc7); o7A <- build(pc7, ADreport = TRUE)
      run_ms(o7, o7A, pc7, list(o7$cpp$par), "7 perturbed data")
      npts <- npts + 1
    }
    if (sp$suit > 0) {
      # Point 7b: other food 1e-9 of its value, so the stomachs' other-prey share sits at
      # its 1e-5 posfun floor
      dat7b <- BS2017MS
      dat7b$other_food <- dat7b$other_food * 1e-9
      pc7b <- msfit(sp, dat = dat7b)
      o7b <- build(pc7b); o7bA <- build(pc7b, ADreport = TRUE)
      set.seed(7)
      run_ms(o7b, o7bA, pc7b, list(o7b$cpp$par, o7b$cpp$par + rnorm(length(o7b$cpp$par), 0, 0.05)),
             "7b stomach floor (other food x 1e-9)")
      npts <- npts + 2
    }
    if (cf == "m1s0") {
      # Point 8: predators 1.5 log units more abundant, so M2 is large
      p8 <- pcg$params; p8$rec_pars[2:3, 1] <- p8$rec_pars[2:3, 1] + 1.5
      pc8 <- modifyList(pcg, list(params = p8))
      o8 <- build(pc8); o8A <- build(pc8, ADreport = TRUE)
      run_ms(o8, o8A, pc8, list(o8$cpp$par), "8 large M2")

      # Point 9: decision 7. Pollock year-1 numbers underflow at init_dev -720, so
      # available food is NaN: the objective is NaN on both sides, and a model recorded
      # there is the start's model back at the start
      at <- function(v, niter = 3) {
        q <- pc0$params; q$init_dev[1, seq_len(BS2017MS$nages[1] - 1)] <- v
        pq <- msfit(sp, inits = q, niter = niter)
        list(pc = pq, o = build(pq))
      }
      pc03 <- msfit(sp, inits = NULL, niter = 3)
      o03 <- build(pc03); o03A <- build(pc03, ADreport = TRUE)
      b9 <- at(-720)
      xq <- b9$o$cpp$par
      f9 <- c(b9$o$cpp$fn(xq), b9$o$rtmb$fn(xq))
      cat(sprintf("  9a init_dev -720: fn C++ %s RTMB %s; at the start, recorded there: C++ %.10f RTMB %.10f (start %.10f)\n",
                  f9[1], f9[2], b9$o$cpp$fn(x0), b9$o$rtmb$fn(x0), o03$cpp$fn(x0)))
      if (!is.nan(f9[1])) fail("reference predates #324: C++ is not NaN at init_dev -720")
      if (!is.nan(f9[2])) fail("9a: RTMB is not NaN at init_dev -720")
      if (abs(b9$o$cpp$fn(x0) - o03$cpp$fn(x0)) > 1e-12 * abs(o03$cpp$fn(x0)) ||
          abs(b9$o$rtmb$fn(x0) - o03$cpp$fn(x0)) > 1e-12 * abs(o03$cpp$fn(x0)))
        fail("9a: a model recorded at init_dev -720 is not the start's model at the start")
      b9b <- at(-720, niter = 1)
      rc9 <- b9b$o$cpp$report(b9b$o$cpp$par); rr9 <- b9b$o$rtmb$report(b9b$o$cpp$par)
      cat(sprintf("  9b the same at niter 1: fn C++ %.10f RTMB %.10f; NaN available food C++ %d RTMB %d, same cells %s\n",
                  b9b$o$cpp$fn(b9b$o$cpp$par), b9b$o$rtmb$fn(b9b$o$cpp$par), sum(is.na(rc9$avail_food)),
                  sum(is.na(rr9$avail_food)), identical(is.na(rc9$avail_food), is.na(rr9$avail_food))))
      if (!any(is.na(rc9$avail_food))) fail("9b: positive control failed, no NaN available food")
      run_ms(b9b$o, build(b9b$pc, ADreport = TRUE), b9b$pc, list(b9b$o$cpp$par), "9b niter 1 at init_dev -720")

      # Point 10: prey numbers exactly 0 (init_dev -800). The C++ decides that branch at
      # taping: its tape recorded at the start must be flagged as a flip there
      b10 <- at(-800)
      x10 <- b10$o$cpp$par
      run_ms(b10$o, build(b10$pc, ADreport = TRUE), b10$pc, list(x10), "10a init_dev -800 (its own tape)")
      cat(sprintf("  10b recorded at -800, at the start: C++ %.10f RTMB %.10f (start %.10f)\n",
                  b10$o$cpp$fn(x0), b10$o$rtmb$fn(x0), o03$cpp$fn(x0)))
      if (abs(b10$o$rtmb$fn(x0) - o03$cpp$fn(x0)) > 1e-12 * abs(o03$cpp$fn(x0)))
        fail("10b: RTMB recorded at -800 does not re-decide at the start")
      nf <- run_ms(o03, o03A, pc03, list(x10), "10c start tape at init_dev -800")
      if (nf != 1) fail("10c: the flip detector missed the exact-zero branch")
      npts <- npts + 6
    }
    if (cf == "m1s4") {
      dat11 <- BS2017MS; dat11$projyr <- dat11$endyr
      pc11 <- msfit(sp, inits = NULL, dat = dat11)
      o11 <- build(pc11); o11A <- build(pc11, ADreport = TRUE)
      set.seed(11)
      run_ms(o11, o11A, pc11, list(o11$cpp$par, o11$cpp$par + rnorm(length(o11$cpp$par), 0, 0.05)),
             "11 no projection years")
      npts <- npts + 2
    }

    cat("\n[derivatives: ADREPORT Jacobian, Hessian, standard errors]\n")
    deriv(o, oA, xg, "4 golden point")
    deriv(o5, o5A, o5$cpp$par, "5 perturbed parameters")
    deriv(o6, o6A, o6$cpp$par, "6 all-free map")
    deriv(o6, o6A, x6j, "6 all-free map, jittered")
    deriv(o, oA, jit[[2]][[length(jit[[2]])]], "2 jittered point b")
    if (sp$suit > 0) deriv(o7b, o7bA, o7b$cpp$par, "7b stomach floor")

    cat("\n[negative controls: each deliberate break must fail where named]\n")
    ctl <- list(
      m2s0 = list(list("Type III square dropped from available food", ".rtmb_predation",
                       "adj <- function(N) if (pm$msmMode == 2)", "adj <- function(N) if (FALSE)",
                       function() detects(build(pcg), xg))),
      m2s4 = list(list("Type III square added to the diet prediction", ".rtmb_diet_hat",
                       "NW <- Nc * pm$W", "NW <- Nc * Nc * pm$W",
                       function() detects(build(pcg), xg))),
      m1s4 = list(
        list("stomach weight exponentiated", ".rtmb_stomach_nll",
             "weighted[sp] <- -sum(diet_comp_weights[sp] * ll[j])",
             "weighted[sp] <- -sum(exp(diet_comp_weights[sp]) * ll[j])",
             function() detects(build(pcg), xg)),
        list("vulnerability transposed (point 5)", ".rtmb_suit_parametric",
             "vul[cbind(ps$sp_r, ps$sp_k)]", "vul[cbind(ps$sp_k, ps$sp_r)]",
             function() detects(build(pc5), o5$cpp$par))),
      m1s0 = list(
        list("suitability window one year late", ".rtmb_prepare_predation",
             "seq(int(d[[\"suit_styr\"]][sp]) - styr + 1L, int(d[[\"suit_endyr\"]][sp]) - styr + 1L)",
             "seq(int(d[[\"suit_styr\"]][sp]) - styr + 2L, int(d[[\"suit_endyr\"]][sp]) - styr + 2L)",
             function() detects(build(pcg), xg)),
        list("other-food floor dropped (cold start)", ".rtmb_suit_empirical",
             ".rtmb_max2(matrix(pm$ofdp, pm$nc, length(u)), 1e-06 * suma)",
             "matrix(pm$ofdp, pm$nc, length(u))", function() detects(build(pc0), x0)),
        list("dynamic B0 predation fed the fished numbers", ".rtmb_objective",
             "do.call(rbind, avgN_dB0)", "avgN_c", function() detects(build(pcg), xg)),
        list("Ftarget dropped from dynamic BF mortality (point 5)", ".rtmb_objective",
             "exp(-MF[, y - 1] - Ft[, y - 1])", "exp(-MF[, y - 1])",
             function() detects(build(pc5), o5$cpp$par)),
        list("one pass more, applying the last M2", ".rtmb_objective",
             "for (iter in seq_len(pd$niter))", "for (iter in seq_len(pd$niter + 1))",
             function() detects(build(pcg), xg)),
        list("initial age structure on M1 + M2", ".rtmb_objective",
             "m1 <- M1[[sp]][, 1]", "m1 <- M[[sp]][, 1]", function() detects(build(pcg), xg)),
        list("M_at_age reported as M1 + the last M2", ".rtmb_objective",
             "M_at_age <- to_array(M, nspp, nages)",
             "M_at_age <- to_array(mapply(function(a, b) a + b, M1, M2, SIMPLIFY = FALSE), nspp, nages)",
             function() detects(build(pcg), xg)),
        list("other food of the wrong species (point 7)", ".rtmb_prepare_predation",
             "other_food_c = other_food[cell_sp]", "other_food_c = rev(other_food)[cell_sp]",
             function() detects(build(pc7), o7$cpp$par)),
        list("prey test dropped (exact-zero prey, point 10a)", ".rtmb_suit_empirical",
             "ok <- 1 - .rtmb_ind_le(prey_denom, 0)", "ok <- 1 + 0 * prey_denom",
             function() detects(build(b10$pc), x10)),
        list("zero-N penalty carried across passes (deep floor)", ".rtmb_objective",
             "zero_N_pen <- numeric(nspp)", "zero_N_pen <- if (iter == 1) numeric(nspp) else zero_N_pen",
             function() detects(build(pcg), x_deep))))
    for (k in ctl[[cf]]) {
      null_ok <- !with_break(k[[2]], k[[3]], k[[3]], k[[5]]())
      caught <- with_break(k[[2]], k[[3]], k[[4]], k[[5]]())
      cat(sprintf("  %-62s %s (unbroken: %s)\n", k[[1]], if (caught) "caught" else "MISSED",
                  if (null_ok) "passes" else "FAILS"))
      if (!caught) fail(paste("negative control missed:", k[[1]]))
      if (!null_ok) fail(paste("negative control's unbroken round trip fails:", k[[1]]))
    }

    cat("\n[coverage]\n")
    cat(sprintf("  points compared: %d (nlminb path %d, jittered %d)\n", cover$n, length(p1), sum(lengths(jit))))
    cat(sprintf("  suit_other cells at the 1e-6 floor per point: %d-%d; max M2 %.3g-%.3g\n",
                min(mc$floor), max(mc$floor), min(mc$maxM2), max(mc$maxM2)))
    cat(sprintf("  points with DynamicSBF != DynamicSB0: %d; with NaN M2: %d; stomachs at the 1e-5 floor: %d-%d\n",
                sum(mc$dsbf), sum(mc$nan > 0), min(mc$stom_floor), max(mc$stom_floor)))
    cat("  REPORTs zero at every point:", paste(sort(names(which(cover$zero))), collapse = ", "), "\n")
    nz <- which(cover$jnz, arr.ind = TRUE)
    rows <- getFromNamespace(".JNLL_ROW_AXIS", "Rceattle")
    cat("  jnll_comp cells non-zero at some point:",
        paste(sprintf("%s[%d]", names(rows)[nz[, 1]], nz[, 2]), collapse = "; "), "\n")

    }
    if (speed != "no") {
    cat("\n[speed, interleaved; reference = the package's own -O2 DLL]\n")
    tape <- matrix(NA, 5, 2)
    for (i in 1:5) {
      tape[i, 1] <- system.time(TMB::MakeADFun(pcg$data, pcg$params, map = pcg$map, DLL = "ceattle", silent = TRUE))[["elapsed"]]
      tape[i, 2] <- system.time(make_rtmb(pcg$data, pcg$params, map = pcg$map))[["elapsed"]]
    }
    os <- build(pcg, dll = "ceattle")
    xa <- xg; xb <- xg + 1e-3
    gr_only <- function(obj, n = 200) system.time(for (i in seq_len(n)) obj$gr(if (i %% 2) xa else xb))[["elapsed"]] / n
    fn_gr <- function(obj, n = 200) system.time(for (i in seq_len(n)) { x <- if (i %% 2) xa else xb; obj$fn(x); obj$gr(x) })[["elapsed"]] / n
    pair <- gro <- matrix(NA, 5, 2)
    for (i in 1:5) { pair[i, ] <- c(fn_gr(os$cpp), fn_gr(os$rtmb)); gro[i, ] <- c(gr_only(os$cpp), gr_only(os$rtmb)) }
    rt <- apply(tape, 2, median); g <- apply(gro, 2, median); pr <- apply(pair, 2, median)
    cat(sprintf("  retape: C++ %.3f s  RTMB %.3f s  ratio %.2f (gate <= 2)\n", rt[1], rt[2], rt[2] / rt[1]))
    cat(sprintf("  fn+gr pair: C++ %.2f ms  RTMB %.2f ms  ratio %.2f (reported)\n", 1e3 * pr[1], 1e3 * pr[2], pr[2] / pr[1]))
    cat(sprintf("  gr alone: C++ %.2f ms  RTMB %.2f ms  ratio %.2f (gate <= 1.2)\n", 1e3 * g[1], 1e3 * g[2], g[2] / g[1]))
    if (rt[2] / rt[1] > 2) fail(sprintf("%s: retape ratio %.2f > 2", cf, rt[2] / rt[1]))
    if (g[2] / g[1] > 1.2) fail(sprintf("%s: gr ratio %.2f > 1.2", cf, g[2] / g[1]))
    }
    if (speed == "only") next
    os <- build(pcg, dll = "ceattle")

    cat("\n[fit: both optimized]\n")
    of <- if (cf == "m1s0") build(pc0, dll = "ceattle") else os
    xs <- if (cf == "m1s0") x0 else xg
    fit_one <- function(obj) suppressWarnings(nlminb(xs, obj$fn, obj$gr, control = list(eval.max = 1e4, iter.max = 1e4)))
    oc <- fit_one(of$cpp); orr <- fit_one(of$rtmb)
    rc <- of$cpp$report(oc$par); rr <- of$rtmb$report(orr$par)
    rel <- function(a, b) max(abs(a - b) / abs(a))
    cat(sprintf("  objective C++ %.10f  RTMB %.10f  (iterations %d, %d)\n", oc$objective, orr$objective,
                oc$iterations, orr$iterations))
    # M2 is scaled by its largest cell: weakly identified vulnerabilities (m1s4's arrowtooth
    # as prey) leave cells of 1e-10 and below that differ by orders of magnitude at equal objectives
    cat(sprintf("  max rel difference at the optima: ssb %.2e  biomass %.2e  R %.2e  M2/max(M2) %.2e; max|g| %.1e, %.1e\n",
                rel(rc$ssb, rr$ssb), rel(rc$biomass, rr$biomass), rel(rc$R, rr$R),
                max(abs(rc$M2_at_age - rr$M2_at_age)) / max(abs(rc$M2_at_age)),
                max(abs(of$cpp$gr(oc$par))), max(abs(of$rtmb$gr(orr$par)))))
    if (oc$convergence == 0) {
      if (orr$convergence != 0) fail(sprintf("%s: the RTMB optimization did not converge", cf))
      if (abs(oc$objective - orr$objective) > 1e-6 * abs(oc$objective)) fail(sprintf("%s: optimized objectives differ by more than 1e-6 relative", cf))
    } else {
      # No optimum to match when the C++ itself does not converge: gate both end points
      cat(sprintf("  NOT GATED: the C++ fit ends in \"%s\"; RTMB ends in \"%s\". Both end points are compared instead.\n",
                  oc$message, orr$message))
      pe <- if (cf == "m1s0") pc0 else pcg
      oe <- build(pe); oeA <- build(pe, ADreport = TRUE)
      run_ms(oe, oeA, pe, list(oc$par, orr$par), "fit end points (C++, RTMB)")
    }
  }

  cat("\n[scope guard: each input the C++ fits differently must be refused]\n")
  g0 <- msfit(specs$m1s0)
  g4 <- msfit(specs$m1s4)
  refusals <- list(
    "msmMode 3" = function(d) { d$msmMode <- 3; d },
    "suitMode 2 (gamma)" = function(d) { d$suitMode[] <- 2; d },
    "suitMode 6 (normal)" = function(d) { d$suitMode[] <- 6; d },
    "suitMode differing between predators" = function(d) { d$suitMode[] <- c(4, 0, 0); d },
    "suitMode 4 at msmMode 0" = function(d) { d$msmMode <- 0; d$suitMode[] <- 4; d },
    "a diet row for one year" = function(d) { d$diet_ctl[1, 7] <- 1990; d },
    "a diet row summed over prey ages" = function(d) { d$diet_ctl[1, 6] <- -1; d },
    "niter 0" = function(d) { d$niter <- 0; d },
    "minimum age 2" = function(d) { d$minage[1] <- 2; d },
    "negative other food" = function(d) { d$other_food[1] <- -1; d },
    "negative weight-at-age" = function(d) { d$weight_obs[d$pop_wt_index[1], 1, 2, 1] <- -1; d })
  for (nm in names(refusals)) {
    e <- try(make_rtmb(refusals[[nm]](g0$data), g0$params, map = g0$map), silent = TRUE)
    ok <- inherits(e, "try-error") && grepl("does not yet support", e)
    cat(sprintf("  %-45s %s\n", nm, if (ok) "refused" else "ACCEPTED"))
    if (!ok) fail(paste("scope guard accepted:", nm))
  }
  d4 <- g4$data; d4$diet_ll_type[] <- 1
  e <- try(make_rtmb(d4, g4$params, map = g4$map), silent = TRUE)
  ok <- inherits(e, "try-error") && grepl("does not yet support", e)
  cat(sprintf("  %-45s %s\n", "Dirichlet-multinomial diet", if (ok) "refused" else "ACCEPTED"))
  if (!ok) fail("scope guard accepted: Dirichlet-multinomial diet")
}

if (fixture == "bs2017ms") {
  run_bs2017ms(configs)
  finish()
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
deriv(o, oA, xg, "4 golden point")
deriv(o5, o5A, o5$cpp$par, "5 perturbed parameters")
deriv(o6, o6A, o6$cpp$par, "6 all-free map")
deriv(o6, o6A, x6j, "6 all-free map, jittered")
deriv(o, oA, jit[[2]][[length(jit[[2]])]], "2 jittered point b")
deriv(o, oA, jit[[1]][[length(jit[[1]])]], "2 jittered point")

# ---- negative controls ------------------------------------------------------------
cat("\n[negative controls: each deliberate break must fail where named]\n")
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

finish()
