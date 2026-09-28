# =============================================================================
# save_config() / load_config() -- a documented, git-diffable YAML round-trip
# for a full Rceattle run configuration.
# =============================================================================
#
# fit_mod() takes the model structure (predation mode, initialization, HCR, and
# the build_*() process specs) plus the estimation knobs (estimateMode,
# random_*, fit_control()) as call-time arguments, so a configuration lives only
# in the fit_mod() call. save_config()/load_config() round-trip that
# configuration to a YAML file with each field's doc string as a comment, so two
# runs diff cleanly and a configuration is archivable alongside an assessment.
# The parameter values (inits/map/bounds) are NOT stored here; only the model +
# estimation structure is. (Cf. SAM's saveConf/loadConf.)
#
# The content is a `run config`: the 12 model_config() fields, the estimation
# controls, and the fit_control() bundle. Everything serializes as plain
# scalars/lists EXCEPT the `formula` objects and Rceattle_prior S3 objects inside
# the build_*() specs' linkages, which get deparse()/{family,p1,p2} handling.

# ---- The run-config object ---------------------------------------------------

# The estimation controls captured alongside model_config (fit_mod() args that
# are NOT part of model_config()), with their fit_mod() defaults.
.RCE_RUN_ESTIMATION_FIELDS <- list(
  estimateMode = 0, random_rec = FALSE, random_q = FALSE, random_sel = FALSE,
  suit_styr = NULL, suit_endyr = NULL
)

# model_config field -> the build_*() constructor that produces / rebuilds it.
.RCE_CONFIG_BUILDERS <- c(
  HCR = "build_hcr", recFun = "build_srr", M1Fun = "build_M1",
  growthFun = "build_growth", qFun = "build_catchability",
  selFun = "build_selectivity", compFun = "build_composition"
)

#' Assemble an Rceattle run-configuration object.
#'
#' @param model_config An `Rceattle_model_config` (from [model_config()]).
#' @param estimateMode,random_rec,random_q,random_sel,suit_styr,suit_endyr
#'   Estimation controls (see [fit_mod()]).
#' @param fit_control An `Rceattle_fit_control` (from [fit_control()]).
#' @return An object of class `Rceattle_run_config`.
#' @keywords internal
#' @noRd
.rce_run_config <- function(mc = model_config(),
                            estimateMode = 0, random_rec = FALSE,
                            random_q = FALSE, random_sel = FALSE,
                            suit_styr = NULL, suit_endyr = NULL,
                            fc = fit_control()) {
  # `mc`/`fc` argument names (not `model_config`/`fit_control`) so the defaults
  # call the constructors rather than shadowing them with the argument promise.
  structure(
    list(model_config = mc, estimateMode = estimateMode,
         random_rec = random_rec, random_q = random_q, random_sel = random_sel,
         suit_styr = suit_styr, suit_endyr = suit_endyr,
         fit_control = fc),
    class = c("Rceattle_run_config", "list")
  )
}

#' @export
print.Rceattle_run_config <- function(x, ...) {
  cat("<Rceattle run_config>\n")
  print(x$model_config)
  cat("  estimateMode:", x$estimateMode, "| random_rec:", x$random_rec,
      "| random_q:", x$random_q, "| random_sel:", x$random_sel, "\n")
  cat("  fit_control : getsd =", x$fit_control$getsd, ", phase =",
      paste(unlist(x$fit_control$phase), collapse = ","), "\n")
  invisible(x)
}

# ---- Serialization primitives ------------------------------------------------

# A one-sided formula <-> its text. deparse1 collapses multi-line formulas to a
# single string; as.formula rebuilds it (dropping the captured environment,
# which the linkage machinery does not use).
.rce_formula_to_chr <- function(f) if (is.null(f)) NULL else paste(deparse(f), collapse = " ")
.rce_chr_to_formula <- function(s) if (is.null(s)) NULL else stats::as.formula(s, env = baseenv())

# An Rceattle_prior <-> {family, p1, p2}. Priors can also be a (possibly nested)
# named list of prior objects, so recurse.
.rce_prior_to_list <- function(p) {
  if (is.null(p)) return(NULL)
  if (inherits(p, "Rceattle_prior"))
    return(list(family = p$family, p1 = p$p1, p2 = p$p2))
  if (is.list(p)) return(lapply(p, .rce_prior_to_list))
  p
}
.rce_prior_from_list <- function(l) {
  if (is.null(l)) return(NULL)
  if (is.list(l) && all(c("family", "p1", "p2") %in% names(l)))
    return(new_prior(l$family, l$p1, l$p2))
  if (is.list(l)) return(lapply(l, .rce_prior_from_list))
  l
}

# A single Rceattle_linkage_spec <-> a plain list. The constructor slots that are
# also constructor arguments are round-tripped; `formula`/`by` become text and
# `priors` become {family,p1,p2}. NULL slots are dropped so a spec re-built with
# the constructor picks up the same defaults.
.RCE_LINKAGE_SPEC_ARGS <- c("formula", "param", "by", "species", "sex", "fleet",
                            "link", "init", "bounds", "priors", "re_group",
                            "est_phase", "observe", "obs_sd", "obs_sd_est",
                            "integrate")

.rce_spec_to_list <- function(spec) {
  # A reference spec for default-omission (formula/by differ and are handled
  # separately). `param` is dropped -- the build_*() validators re-stamp it from
  # the linkage list key on reconstruction.
  ref <- linkage_spec(~ 1)
  out <- list()
  for (nm in .RCE_LINKAGE_SPEC_ARGS) {
    if (nm == "param") next
    v <- spec[[nm]]
    # `by` is handled before the NULL skip. An auto-filled `by` (the user omitted
    # it) is dropped -- reconstructing with the omitted argument re-derives the same
    # per-process base stratum. An *explicit* `by` is always recorded so it survives
    # the round-trip unchanged, including `by = NULL` (a single shared coefficient),
    # via a sentinel. Keying on `by_auto` -- not on equality with `~ species` -- is
    # what keeps an explicit `by = ~ species` on a fleet-keyed process (whose default
    # is now `~ fleet`) from being silently re-resolved on reload.
    if (nm == "by") {
      if (!isTRUE(spec$by_auto)) {
        out$by <- if (is.null(v)) "NULL" else .rce_formula_to_chr(v)
      }
      next
    }
    if (is.null(v)) next
    if (nm == "formula") {
      out[[nm]] <- .rce_formula_to_chr(v)
    } else if (nm == "priors") {
      pv <- .rce_prior_to_list(v)
      if (length(pv) > 0) out$priors <- pv    # drop an empty priors list
    } else if (!identical(v, ref[[nm]])) {     # drop values left at the default
      out[[nm]] <- v
    }
  }
  out
}
.rce_spec_from_list <- function(l) {
  args <- l
  if (!is.null(args$formula)) args$formula <- .rce_chr_to_formula(args$formula)
  if (!is.null(args$priors))  args$priors  <- .rce_prior_from_list(args$priors)
  # `by`: the "NULL" sentinel means shared coefficients (by = NULL); otherwise a
  # deparsed formula string. `args["by"] <- list(NULL)` keeps an explicit NULL
  # element (plain `args$by <- NULL` would delete it and fall back to ~species).
  if (identical(args$by, "NULL")) args["by"] <- list(NULL)
  else if (!is.null(args$by))     args$by <- .rce_chr_to_formula(args$by)
  do.call(linkage_spec, args)
}

# A `linkages` value (NULL, one spec, or a named list of specs / lists of specs)
# <-> plain lists.
.rce_linkages_to_list <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "Rceattle_linkage_spec")) return(.rce_spec_to_list(x))
  if (is.list(x)) return(lapply(x, .rce_linkages_to_list))
  x
}
.rce_linkages_from_list <- function(x) {
  if (is.null(x)) return(NULL)
  # A serialized single spec is a plain list carrying a `formula` field.
  if (is.list(x) && !is.null(x$formula) &&
      (is.character(x$formula) || inherits(x$formula, "formula")))
    return(.rce_spec_from_list(x))
  if (is.list(x)) return(lapply(x, .rce_linkages_from_list))
  x
}

# A build_*() spec object <-> plain list. Only fields that are constructor
# ARGUMENTS are kept (build_growth() adds derived, non-argument fields); a field
# equal to the builder's default is dropped so the file stays minimal and its
# reconstruction re-derives the rest.
.rce_build_to_list <- function(obj, builder_name) {
  builder <- get(builder_name, envir = asNamespace("Rceattle"))
  arg_names <- names(formals(builder))
  default <- builder()
  out <- list()
  for (nm in intersect(names(obj), arg_names)) {
    v <- obj[[nm]]
    if (nm == "linkages") {
      lk <- .rce_linkages_to_list(v)
      if (!is.null(lk) && length(lk) > 0) out$linkages <- lk
      next
    }
    if (!identical(v, default[[nm]])) out[[nm]] <- v
  }
  out
}
.rce_build_from_list <- function(l, builder_name) {
  builder <- get(builder_name, envir = asNamespace("Rceattle"))
  args <- l
  if (!is.null(args$linkages)) args$linkages <- .rce_linkages_from_list(args$linkages)
  do.call(builder, args)
}

# ---- Full run config <-> nested plain list (YAML-ready) ----------------------

#' Convert a run config to a plain nested list ready for YAML.
#' @keywords internal
#' @noRd
.rce_run_config_to_list <- function(rc) {
  mc <- rc$model_config
  mc_def <- model_config()
  # A field the config set is written even at its default, so a saved run
  # imposes it on reload; an unset default field is left out.
  set <- attr(mc, "set") %||% character(0)
  model <- list()
  for (nm in .RCE_MODEL_CONFIG_FIELDS) {
    if (nm %in% names(.RCE_CONFIG_BUILDERS)) {
      b <- .rce_build_to_list(mc[[nm]], unname(.RCE_CONFIG_BUILDERS[nm]))
      if (length(b) > 0 || nm %in% set) model[[nm]] <- b
    } else if (nm %in% set || !identical(mc[[nm]], mc_def[[nm]])) {
      model[nm] <- list(mc[[nm]])   # a set field at NULL is written as null
    }
  }

  est <- list()
  for (nm in names(.RCE_RUN_ESTIMATION_FIELDS))
    if (!identical(rc[[nm]], .RCE_RUN_ESTIMATION_FIELDS[[nm]])) est[[nm]] <- rc[[nm]]

  fc <- list(); fc_def <- fit_control()
  for (nm in names(fc_def))
    if (!identical(rc$fit_control[[nm]], fc_def[[nm]])) fc[[nm]] <- rc$fit_control[[nm]]

  out <- list()
  if (length(model) > 0) out$model       <- model
  if (length(est)   > 0) out$estimation  <- est
  if (length(fc)    > 0) out$fit_control <- fc
  out
}

#' Rebuild a run config from a plain nested list (parsed YAML).
#' @keywords internal
#' @noRd
.rce_run_config_from_list <- function(l) {
  model <- l$model %||% list()
  mc_args <- list()
  for (nm in .RCE_MODEL_CONFIG_FIELDS) {
    if (!nm %in% names(model)) next
    mc_args[nm] <- list(if (nm %in% names(.RCE_CONFIG_BUILDERS))
      .rce_build_from_list(model[[nm]] %||% list(), unname(.RCE_CONFIG_BUILDERS[nm]))
    else model[[nm]])
  }
  mc <- do.call(model_config, mc_args)

  est <- l$estimation %||% list()
  fc  <- do.call(fit_control, l$fit_control %||% list())

  .rce_run_config(
    mc = mc,
    estimateMode = est$estimateMode %||% .RCE_RUN_ESTIMATION_FIELDS$estimateMode,
    random_rec   = est$random_rec   %||% .RCE_RUN_ESTIMATION_FIELDS$random_rec,
    random_q     = est$random_q     %||% .RCE_RUN_ESTIMATION_FIELDS$random_q,
    random_sel   = est$random_sel   %||% .RCE_RUN_ESTIMATION_FIELDS$random_sel,
    suit_styr    = est$suit_styr, suit_endyr = est$suit_endyr,
    fc           = fc
  )
}

# ---- Config-field doc dictionary (for the YAML comments) ----------------------
#
# The column schema (R/0-column_schema.R) documents workbook columns; the
# control-level run-config switches are documented here, with allowed readable
# values pulled from the switch maps in R/0-switches.R.

#' One-line docs (and allowed values) for the run-config fields.
#' @keywords internal
#' @noRd
.rce_config_schema <- function() {
  # Name AND integer code for every value. A saved config writes the switch as
  # the number (`msmMode: 1.0`), so a comment listing only the readable names
  # cannot be used to decode the line it sits above.
  allowed <- function(map) {
    paste0(names(map), " (", unname(map), ")", collapse = " / ")
  }
  d <- function(doc, allowed = NULL) list(doc = doc, allowed = allowed)

  # Model-level switches project from .rce_model_switch_schema(), so the
  # comments written into a saved config list every value the map defines. They
  # were hand-kept here, which is how estimateMode's DebugOptimize came to be a
  # real mode that the YAML comment never mentioned.
  .from_switch_table <- function(nm) {
    row <- .rce_model_switch_schema()[[nm]]
    if (is.null(row)) return(NULL)
    map <- if (!is.na(row$allowed)) {
      tryCatch(get(row$allowed, envir = asNamespace("Rceattle")), error = function(e) NULL)
    } else NULL
    d(row$doc, if (is.null(map)) NULL else allowed(map))
  }
  list(
    # model_config scalars
    msmMode  = .from_switch_table("msmMode"),
    initMode = .from_switch_table("initMode"),
    avgnMode = .from_switch_table("avgnMode"),
    suitMode = .from_switch_table("suitMode"),
    niter    = d("Number of predation iterations"),
    HCR = .from_switch_table("HCR"),
    recFun   = d("Stock-recruit specification (build_srr)"),
    M1Fun    = d("Natural-mortality specification (build_M1)"),
    growthFun= d("Growth specification (build_growth)"),
    qFun     = d("Catchability specification (build_catchability)"),
    selFun   = d("Selectivity specification (build_selectivity)"),
    compFun  = d("Composition specification (build_composition)"),
    # estimation controls
    estimateMode = .from_switch_table("estimateMode"),
    random_rec = d("Estimate recruitment deviations as random effects"),
    random_q   = d("Integrate the Time_varying_q deviations and estimate their sd (linkages integrate either way)"),
    random_sel = d("Integrate the Time_varying_sel deviations and estimate their sd (linkages integrate either way)"),
    suit_styr  = d("First year of the diet/suitability averaging window"),
    suit_endyr = d("Last year of the diet/suitability averaging window"),
    # fit_control knobs (the commonly-tuned ones)
    getsd            = d("Run TMB::sdreport() after optimization"),
    bias.correct     = d("Apply epsilon bias correction in sdreport"),
    getJointPrecision= d("Return the joint fixed+random precision"),
    phase            = d("Phase the optimization (TRUE/FALSE or a phase list)"),
    loopnum          = d("Number of optimizer restarts"),
    newtonsteps      = d("Extra Newton steps after optimization"),
    comp_offset      = d("Composition proportion offset (avoids log(0))"),
    verbose          = d("0 silent / 1 fit updates / 2 fit + TMB progress")
  )
}

# Emit a YAML string with a `# doc` comment above every documented key.
.rce_yaml_with_comments <- function(l) {
  dict <- .rce_config_schema()
  body <- yaml::as.yaml(l, indent = 2)
  out <- character(0)
  for (ln in strsplit(body, "\n", fixed = TRUE)[[1]]) {
    m <- regmatches(ln, regexec("^(\\s*)([A-Za-z0-9_.]+):", ln))[[1]]
    if (length(m) == 3 && !is.null(dict[[m[3]]])) {
      doc <- dict[[m[3]]]$doc
      al  <- dict[[m[3]]]$allowed
      cmt <- if (is.null(al)) doc else paste0(doc, "  (", al, ")")
      out <- c(out, paste0(m[2], "# ", cmt))
    }
    out <- c(out, ln)
  }
  paste(out, collapse = "\n")
}

# ---- Extract a run config from a fit / data list / config object --------------

#' Extract the run configuration from a fit, data list, or config object
#'
#' Returns the [model_config()] structure plus the estimation controls and
#' [fit_control()] bundle as a single `Rceattle_run_config`. Accepts a fitted
#' Rceattle object, a data list holding `$model_config`, an
#' `Rceattle_run_config`, or an `Rceattle_model_config`. Estimation controls and
#' `fit_control` supplied via `...` override any found on the object.
#'
#' @param x A fitted Rceattle object, a data list, an `Rceattle_run_config`, or
#'   an `Rceattle_model_config`.
#' @param ... Estimation controls / `fit_control` to override
#'   (`estimateMode`, `random_rec`, `random_q`, `random_sel`, `suit_styr`,
#'   `suit_endyr`, `fit_control`).
#' @return An `Rceattle_run_config`.
#' @examples
#' # Coerce a model_config() into a full run configuration, overriding the
#' # estimation controls. A fitted model or a data list works the same way.
#' rc <- run_config(model_config(msmMode = 1), estimateMode = "Hindcast")
#' rc$estimateMode
#'
#' # Optimizer and uncertainty settings belong inside fit_control(), not here;
#' # an unrecognized field is an error rather than a silent no-op.
#' rc <- run_config(rc, fit_control = fit_control(getsd = FALSE))
#' @seealso [save_config()], [load_config()], [model_config()].
#' @export
run_config <- function(x, ...) {
  ov <- list(...)
  if (inherits(x, "Rceattle_run_config")) {
    rc <- x
  } else if (inherits(x, "Rceattle_model_config")) {
    rc <- .rce_run_config(mc = x)
  } else {
    # A fit records $run_config (set by fit_mod); otherwise fall back to the
    # data list's $model_config (+ $estimateMode).
    if (!is.null(x$run_config) && inherits(x$run_config, "Rceattle_run_config")) {
      rc <- x$run_config
    } else {
      dl <- if (!is.null(x$data_list)) x$data_list else x
      mc <- dl$model_config %||% model_config()
      rc <- .rce_run_config(mc = mc, estimateMode = dl$estimateMode %||% 0)
    }
  }
  # Apply any explicit overrides. Validate the NAMES against the fixed field set
  # rather than names(rc): assigning NULL drops an element, so a config that has
  # already had a field cleared would otherwise reject a later legitimate
  # override of it. Values are not validated here -- fit_mod() does that.
  valid <- c("model_config", names(.RCE_RUN_ESTIMATION_FIELDS), "fit_control")
  if (length(ov) && (is.null(names(ov)) || any(!nzchar(names(ov))))) {
    stop("All `...` arguments to run_config() must be named.", call. = FALSE)
  }
  bad <- setdiff(names(ov), valid)
  if (length(bad)) {
    stop("Unrecognised run_config() field(s): ", paste(bad, collapse = ", "),
         ".\nValid fields: ", paste(valid, collapse = ", "),
         ".\n(Optimizer and uncertainty settings go inside fit_control().)",
         call. = FALSE)
  }
  for (nm in names(ov)) rc[[nm]] <- ov[[nm]]
  rc
}

# ---- save_config() / load_config() -------------------------------------------

#' Save a model run configuration to a documented YAML file
#'
#' Round-trips a full run configuration, the [model_config()] structure plus
#' the estimation controls and [fit_control()] bundle, to a
#' git-diffable YAML file, with each field's documentation emitted as a comment
#' and a spec-tree + provenance header. Only fields that differ from their
#' defaults are written, so two configurations diff to just their real
#' differences. The parameter values (`inits`/`map`/`bounds`) are NOT stored;
#' pair the config with a saved fit for those.
#'
#' @param x A fitted Rceattle object, a data list holding `$model_config`, an
#'   `Rceattle_run_config`, or an `Rceattle_model_config`.
#' @param file Output path for the `.yaml` file.
#' @param ... Estimation controls / `fit_control` to record (passed to
#'   [run_config()]).
#' @return Invisibly, the `Rceattle_run_config` that was written.
#' @seealso [load_config()], [run_config()], [model_config()], [fit_mod()].
#' @examples
#' cfg <- model_config(msmMode = 1, initMode = "OffsetEquilibrium")
#' f <- file.path(tempdir(), "run.yaml")
#' save_config(cfg, f)
#' identical(load_config(f)$model_config$msmMode, 1)
#' @export
save_config <- function(x, file = "Rceattle_config.yaml", ...) {
  rc <- run_config(x, ...)
  l  <- .rce_run_config_to_list(rc)

  # Provenance + human spec-tree header (comment lines; ignored on read).
  header <- c(
    "# Rceattle run configuration",
    paste0("# saved: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("# Rceattle: ",
           tryCatch(as.character(utils::packageVersion("Rceattle")),
                    error = function(e) "dev")),
    "#",
    paste0("# ", utils::capture.output(print(rc$model_config))),
    "")
  writeLines(c(header, .rce_yaml_with_comments(l)), file)
  invisible(rc)
}

#' Load a model run configuration from a YAML file
#'
#' The inverse of [save_config()]: reads a run-config YAML file and rebuilds the
#' [model_config()] structure, estimation controls, and [fit_control()] bundle
#' (reconstructing linkage formulas and priors). Attach the result to a fit with
#' `fit_mod(data_list, config = load_config("run.yaml"))`, or read
#' `$model_config` off it to attach to a data list.
#'
#' @param file Path to a YAML file written by [save_config()].
#' @return An `Rceattle_run_config`.
#' @examples
#' f <- file.path(tempdir(), "run.yaml")
#' save_config(model_config(msmMode = 1, initMode = "Equilibrium"), f)
#'
#' cfg <- load_config(f)
#' cfg$model_config$msmMode
#' # Apply with fit_mod(data_list, config = cfg).
#' unlink(f)
#' @seealso [save_config()], [fit_mod()].
#' @export
load_config <- function(file) {
  if (!file.exists(file))
    stop("config file not found: '", file, "'", call. = FALSE)
  .rce_run_config_from_list(yaml::read_yaml(file))
}
