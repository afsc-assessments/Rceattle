# =============================================================================
# Rceattle_data S3 class + spec-tree printer
# =============================================================================
#
# build_data() tags its result with the class "Rceattle_data" so that printing a
# data object shows an indented spec tree (dimensions -> fleets -> per-process
# form -> active linkages -> configuration) instead of dumping the ~40-element
# list. The same renderer feeds print.Rceattle() for a fitted model. The class
# is a thin tag: the object is still a plain list to every function that consumes
# it (clean_data(), data_check(), fit_mod(), write_data() all treat it as a
# list), so tagging cannot change a fit.

# Readable labels for the small fixed set of integer-coded fleet columns, used
# only when switch_check() could not normalise a (partial) object to strings.
# Casing matches fleet_map ("Off"/"Fishery"/"Survey") so the fallback prints the
# same tokens as a switch-checked object.
.RCE_FLEET_TYPE_LABEL <- c("0" = "Off", "1" = "Fishery", "2" = "Survey")

# Tag a data list so print()/summary() render the spec tree instead of dumping
# the ~40-element list. Applied by every constructor (build_data, read_data,
# clean_data, combine_data) so the tag survives a write_data()/read_data()
# round-trip. Idempotent, and inert: the object is still a plain list to every
# consumer, and the tag is dropped before anything reaches TMB::MakeADFun(), so
# it cannot change a fit.
.rce_as_data <- function(dl) {
  if (!is.list(dl)) return(dl)
  class(dl) <- unique(c("Rceattle_data", class(dl)))
  dl
}

# Show a switch value by its string alias, so the tree reads the same regardless of
# the underlying integer code (and keeps meaning if a code is ever renumbered). `map`
# is a name(string) -> integer vector; `x` may be strings (kept as-is) or integers
# (reverse-mapped). Unmapped values fall back to their raw form.
.rce_alias <- function(x, map) {
  vapply(x, function(v) {
    if (is.na(v)) return(NA_character_)
    v <- as.character(v)
    if (v %in% names(map)) return(v)                    # already a string alias
    hit <- names(map)[match(v, as.character(map))]      # integer code -> string
    if (!is.na(hit)) hit else v                         # fall back to the raw value
  }, character(1))
}

# Alias a (possibly per-species) switch vector for display: collapse to one token
# when every element is the same, else list them comma-separated.
.rce_alias_show <- function(x, map) {
  s <- .rce_alias(x, map)
  if (length(unique(s)) == 1L) s[1] else paste(s, collapse = ", ")
}

# Every linkage_spec held under one process slot, flattened to a list.
#
# A slot holds either a single spec or a LIST of them: the grammar lets one
# parameter carry several, which is how a shared prior and a fleet-specific
# random walk sit on the same selectivity limb. Reading `$formula` off that
# list gives NULL, so a stacked parameter has to be unpacked before display --
# and stacking is what a real assessment does.
#
# `is.list()` first, and not merely for tidiness: `$` on an atomic vector is an
# error, not NULL, so testing `x$formula` on a stray numeric would throw out of
# the whole tree. Both print methods catch that and fall back to "(spec tree
# unavailable)", which loses every row rather than the one that is malformed.
.rce_linkage_specs <- function(x) {
  if (is.null(x) || !is.list(x) || !length(x)) return(list())
  if (inherits(x, "Rceattle_linkage_spec") || !is.null(x$formula)) return(list(x))
  out <- unlist(lapply(x, .rce_linkage_specs), recursive = FALSE)
  if (is.null(out)) list() else out
}

# One spec's formula as a single line. deparse() splits a long formula across
# elements, so collapse rather than taking the first and silently truncating.
.rce_linkage_formula_text <- function(spec) {
  f <- tryCatch(spec$formula, error = function(e) NULL)
  if (is.null(f)) return("(no formula)")
  txt <- tryCatch(paste(trimws(deparse(f)), collapse = " "), error = function(e) "")
  if (!nzchar(txt)) "(no formula)" else txt
}


#' Render an Rceattle data list as an indented spec tree
#'
#' @param dl A data list (from build_data()/read_data(), or a fitted model's
#'   `$data_list`).
#' @param config Append the attached [model_config()] block, when present.
#' @return A character vector, one element per line (invisibly printable).
#' @keywords internal
#' @noRd
.rce_spec_tree <- function(dl, config = TRUE) {
  `%||%` <- function(a, b) if (is.null(a)) b else a

  # Normalise switches to their canonical string forms so the tree is readable
  # whether it is fed a raw (integer-coded) build_data() object or an already
  # switch-checked fitted data_list. Fall back to the raw object if a partial
  # object cannot be normalised.
  dl2 <- tryCatch(
    suppressWarnings(suppressMessages(Rceattle::switch_check(dl))),
    error = function(e) dl
  )

  L <- character(0)
  add <- function(...) L[[length(L) + 1L]] <<- paste0(...)

  # ---- header + dimensions --------------------------------------------------
  spn <- dl2$spnames %||% (if (!is.null(dl2$nspp)) seq_len(dl2$nspp) else NULL)
  nspp <- dl2$nspp %||% length(spn)
  add("  dimensions")
  add("  \u251c\u2500 species  : ", nspp %||% "?",
      if (!is.null(spn)) paste0("  (", paste(spn, collapse = ", "), ")") else "")
  yr <- paste0(dl2$styr %||% "?", "\u2013", dl2$endyr %||% "?")
  if (!is.null(dl2$projyr) && !identical(dl2$projyr, dl2$endyr)) {
    yr <- paste0(yr, "  (projection \u2192 ", dl2$projyr, ")")
  }
  add("  \u251c\u2500 years    : ", yr)
  agr <- if (!is.null(dl2$minage) && !is.null(dl2$nages)) {
    paste0(dl2$minage, "\u2013", dl2$minage + dl2$nages - 1)
  } else as.character(dl2$nages %||% "?")
  add("  \u2514\u2500 ages     : ", paste(agr, collapse = ", "),
      "  |  lengths ", paste(dl2$nlengths %||% "?", collapse = ", "),
      "  |  sexes ", paste(dl2$nsex %||% "?", collapse = ", "))

  # ---- fleets ---------------------------------------------------------------
  # Coerce to a data.frame: a fitted model$data_list holds the original
  # data.frame, but harden against a rearranged (matrix) fleet_control so a
  # printer can never crash auto-printing a model.
  fc <- tryCatch(as.data.frame(dl2$fleet_control, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (!is.null(fc) && nrow(fc) > 0) {
    add("  fleets (", nrow(fc), ")")
    lab <- function(col, i) {
      v <- if (col %in% colnames(fc)) fc[[col]][i] else NA
      if (is.na(v)) return(NA_character_)
      if (col == "Fleet_type" && v %in% names(.RCE_FLEET_TYPE_LABEL)) {
        return(unname(.RCE_FLEET_TYPE_LABEL[as.character(v)]))
      }
      as.character(v)
    }
    # Which Selectivity_index / Catchability_index values are shared (mirrored)?
    shared <- function(col) {
      if (!col %in% colnames(fc)) return(integer(0))
      t <- table(fc[[col]][!is.na(fc[[col]])])
      as.integer(names(t[t > 1]))
    }
    # Accept the legacy `Q_index` spelling for display of an un-upgraded list.
    q_col <- if ("Catchability_index" %in% colnames(fc)) "Catchability_index" else "Q_index"
    sel_mir <- shared("Selectivity_index"); q_mir <- shared(q_col)
    # Fleets whose catchability block is actually opened: those carrying fitted
    # index observations. A fleet with none gets no q however its Catchability
    # column reads -- a q with no index to inform it is a flat direction -- so
    # say so here rather than leave the column looking authoritative. It can
    # still be estimated by following the lead of a shared q group, which the
    # [shared: q<->N] annotation on the same line shows.
    q_data <- tryCatch(.fleets_with_index(dl2), error = function(e) integer(0))
    n <- nrow(fc)
    # Two fixed left-hand columns -- "[code] name" and fleet type -- are padded to a
    # common width so the " - "-separated form fields line up down the block; the
    # optional sel:/q:/mirror fields trail unpadded.
    id_col   <- vapply(seq_len(n),
                       function(i) paste0("[", fc$Fleet_code[i] %||% i, "] ", fc$Fleet_name[i] %||% ""),
                       character(1))
    type_col <- vapply(seq_len(n),
                       function(i) { t <- lab("Fleet_type", i); if (is.na(t)) "" else t },
                       character(1))
    idw   <- max(nchar(id_col))
    typew <- max(nchar(type_col))
    for (i in seq_len(n)) {
      elbow <- if (i == n) "  \u2514\u2500 " else "  \u251c\u2500 "
      sel   <- lab("Selectivity", i)
      qform <- lab("Catchability", i)
      cells <- c(formatC(id_col[i],   width = -idw),
                 formatC(type_col[i], width = -typew))
      if (!is.na(sel))   cells <- c(cells, paste0("sel: ", sel))
      if (!is.na(qform)) {
        qcell <- paste0("q: ", qform)
        q_would_est <- !(qform %in% c("Fixed", "Analytical", "AnalyticalArith"))
        if (q_would_est && !((fc$Fleet_code[i] %||% i) %in% q_data)) {
          qcell <- paste0(qcell, " (fixed: no index data)")
        }
        cells <- c(cells, qcell)
      }
      # Mirror annotations: fleets that share a selectivity / catchability block.
      mir <- c()
      si <- if ("Selectivity_index" %in% colnames(fc)) fc$Selectivity_index[i] else NA
      qi <- if (q_col %in% colnames(fc)) fc[[q_col]][i] else NA
      if (!is.na(si) && si %in% sel_mir) mir <- c(mir, paste0("sel\u2194", si))
      if (!is.na(qi) && qi %in% q_mir)   mir <- c(mir, paste0("q\u2194", qi))
      if (length(mir)) cells <- c(cells, paste0("[shared: ", paste(mir, collapse = ", "), "]"))
      add(elbow, paste(cells, collapse = " \u2502 "))
    }
  }

  # ---- processes (present on a configured / fitted object) ------------------
  # Show each switch by its string alias (via .rce_alias_show) so the form is legible
  # and stable if the underlying integer code is ever renumbered. A per-species vector
  # that is constant collapses to a single token, else lists one entry per species.
  procs <- c()
  if (!is.null(dl2$srr_fun))      procs <- c(procs, paste0("recruitment  srr_fun      = ", .rce_alias_show(dl2$srr_fun, .SRR_FUNS)))
  if (!is.null(dl2$M1_model))     procs <- c(procs, paste0("M1           M1_model     = ", .rce_alias_show(dl2$M1_model, .M1_MODELS)))
  if (!is.null(dl2$growth_model)) procs <- c(procs, paste0("growth       growth_model = ", .rce_alias_show(dl2$growth_model, .GROWTH_FUN_TO_INT)))
  if (length(procs)) {
    add("  processes")
    for (k in seq_along(procs)) {
      add(if (k == length(procs)) "  \u2514\u2500 " else "  \u251c\u2500 ", procs[k])
    }
  }

  # ---- active linkages ------------------------------------------------------
  link_slots <- c(q = "q_linkages", sel = "sel_linkages", growth = "growth_linkages",
                  M1 = "M1_linkages", srr = "srr_linkages", comp = "comp_linkages")
  # One line per SPEC, not per parameter, so the count is the number of linkages
  # the model actually carries. A parameter holding more than one is tagged
  # [i/n], which says the specs are stacked on the same parameter without
  # implying they are one formula -- they usually name different fleets.
  link_lines <- c()
  for (nm in names(link_slots)) {
    lk <- dl2[[link_slots[[nm]]]]
    if (is.null(lk) || length(lk) == 0) next
    for (j in seq_along(lk)) {
      par_nm <- names(lk)[j] %||% nm
      if (is.na(par_nm) || !nzchar(par_nm)) par_nm <- nm
      specs <- .rce_linkage_specs(lk[[j]])
      if (!length(specs)) next
      for (k in seq_along(specs)) {
        stack <- if (length(specs) > 1L) paste0(" [", k, "/", length(specs), "]") else ""
        link_lines <- c(link_lines, paste0(nm, ": ", par_nm, stack, " ",
                                           .rce_linkage_formula_text(specs[[k]])))
      }
    }
  }
  if (length(link_lines)) {
    add("  linkages (", length(link_lines), ")")
    for (k in seq_along(link_lines)) {
      add(if (k == length(link_lines)) "  \u2514\u2500 " else "  \u251c\u2500 ", link_lines[k])
    }
  }

  # ---- model_config ---------------------------------------------------------
  # msmMode / initMode / HCR are aliased for the same reason as the processes above.
  cfg <- dl2$model_config
  if (isTRUE(config) && !is.null(cfg)) {
    add("  model_config : msmMode ",
        if (!is.null(cfg$msmMode)) .rce_alias_show(cfg$msmMode, msmMode_map) else "?",
        " \u00b7 initMode ",
        if (!is.null(cfg$initMode)) .rce_alias_show(cfg$initMode, initMode_map) else "?",
        if (!is.null(cfg$HCR$HCR)) paste0(" \u00b7 HCR ", .rce_alias_show(cfg$HCR$HCR, hcr_map)) else "")
  }

  L
}

#' Print method for an Rceattle data list
#'
#' Shows the model specification as an indented tree, dimensions, fleets and
#' their selectivity / catchability forms, configured processes, active
#' linkages, and any attached [model_config()], rather than dumping the full
#' data list.
#'
#' `print()` includes the configuration block; `summary()` omits it and reports
#' the data structure alone. Either can be overridden with `config`.
#'
#' @param x An `"Rceattle_data"` object from [build_data()].
#' @param config Show the attached [model_config()] block. Defaults to `TRUE`
#'   for `print()` and `FALSE` for `summary()`. Has no effect on an object that
#'   holds no configuration.
#' @param ... Currently unused.
#' @return `x`, invisibly.
#' @export
print.Rceattle_data <- function(x, config = TRUE, ...) {
  cat("<Rceattle data>\n")
  tree <- tryCatch(.rce_spec_tree(x, config = config),
                   error = function(e) paste0("  (spec tree unavailable: ",
                                              conditionMessage(e), ")"))
  cat(paste(tree, collapse = "\n"), "\n", sep = "")
  invisible(x)
}

#' @rdname print.Rceattle_data
#' @param object An `"Rceattle_data"` object.
#' @export
summary.Rceattle_data <- function(object, config = FALSE, ...) {
  # Forward `config` explicitly: it is matched by this method's own formal, so
  # it never reaches `...` and print() would otherwise fall back to its default.
  print(object, config = config, ...)
  invisible(object)
}
