#' Parse a linkage formula into fixed and random parts
#'
#' @description
#' A linkage formula describes how a process parameter depends on covariates
#' and on time. It has two kinds of term:
#'
#' \describe{
#'   \item{Fixed}{Anything outside a `( | )` bar: `~ temp`, `~ Year` (a linear
#'     trend), `~ factor(Year)` (annual deviates as fixed effects),
#'     `~ cut(Year, breaks = ...)` (time blocks). The RHS is passed to
#'     [stats::model.matrix()], so any function that works there works here.
#'     These have no standard deviation parameter; attach a prior to shrink
#'     them.}
#'   \item{Random}{Anything inside a bar: `(1 | Year)` for independent annual
#'     deviations, or a covariance structure such as `ar1(Year + 0 | fleet)`.
#'     These get a standard deviation (and correlation, where the structure has
#'     one) and are integrated out by the Laplace approximation, unless the
#'     spec sets `integrate = FALSE`, which estimates the deviations as a
#'     penalized fixed effect instead and requires that SD to be fixed.}
#' }
#'
#' Parsing is done by `reformulas`, as in `lme4` and `glmmTMB`, so the syntax
#' behaves as it does there.
#'
#' @param formula a one-sided formula.
#'
#' @return A list with `fixed` (a one-sided formula for `model.matrix()`),
#'   `re_terms` (list of bar expressions), and `re_structures` (character
#'   vector of covariance-structure names, one per term).
#'
#' @keywords internal
#' @noRd
.parse_linkage_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    stop("linkage formula must be one-sided, e.g. ~ temp + (1 | Year)",
         call. = FALSE)
  }

  .check_linkage_specials(formula)

  sf <- reformulas::splitForm(formula, specials = LINKAGE_STRUCTURES)
  list(fixed         = sf$fixedFormula,
       re_terms      = sf$reTrmFormulas,
       re_structures = as.character(sf$reTrmClasses))
}


#' Covariance structures a linkage random effect may use
#'
#' Names follow `glmmTMB` so they read the same way to anyone who has used it.
#' `"us"` is what a plain `(1 | g)` bar resolves to.
#'
#' @keywords internal
#' @noRd
LINKAGE_STRUCTURES <- c(
  "us",      # unstructured (what a plain bar gives)
  "diag",    # independent, structure-free
  "ar1",     # first-order autoregressive
  "rw",      # random walk
  "ou",      # Ornstein-Uhlenbeck
  "exp",     # exponential
  "cs",      # compound symmetry
  "toep"     # Toeplitz
)


#' Reject unrecognized covariance-structure wrappers
#'
#' @description
#' `reformulas::splitForm()` silently reads an unknown wrapper as `"us"`
#' (`bogus(Year + 0 | fleet)` becomes an unstructured term), so any wrapper
#' around a bar that is not in `LINKAGE_STRUCTURES` is an error here.
#'
#' @param formula a one-sided formula.
#' @return `formula` invisibly, on success.
#' @keywords internal
#' @noRd
.check_linkage_specials <- function(formula) {
  # Walk the parse tree for calls of the form f(<something> | <something>).
  wrappers <- character(0)
  walk <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    fn <- as.character(e[[1L]])[1L]
    has_bar <- length(e) >= 2L && any(vapply(
      as.list(e)[-1L],
      function(a) is.call(a) && identical(as.character(a[[1L]])[1L], "|"),
      logical(1)))
    if (has_bar && !fn %in% c("(", "|", "+", "~")) wrappers <<- c(wrappers, fn)
    for (a in as.list(e)[-1L]) walk(a)
    invisible(NULL)
  }
  walk(formula[[2L]])

  unknown <- setdiff(unique(wrappers), LINKAGE_STRUCTURES)
  if (length(unknown) > 0) {
    stop(sprintf(
      paste0("unknown covariance structure(s) in linkage formula: %s.\n",
             "  Available: %s.\n",
             "  A plain `(1 | group)` bar uses \"us\"."),
      paste0(unknown, "()", collapse = ", "),
      paste(LINKAGE_STRUCTURES, collapse = ", ")),
      call. = FALSE)
  }
  invisible(formula)
}


#' env_data columns that enter a model as fixed-effect linkage terms
#'
#' @description
#' `materialize_linkage()` refuses an NA in a fixed-effect linkage covariate,
#' because `model.matrix()` silently drops that row and misaligns the whole
#' design. Anything that wants to withhold covariate values -- a retrospective
#' peel -- therefore has to know which columns it must leave alone.
#'
#' Covariates that appear only inside a bar term (`(1 | Year)`, `ar1(1 | Year)`)
#' are not fixed terms and are not returned: a grouping variable is a level
#' label, and a state-space `observe` covariate is masked in the observation
#' rather than read as a design column.
#'
#' @param data_list a data list carrying the per-process linkage specs
#'   (`srr_linkages`, `M1_linkages`, `growth_linkages`, `q_linkages`,
#'   `sel_linkages`, `comp_linkages`), as `fit_mod()` stores them.
#'
#' @return A character vector of column names, possibly empty.
#'
#' @keywords internal
#' @noRd
.rce_linkage_fixed_covariates <- function(data_list) {
  groups <- c("srr_linkages", "M1_linkages", "growth_linkages",
              "q_linkages", "sel_linkages", "comp_linkages")
  specs <- list()
  collect <- function(x) {
    if (inherits(x, "Rceattle_linkage_spec")) { specs[[length(specs) + 1L]] <<- x }
    else if (is.list(x)) for (e in x) collect(e)
    invisible(NULL)
  }
  for (g in groups) collect(data_list[[g]])
  if (length(specs) == 0L) return(character(0))
  out <- unlist(lapply(specs, function(s) {
    p <- tryCatch(.parse_linkage_formula(s$formula), error = function(e) NULL)
    if (is.null(p)) return(character(0))
    all.vars(p$fixed)
  }), use.names = FALSE)
  unique(out)
}
