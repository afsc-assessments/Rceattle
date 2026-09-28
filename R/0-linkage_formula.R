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
#' The split is the same one `lme4` and `glmmTMB` use, and the parsing is done
#' by `reformulas`, the package both of those rely on, so the syntax behaves
#' the way a user of either would expect.
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
#' `reformulas::splitForm()` silently treats an unrecognized wrapper as the
#' default structure: `bogus(Year + 0 | fleet)` comes back as `"us"` with no
#' warning. A user would then get an unstructured term while believing they
#' had asked for something else, and nothing downstream could tell the
#' difference.
#'
#' So check the wrappers ourselves before handing the formula over. Anything
#' that wraps a bar and is not in `LINKAGE_STRUCTURES` is an error naming both
#' the offending term and the available structures.
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
