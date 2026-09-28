#' Report which data inputs a model configuration requires
#'
#' Introspects the Rceattle data-element catalogue and reports, for a given model
#' configuration, which top-level `data_list` elements are **Required**,
#' **Optional** (used if supplied, otherwise default-filled by [clean_data()]), or
#' **Ignored** (not consulted because the feature that would use them is switched
#' off). It answers "what do I actually need to supply for *this* model?" without
#' having to read the validation code or the switch tables, the same conditions
#' enforced at fit time (they share one declarative table).
#'
#' The configuration can be given either as an existing (possibly partial)
#' `data_list`, its switches are normalized through [clean_data()] /
#' [switch_check()] so the conditions evaluate against filled defaults, or,
#' when no `data_list` is supplied, built from the convenience arguments.
#'
#' Requirements are *conditional*: e.g. `diet_data`, `ration_data` and the
#' bioenergetics scalars are Ignored under single-species (`msmMode = 0`) but
#' Required under multispecies (`msmMode > 0`); `caal_data` is Optional unless
#' `growth_model > 0`; `NByageFixed` is Ignored unless `estDynamics > 0`;
#' `emp_sel` is Required only for fleets with `Selectivity = "Fixed"`;
#' `index_cov` only for `Index_distribution = "MVN"`.
#'
#' @param data_list Optional. An existing Rceattle data list (e.g. from
#'   [read_data()] / [build_data()], or a bundled dataset). When supplied, the
#'   configuration is read from it and the convenience arguments are ignored.
#' @param msmMode Predation mode used to build the configuration when `data_list`
#'   is `NULL`. `0` single-species (default), `> 0` multispecies.
#' @param growth_model Growth mode (scalar or per-species). `> 0` estimates
#'   growth and requires `caal_data`.
#' @param estDynamics Numbers-at-age mode (scalar or per-species). `> 0` fixes
#'   numbers-at-age and requires `NByageFixed`.
#' @param Selectivity Optional character vector of per-fleet selectivity forms
#'   (used only when `data_list` is `NULL`) so the `emp_sel` requirement
#'   (`Selectivity = "Fixed"`) can be evaluated. Named for the `fleet_control`
#'   column it stands in for.
#' @param Index_distribution Optional character vector of per-fleet index
#'   likelihood families (used only when `data_list` is `NULL`) so the
#'   `index_cov` requirement (`Index_distribution = "MVN"`) can be evaluated.
#'   Named for the `fleet_control` column it stands in for.
#' @param Ceq Optional per-species consumption-equation codes (used only when
#'   `data_list` is `NULL`); `> 1` requires an environmental temperature index.
#'
#' @return A `data.frame` with one row per data element and columns:
#'   \describe{
#'     \item{`element`}{the `data_list` element name.}
#'     \item{`category`}{grouping (dimensions / biology / fishery / composition /
#'       predation / environment).}
#'     \item{`status`}{`"Required"`, `"Optional"`, or `"Ignored"`.}
#'     \item{`condition`}{the condition under which the element is required
#'       (`"always"` for the core backbone).}
#'     \item{`default`}{for Optional elements, the [clean_data()] default used
#'       when the element is absent.}
#'   }
#'   Rows are ordered Required, then Optional, then Ignored.
#'
#' @seealso [clean_data()], [build_data()].
#' @examples
#' # Single-species: predation/diet inputs are Ignored, comp_data Optional.
#' data_requirements(msmMode = 0)
#'
#' # Multispecies: diet_data and bioenergetics become Required.
#' data_requirements(msmMode = 1)
#'
#' # From a bundled dataset.
#' data_requirements(BS2017SS)
#' @export
data_requirements <- function(data_list = NULL, msmMode = 0, growth_model = 0,
                              estDynamics = 0, Selectivity = NULL,
                              Index_distribution = NULL, Ceq = NULL) {

  `%||%` <- function(a, b) if (is.null(a)) b else a

  # ---- Normalise the configuration into a switch-populated data list --------
  if (!is.null(data_list)) {
    dl <- tryCatch(
      suppressWarnings(suppressMessages(
        Rceattle::switch_check(Rceattle::clean_data(data_list)))),
      error = function(e) data_list
    )
  } else {
    fc <- NULL
    if (!is.null(Selectivity) || !is.null(Index_distribution)) {
      n <- max(length(Selectivity), length(Index_distribution), 1L)
      fc <- data.frame(
        Fleet_name    = paste0("Fleet_", seq_len(n)),
        Fleet_code    = seq_len(n),
        Selectivity   = if (is.null(Selectivity)) NA_character_
                        else rep_len(as.character(Selectivity), n),
        Index_distribution = if (is.null(Index_distribution)) .rce_column_schema()[["Index_distribution"]]$default
                        else rep_len(as.character(Index_distribution), n),
        stringsAsFactors = FALSE
      )
    }
    dl <- list(
      msmMode = msmMode, growth_model = growth_model,
      estDynamics = estDynamics, Ceq = Ceq,
      fleet_control = fc, nspp = 1L
    )
  }

  # Let an attached model_config() supply the switches, as fit_mod() does, so
  # this report and print(dl) describe the same model.
  dl <- .rce_classify_view(dl)

  # Explicitly-supplied switch arguments override whatever the data list stored,
  # matching fit_mod()'s "arg overrides the data slot" precedence -- so a user
  # can preview a data object under a different mode, e.g.
  # data_requirements(BS2017SS, msmMode = 0). (A bundled data list may store a
  # different msmMode than the fit actually uses.)
  if (!missing(msmMode))     dl$msmMode     <- msmMode
  if (!missing(growth_model)) dl$growth_model <- growth_model
  if (!missing(estDynamics)) dl$estDynamics <- estDynamics
  if (!missing(Ceq))         dl$Ceq         <- Ceq

  # ---- Classify every catalogued element ------------------------------------
  tbl <- .rce_requirement_table()
  parts <- lapply(tbl, function(r) {
    status <- .rce_classify(r, dl)
    data.frame(
      element   = r$element,
      category  = r$category %||% NA_character_,
      status    = status,
      condition = if (isTRUE(r$always_required)) "always"
                  else r$condition_label %||% "",
      default   = if (status == "Optional") (r$default_label %||% "") else "",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL

  # ---- Order: Required, then Optional, then Ignored; then category, element -
  status_rank <- match(out$status, c("Required", "Optional", "Ignored"))
  out <- out[order(status_rank, out$category, out$element), , drop = FALSE]
  rownames(out) <- NULL
  out
}
