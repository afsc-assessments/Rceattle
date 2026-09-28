#' Likelihood profile across one or more parameter cells
#'
#' @description Re-fits an Rceattle model while holding selected cells of a
#'   parameter fixed at user-specified values. Supports profiling a single
#'   cell (e.g. \code{R_log_sd[species = 1]}) and arbitrary N-dimensional
#'   cross-profiles over multiple cells, e.g. \code{log_M1[1, 1, 1]} and
#'   \code{log_M1[1, 2, 1]} jointly, to profile residual M for males against
#'   females. For each grid point the targeted cells are fixed in the TMB
#'   map and the remaining parameters are re-estimated; the result is a
#'   grid of Rceattle models for downstream NLL surfaces.
#'
#' @param fitted an Rceattle model fit using \code{\link{fit_mod}}
#' @param param Name of the parameter to profile. Two ways to specify it:
#'   \describe{
#'     \item{Raw parameter slot}{any name in
#'       \code{Rceattle$estimated_params}; tested for \code{"R_log_sd"},
#'       \code{"rec_pars"}, and \code{"log_M1"}. \code{slots} must index
#'       into the full array and \code{transform} controls the scale.}
#'     \item{Natural-scale alias}{convenience shortcut for the three
#'       documented parameters. Aliases imply \code{transform = "log"}
#'       (values are taken in natural units and log'd before being
#'       substituted) and, for \code{rec_pars}, fill in the column from
#'       the alias name so \code{slots} only needs the species index:
#'       \itemize{
#'         \item \code{"sigmaR"}, \code{"R_sd"} -> \code{R_log_sd}
#'         \item \code{"M1"} -> \code{log_M1}
#'         \item \code{"R0"} -> \code{rec_pars[, 1]}
#'         \item \code{"alpha"} -> \code{rec_pars[, 2]}
#'         \item \code{"beta"} -> \code{rec_pars[, 3]}
#'         \item \code{"q"} -> \code{index_log_q}, base catchability, one entry
#'           per fleet. \code{slots} takes a fleet index or a \strong{fleet
#'           name}. Fleets sharing a \code{Catchability_index} share one q, so
#'           naming any of them profiles the whole group: the rest are added
#'           with a message and \code{joint} becomes \code{"value"}.
#'       }
#'       If \code{transform} is supplied with an alias it is ignored
#'       (with a warning).}
#'   }
#' @param slots A list whose entries are integer index vectors, one entry
#'   per cell to fix. Each entry's length must equal the number of
#'   dimensions of the resolved parameter, 1 for vectors
#'   (\code{R_log_sd}), 2 for matrices (\code{rec_pars}), 3 for 3-D arrays
#'   (\code{log_M1}). When using the \code{"R0"}/\code{"alpha"}/\code{"beta"}
#'   aliases, supply only the species index (length 1); the column is
#'   filled in from the alias. E.g. \code{list(c(1, 2, 1))} fixes
#'   \code{log_M1[1, 2, 1]}; \code{list(c(1, 1, 1), c(1, 2, 1))} fixes both
#'   sex cells for a males-vs-females cross-profile of species 1;
#'   \code{list(1, 2)} with \code{param = "sigmaR"} cross-profiles species
#'   1 and 2. If omitted, defaults to a single species-1 slot shaped to
#'   match the resolved parameter (e.g. \code{list(1)} for
#'   \code{R_log_sd}, \code{list(c(1, 1, 1))} for \code{log_M1},
#'   \code{list(1)} for the \code{rec_pars} aliases) and emits a warning;
#'   pass \code{slots} explicitly to silence the warning. Defaulting
#'   requires \code{length(values) == 1L} (otherwise the user must
#'   explicitly say which cell each grid targets).
#' @param values A list of numeric vectors, one per entry of \code{slots}.
#'   The full grid of fits is \code{expand.grid(values)}, so a single slot
#'   gives a 1-D profile and \emph{k} slots give a \emph{k}-D cross-profile.
#' @param transform How to map user values onto the internal parameter scale
#'   before substituting them into \code{inits}. Either \code{"log"}
#'   (default), \code{"identity"}, or a unary function (e.g.
#'   \code{qlogis}). Applied element-wise to every grid value. Aliases
#'   override this with \code{"log"}.
#' @param joint How the grid moves the cells in \code{slots}. Use a joint mode
#'   to profile a parameter supplied as a vector, such as an empirical
#'   age-based M: under \code{"none"} ten ages over 13 values is \code{13^10}
#'   fits, not 13.
#'   \describe{
#'     \item{\code{"none"}}{(default) each slot gets its own grid and they
#'       cross, giving a \emph{k}-dimensional profile for \emph{k} slots.}
#'     \item{\code{"multiply"}}{one grid, multiplying every cell's current
#'       value. The schedule keeps its shape and moves up or down as a whole.
#'       \code{1} is the fitted model.}
#'     \item{\code{"add"}}{one grid, added to every cell. The schedule keeps
#'       the differences between cells instead. \code{0} is the fitted model.}
#'     \item{\code{"value"}}{one grid; every cell takes that value.}
#'   }
#'   A joint mode needs \code{values} to be one vector.
#'   \code{"multiply"} and \code{"add"} work in natural units, so
#'   \code{transform} must be \code{"log"} or \code{"identity"}. A grid value
#'   that would take a log-scale parameter to zero or below is an error naming
#'   the cell, not a failed fit.
#' @param cores Number of cores to use for parallel fits. Default
#'   \code{NULL} picks \code{parallel::detectCores() - 6}, capped at 2 when
#'   running under \code{R CMD check} (which sets
#'   \code{_R_CHECK_LIMIT_CORES_}). Set to 1 to force sequential execution.
#' @param getsd whether each grid fit runs \code{TMB::sdreport}. The profile
#'   reads only the objective (\code{nll}), so \code{FALSE} is faster with no
#'   effect on the profile. Default \code{NULL} inherits the input model's
#'   setting (\code{TRUE} only if it holds an \code{sdrep}).
#' @param ... Unused; present for consistency with the \code{stats::profile}
#'   generic.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{Rceattle_list}{list of fitted Rceattle models, one per grid
#'       row; entries for non-converged fits are \code{NULL} so positions
#'       stay aligned with \code{grid}.}
#'     \item{grid}{data frame of grid values on the user scale (before
#'       \code{transform}); one column per profiled cell, named
#'       \code{slot_1}, \code{slot_2}, ...}
#'     \item{nll}{numeric vector of joint negative log-likelihoods
#'       (\code{opt$objective}); \code{NA} where the fit did not
#'       converge.}
#'     \item{param}{the profiled parameter name (echoed).}
#'     \item{slots}{the slots list (echoed for downstream plotting).}
#'     \item{alias}{the name you asked for, when it was one of the
#'       natural-scale aliases. Profiling \code{"M1"} returns \code{param =
#'       "log_M1"}, because the model estimates M on the log scale, while
#'       \code{grid} holds the M values you supplied. \code{alias} keeps
#'       \code{"M1"}, so a figure can label the axis in the units profiled
#'       rather than in log units. \code{NA} if you named the parameter slot
#'       directly.}
#'     \item{joint}{the \code{joint} mode used, so a figure and \code{print()}
#'       can say that \code{grid} holds a multiplier or an offset rather than a
#'       parameter value.}
#'   }
#'
#'   Has class \code{"Rceattle_profile"}, so printing it reports whether the
#'   grid brackets the minimum; see \code{\link{print.Rceattle_profile}}. Every
#'   element indexes exactly as before.
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(data_list = BS2017SS,
#'     inits = NULL, file = NULL,
#'     estimateMode = 0, random_rec = FALSE,
#'     msmMode = 0, avgnMode = 0,
#'     phase = FALSE, verbose = 0)
#'
#' # 1-D profile of sigmaR for species 1 (alias form, natural scale)
#' p1 <- profile(ss_run,
#'     param  = "sigmaR",
#'     slots  = list(1),
#'     values = list(seq(0.1, 1.5, by = 0.1)))
#'
#' # Equivalent raw form (log scale, user does the transform)
#' p1_raw <- profile(ss_run,
#'     param     = "R_log_sd",
#'     slots     = list(1),
#'     values    = list(log(seq(0.1, 1.5, by = 0.1))),
#'     transform = "identity")
#'
#' # 2-D cross-profile of M1 across species 1 and 2 (sex 1, age 1).
#' # BS2017SS is single-sex; with a multi-sex model the same form
#' # (e.g. c(1, 1, 1), c(1, 2, 1)) would cross-profile males vs females.
#' p2 <- profile(ss_run,
#'     param  = "M1",
#'     slots  = list(c(1, 1, 1), c(2, 1, 1)),
#'     values = list(seq(0.1, 0.4, length.out = 3),
#'                   seq(0.1, 0.4, length.out = 3)))
#'
#' # 1-D profile of SRR alpha for species 1 (alias drops the rec_pars column)
#' p3 <- profile(ss_run,
#'     param  = "alpha",
#'     slots  = list(1),
#'     values = list(seq(2, 80, length.out = 20)))
#'
#' # An empirical age-based M, moved as a whole. One slot per age bin of
#' # species 1, all on ONE grid of multipliers, so the shape of the schedule is
#' # preserved and 1 is the fitted model. Ages run minage .. minage + nages - 1.
#' sp <- 1
#' p4 <- profile(ss_run,
#'     param  = "M1",
#'     slots  = lapply(seq_len(ss_run$data_list$nages[sp]),
#'                     function(a) c(sp, 1, a)),
#'     values = list(seq(0.6, 1.4, by = 0.05)),
#'     joint  = "multiply")
#'
#' # The same schedule shifted by a constant instead of scaled; 0 is the fit.
#' p5 <- profile(ss_run,
#'     param  = "M1",
#'     slots  = lapply(seq_len(ss_run$data_list$nages[sp]),
#'                     function(a) c(sp, 1, a)),
#'     values = list(seq(-0.1, 0.1, by = 0.02)),
#'     joint  = "add")
#'
#' # Base catchability for a fleet, named rather than counted. If it shares a
#' # Catchability_index the whole group moves with it.
#' p6 <- profile(ss_run,
#'     param  = "q",
#'     slots  = list("BT_Pollock"),
#'     values = list(seq(0.5, 2.0, by = 0.1)))
#' }
#' @importFrom stats profile
#' @method profile Rceattle
#' @export
profile.Rceattle <- function(fitted = NULL,
                          param = NULL,
                          slots = NULL,
                          values = NULL,
                          transform = "log",
                          cores = NULL,
                          getsd = NULL,
                          # Appended rather than placed beside `transform`, so a
                          # positional call keeps meaning what it did.
                          joint = c("none", "multiply", "add", "value"),
                          ...) {

  joint <- match.arg(joint)

  # -- Input validation ----
  if (!inherits(fitted, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }
  # Grid fits inherit the input model's sdreport setting unless overridden;
  # the profile reads only the objective, not sdrep.
  if (is.null(getsd)) getsd <- !is.null(fitted$sdrep)
  if (is.null(param) || !is.character(param) || length(param) != 1L) {
    stop("`param` must be a single character string naming a parameter slot.")
  }
  if (!is.list(values) || length(values) == 0L) {
    stop("`values` must be a non-empty list of numeric grids.")
  }

  # Natural-scale aliases: each maps to a real parameter, implies log()
  # transform, and (for rec_pars aliases) fills in the column index so
  # `slots` only needs the species index.
  alias_table <- list(
    sigmaR = list(param = "R_log_sd", rec_pars_col = NA_integer_),
    R_sd   = list(param = "R_log_sd", rec_pars_col = NA_integer_),
    M1     = list(param = "log_M1",   rec_pars_col = NA_integer_),
    R0     = list(param = "rec_pars", rec_pars_col = 1L),
    alpha  = list(param = "rec_pars", rec_pars_col = 2L),
    beta   = list(param = "rec_pars", rec_pars_col = 3L),
    q      = list(param = "index_log_q", rec_pars_col = NA_integer_)
  )

  alias_name   <- NA_character_
  rec_pars_col <- NA_integer_
  if (param %in% names(alias_table)) {
    alias_name <- param
    a <- alias_table[[alias_name]]

    # Aliases force log transform; warn if user passed something else
    if (!identical(transform, "log")) {
      warning(sprintf(
        "`param = \"%s\"` is a natural-scale alias for `%s`; ignoring the supplied `transform` (aliases imply transform = \"log\").",
        alias_name, a$param
      ))
    }
    transform    <- "log"
    rec_pars_col <- a$rec_pars_col
    param        <- a$param   # resolve to real parameter slot
  }

  if (!param %in% names(fitted$estimated_params)) {
    stop("`param` '", param, "' not found in Rceattle$estimated_params.")
  }

  par_array <- fitted$estimated_params[[param]]
  par_ndim  <- if (is.null(dim(par_array))) 1L else length(dim(par_array))

  # Default `slots` to species 1 (a single profile point shaped to match
  # the resolved parameter). For rec_pars aliases the user slot is just
  # the species index; otherwise it's a 1 for every dimension.
  if (is.null(slots)) {
    user_slot_dim <- par_ndim - if (!is.na(rec_pars_col)) 1L else 0L
    default_user_slot <- rep(1L, user_slot_dim)

    if (length(values) != 1L) {
      stop(sprintf(
        "`slots` not supplied but `values` has %d grids -- the species-1 default only covers one cell. Pass `slots` explicitly to profile multiple cells.",
        length(values)
      ))
    }

    pretty_slot <- if (length(default_user_slot) == 1L) {
      as.character(default_user_slot)
    } else {
      paste0("c(", paste(default_user_slot, collapse = ", "), ")")
    }
    warning(sprintf(
      "`slots` not supplied; defaulting to species 1 (slots = list(%s)). Pass `slots` explicitly to silence this warning.",
      pretty_slot
    ))

    slots <- list(default_user_slot)
  }

  if (!is.list(slots) || length(slots) == 0L) {
    stop("`slots` must be a non-empty list of integer index vectors.")
  }
  # Under a joint mode every slot reads the single grid, so one vector covers
  # any number of cells; only the crossing mode needs one grid per slot.
  if (identical(joint, "none") && length(values) != length(slots)) {
    stop("`values` must be a list with the same length as `slots`.")
  }

  # Append rec_pars column for rec_pars aliases
  if (!is.na(rec_pars_col)) {
    for (k in seq_along(slots)) {
      if (length(slots[[k]]) != 1L) {
        stop(sprintf(
          "Under alias `\"%s\"`, slots[[%d]] should be a single species index (got length %d). The rec_pars column is filled in from the alias name.",
          alias_name, k, length(slots[[k]])
        ))
      }
      slots[[k]] <- c(as.integer(slots[[k]]), rec_pars_col)
    }
  }

  # -- Catchability: name a fleet, and move its whole sharing group ----
  # index_log_q is one entry per fleet, named by Fleet_name. Fleets sharing a
  # Catchability_index share ONE q, so fixing one member alone would leave the
  # rest estimating a common q -- not that fleet's q profiled.
  if (identical(alias_name, "q")) {
    fc <- fitted$data_list$fleet_control
    resolve_fleet <- function(s) {
      if (is.character(s)) {
        hit <- match(s, as.character(fc$Fleet_name))
        if (anyNA(hit)) {
          stop("No fleet named ", paste0("\"", s[is.na(hit)], "\"", collapse = ", "),
               ". Fleets are: ",
               paste(fc$Fleet_name, collapse = ", "), ".", call. = FALSE)
        }
        return(as.integer(hit))
      }
      as.integer(s)
    }
    slots <- lapply(slots, resolve_fleet)

    # Analytical forms solve q from each fleet's own index observations and
    # overwrite index_q, so index_log_q never reaches the likelihood: the grid
    # would move nothing and the profile would come back flat, reading as "the
    # data do not inform q". `.map_switch()` is the one reading of q_map;
    # est_index_q is not carried on a fitted object.
    qform <- suppressWarnings(suppressMessages(
      .map_switch(fc$Catchability, q_map, "Catchability")))
    analytical <- which(qform %in% c(3L, 7L))
    bad <- intersect(unlist(slots), analytical)
    if (length(bad)) {
      stop("Fleet(s) ", paste(fc$Fleet_name[bad], collapse = ", "),
           " use an analytical catchability, which is solved from their own ",
           "index and ignores `index_log_q`. Profiling it would return the ",
           "same fit at every grid point.", call. = FALSE)
    }

    qi <- fc$Catchability_index
    if (!is.null(qi)) {
      asked <- unlist(slots)
      group_of <- function(s) {
        # An NA Catchability_index means no sharing group; `qi == NA` is NA
        # throughout and would select nothing. An analytical fleet solves its
        # own q, so it is not part of the group's shared parameter.
        if (is.na(qi[s])) return(s)
        setdiff(which(!is.na(qi) & qi == qi[s]), analytical)
      }
      groups   <- lapply(asked, group_of)
      expanded <- unique(unlist(groups))
      shared   <- any(vapply(groups, length, integer(1)) > 1L)

      added <- setdiff(expanded, asked)
      if (length(added)) {
        message("Fleet(s) ", paste(fc$Fleet_name[added], collapse = ", "),
                " share a Catchability_index with the fleet(s) profiled, so ",
                "they carry the same q and are moved with it.")
      }
      # Keyed on whether a group is shared, not on the slot count growing: two
      # members of the SAME group named together add nothing, and crossing them
      # would put two values on one q parameter.
      if (shared) {
        slots <- lapply(expanded, function(x) x)
        if (identical(joint, "none")) joint <- "value"
      }
    }
  }

  par_dim <- if (is.null(dim(par_array))) length(par_array) else dim(par_array)

  for (k in seq_along(slots)) {
    if (length(slots[[k]]) != par_ndim) {
      stop(sprintf(
        "slots[[%d]] has length %d but '%s' has %d dimension(s).",
        k, length(slots[[k]]), param, par_ndim
      ))
    }
    if (!all(is.finite(slots[[k]])) || any(slots[[k]] < 1)) {
      stop(sprintf("slots[[%d]] must be a vector of positive integers.", k))
    }
    if (any(slots[[k]] > par_dim)) {
      stop(sprintf(
        "slots[[%d]] = c(%s) is out of bounds for '%s' (dim c(%s)).",
        k,
        paste(slots[[k]], collapse = ", "),
        param,
        paste(par_dim, collapse = ", ")
      ))
    }
  }

  # Build transform fn
  trans_fun <- if (is.function(transform)) {
    transform
  } else if (identical(transform, "log")) {
    log
  } else if (identical(transform, "identity")) {
    function(x) x
  } else {
    stop("`transform` must be \"log\", \"identity\", or a function.")
  }

  # -- Joint modes: one grid moving every named cell together ----
  # Slots otherwise cross, so a 10-age M schedule over 13 values would be 13^10
  # fits rather than 13.
  if (!identical(joint, "none")) {
    if (length(values) != 1L) {
      stop("`joint = \"", joint, "\"` moves every cell in `slots` on ONE grid, ",
           "so `values` must be a single vector; got ", length(values), ".",
           call. = FALSE)
    }
    if (is.function(transform) && joint %in% c("multiply", "add")) {
      # Scaling a stored log is a different model, not a different
      # parameterisation, so the transform has to be invertible.
      stop("`joint = \"", joint, "\"` needs `transform` to be \"log\" or ",
           "\"identity\" -- it acts on the natural scale, and an arbitrary ",
           "transform cannot be inverted.", call. = FALSE)
    }
    names(values) <- "slot_1"
  } else {
    names(values) <- paste0("slot_", seq_along(values))
  }

  # Natural scale <-> stored scale, for the joint modes.
  from_natural <- trans_fun
  to_natural <- if (identical(transform, "log")) exp else function(x) x

  # Build grid (user-scale values; transform applied at fit time)
  grid <- if (identical(joint, "none")) {
    expand.grid(values, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  } else {
    data.frame(slot_1 = values[[1]])
  }
  ngrid <- nrow(grid)

  # Cross-platform parallel via parallel::parLapply on a PSOCK cluster
  # (mirrors jitter()/retrospective()). Respect the CRAN core limit.
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  if (is.null(cores)) {
    cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  } else {
    cores <- max(1L, as.integer(cores))
    if (cran_cap) cores <- min(cores, 2L)
  }
  use_parallel <- ngrid > 1L && cores > 1L

  # Generic [<-]: assign `val` into `arr` at index vector `idx`
  assign_at <- function(arr, idx, val) {
    do.call("[<-", c(list(arr), as.list(idx), list(val)))
  }

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Per-grid-point closure ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  run_one_point <- function(i) {

    inits     <- fitted$estimated_params
    data_list <- fitted$data_list
    map_obj <- fitted$map

    # Substitute the fixed value at each profiled cell. Under "none" slot k
    # reads grid column k; a joint mode gives every slot the one column, and
    # multiply/add act on the cell's current natural-scale value.
    for (k in seq_along(slots)) {
      new_stored <- if (identical(joint, "none")) {
        trans_fun(grid[i, k])
      } else {
        gv  <- grid[i, 1]
        cur <- to_natural(do.call("[", c(list(inits[[param]]),
                                         as.list(slots[[k]]))))
        nat <- switch(joint,
                      multiply = cur * gv,
                      add      = cur + gv,
                      value    = gv)
        if (identical(transform, "log") && (!is.finite(nat) || nat <= 0)) {
          stop(sprintf(
            "`joint = \"%s\"` with grid value %s puts cell c(%s) of '%s' at %s, which is not positive; '%s' is stored on the log scale. Narrow `values`.",
            joint, format(gv), paste(slots[[k]], collapse = ", "), param,
            format(nat), param), call. = FALSE)
        }
        from_natural(nat)
      }
      inits[[param]] <- assign_at(inits[[param]], slots[[k]], new_stored)
    }

    # Force profiled cells to NA
    for (k in seq_along(slots)) {
      map_obj$mapList[[param]] <- assign_at(map_obj$mapList[[param]],
                                            slots[[k]],
                                            NA)
    }
    map_obj$mapFactor <- lapply(map_obj$mapList, factor)

    newmod <-
      suppressMessages(suppressWarnings(
        # Refit with the profiled parameter fixed at its grid value (mapped off
        # in map_obj). estimateMode falls back to 1 -- profile the hindcast fit,
        # not a projection.
        .refit_like(
          data_list        = data_list,
          inits            = inits,
          map              = map_obj,
          estimateMode     = ifelse(data_list$estimateMode < 3, 1, data_list$estimateMode),
          getsd            = getsd,
          srr_mse_switchyr = min(data_list$srr_mse_switchyr, data_list$endyr),
          suit_endyr       = pmin(data_list$suit_endyr, data_list$endyr))
      ))

    if (.refit_converged(newmod)) {
      return(newmod)
    }
    return(NULL)
  } # End run_one_point closure


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Dispatch (parallel via PSOCK or sequential) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if (use_parallel) {
    mod_list <- .parallel_lapply(seq_len(ngrid), run_one_point, min(cores, ngrid), environment())
  } else {
    mod_list <- lapply(seq_len(ngrid), run_one_point)
  }

  # NLL aligned with grid; NA for non-converged
  nll <- vapply(mod_list,
                function(x) if (is.null(x)) NA_real_ else x$opt$objective,
                numeric(1))

  names(mod_list) <- paste0("Fit_", seq_len(ngrid))

  # Still the same list -- $Rceattle_list, $grid, $nll, $param and $slots are
  # unchanged, so `prof$grid$slot_1` and `prof$nll - min(prof$nll)` index as
  # before. The class adds a print method that says whether the grid brackets
  # the minimum, which the bare vectors never did; see print.Rceattle_profile().
  structure(list(
    Rceattle_list = mod_list,
    grid          = grid,
    nll           = nll,
    param         = param,
    slots         = slots,
    alias         = alias_name,
    joint         = joint
  ), class = "Rceattle_profile")
}


#' Print method for a likelihood profile
#'
#' @description Reports whether the grid actually brackets the minimum. A
#' profile whose lowest point is its first or last grid value has not found the
#' optimum, it has run out of grid, and the curve drawn from it understates
#' how far the parameter can move. That is the failure the numbers alone hide,
#' since a partial profile plots as a perfectly ordinary line.
#'
#' @details
#' For a one-dimensional profile the interval reported is the usual
#' profile-likelihood one: the grid values within `cutoff` of the minimum, where
#' the default 1.92 is \eqn{\chi^2_1(0.95)/2}. It is read off the grid, so it is
#' no finer than the spacing of `values`, and it is reported as open on either
#' side the grid does not close. It is also referenced to the best GRID point
#' rather than to the unconstrained MLE, which the object does not hold: the
#' grid minimum sits at or above the MLE, so the interval errs wide. And it is
#' reported as a range, so a profile with a second basin, or a failed point
#' inside the range, is called out as not contiguous rather than left to read
#' as one interval. No interval is given for a cross-profile over
#' two or more cells, the cutoff would be \eqn{\chi^2_k(0.95)/2} and the region
#' is not an interval.
#'
#' Under `random_rec = TRUE` the objective is the Laplace-approximated marginal
#' likelihood, so this is a profile of that, with the usual caveat that the
#' approximation is what is being profiled.
#'
#' @param x A `"Rceattle_profile"` object from [profile.Rceattle()].
#' @param cutoff Objective units above the minimum that bound the reported
#'   interval. Default `1.92`, the 95% profile-likelihood cutoff for one
#'   parameter.
#' @param ... Currently unused.
#' @return `x`, invisibly.
#' @export
print.Rceattle_profile <- function(x, cutoff = 1.92, ...) {
  nll  <- x$nll
  grid <- x$grid
  ok   <- is.finite(nll)
  ngrid <- length(nll)

  if (!any(ok)) {
    .rce_diag_header("profile", "FAIL",
                     paste0("no grid point converged for '", x$param, "'"))
    return(invisible(x))
  }

  best <- which.min(replace(nll, !ok, Inf))
  # An edge minimum on ANY profiled cell: the grid stopped before the surface
  # turned, so the minimum reported is the end of the grid, not the optimum.
  at_edge <- vapply(seq_along(grid), function(j) {
    v <- grid[[j]]
    v[best] == min(v) || v[best] == max(v)
  }, logical(1))

  sev <- if (any(at_edge)) "WARN"
         else if (any(!ok)) "NOTE"
         else "OK"

  # `format(trim = TRUE)`, not formatC("g"): that pads to the width implied by
  # `digits` and would print a grid value as "  0.7".
  num <- function(v) format(v, digits = 4, trim = TRUE)

  # Name the parameter in the units `grid` is in. Under an alias `param` has
  # been resolved to the internal slot, so reporting it would put a minimum of
  # 0.4 next to "log_M1" when 0.4 is an M, not a log M.
  shown <- if (!is.null(x$alias) && !is.na(x$alias)) x$alias else x$param

  # Under a joint mode the grid is a multiplier or an offset applied to every
  # cell, so reporting a minimum "at 1.15" against the bare parameter name would
  # read as the parameter itself being 1.15.
  if (!is.null(x$joint) && !identical(x$joint, "none")) {
    shown <- paste0(shown,
                    switch(x$joint, multiply = " multiplier",
                           add = " offset", value = ""),
                    " on ", length(x$slots), " cell(s)")
  }

  .rce_diag_header(
    "profile", sev,
    paste0(shown, " over ", ngrid, " grid point(s); ",
           if (all(ok)) "all converged"
           else paste0(sum(!ok), " did not converge")))

  at <- vapply(grid, function(v) v[best], numeric(1))
  cat("  minimum -log L : ", formatC(nll[best], format = "f", digits = 4),
      "  at ", paste(names(grid), num(at), sep = " = ", collapse = ", "),
      "\n", sep = "")

  if (length(grid) == 1L) {
    v <- grid[[1]]
    within <- ok & (nll - nll[best] <= cutoff)
    lo <- min(v[within]); hi <- max(v[within])
    cat("  within ", cutoff, " of the minimum : ",
        if (lo == min(v)) "<=" else "", num(lo),
        " to ", if (hi == max(v)) ">=" else "", num(hi),
        "\n", sep = "")
    # min-to-max is an interval only if the grid points between them are all
    # inside it. A profile with a second basin -- or one non-converged point in
    # the middle -- would otherwise be reported as one wide interval that
    # includes values the profile actually rules out.
    span <- which(v >= lo & v <= hi)
    if (!all(within[span])) {
      cat("  not contiguous: ", sum(!within[span]),
          " grid point(s) inside that range are above the cutoff or did not ",
          "converge -- read $nll rather than the range\n", sep = "")
    }
  }

  if (any(at_edge)) {
    cat("  the minimum is at a grid edge (",
        paste(names(grid)[at_edge], collapse = ", "),
        "), so the grid does not bracket it -- widen `values`\n", sep = "")
  }
  invisible(x)
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Likelihood components across a profile ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Which axis each row of `jnll_comp` is accumulated over, i.e. what its columns
# count. The data, selectivity and catchability rows are filled per fleet; the
# priors, penalties and predation rows per species; the linkage random-effect
# row is model-wide and always lands in column 1. Reference-point penalties are
# species rows that also use column 1 for a model-wide term, so a non-zero
# "Reference point penalties" on species 1 may be either.
#
# This is a hand-synced registry: it must stay in lockstep with the `JnllRow`
# enum in src/TMB/ceattle.cpp and the row labels in R/6-rename_output.R. A row
# added there and not here is labelled by component alone, with a warning.
.JNLL_ROW_AXIS <- c(
  "Index data"                 = "fleet",
  "Catch data"                 = "fleet",
  "Composition data"           = "fleet",
  "CAAL data"                  = "fleet",
  "Non-parametric selectivity" = "fleet",
  "Selectivity deviates"       = "fleet",
  "Catchability prior"         = "fleet",
  "Catchability deviates"      = "fleet",
  "Stock-recruit prior"        = "species",
  "Initial abundance deviates" = "species",
  "Recruitment deviates"       = "species",
  "Stock-recruit penalty"      = "species",
  "Reference point penalties"  = "species",
  "Zero n-at-age penalty"      = "species",
  "M prior"                    = "species",
  "M random effects"           = "species",
  "Ration"                     = "species",
  "Ration penalties"           = "species",
  "Stomach content data"       = "species",
  "Linkage-table priors"       = "species",
  "Linkage random effects"     = "model"
)


#' Likelihood components along a profile
#'
#' @description Pulls the per-fleet and per-species negative log-likelihood
#'   components out of every fit in a [profile()] and returns them as one long
#'   data frame, one row per grid point per component.
#'
#'   The total is the least informative curve on a profile: a well-behaved
#'   quadratic total can hide two data sources pulling the parameter in
#'   opposite directions, and their conflict is visible only once the
#'   components are drawn separately. This is the extractor behind
#'   [plot_profile()]; use it directly to tabulate, or to draw the figure
#'   yourself.
#'
#' @details
#' **Which cells are reported.** `jnll_comp` is a component-by-column matrix
#' whose columns mean different things on different rows, fleets on the data,
#' selectivity and catchability rows, species on the priors, penalties and
#' predation rows. Each cell is labelled from the axis its row uses, so a cell
#' becomes e.g. `"Shelikof acoustic: Index data"`. The unit is dropped from the
#' label when the model has only one fleet (or one species) to distinguish.
#' Cells that are zero at every grid point are dropped: they are components the
#' model does not fit.
#'
#' **What `Total` is.** The `"Total"` series is `object$nll`, the objective each
#' grid fit actually minimized. Under `random_rec = TRUE` that is the
#' Laplace-approximated marginal likelihood while the components are the inner
#' joint negative log-likelihood, so the components will not sum to it; the
#' function says so when they differ. Compare the shapes, not the sums.
#'
#' **Non-converged grid points** keep their row with a value of `NA`, so a
#' failed fit leaves a gap in the curve rather than a straight segment drawn
#' across it.
#'
#' @param object An `"Rceattle_profile"` object from [profile.Rceattle()].
#' @param weighted Report the weighted components the optimizer minimized
#'   (`TRUE`, the default, `quantities$jnll_comp`) or the unweighted ones
#'   (`FALSE`, `quantities$unweighted_jnll_comp`). Conflict is normally read off
#'   the weighted components, since those are what moved the fit.
#'   `unweighted_jnll_comp` exists so Francis and McAllister-Ianelli can read a
#'   composition likelihood without its `Comp_weights` multiplier, so only the
#'   rows that hold such a multiplier are filled: composition, CAAL, stomach
#'   content and the two linkage rows. Every other row is zero there and is
#'   dropped as unfitted, so `weighted = FALSE` returns a much smaller set of
#'   series, the index, catch, selectivity, catchability and penalty
#'   components are absent, not flat.
#' @param relative How to place each series on the y axis.
#'   \describe{
#'     \item{`"own"`}{(default) subtract each series' own minimum, so every
#'       curve starts at zero and its minimum marks the value that component
#'       prefers. Objective units, so depth still says how strongly.}
#'     \item{`"scaled"`}{as `"own"`, then divide each series by its own change
#'       over the grid, so every curve runs 0 to 1. Use it when one component
#'       dwarfs the rest and flattens the others onto the axis. It **discards
#'       magnitude**: a component moving 0.02 draws like one moving 40, so
#'       raise `minfraction` with it. A series that does not move stays at
#'       zero.}
#'     \item{`"minimum"`}{subtract each series' value at the grid point where
#'       the total is lowest, so the curves show each component's change away
#'       from the fitted optimum.}
#'     \item{`"none"`}{the raw negative log-likelihoods.}
#'   }
#'   `minfraction` and the series ordering are always computed on the raw
#'   change over the grid, before any of this is applied, so `"scaled"` does
#'   not quietly disable the filter by making every span 1.
#' @param minfraction Drop components whose change over the grid is less than
#'   this fraction of the total's change, as in `r4ss::SSplotProfile()`. Default
#'   `0` keeps everything; `plot_profile()` uses `0.01`. `"Total"` is never
#'   dropped.
#' @param include_total Include the `"Total"` series. Default `TRUE`; `FALSE`
#'   also disables `minfraction`, which is defined against the total.
#'
#' @return A data frame with one row per grid point per retained component:
#'   the profile's `grid` columns (`slot_1`, ...) holding the value profiled
#'   over, then `fit` (grid row index), `component` (the `jnll_comp` row),
#'   `unit` (fleet or species name, `NA` for model-wide rows), `axis`
#'   (`"fleet"`, `"species"` or `"model"`), `series` (the plotting label), and
#'   `value` (the re-zeroed negative log-likelihood). Series are ordered by
#'   decreasing change over the grid, with `"Total"` first. The profile's
#'   `param`, `alias` and the `relative` used are held as attributes.
#'
#' @seealso [plot_profile()] to draw it, [profile.Rceattle()] to produce the
#'   profile.
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(data_list = BS2017SS, inits = NULL, file = NULL,
#'     estimateMode = 0, random_rec = FALSE, msmMode = 0, avgnMode = 0,
#'     phase = FALSE, verbose = 0)
#'
#' prof <- profile(ss_run, param = "M1", slots = list(c(1, 1, 1)),
#'                 values = list(seq(0.2, 0.5, by = 0.05)))
#'
#' comps <- profile_components(prof)
#' head(comps)
#' }
#' @export
profile_components <- function(object,
                               weighted = TRUE,
                               relative = c("own", "scaled", "minimum", "none"),
                               minfraction = 0,
                               include_total = TRUE) {

  # -- Input validation ----
  if (!inherits(object, "Rceattle_profile")) {
    stop("`object` must be an \"Rceattle_profile\" from profile().",
         call. = FALSE)
  }
  relative <- match.arg(relative)
  if (!is.logical(weighted) || length(weighted) != 1L || is.na(weighted)) {
    stop("`weighted` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(include_total) || length(include_total) != 1L ||
      is.na(include_total)) {
    stop("`include_total` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(minfraction) || length(minfraction) != 1L ||
      is.na(minfraction) || minfraction < 0) {
    stop("`minfraction` must be a single non-negative number.", call. = FALSE)
  }

  fits  <- object$Rceattle_list
  grid  <- object$grid
  ngrid <- nrow(grid)
  ok    <- !vapply(fits, is.null, logical(1))
  if (!any(ok)) {
    stop("No grid point in this profile converged, so it carries no ",
         "likelihood components.", call. = FALSE)
  }

  slot <- if (weighted) "jnll_comp" else "unweighted_jnll_comp"
  ref  <- fits[[which(ok)[1]]]
  ref_mat <- ref$quantities[[slot]]
  if (is.null(ref_mat) || !length(ref_mat)) {
    stop("The fits in this profile report no `quantities$", slot, "`.",
         call. = FALSE)
  }
  ref_mat <- as.matrix(ref_mat)

  # -- Label every cell from the axis its row is accumulated over ----
  components <- rownames(ref_mat)
  if (is.null(components)) {
    stop("`quantities$", slot, "` has no row labels; this profile was not ",
         "produced by a current version of fit_mod().", call. = FALSE)
  }
  axis <- unname(.JNLL_ROW_AXIS[components])
  if (anyNA(axis)) {
    warning("Likelihood component(s) ",
            paste0("\"", components[is.na(axis)], "\"", collapse = ", "),
            " are not in the row registry, so their fleet or species is not ",
            "named. Add them to `.JNLL_ROW_AXIS` in R/9-profile.R.",
            call. = FALSE)
    axis[is.na(axis)] <- "model"
  }

  dl <- ref$data_list
  fleet_names <- as.character(dl$fleet_control$Fleet_name)
  sp_names    <- as.character(dl$spnames)
  # A model may carry fewer names than jnll_comp has columns -- the matrix is
  # dimensioned max(n_fleets, nspp) -- so a column past the end is named
  # positionally rather than left blank.
  name_at <- function(nms, j, what) {
    if (j <= length(nms) && !is.na(nms[j]) && nzchar(nms[j])) nms[j]
    else paste(what, j)
  }

  # -- Gather the matrices, aligned with the grid ----
  vals <- array(NA_real_, dim = c(nrow(ref_mat), ncol(ref_mat), ngrid))
  for (i in seq_len(ngrid)) {
    if (!ok[i]) next
    m <- fits[[i]]$quantities[[slot]]
    if (is.null(m) || !identical(dim(as.matrix(m)), dim(ref_mat))) {
      stop("Grid point ", i, " reports a `", slot, "` of a different shape ",
           "than the others; the profile mixes models and cannot be pooled.",
           call. = FALSE)
    }
    vals[, , i] <- as.matrix(m)
  }

  # A cell that is zero at every grid point is a component this model does not
  # fit, not a component that happens to be flat.
  keep <- which(apply(vals, c(1, 2), function(v) any(is.finite(v) & v != 0)),
                arr.ind = TRUE)
  if (!nrow(keep) && !include_total) {
    stop("No likelihood component in this profile is non-zero.", call. = FALSE)
  }

  # The fleet or species name is dropped from the label when there is only one
  # of them: "Recruitment deviates" reads better than "Pollock: Recruitment
  # deviates" on a single-species model, and no ambiguity is created.
  prefix <- c(fleet = length(fleet_names) > 1L,
              species = length(sp_names) > 1L, model = FALSE)

  units  <- character(nrow(keep))
  axes   <- character(nrow(keep))
  labels <- character(nrow(keep))
  for (k in seq_len(nrow(keep))) {
    r <- keep[k, "row"]; cc <- keep[k, "col"]
    axes[k]  <- axis[r]
    units[k] <- switch(axes[k],
                       fleet   = name_at(fleet_names, cc, "Fleet"),
                       species = name_at(sp_names, cc, "Species"),
                       NA_character_)
    labels[k] <- if (!is.na(units[k]) && prefix[[axes[k]]]) {
      paste0(units[k], ": ", components[r])
    } else {
      components[r]
    }
  }

  # Only `Fleet_code` is required to be unique -- it must equal the row number
  # -- so two fleets may share a `Fleet_name`. Left alone, their two series
  # would merge into one line drawn from interleaved values. A repeated label
  # is disambiguated by the fleet or species number, which is unique by
  # construction.
  dup <- labels %in% labels[duplicated(labels)]
  for (k in which(dup)) {
    labels[k] <- if (is.na(units[k])) {
      # A model-wide row with more than one column filled: there is no fleet or
      # species to name it by, so the column is all there is to say.
      paste0(components[keep[k, "row"]], " (column ", keep[k, "col"], ")")
    } else {
      paste0(units[k], " (", axes[k], " ", keep[k, "col"], "): ",
             components[keep[k, "row"]])
    }
  }

  pieces <- lapply(seq_len(nrow(keep)), function(k) {
    data.frame(
      fit       = seq_len(ngrid),
      component = components[keep[k, "row"]],
      unit      = units[k],
      axis      = axes[k],
      series    = labels[k],
      value     = vals[keep[k, "row"], keep[k, "col"], ],
      stringsAsFactors = FALSE)
  })

  # -- The total is the objective, which is not always the components' sum ----
  if (include_total) {
    comp_sum <- apply(vals, 3, sum)
    drift <- suppressWarnings(max(abs(comp_sum[ok] - object$nll[ok])))
    if (is.finite(drift) && drift > 1e-4) {
      warning("The likelihood components sum to the joint negative ",
              "log-likelihood but `Total` is the objective each fit ",
              "minimized, which for a model with random effects is the ",
              "Laplace-approximated marginal; they differ by up to ",
              formatC(drift, format = "g", digits = 3),
              ". Compare the shapes of the curves, not their sums.",
              call. = FALSE)
    }
    pieces <- c(list(data.frame(
      fit = seq_len(ngrid), component = "Total", unit = NA_character_,
      axis = "model", series = "Total", value = object$nll,
      stringsAsFactors = FALSE)), pieces)
  }

  out <- do.call(rbind, pieces)
  out$value[!ok[out$fit]] <- NA_real_

  # How far each series moves, in objective units. Measured HERE because
  # `relative = "scaled"` divides it out: computed after, every span would be 1
  # and `minfraction` would filter nothing -- losing the guard exactly when
  # scaling makes a component moving 0.02 look like one moving 40.
  span <- vapply(split(out$value, out$series),
                 function(v) suppressWarnings(diff(range(v, na.rm = TRUE))),
                 numeric(1))
  span[!is.finite(span)] <- 0

  # -- Re-zero ----
  shift <- switch(
    relative,
    none = stats::setNames(rep(0, length(unique(out$series))),
                           unique(out$series)),
    own  = ,
    scaled = vapply(split(out$value, out$series),
                    function(v) suppressWarnings(min(v, na.rm = TRUE)),
                    numeric(1)),
    minimum = {
      # No usable objective anywhere: which.min() over the all-Inf replacement
      # would silently pick grid point 1 and re-zero everything there, which
      # reads as a fitted optimum that was never found.
      if (!any(is.finite(object$nll))) {
        warning("No grid point reports an objective, so `relative = ",
                "\"minimum\"` has no minimum to re-zero at; values are raw.",
                call. = FALSE)
        stats::setNames(rep(0, length(unique(out$series))), unique(out$series))
      } else {
        best <- which.min(replace(object$nll, !is.finite(object$nll), Inf))
        vapply(split(out, out$series),
               function(d) d$value[match(best, d$fit)], numeric(1))
      }
    })
  shift[!is.finite(shift)] <- 0
  out$value <- out$value - shift[out$series]

  # -- Put every series on its own scale, if asked ----
  # Buys the shape and the position of each minimum; costs magnitude. A series
  # that does not move keeps its zero rather than dividing by it.
  if (identical(relative, "scaled")) {
    denom <- span[out$series]
    out$value <- ifelse(denom > 0, out$value / denom, 0)
  }

  # -- Drop components that barely move, and order by how much they do ----
  if (minfraction > 0 && include_total) {
    if (span[["Total"]] <= 0) {
      warning("`minfraction` filters against the total's change over the ",
              "grid, which is ", formatC(span[["Total"]], format = "g",
                                         digits = 3),
              " here; keeping every component.", call. = FALSE)
    } else {
      drop <- setdiff(names(span)[span < minfraction * span[["Total"]]], "Total")
      out <- out[!out$series %in% drop, , drop = FALSE]
    }
  }

  ord <- names(sort(span[unique(out$series)], decreasing = TRUE))
  if (include_total) ord <- c("Total", setdiff(ord, "Total"))
  out$series <- factor(out$series, levels = ord)

  out <- cbind(grid[out$fit, , drop = FALSE], out)
  rownames(out) <- NULL
  attr(out, "param")    <- object$param
  attr(out, "alias")    <- object$alias
  attr(out, "relative") <- relative
  out
}

