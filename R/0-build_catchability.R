#' Catchability parameters that accept a linkage
#' @keywords internal
#' @noRd
Q_LINKAGE_PARAMS <- c("q")


#' @keywords internal
#' @noRd
.validate_q_linkages <- function(linkages) {
  .validate_process_linkages(linkages, Q_LINKAGE_PARAMS, "q")
}


#' Catchability specification
#'
#' @description
#' Holds environmental linkages on survey/index catchability `q`. The effect of
#' an `env_data` covariate is written as a formula and can hold priors, bounds,
#' and an estimation phase like any other linkage.
#'
#' @param linkages Optional named list of [linkage_spec()] objects keyed by
#'   catchability parameter. The only parameter is `q`. Coefficients are per
#'   fleet by default (`by = ~ fleet`); use the `fleet` argument of
#'   [linkage_spec()] to restrict a spec to particular fleets.
#'
#' @details Catchability is the one process that accepts
#'   `link = "exponential"`, Stock Synthesis's environmental link type 1:
#'   `q^exp(beta * x)`, which MULTIPLIES `log q` where `"log"` shifts it. It
#'   needs an estimated base `q` and a lognormal or t index likelihood, and
#'   cannot share a fleet with a random-effect q linkage; each of those is
#'   refused with the reason. See
#'   `vignette("environmental-linkages-and-priors")`.
#'
#' @return A list of catchability settings for [fit_mod()].
#'
#' @examples
#' \donttest{
#' # One temperature effect on q per fleet. `by = ~ fleet` is the default, and a
#' # fleet without an estimated survey catchability is not linkable, so at fit time
#' # restrict the linkage to the estimated-q fleets (here fleet 7 of BS2017SS).
#' build_catchability(linkages = list(
#'   q = linkage_spec(~ temp, fleet = 7)))
#'
#' # Restrict it to fleets 1 and 3, with a prior on the slope
#' build_catchability(linkages = list(
#'   q = linkage_spec(~ temp, fleet = c(1, 3),
#'                    priors = list(temp = prior_normal(0, 1)))))
#' }
#'
#' @export
build_catchability <- function(linkages = NULL) {
  list(linkages = .validate_q_linkages(linkages))
}


# Catchability forms that do NOT estimate q, so a q linkage must not attach:
# "Fixed" holds index_log_q at its input (a linkage would silently turn a fixed
# q time-varying, contrary to the assessor's Fixed setting), and
# "Analytical"/"AnalyticalArith" solve q from the data (index_log_q is mapped
# out, so a linkage targets a non-free parameter and does nothing). Both are
# rejected up front rather than quietly changing q. See build_map_catchability.
.Q_LINKAGE_UNESTIMATED_FORMS <- c("Fixed", "Analytical", "AnalyticalArith")

# A second, different reason to refuse a q linkage. "Environmental" and "AR1" DO
# estimate q, but they rebuild index_q from their own formula in the cpp,
# assigning over the value that carries q_linkage_offset rather than adding to
# it. The linkage is then accepted, REPORTed in q_linkage_offset as a live
# covariate effect, and never enters the likelihood -- beta_linkage is left free
# with an identically-zero gradient. Refuse it, and point at the linkage
# equivalent: these two forms exist to be replaced by exactly that grammar, so a
# user reaching for both wants the linkage on its own.
.Q_LINKAGE_SELFBUILT_FORMS <- c("Environmental", "AR1")


# Informational one-time message: a catchability/selectivity linkage whose `by` was
# auto-filled (the user omitted it) and that names no fleets applies to EVERY fleet
# of that process -- one coefficient each. Tell the user so the per-fleet expansion
# (and any resulting eligibility error) is not a surprise; naming fleets via
# linkage_spec(fleet = ) restricts it.
.message_auto_fleet_linkages <- function(spec_groups) {
  labs <- c(q = "catchability", sel = "selectivity")
  for (proc in names(labs)) {
    specs <- spec_groups[[proc]]
    if (is.null(specs)) next
    for (nm in names(specs)) {
      val <- specs[[nm]]
      lst <- if (inherits(val, "Rceattle_linkage_spec")) list(val) else val
      for (s in lst) {
        if (isTRUE(s$by_auto) && is.null(s$fleet) && "fleet" %in% all.vars(s$by)) {
          message(sprintf(
            paste0("A %s linkage (`%s`) did not name fleets, so it applies to every ",
                   "eligible fleet -- one coefficient each. Pass ",
                   "linkage_spec(fleet = ) to restrict it to specific fleets."),
            labs[[proc]], nm))
        }
      }
    }
  }
  invisible()
}


# Fleets a set of q linkage rows reaches; NA is the shared sentinel, which the
# cpp expands to every fleet.
.q_linkage_fleets <- function(fleet, fleet_control) {
  f <- unique(fleet)
  if (anyNA(f)) seq_len(nrow(fleet_control)) else as.integer(f)
}


#' Reject q linkages on fleets whose catchability is not estimated
#'
#' @param linkage_table pooled linkage table (may be NULL / empty).
#' @param fleet_control the fleet control table.
#' @return invisibly NULL; errors on a q linkage targeting a fleet whose q is
#'   fixed or analytically solved.
#' @keywords internal
#' @noRd
.check_q_linkage_support <- function(linkage_table, fleet_control) {
  if (is.null(linkage_table) || nrow(linkage_table) == 0L) return(invisible())
  q <- linkage_table[linkage_table$process == "q", , drop = FALSE]
  if (nrow(q) == 0L) return(invisible())

  # The cpp expands the NA sentinel to EVERY fleet, so a shared row has to be
  # checked against all of them, not only the fleets other rows name.
  flts <- .q_linkage_fleets(q[["fleet"]], fleet_control)
  forms <- as.character(fleet_control$Catchability[flts])
  # A fleet does not estimate q if its Catchability holds q fixed / solves it from
  # the data (Fixed / Analytical), or is absent (NA) -- a fleet with no survey index
  # has no catchability to link. A linkage on any of these cannot be estimated.
  self <- flts[!is.na(forms) & forms %in% .Q_LINKAGE_SELFBUILT_FORMS]
  if (length(self) > 0) {
    stop(sprintf(
      paste0("catchability linkage on fleet(s) %s whose Catchability (%s) builds ",
             "q from its own formula: the cpp assigns index_q from that formula, ",
             "overwriting the linkage offset instead of adding to it, so the ",
             "linkage would be reported in q_linkage_offset but never enter the ",
             "likelihood.\n",
             "  Express the whole relationship as the linkage and set ",
             "Catchability to \"Estimated\": a covariate effect is ",
             "build_catchability(linkages = list(q = linkage_spec(~ covariate, ",
             "by = ~ fleet))), an AR1 deviation is linkage_spec(ar1(1 | Year), ",
             "by = ~ fleet); see vignette('environmental-linkages-and-priors')."),
      paste(fleet_control$Fleet_name[self], collapse = ", "),
      paste(unique(as.character(fleet_control$Catchability[self])), collapse = ", ")),
      call. = FALSE)
  }

  bad <- flts[is.na(forms) | forms %in% .Q_LINKAGE_UNESTIMATED_FORMS]
  if (length(bad) > 0) {
    stop(sprintf(
      paste0("catchability linkage on fleet(s) %s whose Catchability (%s) does ",
             "not estimate q: index_log_q is held fixed (Fixed), solved from the ",
             "data (Analytical), or absent (NA -- the fleet has no survey index), ",
             "so a linkage would turn a fixed q time-varying or have no effect.\n",
             "  Give the fleet an estimated survey catchability (\"Estimated\" / ",
             "\"Estimated-with-prior\") with index data, or restrict the linkage to ",
             "estimated-q fleets via linkage_spec(fleet = ...)."),
      paste(fleet_control$Fleet_name[bad], collapse = ", "),
      paste(unique(as.character(fleet_control$Catchability[bad])),
            collapse = ", ")), call. = FALSE)
  }

  ex <- q[q[["link"]] == "exponential", , drop = FALSE]
  if (nrow(ex) > 0L) {
    ex_flts <- intersect(.q_linkage_fleets(ex[["fleet"]], fleet_control), flts)

    # SS3 stores q as a log only for a lognormal/t survey (SS_expval.tpl:413-419),
    # so this is the wrong bridge elsewhere -- but a well-defined model, so warn.
    nat <- .index_fleets_natural_scale(fleet_control)
    bad_fam <- ex_flts[nat[ex_flts]]
    if (length(bad_fam) > 0L) {
      warning(sprintf(paste0(
        "link = \"exponential\" on fleet(s) %s, whose Index_distribution (%s) ",
        "scores the index on the\n  NATURAL scale. Rceattle still holds q on ",
        "the log scale, so the model is well defined,\n  but it is NOT Stock ",
        "Synthesis's environmental link type 1 for such a fleet: SS3 reads the ",
        "q\n  slot arithmetically under a natural-scale survey ",
        "(SS_expval.tpl:413-419), where its type 1 is\n  q * exp(beta * x) -- ",
        "Rceattle's link = \"log\". If you are bridging an SS3 model, use ",
        "\"log\"."),
        paste(fleet_control$Fleet_name[bad_fam], collapse = ", "),
        paste(unique(as.character(
          fleet_control[["Index_distribution"]][bad_fam])), collapse = ", ")),
        call. = FALSE)
    }

    # The link multiplies log q, so it needs a base map_linkage_adjuster() leaves
    # free: a slope-only group or a fixed intercept freezes it. NA %in% NA is TRUE,
    # so the shared sentinel frees a shared slope.
    is_icept <- !is.na(q[["design_col"]]) & q[["design_col"]] == "(Intercept)"
    free <- q[["fleet"]][is_icept & as.integer(q[["est_phase"]]) != 0L]
    frozen <- unique(ex[["fleet"]][!(ex[["fleet"]] %in% free)])
    if (length(frozen)) {
      bad_base <- if (anyNA(frozen)) ex_flts else
        intersect(as.integer(frozen), ex_flts)
      stop(sprintf(paste0(
        "link = \"exponential\" needs an estimated base catchability to ",
        "multiply, but the linkage on\n  fleet(s) %s holds index_log_q fixed: ",
        "the formula has no intercept (`~ 0 + x` / `~ x - 1`),\n  or its ",
        "intercept is fixed at 0 (est_phase 0), and either makes ",
        "map_linkage_adjuster()\n  mask index_log_q.\n",
        "  Give the formula an estimated intercept, or use link = \"log\"."),
        paste(fleet_control$Fleet_name[bad_base], collapse = ", ")),
        call. = FALSE)
    }

  }

  # An intercept prior and a Catchability_index group's own q prior (on its lead fleet,
  # as in the template) both penalize the shared log q. A row with no fleet targets fleet 1.
  pri <- q[q$design_col == "(Intercept)" & !is.na(q$prior_family) &
             q$prior_family != "none", , drop = FALSE]
  if (nrow(pri) > 0L) {
    tgt  <- unique(ifelse(is.na(pri$fleet), 1L, as.integer(pri$fleet)))
    grp  <- fleet_control$Catchability_index
    lead <- .group_lead(grp, fleet_control$Fleet_type %in% c("Off", 0, "0")) == 1L
    tgt_lead <- vapply(tgt, function(f) {
      if (is.na(grp[f])) return(as.integer(f))
      which(lead & grp %in% grp[f])[1]
    }, integer(1))
    both <- tgt[as.character(fleet_control$Catchability[tgt_lead]) %in%
                  c("2", "Estimated-with-prior")]
    if (length(both) > 0L) {
      stop(sprintf(paste0(
        "fleet(s) %s carry a prior on the q linkage's intercept, and their q (led by %s) ",
        "has Catchability = \"Estimated-with-prior\". Both are priors on the same q, so it ",
        "would be counted twice. Keep one: set the lead's Catchability to \"Estimated\", or ",
        "drop the intercept prior."),
        paste(fleet_control$Fleet_name[both], collapse = ", "),
        paste(unique(fleet_control$Fleet_name[tgt_lead[tgt %in% both]]), collapse = ", ")),
        call. = FALSE)
    }
  }

  # Fleets sharing a Catchability_index estimate ONE coefficient per design
  # column, and the fixed-beta prior loop runs over every ROW of the linkage
  # table with no lead gate -- so a prior named on two members of the group is
  # evaluated twice on the one coefficient they share. Doubling a Gaussian log
  # prior divides its variance by two: a stated SD of 0.1 is enforced as 0.0707.
  # Same reasoning and same wording as the selectivity guard in
  # `.check_sel_linkage_support()`; this is the catchability half of it. Keep
  # the prior on one row and give the other members a plain row, which is what
  # the coverage rule asks for.
  qpri <- q[!is.na(q[["prior_family"]]) & q[["prior_family"]] != "none" &
              !is.na(q[["fleet"]]), , drop = FALSE]
  if (nrow(qpri) > 0L) {
    qgrp <- fleet_control[["Catchability_index"]]
    for (g in unique(stats::na.omit(qgrp))) {
      mem <- which(!is.na(qgrp) & qgrp == g)
      if (length(mem) < 2L) next
      key <- paste(qpri[["param"]], qpri[["design_col"]], sep = "\r")
      for (k in unique(key)) {
        hit <- unique(qpri[["fleet"]][key == k & qpri[["fleet"]] %in% mem])
        if (length(hit) > 1L) {
          stop(sprintf(paste0(
            "a catchability prior on `%s` names fleet(s) %s, which share ",
            "Catchability_index %d and so estimate ONE coefficient: the prior ",
            "would be counted once per sharing fleet, which divides its stated ",
            "variance by that count. Put the prior on one of them and give the ",
            "rest a plain linkage row for the same column."),
            qpri[["param"]][key == k][1],
            paste(sprintf("'%s'", fleet_control[["Fleet_name"]][hit]),
                  collapse = ", "), g), call. = FALSE)
        }
      }
    }
  }

  # Fleets sharing a Catchability_index estimate one q block, so a linkage on
  # the group has to reach all of it and mean one thing.
  .stop_if_mirrored_block_linkage(linkage_table, fleet_control, "q")
  invisible()
}
