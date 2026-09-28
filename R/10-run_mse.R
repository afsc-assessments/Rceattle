# Derived quantities kept when an MSE operating/estimation model object is
# slimmed for storage; every other quantity is dropped to keep saved runs small.
# Shared by run_mse() and mse_summary() so the retained set stays in one place.
.mse_keep_quantities <- c(
  "catch_hat", "catch_sd", "index_hat", "index_sd",
  "log_catch_sd", "log_index_sd",          # deprecated spellings, dropped next release
  "ssb_depletion", "biomass_depletion", "biomass", "ssb",
  "BO", "SB0", "SBF", "F_spp", "R",
  "M1_at_age", "M_at_age", "avg_rec",
  "DynamicB0", "DynamicSB0", "DynamicSBF",
  "SPR0", "SPRlimit", "SPRtarget", "Ftarget",
  "B_eaten", "B_eaten_as_prey", "Flimit"
)

# Assessment years ------------------------------------------------------------
#
# `assessment_period` is read two ways. A single number is a fixed cycle: an
# assessment every n years, counted from the operating model's terminal year to
# the end of the MSE. A vector is the schedule itself -- the exact years an
# assessment is completed. The two are different management questions: one
# assessment missed inside an otherwise biennial cycle is a single episode of
# stale advice, while a triennial period is a permanent policy, and the cost of
# one is not the cost of the other.
#
# `max_yr` is the last year an assessment can be run in: the earlier of the two
# models' projection horizons, and of `endyr` where the caller set one.
# `proj_last` is the horizon mse_summary() reads, which is the two models' own
# projyr and takes no notice of `endyr` -- so the two differ whenever the caller
# shortens the MSE, and the short-schedule warning has to test against the
# second.
.mse_assess_years <- function(assessment_period, om_endyr, max_yr,
                              proj_last = max_yr) {
  if (!is.numeric(assessment_period) || length(assessment_period) == 0 ||
      anyNA(assessment_period)) {
    stop("`assessment_period` must be a number of years, or a vector of the ",
         "years an assessment is completed.", call. = FALSE)
  }

  fractional <- assessment_period[assessment_period != round(assessment_period)]
  if (length(fractional)) {
    stop("`assessment_period` must be whole years; got ",
         paste(fractional, collapse = ", "), ".", call. = FALSE)
  }

  # * A period ----
  if (length(assessment_period) == 1) {
    if (assessment_period < 1) {
      stop("`assessment_period` must be at least 1 year; got ",
           assessment_period, ".", call. = FALSE)
    }

    # A schedule of one year cannot be told from a period, and the two readings
    # of, say, 2029 are a world apart. Refused rather than guessed at. A period
    # simply longer than the projection lands here too, and is named as such
    # rather than surfacing downstream as a `seq()` error.
    if (om_endyr + assessment_period > max_yr) {
      stop("A single `assessment_period` is the number of years BETWEEN ",
           "assessments, not the year of one: ", om_endyr, " + ",
           assessment_period, " is past the last year an assessment can be ",
           "run in (", max_yr, "), so no assessment would run. Give the ",
           "period, or list two or more assessment years.", call. = FALSE)
    }

    return(.mse_warn_short_schedule(
      seq(from = om_endyr + assessment_period, to = max_yr,
          by = assessment_period),
      proj_last))
  }

  # * An explicit schedule ----
  # round() rather than as.integer(): the years are already whole, and keeping
  # them the same storage mode `seq()` returns keeps the two paths
  # indistinguishable to everything downstream.
  assess_yrs <- sort(unique(round(assessment_period)))

  # Named and refused rather than dropped. A year outside the window cannot host
  # an assessment, and quietly removing it would run a schedule the caller did
  # not ask for and could not see anywhere in the result.
  early <- assess_yrs[assess_yrs <= om_endyr]
  if (length(early)) {
    stop("`assessment_period` years must be after the operating model's ",
         "terminal year (", om_endyr, "); got ",
         paste(early, collapse = ", "), ".", call. = FALSE)
  }

  late <- assess_yrs[assess_yrs > max_yr]
  if (length(late)) {
    stop("`assessment_period` years must be within the projection horizon ",
         "(through ", max_yr, "); got ", paste(late, collapse = ", "), ".",
         call. = FALSE)
  }

  .mse_warn_short_schedule(assess_yrs, proj_last)
}

# A schedule that stops short of the projection horizon leaves the trailing
# years with no catch at all: run_mse() only ever fills catch up to the last
# assessment, so those years keep the NA the projection rows were created with.
#
# mse_summary() then counts them anyway. `projyrs` there runs to the models'
# own projyr, so the NA years sit in the denominator of P(Closed) but not its
# numerator, in both halves of Catch IAV, and outside the na.rm mean that is
# Average Catch. All three come out low, with nothing in the table saying which
# years were summarised. A biennial cycle over an odd number of projection
# years lands here every time -- 2027(2)2049 against a 2050 horizon.
#
# Warned rather than refused: a schedule that stops short is a legitimate
# design, it just has to be summarised over the years it covers.
#
# Tested against `proj_last`, the models' own projyr, NOT the last year an
# assessment may run in. `endyr` shortens the second and not the first, so a
# caller who ends the MSE early leaves every remaining projection year unfished
# -- the same understatement, over more years.
.mse_warn_short_schedule <- function(assess_yrs, proj_last) {
  last <- max(assess_yrs)
  if (last < proj_last) {
    unfished <- if (last + 1 == proj_last) {
      as.character(proj_last)
    } else {
      paste0(last + 1, "-", proj_last)
    }
    warning("The last assessment is in ", last, " but the models project to ",
            proj_last, ", so catch in ", unfished, " is never set. Those years ",
            "carry NA, and mse_summary() understates Average Catch, Catch IAV ",
            "and P(Closed) because it summarises over the whole projection. ",
            "Either set projyr to ", last, " in both models, or run ",
            "assessments through ", proj_last, ".", call. = FALSE)
  }
  assess_yrs
}

# Catch multiplier ------------------------------------------------------------
#
# `catch_mult` scales the catch the estimation model's control rule recommends,
# before the cap and before the exploitable-biomass limit. A single number or a
# vector of length nspp applies in every projection year, which is a permanent
# change in harvest policy. A data.frame of Year / Species / mult applies only
# in the pairs it lists -- a buffer held for the years advice is stale, say --
# and every pair it does not list is multiplied by 1.
#
# `Species` is the species number, the same code the catch data carries.
# Validated up front so a mistyped year or species number is reported once at
# the call, not silently ignored inside every simulation as a multiplier of 1.
.mse_check_catch_mult <- function(catch_mult, nspp, proj_first, proj_last) {
  if (is.null(catch_mult)) return(NULL)

  # * Year- and species-indexed ----
  if (is.data.frame(catch_mult)) {
    absent <- setdiff(c("Year", "Species", "mult"), names(catch_mult))
    if (length(absent)) {
      stop("A `catch_mult` data.frame needs columns Year, Species and mult; ",
           "missing ", paste(absent, collapse = ", "), ".", call. = FALSE)
    }
    if (nrow(catch_mult) == 0) {
      stop("`catch_mult` has no rows. Pass NULL for no multiplier.",
           call. = FALSE)
    }
    if (!is.numeric(catch_mult$mult) || anyNA(catch_mult$mult) ||
        any(!is.finite(catch_mult$mult)) || any(catch_mult$mult < 0)) {
      stop("`catch_mult$mult` must be finite and non-negative.", call. = FALSE)
    }
    if (!is.numeric(catch_mult$Species) || anyNA(catch_mult$Species) ||
        any(!catch_mult$Species %in% seq_len(nspp))) {
      stop("`catch_mult$Species` must be a species number in 1:", nspp, ".",
           call. = FALSE)
    }
    if (!is.numeric(catch_mult$Year) || anyNA(catch_mult$Year) ||
        any(catch_mult$Year != round(catch_mult$Year))) {
      stop("`catch_mult$Year` must be a whole year.", call. = FALSE)
    }

    off <- catch_mult$Year[catch_mult$Year < proj_first |
                             catch_mult$Year > proj_last]
    if (length(off)) {
      stop("`catch_mult$Year` must be a year the assessment schedule covers (",
           proj_first, ":", proj_last, "); got ",
           paste(sort(unique(off)), collapse = ", "), ".", call. = FALSE)
    }

    key <- paste(catch_mult$Year, catch_mult$Species)
    if (anyDuplicated(key)) {
      stop("`catch_mult` gives more than one multiplier for the same year and ",
           "species: ", paste(unique(key[duplicated(key)]), collapse = "; "),
           ".", call. = FALSE)
    }

    return(catch_mult)
  }

  # * One multiplier per species, every year ----
  if (!is.numeric(catch_mult) || anyNA(catch_mult) ||
      any(!is.finite(catch_mult)) || any(catch_mult < 0)) {
    stop("`catch_mult` must be finite and non-negative.", call. = FALSE)
  }
  if (length(catch_mult) == 1) {
    catch_mult <- rep(catch_mult, nspp)
  }
  if (length(catch_mult) != nspp) {
    stop("catch_mult is not length 1 or length nspp")
  }
  catch_mult
}

# Apply the multiplier to the catch rows filled for this assessment interval.
# `year` and `species` are those rows' own Year and Species columns.
.mse_apply_catch_mult <- function(catch, year, species, catch_mult) {
  if (is.data.frame(catch_mult)) {
    idx <- match(paste(year, species),
                 paste(catch_mult$Year, catch_mult$Species))
    mult <- catch_mult$mult[idx]
    mult[is.na(idx)] <- 1   # a pair the caller did not list is left alone
    return(catch * mult)
  }
  catch * catch_mult[species]
}

# Projection catch -------------------------------------------------------------
#
# Catch recorded past a model's terminal year is arbitrary "future" data, the
# same as the index and composition rows the MSE setup filters out. It cannot be
# dropped the same way: those rows ARE the projection, and each assessment fills
# them with the catch its control rule recommends. They are blanked to NA
# instead, which is the state `clean_data()` creates a projection row in.
#
# Nothing the model uses is lost. The likelihood scores only `Year <= endyr`, so
# these rows never entered the fit, and no harvest control rule reads observed
# catch -- HCR 1 (constant catch) sums `catch_hat`, not the data. Section 3 of
# the assessment loop already overwrites every year in the interval with the
# operating model's realized catch, so a value left here survives one assessment
# at most.
#
# What it costs to leave them is the fill in section 1, which selects rows on
# `is.na(Catch)`. A projection year arriving with a number is skipped: no catch
# is written for it, the interval total comes out 0, and the operating model is
# REBUILT rather than advanced -- the no-catch path, `estimateMode = 3` with the
# map dropped. A workbook whose terminal year lags its catch series (catch
# through 2023, `endyr` still 2019) lands there for every assessment before the
# series ends.
#
# Warned rather than blanked in silence. Catch recorded for a year the model
# does not fit usually means `endyr` is out of date, and conditioning on those
# years is a different MSE than projecting over them.
.mse_blank_proj_catch <- function(catch_data, endyr, model) {
  if (is.null(catch_data) || !nrow(catch_data)) return(catch_data)

  filled <- which(catch_data$Year > endyr & !is.na(catch_data$Catch))
  if (!length(filled)) return(catch_data)

  yrs <- sort(unique(catch_data$Year[filled]))
  warning("The ", model, " model records catch in ",
          paste(yrs, collapse = ", "), ", past its terminal year (", endyr,
          "). Those are projection years, whose catch the MSE sets from the ",
          "control rule, so the recorded values are dropped.", call. = FALSE)

  catch_data$Catch[filled] <- NA
  catch_data
}

# Projection rows for the fixed biological inputs ------------------------------
#
# `weight` (weight-at-age, kg) and `ration_data` (annual foraging days) are
# carried into the projection by repeating each series' terminal year, one row
# per projection year -- the MSE holds both at the operating model's last
# hindcast year.
#
# A series supplied for Year 0 is exempt. That is the time-invariant convention,
# and `rearrange_data()` already fills every hindcast year from the single row,
# so there is nothing for the projection to extend. Expanding one replaces a
# legal one-row series with one naming the projection years and no year before
# them; `data_check()` reads any index carrying more than one year as
# time-varying and requires it to span `styr..endyr`, so the first assessment --
# which advances `endyr` -- fails with "Weight data for index = 4 & sex = 1 does
# not span all hindcast years". The check is right; the expansion was giving it
# something to reject.
#
# Grouped on the series id and `Sex`, and carried forward from the group's last
# row at or before `endyr` -- the terminal HINDCAST year, which is the value
# `run_mse()` documents holding through the projection. A workbook whose series
# runs past its own `endyr` therefore does not project its post-terminal rows
# forward; those rows are kept where they are, and a projection year that
# already has a row of its own keeps it rather than gaining a second. Two rows
# for one year is not inert: `rearrange_data()` assigns row by row into the
# weight array, so the later row wins, and the carried-forward value would
# silently replace an observation. `clean_data()` skips an already-present
# projection year for catch, and the `NByageFixed` expansion below does the
# same.
.mse_expand_fixed_input <- function(dat, id_col, endyr, proj_yrs) {
  if (is.null(dat) || !nrow(dat)) return(dat)

  key <- paste(dat[[id_col]], dat$Sex)
  out <- dat
  for (k in unique(key)) {
    rows <- which(key == k)
    if (all(dat$Year[rows] == 0)) next   # time-invariant; already spans the projection

    add <- setdiff(proj_yrs, dat$Year[rows])
    if (!length(add)) next

    # Year 0 alongside dated rows is a mixed series data_check() rejects; take
    # the LATEST hindcast row, falling back to the group's latest row so a
    # series ending before styr still carries something forward.
    #
    # By year, not by table position: nothing sorts `weight` on read, so the
    # last row of a group need not be its last year. Taking it by position
    # carries the wrong weight-at-age through the whole projection.
    hind <- rows[dat$Year[rows] > 0 & dat$Year[rows] <= endyr]
    src  <- if (length(hind)) hind[which.max(dat$Year[hind])]
            else rows[which.max(dat$Year[rows])]

    proj <- dat[rep(src, length(add)), , drop = FALSE]
    proj$Year <- add
    out <- rbind(out, proj)
  }

  out <- out[order(out[[id_col]], out$Year), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Run a management strategy evaluation
#'
#' @description Runs a forward-projecting management strategy evaluation (MSE). Projected selectivity, catchability, foraging days, and weight-at-age are held at the operating model's terminal hindcast year. Survey SD is set to the average over the historical time series, and composition sample size is held at the last year. There is no implementation error and no observation error on catch.
#'
#' @param om CEATTLE model object exported from \code{Rceattle}
#' @param em CEATTLE model object exported from \code{Rceattle}
#' @param nsim Number of simulations to run (default 10)
#' @param start_sim First simulation number to start at. Useful if the code stops at specific seed/sim (default = 1).
#' @param assessment_period Assessment schedule. A single number is the period in years between assessments, counted from the operating model's terminal year. A vector is the explicit set of years an assessment is completed. Default = 1.
#' @param sampling_period Period of years data sampling is conducted. Single value or vector the same length as the number of fleets.
#' @param simulate_data Include simulated random error proportional to that estimated/provided for the data from the OM.
#' @param regenerate_past Refits the EM to historical/conditioning data prior to the MSE, where the data are generated from the OM with \code{simulate_data = TRUE} or without \code{simulate_data = FALSE} sampling error.
#' @param sample_rec Include resampled recruitment deviations from the hindcast in the OM projection. Resampled deviations are used rather than drawing from N(0, sigmaR) because the initial deviations bias R0 low. If FALSE, uses the single deviation described in [sample_rec()].
#' @param rec_trend Linear increase or decrease in mean recruitment from \code{endyr} to \code{projyr}. This is the terminal multiplier \code{mean rec * (1 + (rec_trend/projection years) * 1:projection years)}. Can be of length 1 or of length nspp. If length 1, all species get the same trend.
#' @param fut_sample future sampling effort relative to last year.  \code{ Log_sd * 1 / fut_sample} for index and \code{ Sample_size * fut_sample} for comps
#' @param cap A cap on the catch in the projection. Can be a single number applied to all species (proportional to recommended catch) or vector of length \code{nspp} applied to each species. Default = NULL
#' @param catch_mult A multiplier applied to the catch the control rule recommends. A single number or a vector of length \code{nspp} applies in every projection year; a \code{data.frame} with columns \code{Year}, \code{Species} and \code{mult} applies only in the year and species pairs it lists, where \code{Year} is a year the assessment schedule covers. Default = NULL
#' @param loopnum number of times to re-start optimization (where \code{loopnum=3} sometimes achieves a lower final gradient than \code{loopnum=1})
#' @param file (Optional) Filename where each OM simulation with EMs will be saved. If NULL, no files are saved.
#' @param dir (Optional) Directory where each OM simulation is saved
#' @param seed seed for the simulation
#' @param regenerate_seed seed for regenerating data
#' @param timeout length of time (minutes) estimation will run before stopping a sim (default 999 minutes)
#' @param endyr Terminal year of the MSE projection. Default = NA uses \code{projyr} from the operating model.
#' @param cores Number of cores to use for parallel simulations. Default
#'   \code{NULL} picks \code{parallel::detectCores() - 6}, capped at 2 when
#'   running under \code{R CMD check} (which sets
#'   \code{_R_CHECK_LIMIT_CORES_}). Set to 1 to force sequential execution.
#'
#' @details
#' # Assessment schedule
#'
#' A single `assessment_period` is a fixed cycle. A vector is the schedule
#' itself, for a design whose years are not evenly spaced, one assessment
#' missed inside an otherwise biennial cycle, for instance:
#'
#' ```
#' biennial <- seq(om$data_list$endyr + 2, om$data_list$projyr, by = 2)
#' run_mse(om, em, assessment_period = setdiff(biennial, 2031))
#' ```
#'
#' Every year must be after the operating model's terminal year and within the
#' projection horizon; a year outside that window is an error rather than being
#' dropped. One missed assessment and a permanently longer cycle are different
#' questions, and the second does not answer the first.
#'
#' A schedule must name two or more years. One year on its own cannot be told
#' from a period, and the two readings are a world apart, so it is refused
#' rather than guessed at.
#'
#' Two schedules run on the same `seed` share their recruitment deviations,
#' which are drawn once per replicate before the first assessment, and each
#' assessment's observation draw is seeded on that assessment's own year rather
#' than on wherever earlier assessments left the stream. So an assessment in
#' year `Y` starts from the same place in every schedule that assesses `Y`, and
#' a divergence at one assessment no longer displaces every assessment after it.
#'
#' Common random numbers are not complete, and the gap is worth knowing before
#' designing a comparison. One `sim_mod()` call draws every year in the
#' assessment interval, under that assessment's seed --- so a year sitting
#' *inside* a longer interval is drawn under a different seed than the same year
#' in a schedule that assessed it directly. Two schedules that differ in which
#' years they assess therefore realize different observation error in the years
#' between, even where the stock is in the same state. Per-observation-year
#' streams would close it; `inst/dev/TODO-mse-horizon.md` records what that
#' would take.
#'
#' Past the point where two schedules genuinely diverge --- different catch
#' taken, different years surveyed --- the draws necessarily diverge too. That
#' is a real difference in the runs, not an artefact.
#'
#' Make the last assessment year reach the projection horizon. Catch is only
#' ever filled up to the last assessment, so a schedule that stops short leaves
#' the trailing years at `NA`, and `mse_summary()` summarises over the whole
#' projection --- understating Average Catch, Catch IAV and P(Closed) with
#' nothing in the table to say which years it covered. A biennial cycle over an
#' odd number of projection years lands here every time, so this is warned
#' about, for a period and a schedule alike.
#'
#' # Catch multiplier
#'
#' `catch_mult` multiplies the catch the estimation model's control rule
#' recommends, before `cap` and before the exploitable-biomass limit. Given as
#' a `data.frame` it applies only in the years and species it lists, which is
#' how a buffer held for the years advice is stale is expressed.
#'
#' Mind which years those are. The assessment in year `Y` sets catch for `Y+1`
#' onward, so a missed assessment in 2031 leaves **2032 and 2033** on 2029's
#' advice --- 2031 itself is set by the 2029 assessment either way, and cutting
#' it changes a year the missed assessment never touched:
#'
#' ```
#' buffer <- expand.grid(Year = 2032:2033, Species = seq_len(om$data_list$nspp))
#' buffer$mult <- 0.90
#' run_mse(om, em, assessment_period = setdiff(biennial, 2031),
#'         catch_mult = buffer)
#' ```
#'
#' `Species` is the species number, matching the catch data's own column, and
#' any (year, species) pair the table omits is multiplied by 1. `Year` must be
#' a year the assessment schedule actually covers --- from the operating
#' model's terminal year to the last assessment, which is where catch is
#' filled, not the whole projection.
#'
#' Note that this reduces catch, not ABC. Where realized catch sits well below
#' ABC, GOA arrowtooth flounder, for one, reducing ABC changes removals only
#' to the extent the fishery attains it, while reducing catch changes them in
#' full. Either scale the multiplier by recent attainment,
#' `1 - (1 - mult) * attainment`, or report the unscaled result as an upper
#' bound on the effect of the reduction.
#'
#' # Projection catch
#'
#' Catch recorded past a model's terminal year is blanked to `NA` at setup, with
#' a warning naming the years. Those years are the projection, and the MSE sets
#' their catch from the control rule; the likelihood never scored them either,
#' since it fits only `Year <= endyr`. This is what a workbook looks like when
#' `endyr` has fallen behind the catch series, catch through 2023 with `endyr`
#' still 2019, and it is worth resolving before running the MSE, because
#' conditioning the assessment on those years is a different question from
#' projecting over them.
#'
#' @return A list of operating models (differ by simulated recruitment determined by \code{nsim}) and estimation models fit to each operating model (differ by terminal year).
#' @examples
#' \dontrun{
#' data(BS2017SS)
#' om <- fit_mod(BS2017SS, estimateMode = "Estimate", HCR = build_hcr("NPFMC"))
#' # Closed loop: assess every 2 years, survey every year, 10 simulations.
#' mse <- run_mse(om = om, em = om, nsim = 10,
#'                assessment_period = 2, sampling_period = 1)
#' mse_summary(mse)$species
#' }
#' @export
run_mse <- function(om, em, nsim = 10, start_sim = 1, assessment_period = 1, sampling_period = 1, simulate_data = TRUE, regenerate_past = FALSE, sample_rec = TRUE, rec_trend = 0, fut_sample = 1, cap = NULL, catch_mult = NULL, seed = 666, regenerate_seed = seed, loopnum = 1, file = NULL, dir = NULL, timeout = 999, endyr = NA, cores = NULL){

  # om = om; em = em; nsim = 1; start_sim = 1; assessment_period = 1; sampling_period = 1; simulate_data = TRUE; regenerate_past = FALSE; sample_rec = FALSE; rec_trend = 0; fut_sample = 1; cap = NULL; catch_mult = NULL; seed = 666; regenerate_seed = seed; loopnum = 1; file = NULL; dir = NULL; endyr = NA; timeout = 999

  if (missing(om) || is.null(om) || !inherits(om, "Rceattle")) {
    stop("`om` must be an Rceattle model object (see ?fit_mod).")
  }
  if (missing(em) || is.null(em) || !inherits(em, "Rceattle")) {
    stop("`em` must be an Rceattle model object (see ?fit_mod).")
  }

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # MSE SETUP ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  '%!in%' <- function(x,y)!('%in%'(x,y))
  set.seed(regenerate_seed)

  Rceattle_OM_list <- list()
  Rceattle_EM_list <- list()

  # * Input checks ----
  # - Set om to project from R0
  om$data_list$proj_mean_rec = TRUE # - Sample rec devs assuming this down the line

  # - Adjust cap
  if(!is.null(cap)){
    if(!length(cap) %in% c(1, om$data_list$nspp)){
      stop("cap is not length 1 or length nspp")
    }
  }

  # - The two models must share a hindcast and a projection horizon
  # The assessment loop reads the row positions of catch rows to fill off the
  # ESTIMATION model's catch table and indexes the OPERATING model's
  # `max_catch_hat` with them, so the two tables have to describe the same rows
  # in the same order. Terminal and projection years that disagree break that
  # silently -- the exploitable-biomass limit is then read off the wrong years,
  # or off NA -- so it is refused here rather than left to surface as catch
  # advice that is wrong without being obviously wrong.
  if(om$data_list$endyr != em$data_list$endyr ||
     om$data_list$projyr != em$data_list$projyr){
    stop("The operating and estimation models must share a terminal year and a ",
         "projection horizon; got endyr ", om$data_list$endyr, " / ",
         em$data_list$endyr, " and projyr ", om$data_list$projyr, " / ",
         em$data_list$projyr, ". The MSE fills catch by row position across ",
         "both models, which requires the same years in the same order.",
         call. = FALSE)
  }

  # - Blank catch recorded past either model's terminal year
  # Ahead of the Proj_F_proportion check below, which compares recorded catch
  # against projected exploitable biomass: catch in a projection year the
  # operating model fishes at F = 0 fails that comparison and forces a rebuild
  # the model does not need.
  om$data_list$catch_data <- .mse_blank_proj_catch(
    om$data_list$catch_data, om$data_list$endyr, "operating")
  em$data_list$catch_data <- .mse_blank_proj_catch(
    em$data_list$catch_data, em$data_list$endyr, "estimation")

  # na.rm: Proj_F_proportion is NA for fleets that never take catch (surveys),
  # which is a legitimate workbook value. Without it the sum is NA and the `if`
  # errors with "missing value where TRUE/FALSE needed" instead of reporting
  # whether any fleet is set to take projected F.
  if(sum(om$data_list$fleet_control$Proj_F_proportion, na.rm = TRUE) == 0){
    stop("F prop per fleet 'Proj_F_proportion' is zero")
  }

  # ** Refit OM if Proj_F_proportion was not activated ----
  if(sum((om$data_list$catch_data$Catch > 0) - (om$quantities$max_catch_hat > 0), na.rm = TRUE) > 0){
    # -- Set estimate mode back to original
    estimate_mode_base <- om$data_list$estimateMode

    # Rerun OM in debug mode to make sure F-prop is set correctly. Reuses the
    # OM's own configuration unchanged (estimateMode = 3 builds without
    # optimizing, so loopnum is inert here).
    om <- .refit_like(
      data_list    = om$data_list,
      inits        = om$estimated_params,
      map          = om$map,
      estimateMode = 3)

    # Adjust back
    om$data_list$estimateMode <- estimate_mode_base
  }

  # - Years for simulations
  hind_yrs <- (em$data_list$styr) : em$data_list$endyr
  hind_nyrs <- length(hind_yrs)
  om_proj_yrs <- (om$data_list$endyr + 1) : om$data_list$projyr
  om_proj_nyrs <- length(om_proj_yrs)

  em_proj_yrs <- (em$data_list$endyr + 1) : em$data_list$projyr
  em_proj_nyrs <- length(em_proj_yrs)
  nflts = nrow(om$data_list$fleet_control)

  # - N sel ages for sel coff dev
  if(all(is.na(om$data_list$fleet_control$N_sel_bins))){
    n_sel_bins_om = dim(om$estimated_params$sel_coff_dev)[3]
  } else {
    n_sel_bins_om <- max(om$data_list$fleet_control$N_sel_bins, na.rm = TRUE)
  }

  if(all(is.na(em$data_list$fleet_control$N_sel_bins))){
    n_sel_bins_em = dim(em$estimated_params$sel_coff_dev)[3]
  } else {
    n_sel_bins_em <- max(em$data_list$fleet_control$N_sel_bins, na.rm = TRUE)
  }

  # - Assessment period, or an explicit schedule of assessment years
  assess_yrs <- .mse_assess_years(
    assessment_period,
    om_endyr  = om$data_list$endyr,
    max_yr    = min(c(om$data_list$projyr, em$data_list$projyr, endyr),
                    na.rm = TRUE),
    proj_last = min(om$data_list$projyr, em$data_list$projyr))

  # - Adjust catch multiplier
  # Checked against the years catch is actually FILLED for, which is the
  # assessment schedule's own span -- not the projection horizon. The loop only
  # writes catch up to the last assessment, so a multiplier named for a later
  # year applies nowhere, which is the silent no-op this validation exists to
  # catch.
  catch_mult <- .mse_check_catch_mult(catch_mult,
                                      nspp       = om$data_list$nspp,
                                      proj_first = om$data_list$endyr + 1,
                                      proj_last  = max(assess_yrs))

  # - Data sampling period
  if(length(sampling_period)==1){
    sampling_period = rep(sampling_period, nflts)

  }

  if(nflts != nrow(em$data_list$fleet_control)){
    stop("OM and EM fleets do not match or sampling period length is mispecified")
  }
  if(nflts != length(sampling_period)){
    stop("Sampling period length is mispecified, does not match number of fleets")
  }

  # - Set up years of data we are sampling for each fleet
  sample_yrs <- lapply(sampling_period, function(x) seq(from = em$data_list$endyr + x, to = em$data_list$projyr,  by = x))
  fleet_id <- sample_yrs
  # sampling_period is given per fleet in fleet_control row order, but the table
  # built here is matched against the data's Fleet_code column. Carry the fleet's
  # own code rather than its row position: data_check() does require the two to
  # agree, but nothing here should depend on that silently.
  fleet_codes <- em$data_list$fleet_control$Fleet_code
  for(i in 1:length(sample_yrs)){
    fleet_id[[i]] <- replace(fleet_id[[i]], values = fleet_codes[i])
  }
  sample_yrs = data.frame(Fleet_code = unlist(fleet_id), Year = unlist(sample_yrs))


  # * Filter arbitrary "future" data ----
  # -- index_data
  om$data_list$index_data <- om$data_list$index_data |>
    dplyr::filter(abs(Year) <= om$data_list$endyr)
  em$data_list$index_data <- em$data_list$index_data |>
    dplyr::filter(abs(Year) <= em$data_list$endyr)

  # -- comp_data
  om$data_list$comp_data <- om$data_list$comp_data |>
    dplyr::filter(abs(Year) <= om$data_list$endyr)
  em$data_list$comp_data <- em$data_list$comp_data |>
    dplyr::filter(abs(Year) <= em$data_list$endyr)

  # -- caal_data
  om$data_list$caal_data <- om$data_list$caal_data |>
    dplyr::filter(abs(Year) <= om$data_list$endyr)
  em$data_list$caal_data <- em$data_list$caal_data |>
    dplyr::filter(abs(Year) <= em$data_list$endyr)


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Regenerate past data from OM and refit EM ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if(regenerate_past){

    # - Simulate index and comp data and updatae EM
    sim_dat <- sim_mod(om, simulate = FALSE)

    em$data_list$index_data <- sim_dat$index_data
    em$data_list$comp_data <- sim_dat$comp_data
    em$data_list$caal_data <- sim_dat$caal_data

    # Re-estimate
    em <- .refit_like(
      data_list    = em$data_list,
      inits        = em$estimated_params,
      estimateMode = ifelse(em$data_list$estimateMode < 3, 0, em$data_list$estimateMode),
      loopnum      = loopnum)

    # Update avg F given model fit to regenerated data
    if(em$data_list$HCR == 2){

      # - Get avg F
      avg_F <- exp(em$estimated_params$log_F) # Average F from last 5 years
      avg_F <- rowMeans(avg_F[,(ncol(avg_F)-4) : ncol(avg_F)])
      avg_F <- data.frame(avg_F = avg_F, spp = em$data_list$fleet_control$Species)
      avg_F <- avg_F |>
        dplyr::group_by(spp) |>
        dplyr::summarise(avg_F = sum(avg_F)) |>
        dplyr::arrange(spp)

      # - Update model: project on the recomputed average F (input-F HCR).
      em <- .refit_like(
        data_list    = em$data_list,
        inits        = em$estimated_params,
        estimateMode = 2,   # Don't estimate
        HCR          = build_hcr(HCR     = 2,   # Input F
                                 Ftarget = avg_F$avg_F,
                                 Ptarget = em$data_list$Ptarget,
                                 Plimit  = em$data_list$Plimit),
        loopnum      = loopnum)
    }
  }

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Expand OM data-dim ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # -- index_data
  proj_srv <- om$data_list$index_data |>
    dplyr::group_by(Fleet_code) |>
    dplyr::slice(rep(dplyr::n(),  om_proj_nyrs)) |>
    dplyr::mutate(Year = -om_proj_yrs)
  proj_srv$Log_sd <- proj_srv$Log_sd * 1/fut_sample
  proj_srv$Observation <- NA
  om$data_list$index_data  <- rbind(om$data_list$index_data, proj_srv)
  om$data_list$index_data <- dplyr::arrange(om$data_list$index_data, Fleet_code, abs(Year))

  # -- Nbyage
  if(nrow(om$data_list$NByageFixed) > 0){
    proj_nbyage <- om$data_list$NByageFixed |>
      dplyr::group_by(Species, Sex) |>
      dplyr::slice(rep(dplyr::n(),  om_proj_nyrs)) |>
      dplyr::mutate(Year = om_proj_yrs)
    proj_nbyage <- proj_nbyage[which(om_proj_yrs %!in% om$data_list$NByageFixed$Year),] # Subset rows already forcasted
    om$data_list$NByageFixed  <- rbind(om$data_list$NByageFixed, proj_nbyage)
    om$data_list$NByageFixed <- dplyr::arrange(om$data_list$NByageFixed, Species, Year)
  }

  # -- comp_data
  if(nrow(om$data_list$comp_data) > 0){
    proj_comp <- om$data_list$comp_data |>
      dplyr::group_by(Fleet_code, Sex) |>
      dplyr::slice(rep(dplyr::n(),  om_proj_nyrs)) |>
      dplyr::mutate(Year = -om_proj_yrs)
    proj_comp$Sample_size <- proj_comp$Sample_size * fut_sample # Adjust future sampling effort
    proj_comp <- proj_comp |>
      dplyr::mutate_at(vars(matches("Comp_")), ~ 1)
    om$data_list$comp_data  <- rbind(om$data_list$comp_data, proj_comp)
    om$data_list$comp_data <- dplyr::arrange(om$data_list$comp_data, Fleet_code, abs(Year))
  }

  # -- caal_data
  if(nrow(om$data_list$caal_data) > 0){
    proj_caal <- om$data_list$caal_data |>
      dplyr::group_by(Fleet_code, Sex) |>
      dplyr::slice(rep(dplyr::n(),  om_proj_nyrs)) |>
      dplyr::mutate(Year = -om_proj_yrs)
    proj_caal$Sample_size <- proj_caal$Sample_size * fut_sample # Adjust future sampling effort
    proj_caal <- proj_caal |>
      dplyr::mutate_at(vars(matches("CAAL_")), ~ 1)
    om$data_list$caal_data  <- rbind(om$data_list$caal_data, proj_caal)
    om$data_list$caal_data <- dplyr::arrange(om$data_list$caal_data, Fleet_code, abs(Year))
  }

  # -- emp_sel - Use terminal year
  if(nrow(om$data_list$emp_sel) > 0){
    proj_emp_sel <- om$data_list$emp_sel |>
      dplyr::group_by(Fleet_code, Sex) |>
      dplyr::slice(rep(dplyr::n(),  om_proj_nyrs)) |>
      dplyr::mutate(Year = om_proj_yrs)
    om$data_list$emp_sel  <- rbind(om$data_list$emp_sel, proj_emp_sel)
    om$data_list$emp_sel <- dplyr::arrange(om$data_list$emp_sel, Fleet_code, Year)
  }

  # -- weight
  #FIXME ignores forecasted growth
  om$data_list$weight <- .mse_expand_fixed_input(
    om$data_list$weight, "Wt_index", om$data_list$endyr, om_proj_yrs)

  # -- ration_data
  om$data_list$ration_data <- .mse_expand_fixed_input(
    om$data_list$ration_data, "Species", om$data_list$endyr, om_proj_yrs)

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Expand EM data-dim ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

  #FIXME - assuming same as terminal year of hindcast
  # -- EM emp_sel - Use terminal year
  if(nrow(em$data_list$emp_sel) > 0){
    proj_emp_sel <- em$data_list$emp_sel |>
      dplyr::group_by(Fleet_code, Sex) |>
      dplyr::slice(rep(dplyr::n(),  em_proj_nyrs)) |>
      dplyr::mutate(Year = em_proj_yrs)
    em$data_list$emp_sel  <- rbind(em$data_list$emp_sel, proj_emp_sel)
    em$data_list$emp_sel <- dplyr::arrange(em$data_list$emp_sel, Fleet_code, Year)
  }

  # -- EM weight
  em$data_list$weight <- .mse_expand_fixed_input(
    em$data_list$weight, "Wt_index", em$data_list$endyr, em_proj_yrs)

  # -- EM ration_data
  em$data_list$ration_data <- .mse_expand_fixed_input(
    em$data_list$ration_data, "Species", em$data_list$endyr, em_proj_yrs)


  # Cross-platform parallel via parallel::parLapply on a PSOCK cluster.
  # We do *not* use foreach::%dopar% here: under nested test_that
  # backtraces it triggered an 'evaluation nested too deeply: infinite
  # recursion' abort inside rlang's expression deparser (foreach
  # captures call frames that recurse during error formatting). PSOCK
  # clusters work identically on Windows and Unix and avoid that.
  # Respect the CRAN core limit ('_R_CHECK_LIMIT_CORES_' is set during
  # R CMD check; parallel::makeCluster errors if we exceed 2 cores then).
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  if (is.null(cores)) {
    cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  } else {
    cores <- max(1L, as.integer(cores))
    if (cran_cap) cores <- min(cores, 2L)
  }
  use_parallel <- nsim > 1L && cores > 1L

  # TODO: extract run_one_sim as a top-level internal helper with
  # explicit args (om, em, seed, assess_yrs, ...) for testability.
  # Inline closure for now to ship the foreach -> parLapply migration.
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Per-simulation closure (run_one_sim) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  run_one_sim <- function(sim) {

    set.seed(seed = seed + sim) # setting unique seed for each simulation
    kill_sim <- list(kill_sim = FALSE, failure = NA)

    # Set models objects
    sim_list <- list(EM = list())# , OM = list())
    sim_list$EM[[1]] <- em
    # sim_list$OM[[1]] <- om

    em_use <- em
    om_use <- om

    # Sample recruitment
    om_use <- Rceattle::sample_rec(om_use, sample_rec = sample_rec, update_model = FALSE, rec_trend = rec_trend)

    # A seed per projection year for the observation draws, keyed by the YEAR
    # rather than by an assessment's position in the schedule, and drawn over a
    # year list that does not depend on the schedule at all. Each assessment
    # then starts its draw from the same place whatever schedule the run is on,
    # so a divergence in one assessment cannot displace the stream every later
    # assessment draws from. That compounding is what stopped two schedules
    # being comparable replicate by replicate.
    #
    # Drawn AFTER sample_rec(), so recruitment deviations consume the
    # per-simulation stream exactly as they did before and are unchanged.
    draw_seeds <- stats::setNames(
      sample.int(.Machine$integer.max, length(om_proj_yrs)),
      as.character(om_proj_yrs))

    # Run through assessment years
    for(k in 1:length(assess_yrs)){

      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 1. Get recommended catch from the EM-HCR ----
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      new_years <- om_proj_yrs[which(om_proj_yrs <= assess_yrs[k] & om_proj_yrs > om_use$data_list$endyr)]

      # - Get projected catch data from EM
      new_catch_data <- em_use$data_list$catch_data
      dat_fill_ind <- which(new_catch_data$Year %in% new_years & is.na(new_catch_data$Catch))
      new_catch_data$Catch[dat_fill_ind] <- em_use$quantities$catch_hat[dat_fill_ind]

      # * Catch multiplier ----
      if(!is.null(catch_mult)){
        new_catch_data$Catch[dat_fill_ind] <- .mse_apply_catch_mult(
          catch      = new_catch_data$Catch[dat_fill_ind],
          year       = new_catch_data$Year[dat_fill_ind],
          species    = new_catch_data$Species[dat_fill_ind],
          catch_mult = catch_mult)
      }

      # * Apply cap ----
      if(!is.null(cap)){
        # Applied across species
        if(length(cap) == 1){
          # The cap is an ANNUAL ceiling on total removals across species, so it
          # is applied within each projection year. `dat_fill_ind` spans the whole
          # assessment interval, which is more than one year whenever
          # assessment_period > 1; summing over it would hold a multi-year total
          # to one year's cap. At assessment_period = 1 each group is a single
          # year and this is the previous calculation exactly.
          for(cap_yr in unique(new_catch_data$Year[dat_fill_ind])){
            yr_ind <- dat_fill_ind[new_catch_data$Year[dat_fill_ind] == cap_yr]
            yr_tot <- sum(new_catch_data$Catch[yr_ind])
            if(yr_tot > cap){
              new_catch_data$Catch[yr_ind] <- cap * new_catch_data$Catch[yr_ind] / yr_tot
            }
          }
        } else { # Species-specific
          new_catch_data$Catch[dat_fill_ind] <- ifelse(new_catch_data$Catch[dat_fill_ind] > cap[new_catch_data$Species[dat_fill_ind]], cap[new_catch_data$Species[dat_fill_ind]], new_catch_data$Catch[dat_fill_ind])
        }
      }

      # * Exploitable biomass limit ----
      # - If projected catch > exploitable biomass in OM, reduce to exploitable biomass
      exploitable_biomass_data <- om_use$data_list$catch_data
      exploitable_biomass_data$Catch[dat_fill_ind] <- om_use$quantities$max_catch_hat[dat_fill_ind]

      new_catch_data$Catch[dat_fill_ind] <- ifelse(new_catch_data$Catch[dat_fill_ind] > exploitable_biomass_data$Catch[dat_fill_ind],
                                                   exploitable_biomass_data$Catch[dat_fill_ind],
                                                   new_catch_data$Catch[dat_fill_ind])

      new_catch_switch <- sum(new_catch_data$Catch[dat_fill_ind]) #Switch to turn off re-running OM if new catch = 0

      # - Update catch data in OM and EM
      om_use$data_list$catch_data <- new_catch_data
      em_use$data_list$catch_data <- new_catch_data


      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 2. Update the OM ----
      # - Estimate Fdev and update dynamics
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # - Update endyr of OM
      nyrs_hind <- om_use$data_list$endyr - om_use$data_list$styr + 1
      om_use$data_list$endyr <- assess_yrs[k]

      # The operating model is refit over its WHOLE projection every
      # assessment, not just as far as the next one. `sim_mod()` draws once per
      # observation row, so a shorter horizon would carry fewer rows, advance
      # the random stream less far, and make the draws depend on when the next
      # assessment falls -- which is the very thing a comparison of two
      # assessment schedules is trying to hold fixed.
      # inst/dev/TODO-mse-horizon.md has the design for recovering the runtime
      # saving without that dependence.

      # * Update parameters ----
      # -- log_F
      om_use$estimated_params$log_F <- cbind(om_use$estimated_params$log_F, matrix(0, nrow= nrow(om_use$estimated_params$log_F), ncol = length(new_years)))

      # -- M1_dev
      #FIXME - simulate
      # om_use$estimated_params$log_M1_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- om_use$estimated_params$log_M1_dev[,,,nyrs_hind]

      # -- Time-varying survey catchability - assume last year, filled by columns
      om_use$estimated_params$index_q_dev <- cbind(om_use$estimated_params$index_q_dev, matrix(om_use$estimated_params$index_q_dev[,ncol(om_use$estimated_params$index_q_dev)], nrow= nrow(om_use$estimated_params$index_q_dev), ncol = length(new_years)))

      # -- Time-varying selectivity - assume last year, filled by columns.
      # Use the fitted arrays' own sex extent (max_sex), not a hardcoded 2: a
      # single-sex model has sex-dim 1, and forcing 2 silently recycles the
      # fitted values into a phantom second sex (and, since build_params sizes
      # these by max_sex, mismatches the parameter template on the refit).
      n_sex_om <- dim(om_use$estimated_params$log_sel_slp_dev)[3]
      log_sel_slp_dev = array(0, dim = c(2, nflts, n_sex_om, nyrs_hind + length(new_years)))  # selectivity deviation parameters for logistic
      sel_inf_dev = array(0, dim = c(2, nflts, n_sex_om, nyrs_hind + length(new_years)))  # selectivity deviation parameters for logistic
      sel_coff_dev = array(0, dim = c(nflts, n_sex_om, n_sel_bins_om, nyrs_hind + length(new_years)))  # selectivity deviation parameters for non-parametric

      log_sel_slp_dev[,,,1:nyrs_hind] <- om_use$estimated_params$log_sel_slp_dev
      sel_inf_dev[,,,1:nyrs_hind] <- om_use$estimated_params$sel_inf_dev
      sel_coff_dev[,,,1:nyrs_hind] <- om_use$estimated_params$sel_coff_dev

      log_sel_slp_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- log_sel_slp_dev[,,,nyrs_hind]
      sel_inf_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_inf_dev[,,,nyrs_hind]
      sel_coff_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_coff_dev[,,,nyrs_hind]

      om_use$estimated_params$log_sel_slp_dev <- log_sel_slp_dev
      om_use$estimated_params$sel_inf_dev <- sel_inf_dev
      om_use$estimated_params$sel_coff_dev <- sel_coff_dev


      # * Update map ----
      # -(Only new parameter we are estimating in OM is the log_F of the new years)
      om_use$map <- build_map(
        data_list = om_use$data_list,
        params = om_use$estimated_params,
        debug = TRUE,
        random_rec = om_use$data_list$random_rec)
      om_use$map$mapFactor$dummy <- as.factor(NA); om_use$map$mapList$dummy <- NA


      # -- Estimate terminal F for catch
      new_f_yrs <- (ncol(om_use$map$mapList$log_F) - length(new_years) + 1) : ncol(om_use$map$mapList$log_F) # - Years of new F
      f_fleets <- om_use$data_list$fleet_control$Fleet_code[which(om_use$data_list$fleet_control$Fleet_type == "Fishery")] # Fleet rows for F
      om_use$map$mapList$log_F[f_fleets,new_f_yrs] <- replace(om_use$map$mapList$log_F[f_fleets,new_f_yrs], values = 1:length(om_use$map$mapList$log_F[f_fleets,new_f_yrs]))

      # -- Map out Fdev for years with 0 catch to very low number
      zero_catch <- om_use$data_list$catch_data |>
        dplyr::filter(Year <= om_use$data_list$endyr &
                        Catch == 0) |>
        dplyr::mutate(Year = Year - om_use$data_list$styr + 1) |>
        dplyr::select(Fleet_code, Year) |>
        as.matrix()
      om_use$estimated_params$log_F[zero_catch] <- -999
      om_use$map$mapList$log_F[zero_catch] <- NA
      om_use$map$mapFactor$log_F <- factor(om_use$map$mapList$log_F)
      rm(zero_catch)

      # -- Set estimate mode
      estimate_mode_base <- om_use$data_list$estimateMode
      estimate_mode_use <- ifelse(
        new_catch_switch == 0, 3, # Run in debug mode if catch is 0 for all species
        ifelse(
          estimate_mode_base < 3, 1, # Estimate hindcast only if estimating
          estimate_mode_base)
      )

      if(new_catch_switch == 0){
        om_use$map = NULL
      }

      # * Fit OM with new catch data ----
      kill_sim <- tryCatch({
        R.utils::withTimeout({
          suppressMessages(
            # Advance the OM with the new catch. The stock-recruit reference
            # period (srr_*) and the suitability window (suit_*) are pinned to
            # the PRISTINE om$ (not the advancing om_use$) so both stay fixed
            # through the projection -- critical for multispecies models, whose
            # predation suitability must not drift over the MSE. The OM projects
            # on mean recruitment across sim iterations (proj_mean_rec = TRUE).
            om_use <- .refit_like(
              data_list        = om_use$data_list,
              inits            = om_use$estimated_params,
              map              = om_use$map,
              estimateMode     = estimate_mode_use,
              loopnum          = loopnum,
              proj_mean_rec    = TRUE,
              srr_mse_switchyr = om$data_list$srr_mse_switchyr,
              srr_hat_styr     = om$data_list$srr_hat_styr,
              srr_hat_endyr    = om$data_list$srr_hat_endyr,
              suit_styr        = om$data_list$suit_styr,
              suit_endyr       = om$data_list$suit_endyr)
          )
          return(list(kill_sim = FALSE, failure = NA))
        },
        timeout = 60*timeout)
      },
      error = function(e){
        return(list(kill_sim = TRUE, failure = "OM"))
      },
      TimeoutException = function(e){
        return(list(kill_sim = TRUE, failure = "OM"))
      })

      if(kill_sim$kill_sim){
        # Nothing to put back: a killed simulation returns only its failure
        # marker, so the partially advanced operating model is discarded.
        break()
      }

      # -- Set estimate mode back to original
      om_use$data_list$estimateMode <- estimate_mode_base

      # The operating model AS FITTED, kept because section 3 below overwrites
      # om_use's catch_data with the realized catch. sim_mod() reads the
      # quantities positionally against the data frames beside them, so it must
      # see the tables the fit was taped on, not those tables edited afterwards.
      om_fit <- om_use


      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 3. Get actual catch from OM ----
      # - Maybe the OM can't support the TAC
      # (should be minimal with the < exploitable biomass check)
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

      # - Get realized catch data from OM
      new_catch_data <- om_use$data_list$catch_data
      dat_fill_ind <- which(new_catch_data$Year %in% new_years)
      new_catch_data$Catch[dat_fill_ind] <- om_use$quantities$catch_hat[dat_fill_ind] # Catch from OM

      # - Update catch data in OM and EM
      om_use$data_list$catch_data <- new_catch_data
      em_use$data_list$catch_data <- new_catch_data


      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 4. Simulate data from OM ----
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # - Simulate new survey and comp data
      # Seeded on this assessment's year, not on wherever the stream happened
      # to be left by the assessments before it.
      set.seed(draw_seeds[[as.character(assess_yrs[k])]])
      sim_dat <- Rceattle::sim_mod(om_fit, simulate = simulate_data)

      years_include <- sample_yrs[which(sample_yrs$Year > em_use$data_list$endyr & sample_yrs$Year <= assess_yrs[k]),]

      # sample_yrs pairs each fleet with the years THAT fleet is sampled, so the
      # rows to add are the (fleet, year) pairs it lists -- not every fleet in
      # every year. Matching the year set and the fleet set separately is the
      # same thing only when one year advances per assessment; with a longer
      # assessment period it also admits an every-other-year fleet in its off
      # years, giving the estimation model survey and composition data the
      # sampling design says were never collected.
      sampled_key <- paste(years_include$Fleet_code, years_include$Year)

      # -- Add newly simulated survey data to EM and OM
      # - Get simulated survey data
      new_index_data <- sim_dat$index_data |>
        dplyr::filter(paste(Fleet_code, abs(Year)) %in% sampled_key) |>
        dplyr::mutate(Year = -Year)

      # - Add to EM and OM
      om_use$data_list$index_data <- om_use$data_list$index_data |>
        dplyr::filter(!(paste(Fleet_code, abs(Year)) %in% sampled_key)) |>
        rbind(new_index_data |>
                dplyr::mutate(Year = -abs(Year))) |>
        dplyr::arrange(Fleet_code, abs(Year))

      em_use$data_list$index_data <- em_use$data_list$index_data |>
        rbind(new_index_data) |>
        dplyr::arrange(Fleet_code, abs(Year))


      # -- Add newly simulated comp data to EM & OM
      # - Simulated comp data
      new_comp_data <- sim_dat$comp_data |>
        dplyr::filter(paste(Fleet_code, abs(Year)) %in% sampled_key) |>
        dplyr::mutate(Year = -Year)

      new_comp_data$Sample_size <- new_comp_data$Sample_size * as.numeric(rowSums(dplyr::select(new_comp_data, dplyr::contains("Comp_"))) > 0) # Set sample size to 0 if catch is 0
      new_comp_data <- new_comp_data |>
        dplyr::mutate_at(dplyr::vars(dplyr::contains("Comp_")), ~ .x + 1 * (Sample_size == 0)) # Set all values to 1 if catch is 0

      # - Add to EM and OM
      om_use$data_list$comp_data <- om_use$data_list$comp_data |>
        dplyr::filter(!(paste(Fleet_code, abs(Year)) %in% sampled_key)) |>
        rbind(new_comp_data |>
                dplyr::mutate(Year = -abs(Year))) |>
        dplyr::arrange(Fleet_code, abs(Year))

      em_use$data_list$comp_data <- em_use$data_list$comp_data |>
        rbind(new_comp_data) |>
        dplyr::arrange(Fleet_code, abs(Year))

      # -- Add newly simulated CAAL to EM & OM
      # - Simulated caal data
      new_caal_data <- sim_dat$caal_data |>
        dplyr::filter(paste(Fleet_code, abs(Year)) %in% sampled_key) |>
        dplyr::mutate(Year = -Year)

      new_caal_data$Sample_size <- new_caal_data$Sample_size * as.numeric(rowSums(dplyr::select(new_caal_data, dplyr::contains("CAAL_"))) > 0) # Set sample size to 0 if catch is 0
      new_caal_data <- new_caal_data |>
        dplyr::mutate_at(dplyr::vars(dplyr::contains("CAAL_")), ~ .x + 1 * (Sample_size == 0)) # Set all values to 1 if catch is 0

      # - Add to EM and OM
      om_use$data_list$caal_data <- om_use$data_list$caal_data |>
        dplyr::filter(!(paste(Fleet_code, abs(Year)) %in% sampled_key)) |>
        rbind(new_caal_data |>
                dplyr::mutate(Year = -abs(Year))) |>
        dplyr::arrange(Fleet_code, abs(Year))

      em_use$data_list$caal_data <- em_use$data_list$caal_data |>
        rbind(new_caal_data) |>
        dplyr::arrange(Fleet_code, abs(Year))

      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 5. Update EM and HCR ----
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # Update end year and re-estimate
      em_use$data_list$endyr <- assess_yrs[k]

      # Update parameter size and use previous estimates
      # -- log_F
      em_use$estimated_params$log_F <- cbind(em_use$estimated_params$log_F, matrix(0, nrow= nrow(em_use$estimated_params$log_F), ncol = length(new_years)))

      # # -- log_M1_dev
      # em_use$estimated_params$log_M1_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- em_use$estimated_params$log_M1_dev[,,,nyrs_hind]

      # -- Time-varying survey catchability - assume last year, filled by columns
      em_use$estimated_params$index_q_dev <- cbind(em_use$estimated_params$index_q_dev, matrix(em_use$estimated_params$index_q_dev[,ncol(em_use$estimated_params$index_q_dev)], nrow= nrow(em_use$estimated_params$index_q_dev), ncol = length(new_years)))

      # -- Time-varying selectivity - assume last year, filled by columns.
      # Sex extent from the fitted arrays (max_sex), not a hardcoded 2 -- see
      # the OM extension above.
      n_sex_em <- dim(em_use$estimated_params$log_sel_slp_dev)[3]
      log_sel_slp_dev = array(0, dim = c(2, nflts, n_sex_em, nyrs_hind + length(new_years)))  # selectivity deviation parameters for logistic
      sel_inf_dev = array(0, dim = c(2, nflts, n_sex_em, nyrs_hind + length(new_years)))  # selectivity deviation parameters for logistic
      sel_coff_dev = array(0, dim = c(nflts, n_sex_em, n_sel_bins_em, nyrs_hind + length(new_years)))  # selectivity deviation parameters for non-parametric

      log_sel_slp_dev[,,,1:nyrs_hind] <- em_use$estimated_params$log_sel_slp_dev
      sel_inf_dev[,,,1:nyrs_hind] <- em_use$estimated_params$sel_inf_dev
      sel_coff_dev[,,,1:nyrs_hind] <- em_use$estimated_params$sel_coff_dev

      # - Initialize new years with last year
      log_sel_slp_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- log_sel_slp_dev[,,,nyrs_hind]
      sel_inf_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_inf_dev[,,,nyrs_hind]
      sel_coff_dev[,,,(nyrs_hind + 1):(nyrs_hind + length(new_years))] <- sel_coff_dev[,,,nyrs_hind]

      em_use$estimated_params$log_sel_slp_dev <- log_sel_slp_dev
      em_use$estimated_params$sel_inf_dev <- sel_inf_dev
      em_use$estimated_params$sel_coff_dev <- sel_coff_dev


      # Re-estimate
      kill_sim <- tryCatch({
        R.utils::withTimeout({
          suppressMessages(
            # Re-assess with the newly sampled data. srr_mse_switchyr switches
            # the stock-recruit relationship at the EM's CURRENT assessment
            # endyr (which advances each iteration), not the stored value. The
            # HCR is reconstructed from em_use$, whose HCRorder equals the
            # pristine em$ (copied at `em_use <- em`, never modified in the loop).
            em_use <- .refit_like(
              data_list        = em_use$data_list,
              inits            = em_use$estimated_params,
              estimateMode     = ifelse(em_use$data_list$estimateMode < 3, 0, em_use$data_list$estimateMode),
              loopnum          = loopnum,
              srr_mse_switchyr = em_use$data_list$endyr)
          )
          return(list(kill_sim = FALSE, failure = NA))
        },
        timeout = 60*timeout)
      },
      error = function(ex) {
        return(list(kill_sim = TRUE, failure = "EM"))
      },
      TimeoutException = function(ex) {
        return(list(kill_sim = TRUE, failure = "EM"))
      })

      # A failed re-assessment leaves em_use holding the PREVIOUS year's model,
      # so the simulation cannot carry on: it would store that stale assessment
      # under this year's name and hand its catch advice to the next iteration.
      # Stop the simulation the same way a failed operating-model refit does.
      if(kill_sim$kill_sim){
        break()
      }
      # plot_biomass(list(em_use, om_use), model_names = c("EM", "OM"))
      # End year of assessment

      # - Remove unneeded bits for memory
      em_use$initial_params <- NULL
      em_use$bounds <- NULL
      em_use$map <- NULL
      em_use$phase_params <- NULL
      em_use$obj <- NULL
      em_use$opt <- NULL
      em_use$sdrep <- NULL
      em_use$quantities[names(em_use$quantities) %!in% .mse_keep_quantities] <- NULL

      sim_list$EM[[k+1]] <- em_use
      message(paste0("Sim ",sim, " - EM Year ", assess_yrs[k], " COMPLETE"))
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # 6. End year loop ----
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    }


    # - Rename models
    if (kill_sim$kill_sim) {
      # A simulation that broke off part way through never reached the terminal
      # assessment year, so its operating and estimation models describe a state
      # the MSE did not arrive at. Return only the marker: there is nothing here
      # a summary should average over, and keeping the partial models invites
      # them being read as if the simulation had run to completion.
      sim_list <- list(use_sim = FALSE, failure = kill_sim$failure)
    } else {
      sim_list$use_sim <- TRUE
      sim_list$failure <- NA
      sim_list$OM <- om_use # OM
      # The unfished reference run is a separate numerical problem and can fail
      # on its own (it sdreports a model with F pinned off across the
      # projection). The simulation itself is complete at this point, and its
      # assessments are the catch-advice record, so a failure here must not
      # discard it: dropping these would remove simulations in a
      # stock-state-dependent way and bias every performance metric.
      # Assigned through `[` rather than `$` so a failure leaves OM_no_F present
      # and NULL: `sim_list$OM_no_F <- NULL` would DELETE the element, and the
      # simulation would then be indistinguishable by name from one that never
      # attempted the unfished run.
      # No fishing after the original OM's terminal year; the advanced om_use ends at the last assessment.
      sim_list["OM_no_F"] <- list(tryCatch(remove_F(om_use, styr = om$data_list$endyr + 1), error = function(e) {
        # Recorded on the object, not just warned about: this runs in a parallel
        # worker, whose warnings are discarded, and the simulation is still
        # usable for everything that does not compare against the unfished run.
        # mse_summary() drops these from the no-F metrics only.
        sim_list$failure <<- paste0("OM_no_F: ", conditionMessage(e))
        NULL
      }))
      names(sim_list$EM) <- c("EM", paste0("OM_Sim_",sim,". EM_yr_", assess_yrs))
    }

    # - Save
    if(!is.null(dir)){
      dir.create(dir, showWarnings = FALSE, recursive = TRUE)
      saveRDS(sim_list, file = paste0(dir, "/", file, "EMs_from_OM_Sim_",sim, ".rds"))
      # Hand back only the outcome, not the models: run_mse() returns NULL in
      # this mode, and the dispatcher needs to report attrition without reading
      # every saved simulation back off disk.
      list(use_sim = sim_list$use_sim, failure = sim_list$failure)
    } else{
      sim_list # Return simlist
    }
  } # End run_one_sim closure


  # Contain a simulation's failure to that simulation. The refits are already
  # guarded, but the surrounding data reshaping is not, and neither the
  # sequential nor the parallel dispatch has a per-item handler -- so one
  # unguarded error would discard every other simulation of the run. Anything
  # that escapes is reported as a failed simulation, keeping the message so the
  # cause is still visible.
  run_one_sim_guarded <- function(sim) {
    tryCatch(run_one_sim(sim), error = function(e) {
      # Tagged so this stays distinguishable from the enumerated "OM"/"EM"
      # failures the assessment loop reports: those are expected outcomes, this
      # is a bug or a bad input.
      marker <- list(use_sim = FALSE,
                     failure = paste0("unexpected: ", conditionMessage(e)))
      message("Sim ", sim, " FAILED: ", conditionMessage(e))
      if (!is.null(dir)) {
        # Writing the marker must not itself escape, or the containment this
        # handler exists for is lost for every other simulation.
        tryCatch({
          dir.create(dir, showWarnings = FALSE, recursive = TRUE)
          saveRDS(marker,
                  file = paste0(dir, "/", file, "EMs_from_OM_Sim_", sim, ".rds"))
        }, error = function(e2) {
          warning("Sim ", sim, ": could not record the failure (",
                  conditionMessage(e2), ").", call. = FALSE)
        })
      }
      marker
    })
  }

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Dispatch sims (parallel via PSOCK or sequential) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if (use_parallel) {
    # FORK where possible (inherits the loaded package + the large OM/EM objects
    # via copy-on-write); PSOCK fallback exports them. See .parallel_lapply().
    sim_list <- .parallel_lapply(start_sim:nsim, run_one_sim_guarded, min(cores, nsim), environment())
  } else {
    sim_list <- lapply(start_sim:nsim, run_one_sim_guarded)
  }

  names(sim_list) <- paste0("Sim_", start_sim:nsim)

  # Report attrition here rather than leaving it to whatever reads the results.
  # With `dir` set the return value is NULL, so a run in which every simulation
  # failed would otherwise be indistinguishable from one in which every
  # simulation succeeded -- the same number of files either way.
  n_failed <- sum(vapply(sim_list, function(x) isFALSE(x$use_sim), logical(1)))
  if (n_failed) {
    warning(n_failed, " of ", length(sim_list),
            " simulations did not complete; they carry use_sim = FALSE and no ",
            "models. Filter on use_sim before summarising.", call. = FALSE)
  }
  # Worker warnings are discarded by the parallel cluster, so a completed
  # simulation whose unfished reference run failed is reported from here.
  n_no_f <- sum(vapply(sim_list, function(x) {
    isTRUE(x$use_sim) && grepl("^OM_no_F: ", paste(x$failure))
  }, logical(1)))
  if (n_no_f) {
    warning(n_no_f, " of ", length(sim_list),
            " simulations completed but their unfished reference run failed; ",
            "these carry use_sim = TRUE and OM_no_F = NULL, and are excluded ",
            "from the no-F metrics only.", call. = FALSE)
  }

  if(is.null(dir)){
    return(sim_list)
  }
}
