#' Plot composition fits and residuals
#'
#' @description
#' Diagnostic plots for age / length composition data from a fitted [Rceattle]
#' model, drawn with ggplot2 for a consistent look with [plot.rceattle_osa()].
#' With the default `residual_type = "pearson"` three figures are produced per
#' fleet x composition type:
#' \itemize{
#'   \item **Pearson residual bubbles** by year and bin, faceted by fleet (and,
#'     for joint-sex data, by sex); red = positive, blue = negative, sized by
#'     magnitude. The Pearson residual is
#'     the same form used by [residuals.Rceattle()], standardized by the
#'     variance the fleet's own likelihood assumes, so the weight column and any
#'     Dirichlet-multinomial overdispersion are already in the denominator.
#'   \item **Annual composition**: observed (shaded area) vs fitted (line)
#'     proportion at age / length, one panel per year. Joint-sex data are
#'     mirrored (females up, males down).
#'   \item **Aggregated composition**: observed vs fitted, counts pooled over
#'     the fitted hindcast years. The interval holds 95% of the data the model
#'     predicts (not a confidence interval on the mean), so observations outside
#'     it indicate misfit. Poor where the expected count is below about 10, and
#'     slightly narrow throughout, because it treats the fitted proportions as
#'     known rather than estimated from these same data.
#' }
#' The shaded area and fitted line span only the observed bins (they do not
#' extend past the first/last bin), and bins with zero observed proportion are
#' retained (only `NA` bins are dropped), so the curves are not interpolated
#' across gaps.
#'
#' @param Rceattle A single fitted `Rceattle` model.
#' @param file Optional filename stem; when supplied, each figure is also saved
#'   as a PNG.
#' @param model_names Unused (kept for back-compatibility).
#' @param species Optional species code(s) to plot (matched against the
#'   `Species` column); `NULL` (default) plots all. Mirrors the `species`
#'   argument of [residuals.Rceattle()] and [plot.rceattle_osa()].
#' @param cex Unused (kept for back-compatibility).
#' @param lwd Width of the fitted-composition line. Default `3`.
#' @param right_adj Unused (kept for back-compatibility).
#' @param residual_type `"pearson"` (default) for the ggplot2 Pearson-residual
#'   and composition-fit figures drawn here, or `"osa"` to instead draw the
#'   one-step-ahead residual diagnostics via [osa_residuals()] and
#'   [plot.rceattle_osa()], a Q-Q plot (with SDNR / tail annotation) alongside
#'   signed OSA- and Pearson-residual bubbles. The `"osa"` path builds its
#'   observation data on demand, so it works with any fit.
#' @param add_agg_ci Logical. Draw the 95% prediction interval on the aggregated
#'   composition figure. Default `TRUE`.
#' @param add_agg_n Logical. Annotate the aggregated composition figure with the
#'   input sample size, the effective sample size the likelihood assumed, and
#'   the McAllister-Ianelli effective sample size this fit implies. Default
#'   `TRUE`.
#'
#' @return Invisibly, a named list of the `ggplot` objects. Called for its side
#'   effect of drawing (and optionally saving) the figures.
#' @importFrom rlang .data
#' @export
plot_comp <- function(Rceattle, file = NULL, model_names = NULL, species = NULL,
                      cex = 3, lwd = 3, right_adj = 0,
                      residual_type = c("pearson", "osa"),
                      add_agg_ci = TRUE, add_agg_n = TRUE) {

  if (!inherits(Rceattle, "Rceattle")) {
    stop("Please only use one Rceattle model")
  }
  residual_type <- match.arg(residual_type)

  # One-step-ahead residual diagnostics: delegate to the rceattle_osa plot.
  if (residual_type == "osa") {
    osa <- osa_residuals(Rceattle, source = "comp")
    return(invisible(plot(osa, species = species)))
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to plot composition data.")
  }

  long <- .comp_resid_long(Rceattle, species = species)
  if (is.null(long) || nrow(long) == 0) {
    warning("No composition data to plot for the requested species.")
    return(invisible(NULL))
  }

  sex_cols  <- c(combined = "black", female = "#d7301f", male = "#2c7fb8")
  sign_cols <- c(positive = "#d7301f", negative = "#2c7fb8")
  save_png  <- function(p, tag, w = 9, h = 6.5) {
    if (!is.null(file)) {
      ggplot2::ggsave(paste0(file, "_comp_", tag, ".png"), p,
                      width = w, height = h, units = "in", dpi = 300)
    }
  }
  plots <- list()

  # ---- Pearson residual bubbles (faceted by fleet x type [x sex]) ----
  pear <- long[is.finite(long$pearson), , drop = FALSE]
  if (nrow(pear) > 0) {
    # Same fixed size scale as plot.rceattle_osa()'s bubbles: the same residuals
    # appear in both figures, so a free scale would draw them at different sizes.
    pear$pearson <- .rce_truncate_resid(pear$pearson,
                                        "Composition Pearson residuals")
    pear$sign <- ifelse(pear$pearson >= 0, "positive", "negative")
    g <- ggplot2::ggplot(pear, ggplot2::aes(.data$Year, .data$bin)) +
      ggplot2::geom_point(ggplot2::aes(size = abs(.data$pearson),
                                       colour = .data$sign), alpha = 0.8) +
      ggplot2::scale_colour_manual(values = sign_cols, guide = "none") +
      ggplot2::scale_size_continuous(
        breaks = seq(0, .RCE_BUBBLE_MAX, length.out = 4L),
        limits = c(0, .RCE_BUBBLE_MAX),
        range = c(0.1, 3), name = "|Pearson|") +
      ggplot2::facet_wrap(~ source, scales = "free_y") +
      ggplot2::labs(x = "Year", y = "Age / length bin",
                    title = "Composition Pearson residuals") +
      ggplot2::theme_bw(base_size = 10)
    plots$pearson <- g; print(g); save_png(g, "pearson")
  }

  # ---- Composition fit: observed (area) vs fitted (line) ----
  # Mirror joint-sex (females up, males down). geom_area spans only the observed
  # bins, so it cannot extend past the data.
  long$y_obs <- ifelse(long$sex_grp == "male", -long$obs, long$obs)
  long$y_hat <- ifelse(long$sex_grp == "male", -long$hat, long$hat)

  # Observed proportion as a shaded area, fitted as a line; joint-sex males are
  # mirrored below zero (y_obs / y_hat). geom_area spans only the observed bins.
  comp_area_plot <- function(df, ylab, title) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data$bin)) +
      ggplot2::geom_area(ggplot2::aes(y = .data$y_obs, group = .data$sex_grp),
                         position = "identity", fill = "grey80") +
      ggplot2::geom_line(ggplot2::aes(y = .data$y_hat, group = .data$sex_grp,
                                      colour = .data$sex_grp),
                         linewidth = lwd / 3) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
      ggplot2::scale_colour_manual(values = sex_cols, guide = "none") +
      ggplot2::labs(x = df$bin_lab[1], y = ylab, title = title) +
      ggplot2::theme_bw(base_size = 9)
  }

  for (nm in unique(long$panel)) {
    d <- long[long$panel == nm, , drop = FALSE]

    # Annual (one panel per year)
    g_a <- comp_area_plot(d, "Proportion", nm) + ggplot2::facet_wrap(~ Year)
    plots[[paste0("annual_", nm)]] <- g_a; print(g_a)
    save_png(g_a, paste0("annual_fleet", d$Fleet[1], "_", d$type_lab[1]))

    # Aggregated across hindcast years, pooling counts
    agg <- .comp_aggregate(d, endyr = Rceattle$data_list$endyr)
    if (nrow(agg) > 0) {
      g_g <- comp_area_plot(agg, "Pooled proportion", paste(nm, "(aggregated)"))
      if (add_agg_ci && any(is.finite(agg$y_lwr))) {
        g_g <- g_g + ggplot2::geom_linerange(
          data = agg, ggplot2::aes(x = .data$bin, ymin = .data$y_lwr,
                                   ymax = .data$y_upr, colour = .data$sex_grp),
          alpha = 0.55, linewidth = 0.4, show.legend = FALSE)
      }
      if (add_agg_n) {
        g_g <- g_g + ggplot2::annotate(
          "text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3, size = 2.8,
          label = .comp_n_label(Rceattle, d$Fleet[1], agg$ISS[1], agg$ESS[1]))
      }
      plots[[paste0("aggregated_", nm)]] <- g_g; print(g_g)
      save_png(g_g, paste0("aggregated_fleet", d$Fleet[1], "_", d$type_lab[1]))
    }
  }

  invisible(plots)
}


#' The effective sample size the composition likelihood assumed
#'
#' Recovered from the same assumed sd the prediction band is built on, rather
#' than from the weight column: `sd^2 = p (1 - p) / N_eff` holds for every
#' family, so `p (1 - p) / sd^2` returns `N_eff` whether the weight entered as a
#' multinomial multiplier on the input N or as a Dirichlet-multinomial
#' concentration. That keeps one number behind the band and the annotation, and
#' needs no family branching here.
#'
#' `N_eff` belongs to the OBSERVATION, so the rows are grouped on the key the
#' likelihood scores one density over, fleet, species, sex and year, the same
#' key `data_check()` requires to be unique. Year alone is not that key: a fleet
#' may hold a female-only and a male-only row in one year, which are two
#' observations, and pooling them would return half the effective sample size
#' beside an `ISS` that correctly counts both. Joint-sex rows (`Sex == 3`) stay
#' one observation, because one multinomial spans both sexes.
#'
#' Within an observation `N_eff` is constant across bins on the offset-inflated
#' scale the `sd` was built on, but this reads the RAW fitted proportion, where
#' `comp_offset` dominates the smallest bins, at `hat = 1e-6` the recovery is
#' an order of magnitude out. The largest fitted bin is therefore load-bearing,
#' not merely tidier; there the error is below ~1e-3 relative, growing with the
#' number of bins in the observation.
#'
#' Summed over observations, because the pooled count's variances add.
#'
#' @param d One panel's rows, already filtered to the fitted hindcast years.
#' @return The summed effective sample size, or `NA` if any observation's
#'   assumed sd is unusable, a partial sum would read as a downweighted fleet.
#' @noRd
.comp_assumed_ess <- function(d) {
  if (is.null(d$Sd)) return(NA_real_)
  obs <- split(seq_len(nrow(d)), paste(d$Fleet, d$Species, d$Sex, d$Year))
  n_eff <- vapply(obs, function(i) {
    # An infinite sd is Comp_weights = 0, i.e. the fleet's composition dropped
    # from the likelihood: an effective sample size of zero, not an unknown one.
    ok <- !is.na(d$Sd[i]) & d$Sd[i] > 0
    if (!any(ok)) return(NA_real_)
    # Largest fitted bin among the usable ones, so one unusable bin costs the
    # observation precision rather than the whole panel its annotation.
    i <- i[ok]
    j <- i[which.max(d$hat[i])]
    d$hat[j] * (1 - d$hat[j]) / d$Sd[j]^2
  }, numeric(1))
  if (anyNA(n_eff)) NA_real_ else sum(n_eff)
}


#' Input and effective sample size label for the aggregated composition
#'
#' Three numbers an analyst weighing a reweighting decision needs side by side:
#' the input sample size summed over the fitted hindcast years, the effective
#' sample size the likelihood actually assumed (`Comp_weights` and any
#' Dirichlet-multinomial overdispersion already in it), and the effective sample
#' size this fit's own residuals imply, the McAllister-Ianelli tuning target
#' `fit_mod()` computes as a harmonic mean across years, which is the unbiased
#' scale to average a ratio estimator on. The gap between the second and third
#' is the reweighting decision, so both are named rather than either standing in
#' for "ESS" alone.
#'
#' The McAllister-Ianelli line is drawn only on a multinomial fleet. A
#' Dirichlet-multinomial estimates its own weight inside the likelihood, so
#' [reweight_comps()] names and skips those fleets; printing an external tuning
#' target beside one would invite the adjustment the package refuses to make.
#'
#' @param Rceattle The fitted model.
#' @param fleet The fleet code.
#' @param iss The summed input sample size.
#' @param ess The summed effective sample size the likelihood assumed.
#' @return A one- to three-line label.
#' @noRd
.comp_n_label <- function(Rceattle, fleet, iss, ess = NA_real_) {
  n_fmt <- function(x) format(round(x), big.mark = ",")
  lab <- sprintf("ISS = %s", n_fmt(iss))
  if (is.finite(ess)) {
    lab <- sprintf("%s\nESS = %s (likelihood)", lab, n_fmt(ess))
  }
  fc  <- Rceattle$data_list$fleet_control
  # The family only, not through .rce_comp_family(): that resolver warns about a
  # missing Dirichlet-multinomial weight, which residuals() has already said
  # once by the time a label is being drawn.
  fam <- .rce_family_code(.rce_switch_column(fc, "Comp_distribution"),
                          "Comp_distribution")
  i   <- match(fleet, fc$Fleet_code)
  is_dm <- !is.na(i) && length(fam) >= i && isTRUE(fam[i] == 1L)

  w <- fc$Comp_weights_mcallister
  if (!is.null(w) && !is_dm) {
    if (!is.na(i) && is.finite(w[i])) {
      lab <- sprintf("%s\nESS = %s (McAllister-Ianelli)", lab,
                     n_fmt(w[i] * iss))
    }
  }
  lab
}


#' Tidy long-format composition observed / fitted / Pearson table
#'
#' Reuses [residuals.Rceattle()] (`type = "pearson"`, `source = "comp"`) for the
#' observed and fitted proportions and the Pearson residual, the single source
#' of truth, then adds the plotting-only columns: joint-sex bins re-based onto
#' a single age/length axis (males stored as bins `nbin+1 .. 2*nbin` are mapped
#' to `1 .. nbin` and tagged `sex_grp = "male"`) and the facet labels. Zero
#' observed proportions are kept; only `NA` observed/fitted are dropped (by
#' `residuals()`).
#' @param Rceattle A fitted `Rceattle` model.
#' @param species Optional species code(s) to keep.
#' @return A data frame, or `NULL` when there are no composition residuals.
#' @keywords internal
.comp_resid_long <- function(Rceattle, species = NULL) {
  r <- stats::residuals(Rceattle, type = "pearson", source = "comp",
                        species = species)
  if (is.null(r) || nrow(r) == 0) return(NULL)
  nages    <- Rceattle$data_list$nages
  nlengths <- Rceattle$data_list$nlengths

  long <- data.frame(
    Fleet      = r$Fleet_code,
    Fleet_name = r$Fleet_name,
    Species    = r$Species,
    Sex        = r$Sex,
    comp_type  = r$Age0_Length1,
    Year       = r$Year,
    N          = r$Sample_size,
    bin        = r$Bin,
    obs        = r$Observed,      # observed proportion
    hat        = r$Fitted,        # fitted proportion
    pearson    = r$Residual,      # (obs - hat) / Sd
    Sd         = r$Sd,            # sd the fleet's own likelihood assumes
    stringsAsFactors = FALSE)

  # Joint-sex (Sex == 3) stacks females in bins 1..nbin and males in
  # nbin+1..2*nbin; re-base the male bins onto the same 1..nbin axis and label
  # each row's sex so the two can be drawn together (males mirrored below 0).
  bins_per_sex <- ifelse(long$comp_type == 1, nlengths[long$Species],
                         nages[long$Species])
  is_male <- long$Sex == 3 & long$bin > bins_per_sex
  long$bin[is_male] <- long$bin[is_male] - bins_per_sex[is_male]
  long$sex_grp <- ifelse(long$Sex == 1, "female",
                  ifelse(long$Sex == 2, "male",
                  ifelse(long$Sex == 3, ifelse(is_male, "male", "female"),
                         "combined")))

  long$type_lab <- ifelse(long$comp_type == 1, "length", "age")
  long$bin_lab  <- ifelse(long$comp_type == 1, "Length bin", "Age bin")
  long$panel    <- paste0(long$Fleet_name, " - ", long$type_lab, " comp")
  long$source   <- paste0(long$Fleet_name, "\n", long$type_lab, " comp",
                          ifelse(long$Sex == 3, paste0(" - ", long$sex_grp), ""))
  long
}


#' Pool composition counts across the fitted hindcast years
#'
#' Multiplies each year's proportions by its input sample size, sums the counts
#' over years (for one fleet x type panel), and puts the total back on the
#' proportion scale; joint-sex groups keep their shared normalization (females +
#' males sum to 1). Pooling counts rather than averaging proportions is what
#' stops a year with 20 otoliths counting as much as one with 2000.
#'
#' The pooled count is a sum of independent draws, so its variance is the exact
#' sum of the per-year variances the fleet's own likelihood assumes (the `Sd`
#' column). The resulting interval treats the fitted proportions as known, so it
#' is slightly narrower than one that held the estimation uncertainty in
#' `p_hat` as well.
#'
#' @param d One panel's rows from [.comp_resid_long()].
#' @param endyr Last year the composition likelihood fits, or `NULL` to keep
#'   every hindcast row. Default `NULL`.
#' @return A data frame with `bin`, `sex_grp`, the pooled `obs` / `hat`
#'   proportions and their counts `o_n` / `e_n`, the pooled count variance
#'   `v_n`, the 95% prediction limits `lwr` / `upr`, the summed input and
#'   assumed effective sample sizes `ISS` / `ESS`, and the sex-mirrored drawing
#'   columns `y_obs`, `y_hat`, `y_lwr`, `y_upr`.
#' @keywords internal
.comp_aggregate <- function(d, endyr = NULL) {
  d <- d[d$Year > 0, , drop = FALSE]              # hindcast only
  # The C++ composition likelihood gates on Year <= endyr, so a projection row
  # carries no fit and must not enter the aggregate or its input sample size.
  if (!is.null(endyr)) d <- d[d$Year <= endyr, , drop = FALSE]
  if (nrow(d) == 0) return(d[, c("bin", "sex_grp"), drop = FALSE])

  # Pool counts, not proportions: averaging proportions weights a year with 20
  # otoliths like one with 2000. The pooled count is a sum of independent draws,
  # so its variance is exactly the sum of the per-year variances.
  d$o_n <- d$obs * d$N
  d$e_n <- d$hat * d$N
  # Var of the pooled count in a bin = sum_y (N_y * sd_yb)^2, with sd_yb the sd
  # the fleet's own likelihood assumes for that year's proportion.
  d$v_n <- if (is.null(d$Sd)) NA_real_ else (d$N * d$Sd)^2

  # na.action = na.pass: the formula interface drops whole rows on any NA, so a
  # bin with no assumed sd would vanish from the composition rather than just
  # lose its interval. A variance missing any year is NA rather than a partial
  # sum, which would understate the spread and draw a band that is too narrow.
  agg <- stats::aggregate(cbind(o_n, e_n, v_n) ~ bin + sex_grp, data = d,
                          FUN = function(x) if (anyNA(x)) NA_real_ else sum(x),
                          na.action = stats::na.pass)

  # Displayed on the proportion scale, as before, so the figure keeps its axis.
  tot_o <- sum(agg$o_n); tot_e <- sum(agg$e_n)
  agg$obs <- agg$o_n / tot_o
  agg$hat <- agg$e_n / tot_e
  # 95% range of the data the model predicts, not a CI on the mean. Normal
  # approximation on an exact variance; poor below ~10 expected counts.
  agg$lwr <- (agg$e_n - 1.96 * sqrt(agg$v_n)) / tot_e
  agg$upr <- (agg$e_n + 1.96 * sqrt(agg$v_n)) / tot_e

  agg$bin_lab <- d$bin_lab[1]
  # Proportions sum to one across a row's bins (both sexes when joint), so the
  # pooled observed count is the summed input sample size.
  agg$ISS <- tot_o
  agg$ESS <- .comp_assumed_ess(d)
  mir <- ifelse(agg$sex_grp == "male", -1, 1)
  agg$y_obs <- mir * agg$obs
  agg$y_hat <- mir * agg$hat
  agg$y_lwr <- mir * agg$lwr
  agg$y_upr <- mir * agg$upr
  agg
}




#' Plot diet composition fits
#'
#' @description
#' If year == 0, diet data are averaged from suit_styr to suit_endyr
#' If prey_age >= 0 diet data are diet proportion of prey-at-age in predator-at-age
#' If prey_age < 0 diet data are diet proportion of prey-spp in predator-at-age (sum across prey ages)
#' If prey_age < 0 and pred_age < 0, diet data are mean diet proportion of prey-spp in predator-spp (sum across prey ages and take mean across predator ages)
#' If prey_age < 0 and pred_age < -500, diet data are weighted mean diet proportion of prey-spp in predator-spp (sum across prey ages and take weighted mean across predator ages)
#'
#'
#' @param Rceattle Single or list of Rceattle model objects exported from \code{Rceattle}
#' @param file name of a file to identified the files exported by the function.
#' @param species Species names for legend
#'
#' @return Returns and saves a figure
#' @importFrom rlang .data
#' @importFrom ggplot2 theme element_blank aes ggplot geom_point theme_classic ylim ylab xlim xlab ggtitle ggsave
#' @export
plot_diet_comp <-
  function(Rceattle,
           file = NULL,
           species = NULL) {

    # Make sure we are using only one model
    if(!inherits(Rceattle, "Rceattle")){
      stop("Please only use one Rceattle model")
    }
    data_list <- Rceattle$data_list

    # Species names
    if(is.null(species)){
      species =  Rceattle$data_list$spnames
    }

    # Get colors
    colvec=c("red", "blue", "black")
    plots <- list()  # collected per predator/prey/year cowplot grids

    # * Extract data objects ----
    # Diet observed proportions, fitted, and the Pearson residual all come from
    # residuals(type = "pearson", source = "diet") -- the single source of truth
    # for residual computation -- then aliased to the column names this plotting
    # loop expects.
    comp_data <- stats::residuals(Rceattle, type = "pearson", source = "diet")
    comp_data$Pred <- comp_data$Species
    comp_data$Stomach_proportion_by_weight <- comp_data$Observed
    comp_data$Est <- comp_data$Fitted
    comp_data$pearson <- comp_data$Residual

    # Loop around predators ----
    for(pred in 1:data_list$nspp) {
      for(pred_sex_idx in 1:data_list$nsex[pred]){
        for(prey in 1:data_list$nspp) {

          # * Get sex for legend ----
          if(data_list$nsex[pred] > 1){
            # Sex-specific predator label for the legend (1 = female, 2 = male)
            pred_legend <- paste("Pred-", species[pred], ifelse(pred_sex_idx == 1, "female", "male"))
            current_pred_sex = pred_sex_idx - 1
          } else{
            pred_legend <- paste("Pred-", species[pred])
            current_pred_sex = 0
          }

          # * Extract comps ----
          comp_tmp <- comp_data |>
            dplyr::filter(.data$Pred == pred & .data$Pred_sex == current_pred_sex & .data$Prey == prey)

          yrs <- sort(unique(comp_tmp$Year))
          nyrs <- length(yrs)

          # * Plot annual comps ----
          for(yr in 1:nyrs){
            comp_tmp_yr <- comp_tmp |>
              dplyr::filter(.data$Year == yrs[yr] ) |>
              dplyr::mutate(Prey_age = ifelse(.data$Prey_sex == 2, -.data$Prey_age, .data$Prey_age))

            plot_obs <- comp_tmp_yr |>
              dplyr::filter(.data$Stomach_proportion_by_weight > 0) |>
              ggplot2::ggplot(ggplot2::aes(x = .data$Pred_age, y = .data$Prey_age, size = .data$Stomach_proportion_by_weight)) +
              ggplot2::geom_point(alpha = 1) +
              ggplot2::theme_classic() +
              ggplot2::ylim(range(comp_tmp_yr$Prey_age)) +
              ggplot2::ylab(paste(species[prey], "age")) +
              ggplot2::xlim(range(comp_tmp_yr$Pred_age)) +
              ggplot2::xlab(paste(pred_legend, "age")) +
              ggplot2::ggtitle(paste("Observed diet: year", yrs[yr])) +
              ggplot2::theme(legend.position = c(0.25, 0.7),
                             legend.title = ggplot2::element_blank())

            plot_est <- comp_tmp_yr |>
              dplyr::filter(.data$Stomach_proportion_by_weight > 0) |>
              ggplot2::ggplot(ggplot2::aes(x = .data$Pred_age, y = .data$Prey_age, size = .data$Est)) +
              ggplot2::geom_point(alpha = 1) +
              ggplot2::theme_classic() +
              ggplot2::ylim(range(comp_tmp_yr$Prey_age)) +
              ggplot2::ylab(paste(species[prey], "age")) +
              ggplot2::xlim(range(comp_tmp_yr$Pred_age)) +
              ggplot2::xlab(paste(pred_legend, "age")) +
              ggplot2::ggtitle(paste("Estimated diet: year", yrs[yr])) +
              ggplot2::theme(legend.position = "none")

            plot_pear <- comp_tmp_yr |>
              dplyr::filter(.data$Stomach_proportion_by_weight > 0) |>
              ggplot2::ggplot(ggplot2::aes(x = .data$Pred_age, y = .data$Prey_age,
                                           size = abs(.data$pearson),
                                           color = .data$pearson < 0)) +
              ggplot2::geom_point(alpha = 1) +
              ggplot2::theme_classic() +
              ggplot2::ylim(range(comp_tmp_yr$Prey_age)) +
              ggplot2::ylab(paste(species[prey], "age")) +
              ggplot2::xlim(range(comp_tmp_yr$Pred_age)) +
              ggplot2::xlab(paste(pred_legend, "age")) +
              ggplot2::ggtitle(paste("Pearson residual: year", yrs[yr])) +
              ggplot2::theme(legend.position = c(0.25, 0.8))

            p1 <- cowplot::plot_grid(plot_obs, plot_est, plot_pear, nrow = 1)
            plots[[length(plots) + 1L]] <- p1
            print(p1)

            if (!is.null(file)) {
              filename <- paste0(file, "_aggregated_diet_comps_year", yr,"_", pred_legend, "_prey_", species[prey],".png")
              ggplot2::ggsave(filename = filename, plot = p1, width = 10, height = 6.5, units = "in", dpi = 300)
            }
          }
        }
      }
    }
    invisible(plots)
  }


#' Plot diet composition fits (bubble/grid diagnostic)
#'
#' @description
#' `plot_diet_comp1()` is an alias of [plot_diet_comp()]: it produces the
#' predator-age x prey-age bubble grids (observed / estimated / Pearson residual)
#' for each predator-prey-year interaction. It is provided as a named entry point
#' alongside [plot_diet_comp2()] (the aggregation-aware line / bar / bubble plots)
#' for scripts that call both by name.
#'
#' @param Rceattle Single Rceattle model object exported from \code{Rceattle}.
#' @param file Optional file-name prefix for saved figures.
#' @param species Optional species names for the legend.
#'
#' @return Invisibly returns a list of the printed plot grids (see
#'   [plot_diet_comp()]).
#' @examples
#' \dontrun{
#' plot_diet_comp1(my_msm_fit)
#' }
#' @seealso [plot_diet_comp()], [plot_diet_comp2()]
#' @export
plot_diet_comp1 <- plot_diet_comp


#' Assemble diet observed-vs-fitted proportions for plotting
#'
#' @description
#' Internal helper that builds the tidy data frame consumed by
#' [plot_diet_comp2()]. Observed proportions, fitted proportions and the Pearson
#' residual all come from `residuals(type = "pearson", source = "diet")` (the
#' single source of truth); the fitted value is `quantities$diet_hat[, 2]`,
#' matched row-for-row to `data_list$diet_data`, so any prey-age / year
#' aggregation is already resolved in the C++ model.
#'
#' Observed 95% intervals (`lower_95` / `upper_95`) are a normal approximation to
#' the binomial proportion, \eqn{\hat p \pm 1.96\sqrt{p(1-p)/N}}, using the same
#' Sample_size as the Pearson denominator (`diet_data` holds no stored CI
#' columns). Estimated 95% intervals (`Est_Lower` / `Est_Upper`) are added only
#' when the `sdrep` exposes a `diet_hat` standard error; the C++ template
#' `REPORT()`s but does not `ADREPORT()` `diet_hat`, so these are unavailable
#' unless the model is changed to `ADREPORT(diet_hat)`.
#'
#' @param Rceattle A single Rceattle model object.
#' @return A data frame with `Pred`, `Prey`, `Pred_sex`, `Prey_sex`, `Pred_age`,
#'   `Prey_age`, `Year`, `Sample_size`, `Observed`, `Est`, `Pearson`,
#'   `AbsPearson`, `lower_95`, `upper_95`, and (when available) `Est_Lower` /
#'   `Est_Upper`; or `NULL` when the model has no diet data.
#' @noRd
.diet_plot_data <- function(Rceattle) {

  dd <- Rceattle$data_list$diet_data
  if (is.null(dd) || nrow(dd) == 0 || is.null(Rceattle$quantities$diet_hat)) {
    return(NULL)
  }

  # residuals(source = "diet") is the single source of truth for observed,
  # fitted (diet_hat[, 2]) and Pearson; rows align 1:1 with diet_data.
  r <- stats::residuals(Rceattle, type = "pearson", source = "diet")

  plot_data <- data.frame(
    Pred        = r$Species,
    Prey        = r$Prey,
    Pred_sex    = r$Pred_sex,
    Prey_sex    = r$Prey_sex,
    Pred_age    = r$Pred_age,
    Prey_age    = r$Prey_age,
    Year        = r$Year,
    Sample_size = r$Sample_size,
    Observed    = r$Observed,
    Est         = r$Fitted,
    Pearson     = r$Residual,
    stringsAsFactors = FALSE
  )
  plot_data$AbsPearson <- abs(plot_data$Pearson)

  # Observed 95% CI: normal approximation to the binomial proportion, matching
  # the Sample_size used in dev's Pearson denominator (no stored CI columns).
  se_obs <- sqrt(plot_data$Observed * (1 - plot_data$Observed) / plot_data$Sample_size)
  plot_data$lower_95 <- pmax(0, plot_data$Observed - 1.96 * se_obs)
  plot_data$upper_95 <- pmin(1, plot_data$Observed + 1.96 * se_obs)

  # Estimated 95% CI from the sdreport, if diet_hat was ADREPORT'd. The model
  # only REPORT()s diet_hat, so this branch is normally skipped. Both shapes are
  # handled: a vector gives n SDs, an n x 2 matrix gives 2n, of which the second
  # column (the fitted proportion) is used.
  if (!is.null(Rceattle$sdrep)) {
    sd_indices <- which(names(Rceattle$sdrep$value) == "diet_hat")
    if (length(sd_indices) > 0) {
      all_sd_values <- Rceattle$sdrep$sd[sd_indices]
      n_rows_data   <- nrow(plot_data)
      sd_values_to_use <- NULL
      if (length(all_sd_values) == n_rows_data) {
        sd_values_to_use <- all_sd_values
      } else if (length(all_sd_values) == n_rows_data * 2) {
        sd_values_to_use <- all_sd_values[(n_rows_data + 1):(n_rows_data * 2)]
      }
      if (!is.null(sd_values_to_use)) {
        plot_data$Est_sd    <- sd_values_to_use
        plot_data$Est_Lower <- pmax(0, plot_data$Est - 1.96 * plot_data$Est_sd)
        plot_data$Est_Upper <- pmin(1, plot_data$Est + 1.96 * plot_data$Est_sd)
      }
    }
  }

  plot_data
}


#' Plot diet composition fits (aggregation-aware)
#'
#' @description
#' Diagnostic plots for diet-composition fits that adapt to how each
#' predator-prey interaction is aggregated (see [plot_diet_comp()] for the
#' aggregation conventions):
#' \itemize{
#'   \item prey-age aggregated (predator age resolved): line plot of observed vs
#'     estimated diet proportion against predator age, with 95% CI ribbons;
#'   \item predator aggregated (prey age resolved): line plot against prey age;
#'   \item both aggregated: dodged bar plot of observed vs estimated proportion;
#'   \item fully disaggregated: predator-age x prey-age bubble grids (observed /
#'     estimated / Pearson residual).
#' }
#'
#' The observed proportions, fitted proportions and Pearson residuals come from
#' the diet data path (`quantities$diet_hat` / `residuals(source = "diet")` /
#' `data_list$diet_data`) via the internal `.diet_plot_data()` helper. Observed
#' 95% CI ribbons are a normal approximation to the binomial proportion from
#' `Sample_size`. Estimated 95% CI ribbons are drawn only when the `sdreport`
#' exposes a `diet_hat` standard error; the C++ template `REPORT()`s but does
#' not `ADREPORT()` `diet_hat`, so the estimated ribbon is unavailable unless the
#' template is changed to `ADREPORT(diet_hat)` (the code path is retained so it
#' activates automatically if that ever happens, or when `sdrep` is `NULL`).
#'
#' @param Rceattle A single Rceattle model object.
#' @param file Optional file-name prefix for saved figures.
#' @param species Optional species names for the legend.
#'
#' @return Invisibly returns a list of the printed plot objects.
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#' plot_diet_comp2(my_msm_fit)
#' }
#' @seealso [plot_diet_comp()], [plot_diet_comp1()]
#' @export
plot_diet_comp2 <- function(Rceattle, file = NULL, species = NULL) {

  # 1. SETUP & DATA PREPARATION ----
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for this function.")
  if (!requireNamespace("tidyr", quietly = TRUE))   stop("tidyr is required for this function.")
  if (!requireNamespace("dplyr", quietly = TRUE))   stop("dplyr is required for this function.")
  if (!requireNamespace("cowplot", quietly = TRUE)) stop("cowplot is required for this function.")
  if (!inherits(Rceattle, "Rceattle")) stop("Input 'Rceattle' must be a single object of class Rceattle.")
  if (is.null(species)) species <- Rceattle$data_list$spnames

  plot_data <- .diet_plot_data(Rceattle)
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    message("No diet data to plot.")
    return(invisible(NULL))
  }

  # 2. PLOTTING LOGIC ----
  plot_list <- list()

  for (pred_ind in 1:Rceattle$data_list$nspp) {
    for (prey_ind in 1:Rceattle$data_list$nspp) {

      subset_data <- plot_data |>
        dplyr::filter(.data$Pred == pred_ind, .data$Prey == prey_ind)
      if (nrow(subset_data) == 0) next

      # Detect aggregation level for this predator-prey subset only. A negative
      # age flags an aggregated (summed / averaged) dimension in diet_data.
      is_prey_age_agg <- any(subset_data$Prey_age < 0)
      is_pred_age_agg <- any(subset_data$Pred_age < 0)

      pred_legend <- paste("Pred-", species[pred_ind])

      # CASE 1: PREY-AGE AGGREGATED (predator age on x-axis) ----
      if (is_prey_age_agg && !is_pred_age_agg) {

        message(paste("Generating line plot (Pred Age) for Pred:", species[pred_ind], "- Prey:", species[prey_ind]))

        p <- ggplot2::ggplot(subset_data, ggplot2::aes(x = .data$Pred_age)) +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower_95, ymax = .data$upper_95, fill = "Observed"), alpha = 0.3) +
          {if ("Est_Lower" %in% names(subset_data)) {
            ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$Est_Lower, ymax = .data$Est_Upper, fill = "Estimated"), alpha = 0.3)
          }} +
          ggplot2::geom_line(ggplot2::aes(y = .data$Est, color = "Estimated", linetype = "Estimated"), linewidth = 1, alpha = 0.7) +
          ggplot2::geom_point(ggplot2::aes(y = .data$Est, color = "Estimated"), size = 2.5, alpha = 0.7) +
          ggplot2::geom_line(ggplot2::aes(y = .data$Observed, color = "Observed", linetype = "Observed"), linewidth = 1) +
          ggplot2::geom_point(ggplot2::aes(y = .data$Observed, color = "Observed"), size = 2.5) +
          ggplot2::facet_wrap(~ Year, scales = "free_y", labeller = ggplot2::labeller(Year = ~paste("Year:", .))) +
          ggplot2::facet_grid(~ Pred_sex, scales = "free_y", labeller = ggplot2::labeller(Pred_sex = ~paste("Sex:", .))) +
          ggplot2::scale_color_manual(name = "Source", values = c("Observed" = "black", "Estimated" = "darkred")) +
          ggplot2::scale_linetype_manual(name = "Source", values = c("Observed" = "dashed", "Estimated" = "solid")) +
          ggplot2::scale_fill_manual(name = "95% CI", values = c("Observed" = "grey50", "Estimated" = "darkred")) +
          ggplot2::labs(x = "Predator Age", y = "Diet Proportion", title = paste("Diet of", species[pred_ind], "on", species[prey_ind])) +
          ggplot2::theme_bw()

        print(p)
        plot_list[[length(plot_list) + 1]] <- p
        if (!is.null(file)) ggplot2::ggsave(paste0(file, "_diet_line_Pred", pred_ind, "_Prey", prey_ind, ".png"), p, width = 7, height = 5)

        # CASE 2: PREDATOR AGGREGATED (prey age on x-axis) ----
      } else if (!is_prey_age_agg && is_pred_age_agg) {

        message(paste("Generating line plot (Prey Age) for Pred:", species[pred_ind], "- Prey:", species[prey_ind]))

        p <- ggplot2::ggplot(subset_data, ggplot2::aes(x = .data$Prey_age)) +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower_95, ymax = .data$upper_95, fill = "Observed"), alpha = 0.3) +
          {if ("Est_Lower" %in% names(subset_data)) {
            ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$Est_Lower, ymax = .data$Est_Upper, fill = "Estimated"), alpha = 0.3)
          }} +
          ggplot2::geom_line(ggplot2::aes(y = .data$Est, color = "Estimated", linetype = "Estimated"), linewidth = 1, alpha = 0.7) +
          ggplot2::geom_point(ggplot2::aes(y = .data$Est, color = "Estimated"), size = 2.5, alpha = 0.7) +
          ggplot2::geom_line(ggplot2::aes(y = .data$Observed, color = "Observed", linetype = "Observed"), linewidth = 1) +
          ggplot2::geom_point(ggplot2::aes(y = .data$Observed, color = "Observed"), size = 2.5) +
          ggplot2::facet_wrap(~ Year + Pred_sex, scales = "free_y",
                              labeller = ggplot2::labeller(Year = ~paste("Year:", .), Pred_sex = ~paste("Sex:", .))) +
          ggplot2::scale_color_manual(name = "Source", values = c("Observed" = "black", "Estimated" = "darkblue")) +
          ggplot2::scale_linetype_manual(name = "Source", values = c("Observed" = "dashed", "Estimated" = "solid")) +
          ggplot2::scale_fill_manual(name = "95% CI", values = c("Observed" = "grey50", "Estimated" = "darkblue")) +
          ggplot2::labs(x = "Prey Age", y = "Diet Proportion",
                        title = paste("Prey Age Composition in Diet\nPred:", species[pred_ind], "eating", species[prey_ind])) +
          ggplot2::theme_bw()

        print(p)
        plot_list[[length(plot_list) + 1]] <- p
        if (!is.null(file)) ggplot2::ggsave(paste0(file, "_diet_preyage_Pred", pred_ind, "_Prey", prey_ind, ".png"), p, width = 7, height = 5)

        # CASE 3: BOTH AGGREGATED (bar plot) ----
      } else if (is_prey_age_agg && is_pred_age_agg) {

        message(paste("Generating bar plot for Pred:", species[pred_ind], "- Prey:", species[prey_ind]))

        plot_data_long <- subset_data |>
          tidyr::pivot_longer(cols = c("Observed", "Est"), names_to = "Source", values_to = "Proportion") |>
          dplyr::mutate(Prey_name = species[.data$Prey])

        p_fit <- ggplot2::ggplot(plot_data_long, ggplot2::aes(x = .data$Prey_name, y = .data$Proportion, fill = .data$Source)) +
          ggplot2::geom_bar(stat = "identity", position = "dodge") +
          ggplot2::facet_wrap(~ Year, scales = "free_y", labeller = ggplot2::labeller(Year = ~paste("Year:", .))) +
          ggplot2::scale_fill_manual(name = "Source", values = c("Observed" = "grey50", "Est" = "red")) +
          ggplot2::labs(x = "Prey Species", y = "Diet Proportion", title = paste("Fit to Aggregated Diet for Predator:", species[pred_ind])) +
          ggplot2::theme_bw()

        print(p_fit)
        plot_list[[length(plot_list) + 1]] <- p_fit
        if (!is.null(file)) ggplot2::ggsave(paste0(file, "_diet_barplot_Pred", pred_ind, ".png"), p_fit, width = 7, height = 5)

        # CASE 4: FULLY DISAGGREGATED (bubble plots) ----
      } else {

        message(paste("Generating bubble plots for Pred:", species[pred_ind], "- Prey:", species[prey_ind]))
        yrs <- sort(unique(subset_data$Year))

        for (i in seq_along(yrs)) {
          current_yr  <- yrs[i]
          comp_tmp_yr <- subset_data |> dplyr::filter(.data$Year == current_yr)
          if (sum(comp_tmp_yr$Observed, na.rm = TRUE) == 0) next

          title <- paste(pred_legend, "- Prey:", species[prey_ind], "- Year:", current_yr)
          if (current_yr == 0) title <- paste(pred_legend, "- Prey:", species[prey_ind], "(Avg over Years)")

          p_obs <- ggplot2::ggplot(comp_tmp_yr, ggplot2::aes(x = .data$Pred_age, y = .data$Prey_age, size = .data$Observed)) +
            ggplot2::geom_point(alpha = 0.7) + ggplot2::theme_classic() +
            ggplot2::labs(x = "Predator Age", y = "Prey Age", title = "Observed", size = "Prop.")

          p_est <- ggplot2::ggplot(comp_tmp_yr, ggplot2::aes(x = .data$Pred_age, y = .data$Prey_age, size = .data$Est)) +
            ggplot2::geom_point(alpha = 0.7) + ggplot2::theme_classic() +
            ggplot2::labs(x = "Predator Age", y = "Prey Age", title = "Estimated", size = "Prop.")

          p_pear <- ggplot2::ggplot(comp_tmp_yr, ggplot2::aes(x = .data$Pred_age, y = .data$Prey_age, size = .data$AbsPearson, color = .data$Pearson < 0)) +
            ggplot2::geom_point(alpha = 0.7) + ggplot2::theme_classic() +
            ggplot2::labs(x = "Predator Age", y = "Prey Age", title = "Pearson Residuals", size = "Abs(Resid)")

          p1 <- cowplot::plot_grid(p_obs, p_est, p_pear, nrow = 1)
          p1 <- cowplot::ggdraw(p1) + cowplot::draw_label(title, x = 0.5, y = 0.98)
          print(p1)
          plot_list[[length(plot_list) + 1]] <- p1

          if (!is.null(file)) {
            ggplot2::ggsave(paste0(file, "_diet_bubble_Pred", pred_ind, "_Prey", prey_ind, "_Yr", current_yr, ".png"), p1, width = 10, height = 6)
          }
        }
      }
    }
  }
  return(invisible(plot_list))
}
