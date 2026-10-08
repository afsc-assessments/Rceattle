# The plotters' SAVE path had no coverage at all: on 5.54.1 no test passed a
# non-NULL `file=` to any plotter, and the string "png" appeared nowhere under
# `tests/` or `tools/`. Three different implementations write the files --
# `.save_ggplot()` (19 call sites), `plot_comp()`'s own `save_png()` closure
# (1 definition, 3 uses), and 5 inline `ggplot2::ggsave()` calls in
# `R/7-plot_comp.R`, each with its own name and figure size.
#
# The filenames are an interface: 215 call sites across the sibling assessment
# repos pass `file =` to a plotter (169 in Rceattle-models, 34 in GOA-ATF-ESP,
# 12 in GOA_circlulation_study). None of them globs the result, so a rename
# would not error anywhere -- it would just leave a differently-named figure.
# `ggplot_build()` never touches saving, so nothing else can see this.
#
# So this records what every plotter writes today, keyed off the `file` prefix
# it was given. It is a characterisation test, not a statement that these are
# good names -- two of them are recorded below as defects.
#
# 29 of the 31 exported plotters write a file on the single-species BS2017SS
# fixture at `estimateMode = 3`. Only `plot_form` (a stub with no `file`
# formal) and `plot_profile` (needs an `Rceattle_profile`) do not, and they are
# named as exclusions. An earlier version of this file covered 10 and said the
# other 19 "need a multispecies fit, diet data, a profile object or comp/OSA
# inputs" -- that was wrong for 18 of the 19, including every diet plotter,
# which this fixture drives.
#
# WHERE THIS RUNS: nightly, not on a pull request. The file-level
# `skip_on_cran()` is the convention for the whole `test-plot-*.R` family (all
# 14 carry one) because each builds a real `fit_mod()` object, and
# `R-CMD-check.yaml` sets `NOT_CRAN: "false"`. `test-coverage` runs it but
# passes `stop_on_failure = FALSE`. So `deep-checks.yaml` is the only job where
# a failure here is fatal -- run this by hand in the PR that changes a save
# path rather than trusting CI to catch it.

testthat::skip_on_cran()

saved_files <- function(fit, fn, ...) {
  # tempfile() for uniqueness rather than runif(): this file would otherwise
  # advance the global RNG stream 29 times, and later unseeded test files draw
  # from it.
  dir <- tempfile(paste0("rcesave-", fn, "-"))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  prefix <- file.path(dir, "p")
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  suppressWarnings(suppressMessages(
    do.call(getExportedValue("Rceattle", fn), list(fit, file = prefix, ...))))
  list(names = sort(basename(list.files(dir, pattern = "[.]png$"))),
       sizes = file.size(sort(list.files(dir, pattern = "[.]png$",
                                         full.names = TRUE))))
}

# One figure, one exact filename. Measured, not guessed: the suffix is
# per-plotter and not derived from the function name -- plot_recruitment writes
# "_R_trajectory", plot_catch "_fishery_catch", plot_mortality
# "_mortality_at_age". That irregularity is why a unification needs this first.
.SAVE_ONE <- list(
  plot_b_eaten                 = "p_biomass_eaten.png",
  plot_b_eaten_prop            = "p_biomass_eaten_by_predator.png",
  plot_biomass                 = "p_biomass_trajectory.png",
  plot_catch                   = "p_fishery_catch.png",
  plot_catchability            = "p_catchability.png",
  plot_data                    = "p_data_plot.png",
  plot_depletion               = "p_biomass_depletion_trajectory.png",
  plot_depletionSSB            = "p_ssb_depletion_trajectory.png",
  plot_exploitable_biomass     = "p_exploitable_biomass_trajectory.png",
  plot_f                       = "p_f_trajectory.png",
  plot_index                   = "p_survey_indices.png",
  plot_indexresidual           = "p_index_residuals.png",
  plot_logindex                = "p_survey_indices.png",
  plot_m_at_age                = "p_m_at_age1.png",
  plot_m2_at_age_prop          = "p_m2_at_age_prop1.png",
  plot_maturity                = "p_maturity.png",
  plot_mortality               = "p_mortality_at_age.png",
  plot_ration                  = "p_ration1plus.png",
  plot_recruitment             = "p_R_trajectory.png",
  plot_selectivity             = "p_selectivity.png",
  plot_selectivity_vs_maturity = "p_selectivity_vs_maturity.png",
  plot_ssb                     = "p_ssb_trajectory.png",
  plot_ssb_depletion           = "p_ssb_depletion_trajectory.png",
  plot_stock_recruit           = "p_stock_recruit.png",
  plot_timeseries              = "p_biomass_trajectory.png")

# One figure per fleet, or per predator x prey. The count and the naming
# PATTERN are pinned rather than 42 literal names, because the diet names embed
# `spnames` from the dataset and a species rename is not what this test is
# about. The pattern is what an assessment script reads.
.SAVE_MANY <- list(
  plot_comp = list(
    n = 15L,
    re = paste0("^p_comp_(pearson|annual_fleet[0-9]+_(age|length)",
                "|aggregated_fleet[0-9]+_(age|length))[.]png$")),
  plot_diet_comp = list(
    n = 9L, re = "^p_aggregated_diet_comps_year[0-9]+_Pred- .+_prey_.+[.]png$"),
  plot_diet_comp1 = list(
    n = 9L, re = "^p_aggregated_diet_comps_year[0-9]+_Pred- .+_prey_.+[.]png$"),
  plot_diet_comp2 = list(
    n = 9L, re = "^p_diet_bubble_Pred[0-9]+_Prey[0-9]+_Yr[0-9]+[.]png$"))

# Each needs something this fixture cannot supply. Anything exported and not
# here must write a file.
.SAVE_NONE <- c(
  plot_form    = "a stub that refuses the Kinzey & Punt forms; takes no `file`",
  plot_profile = "takes an Rceattle_profile from profile(), not a fit")


.save_fixture <- function() {
  data("BS2017SS", package = "Rceattle", envir = environment())
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017SS, estimateMode = 3, msmMode = 0,
    fit_control = Rceattle::fit_control(verbose = 0))))
}


testthat::test_that("every plotter writes the file names it writes today", {
  testthat::skip_if_not_installed("TMB")
  ss <- .save_fixture()

  for (fn in names(.SAVE_ONE)) {
    got <- saved_files(ss, fn)
    testthat::expect_equal(got$names, .SAVE_ONE[[fn]], info = fn)
    # A name is not a figure: list.files() would be satisfied by a 0-byte file,
    # so assert the PNG has content.
    testthat::expect_gt(got$sizes[1], 1000, label = paste(fn, "file size"))
  }

  for (fn in names(.SAVE_MANY)) {
    spec <- .SAVE_MANY[[fn]]
    got <- saved_files(ss, fn)
    testthat::expect_length(got$names, spec$n)
    testthat::expect_equal(got$names[!grepl(spec$re, got$names)],
                           character(0), info = fn)
    testthat::expect_true(all(got$sizes > 1000), info = paste(fn, "sizes"))
  }
})


testthat::test_that("every exported plotter is covered or named as excluded", {
  # Keeps this file honest about its own coverage, the way test-plot-smoke.R
  # does: the uncovered ones are NAMED with a reason, not counted. A count is
  # not a set -- an earlier version asserted `length(uncovered) == 19L`, which
  # is implied by the other two assertions and stays green through a rename.
  exported <- sort(grep("^plot_", getNamespaceExports("Rceattle"),
                        value = TRUE))
  claimed <- c(names(.SAVE_ONE), names(.SAVE_MANY), names(.SAVE_NONE))

  testthat::expect_equal(sort(claimed), exported)
  testthat::expect_equal(anyDuplicated(claimed), 0L)
  # An exclusion for something no longer exported is a stale exclusion.
  testthat::expect_equal(setdiff(names(.SAVE_NONE), exported), character(0))

  # `file` is not always a formal: plot_logindex is `function(...)` forwarding
  # to plot_index(), so a formals test would drop it -- and
  # GOA-ATF-ESP/R/Run_2025_ceattle.R:274 calls it with `file =`. The tables
  # above are keyed on what a plotter WRITES, which is why they catch it.
  takes_file <- Filter(function(f) {
    "file" %in% names(formals(get(f, envir = asNamespace("Rceattle"))))
  }, exported)
  testthat::expect_false("plot_logindex" %in% takes_file)
  testthat::expect_true("plot_logindex" %in% names(.SAVE_ONE))
})


testthat::test_that("two save paths collide, and one filename is unusable", {
  # Recorded as current behaviour so a fix changes this block deliberately
  # rather than silently. Both are defects, not conventions.
  #
  # 1. plot_index()'s suffix is unconditional (R/7-plot_diagnostics.R,
  #    `suffix = "survey_indices"`), so the natural-scale and the log-scale
  #    index figure write the SAME file. plot_logindex() forwards with
  #    log = TRUE, so a script calling both keeps only whichever ran last.
  testthat::expect_equal(.SAVE_ONE[["plot_index"]],
                         .SAVE_ONE[["plot_logindex"]])

  # 2. plot_diet_comp() builds its name from `paste("Pred-", spnames[i])`, so
  #    every diet filename carries a space after the dash and the spaces inside
  #    a species name -- "p_aggregated_diet_comps_year1_Pred- Arrowtooth
  #    flounder_prey_Cod.png". plot_timeseries() already sanitises a species
  #    name for its CSV with gsub("[^A-Za-z0-9]+", "_", ...); the diet plotters
  #    do not use it.
  testthat::expect_match(.SAVE_MANY[["plot_diet_comp"]]$re, "Pred- ",
                         fixed = TRUE)

  # plot_biomass / plot_timeseries and plot_depletionSSB / plot_ssb_depletion
  # also share a filename, but each pair draws the SAME figure (one is an alias
  # or a default of the other), so those are correct.
  testthat::expect_equal(.SAVE_ONE[["plot_biomass"]],
                         .SAVE_ONE[["plot_timeseries"]])
  testthat::expect_equal(.SAVE_ONE[["plot_depletionSSB"]],
                         .SAVE_ONE[["plot_ssb_depletion"]])
})
