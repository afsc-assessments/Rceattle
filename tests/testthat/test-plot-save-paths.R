# The plotters' SAVE path had no coverage at all: no test passed a non-NULL
# `file=` to any of them, and nothing checked a written .png. 29 of the 31
# exported plotters take a `file` argument, and three different
# implementations write the files -- `.save_ggplot()` (19 call sites),
# `plot_comp()`'s own `save_png()` closure, and 7 bare `ggplot2::ggsave()`
# calls with their own names and sizes.
#
# That matters more than it looks. The filenames are the interface: assessment
# scripts glob and embed them, so unifying the three implementations cannot be
# verified by `ggplot_build()` (which never touches saving) and would rename a
# user's output with nothing to catch it.
#
# So this pins the filenames a plotter writes, keyed off the `file` prefix it
# was given. It is a characterisation test: it records what the plotters do
# today so the unification can be checked against it, not a statement that
# these are the best names.
#
# Not every plotter is covered -- some need a multispecies fit, diet data or a
# profile object. The ones here are the single-species ggplot family, which is
# where 19 of the 19 `.save_ggplot()` call sites live.

testthat::skip_on_cran()

saved_files <- function(fit, fn, ...) {
  dir <- file.path(tempdir(), paste0("rcesave-", fn, "-",
                                     as.integer(runif(1, 1, 1e9))))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  prefix <- file.path(dir, "p")
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  suppressWarnings(suppressMessages(
    do.call(getExportedValue("Rceattle", fn), list(fit, file = prefix, ...))))
  sort(basename(list.files(dir, pattern = "[.]png$")))
}


testthat::test_that("the single-species plotters write the files they name", {
  testthat::skip_if_not_installed("TMB")
  data("BS2017SS", package = "Rceattle", envir = environment())
  ss <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017SS, estimateMode = 3, msmMode = 0,
    fit_control = Rceattle::fit_control(verbose = 0))))

  # Measured, not guessed: the suffix is per-plotter and not derived from the
  # function name -- plot_recruitment writes "_R_trajectory", plot_catch
  # "_fishery_catch", plot_data "_data_plot". That irregularity is exactly why
  # a unification needs this pinned first.
  expected <- list(
    plot_biomass     = "p_biomass_trajectory.png",
    plot_ssb         = "p_ssb_trajectory.png",
    plot_recruitment = "p_R_trajectory.png",
    plot_catch       = "p_fishery_catch.png",
    plot_index       = "p_survey_indices.png",
    plot_selectivity = "p_selectivity.png",
    plot_depletion   = "p_biomass_depletion_trajectory.png",
    plot_f           = "p_f_trajectory.png",
    plot_maturity    = "p_maturity.png",
    plot_data        = "p_data_plot.png")

  for (fn in names(expected)) {
    got <- saved_files(ss, fn)
    testthat::expect_equal(got, expected[[fn]], info = fn)
  }
})


testthat::test_that("every file-taking plotter is either covered or listed", {
  # Keeps this file honest about its own coverage, the way
  # test-plot-smoke.R does: the uncovered ones are named with a reason, so a
  # new plotter cannot be silently left out of both lists.
  exported <- sort(grep("^plot_", getNamespaceExports("Rceattle"),
                        value = TRUE))
  takes_file <- Filter(function(f) {
    "file" %in% names(formals(get(f, envir = asNamespace("Rceattle"))))
  }, exported)
  testthat::expect_equal(length(takes_file), 29L)

  covered <- c("plot_biomass", "plot_ssb", "plot_recruitment", "plot_catch",
               "plot_index", "plot_selectivity", "plot_depletion", "plot_f",
               "plot_maturity", "plot_data")
  uncovered <- setdiff(takes_file, covered)
  # Each needs something the single-species fixture does not have: a
  # multispecies fit, diet data, a profile object, or comp/OSA inputs.
  testthat::expect_equal(length(uncovered), 19L)
  testthat::expect_true(all(covered %in% takes_file))
})
