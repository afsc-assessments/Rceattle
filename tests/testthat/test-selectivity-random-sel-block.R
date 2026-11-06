# Provenance: the developer backlog, Tier 0 -- the R/3-build_map.R marker
# "will fail if random_sel = TRUE?" on the Time_varying_sel = "Block" arm.
#
# Confirmed 2026-08-23 against dev. build_map() maps the block parameters into
# log_sel_slp_dev / sel_inf_dev and turns the main log_sel_slp / sel_inf off,
# and fit_mod() used to add those arrays to TMB's `random` unconditionally. But
# the template's densities on them are gated on flt_varying_sel in {1,2}
# (IID/AR1) and {4,5} (RandomWalk*); Block is code 3, documented as "time blocks
# with no penalty". Measured on a Logistic 2-block fixture at estimateMode = 3:
# 8 parameters declared random, jnll_comp(JNLL_SEL_DEV, ) identically 0, and a
# NaN objective; a real fit died with "NA/NaN gradient evaluation". No bundled
# dataset reaches it -- BS2017SS/BS2017MS are Time_varying_sel 0, GOA2018SS is 0
# and 5 -- which is why it survived.

# obj$report() returns jnll_comp unnamed -- the display names in
# R/6-rename_output.R are put on $quantities, not on the report -- so read the
# row from the JnllRow enum rather than hardcoding it. CLAUDE.md: refer to a row
# by its constant.
jnll_sel_dev_row <- function() {
  dir <- c("src/TMB", testthat::test_path("..", "..", "src", "TMB"))
  dir <- dir[dir.exists(dir)]
  testthat::skip_if(length(dir) == 0, "src/TMB not available")
  src <- readLines(file.path(dir[1], "ceattle.cpp"), warn = FALSE)
  hit <- grep("JNLL_SEL_DEV\\s*=\\s*[0-9]+", src, value = TRUE)
  testthat::expect_length(hit, 1)
  as.integer(sub(".*JNLL_SEL_DEV\\s*=\\s*([0-9]+).*", "\\1", hit)) + 1L
}

blocked_data <- function(tv_sel, nyrs = 10, nages = 5) {
  d <- make_test_data(nyrs = nyrs, nages = nages, seed = 42)
  d$fleet_control$Selectivity        <- "Logistic"
  d$fleet_control$Time_varying_sel   <- tv_sel
  d$fleet_control$Selectivity_index  <- seq_len(nrow(d$fleet_control))
  d$fleet_control$Age_first_selected <- 1
  if (tv_sel == "Block") {
    half <- floor(nyrs / 2)
    d$catch_data$Selectivity_block <-
      ifelse(d$catch_data$Year - d$styr + 1 <= half, 1L, 2L)
    d$index_data$Selectivity_block <-
      ifelse(d$index_data$Year - d$styr + 1 <= half, 1L, 2L)
  }
  d
}

build_only <- function(d, random_sel) {
  suppressMessages(suppressWarnings(
    fit_mod(data_list = d, inits = NULL, estimateMode = 3, msmMode = 0,
            random_rec = FALSE, random_sel = random_sel,
            fit_control = fit_control(verbose = 0))
  ))
}


test_that("random_sel = TRUE is refused on a fleet with selectivity blocks", {
  testthat::skip_on_cran()

  expect_error(
    build_only(blocked_data("Block"), random_sel = TRUE),
    "Time_varying_sel"
  )
  # The message has to name the way out, not just the problem -- the failure it
  # replaces was TMB's "NA/NaN gradient evaluation".
  expect_error(build_only(blocked_data("Block"), random_sel = TRUE), "RandomWalk")
})


test_that("selectivity blocks still fit as the fixed effects they are", {
  testthat::skip_on_cran()

  m <- build_only(blocked_data("Block"), random_sel = FALSE)

  expect_length(m$obj$env$random, 0)
  expect_true(is.finite(m$obj$fn(m$obj$par)))

  # Blocks carry no penalty, so the deviate row stays empty. If this ever goes
  # non-zero, Block has acquired a density and the refusal above should be
  # revisited rather than kept.
  expect_true(all(m$obj$report(m$obj$env$last.par)$jnll_comp[jnll_sel_dev_row(), ] == 0))
})


test_that("the time-varying modes that DO carry a density still integrate out", {
  testthat::skip_on_cran()

  for (tv in c("IID", "RandomWalk")) {
    m <- build_only(blocked_data(tv), random_sel = TRUE)
    expect_gt(length(m$obj$env$random), 0)
    expect_true(is.finite(m$obj$fn(m$obj$par)))
    expect_true(any(m$obj$report(m$obj$env$last.par)$jnll_comp[jnll_sel_dev_row(), ] > 0))
  }

  # Time_varying_sel = "Off" leaves every deviate mapped out. TMB drops the
  # fully-mapped names itself, so the fit is unchanged -- but the names stay in
  # `random` because TMBphase() reads its length (R/6-phaser.R).
  m_off <- build_only(blocked_data("Off"), random_sel = TRUE)
  expect_length(m_off$obj$env$random, 0)
  expect_true(is.finite(m_off$obj$fn(m_off$obj$par)))
})


test_that("Block is the ONLY time-varying mode the template leaves unscored", {
  # The drift guard, and the reason the refusal above is scoped to Block alone.
  # Every JNLL_SEL_DEV gate keyed on the mode is a `flt_varying_sel(flt) == N`
  # comparison, so the modes with no density are exactly the tv_sel_map codes
  # that never appear in one. "Off" estimates no deviate and needs none; "Block"
  # estimates one and gets none. A mode added to tv_sel_map without a density --
  # or a density added for Block -- turns this red.
  dir <- c("src/TMB", testthat::test_path("..", "..", "src", "TMB"))
  dir <- dir[dir.exists(dir)]
  testthat::skip_if(length(dir) == 0, "src/TMB not available")

  src <- paste(readLines(file.path(dir[1], "ceattle.cpp"), warn = FALSE),
               collapse = "\n")
  src <- gsub("(?s)/\\*.*?\\*/", "", src, perl = TRUE)
  src <- gsub("//[^\n]*", "", src, perl = TRUE)

  scored <- regmatches(src, gregexpr("flt_varying_sel\\(flt\\)\\s*==\\s*[0-9]+",
                                     src, perl = TRUE))[[1]]
  scored <- sort(unique(as.integer(sub(".*==\\s*", "", scored))))

  expect_setequal(scored, c(1L, 2L, 4L, 5L))
  expect_setequal(setdiff(unname(tv_sel_map), scored),
                  unname(tv_sel_map[c("Off", "Block")]))
})
