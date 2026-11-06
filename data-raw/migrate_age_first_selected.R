# =============================================================================
# Move the bundled datasets' age-based floors to Age_first_selected
# =============================================================================
#
# On an age-based fleet the first selected age is `Age_first_selected`, an
# absolute age, and `Bin_first_selected` there is its deprecated spelling. The
# bundled datasets carried the old spelling, so fitting one printed the
# deprecation message. This moves each age-based fleet's floor to
# `Age_first_selected` (age = minage + Bin_first_selected - 1), leaves a bin of
# 1 blank, since blank is selected from minage whatever minage is, and blanks
# `Bin_first_selected` there; length-based fleets are untouched. The
# verification asserts the first selected bin the model reads is unchanged.
#
# Re-run with:  Rscript data-raw/migrate_age_first_selected.R
# =============================================================================

suppressMessages(pkgload::load_all(".", quiet = TRUE))

first_selected <- Rceattle:::.rce_first_selected

datasets <- sub("[.]rda$", "", list.files("data", pattern = "[.]rda$"))

for (ds in datasets) {
  f <- file.path("data", paste0(ds, ".rda"))
  e <- new.env()
  load(f, envir = e)
  if (!ds %in% ls(e)) next
  before <- get(ds, envir = e)
  if (!is.list(before) || is.null(before$fleet_control)) next
  fc <- before$fleet_control
  if (!"Bin_first_selected" %in% names(fc)) next

  res  <- first_selected(before)
  move <- res$from_bin
  if (!any(move)) next

  after <- before
  afs <- if ("Age_first_selected" %in% names(fc)) fc$Age_first_selected else rep(NA_integer_, nrow(fc))
  # Bin 1 is the first age whatever minage is, so it becomes a blank, not an age.
  floor <- move & fc$Bin_first_selected > 1
  afs[floor] <- res$age[floor]
  bfs <- fc$Bin_first_selected
  bfs[move] <- NA
  fc$Bin_first_selected <- as.integer(bfs)
  fc$Age_first_selected <- as.integer(afs)
  # Schema order: Age_first_selected directly after Bin_first_selected.
  nms <- setdiff(names(fc), "Age_first_selected")
  at  <- match("Bin_first_selected", nms)
  fc  <- fc[, append(nms, "Age_first_selected", after = at), drop = FALSE]
  after$fleet_control <- fc

  # The model must read the same first selected bin; a blank is bin 1.
  res2 <- first_selected(after)
  bin1 <- function(b) ifelse(is.na(b), 1, b)
  stopifnot(identical(bin1(res$bin), bin1(res2$bin)), !any(res2$from_bin))

  assign(ds, after, envir = e)
  save(list = ds, file = f, envir = e, compress = "xz")
  cat(sprintf("migrated %-22s (%d age-based fleets)\n", ds, sum(move)))
}
