# A figure's theme sets whether it has gridlines, and at what font size it is
# legible. Two plotters each used two of them INSIDE ONE FUNCTION, so a single
# call drew panels a reader is meant to compare at different sizes, or with and
# without gridlines:
#
#   plot_comp         theme_bw(base_size = 10) on the Pearson residual bubbles,
#                     theme_bw(base_size = 9) on the composition fit panels
#   plot_diet_comp2   theme_bw() on the line and bar fits (cases 1-3),
#                     theme_classic() on the bubble panels (case 4)
#
# Both are now one look per function, each resolved on the arm the rest of the
# package already agreed with:
#
#   plot_comp        base_size 10. Its bubbles are pinned to
#                    plot.rceattle_osa()'s size scale -- "the same residuals
#                    appear in both figures" (7-plot_comp.R) -- and plot_osa
#                    draws them at 10, so 9 was the odd one.
#   plot_diet_comp2  theme_classic. Its case-4 panels are the same
#                    observed / estimated / Pearson bubble triptych
#                    plot_diet_comp draws, and that function is theme_classic
#                    throughout, so cases 1-3 were the odd ones.
#
# The package keeps three looks, which is why this guard asserts consistency
# per function rather than one theme everywhere:
#
#   .rceattle_theme()        theme_classic, no gridlines -- trajectory figures
#                            (20 of the 24 plotters that set a theme)
#   theme_bw(base_size = 10) gridlines -- the residual and composition
#                            diagnostics, plot_osa and plot_comp
#   theme_classic()          the diet family, plot_diet_comp[2]
#
# Unifying those is a visual decision, not a drift fix: it would take the
# gridlines off every OSA and composition panel. Nothing here blocks it; the
# second block records the current convention so the change is a conscious one.
#
# Parsed from the call objects, not from the source lines, so a reformatted or
# line-wrapped call cannot blind it.

.theme_calls_by_function <- function(path) {
  theme_fns <- c("theme_bw", "theme_classic", "theme_minimal", "theme_grey",
                 "theme_gray", "theme_light", "theme_dark", "theme_linedraw",
                 "theme_void", ".rceattle_theme")

  # `ggplot2::theme_bw` and a bare `theme_bw` are the same call to a reader.
  # The operator has to be matched as a NAME: as.character() on a call flattens
  # it to operator-first, so `(grDevices::colorRampPalette(x))(n)` -- a real
  # call in 7-plot_helpers.R -- reads as "::" one level too deep and then has
  # no third element.
  call_name <- function(head) {
    if (is.name(head)) return(as.character(head))
    if (is.call(head) && length(head) == 3L && is.name(head[[1]]) &&
        as.character(head[[1]]) %in% c("::", ":::")) {
      return(as.character(head[[3]]))
    }
    NA_character_
  }

  # base_size is every theme_*()'s first formal, so a positional first
  # argument is a base_size. .rceattle_theme() shares that signature.
  base_size_of <- function(e) {
    if (length(e) < 2L) return("default")
    args <- as.list(e)[-1L]
    nms <- names(args)
    if (!is.null(nms) && "base_size" %in% nms) {
      return(paste(deparse(args[["base_size"]]), collapse = ""))
    }
    if (is.null(nms) || !nzchar(nms[1])) {
      return(paste(deparse(args[[1]]), collapse = ""))
    }
    "default"
  }

  found <- character()
  walk <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    nm <- call_name(e[[1]])
    if (!is.na(nm) && nm %in% theme_fns) {
      found <<- c(found, paste0(nm, "(", base_size_of(e), ")"))
    }
    for (k in seq_along(e)) {
      if (!identical(e[[k]], quote(expr = ))) walk(e[[k]])
    }
    invisible(NULL)
  }

  exprs <- parse(path, keep.source = FALSE)
  out <- list()
  for (i in seq_along(exprs)) {
    e <- exprs[[i]]
    if (!(is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=") &&
          is.call(e[[3]]) && as.character(e[[3]][[1]])[1] == "function")) next
    found <- character()
    walk(e[[3]])
    if (length(found)) out[[as.character(e[[2]])]] <- unique(found)
  }
  out
}


testthat::test_that("no plotter draws one figure family at two looks", {
  dir <- c("R", testthat::test_path("..", "..", "R"))
  dir <- dir[dir.exists(dir)]
  testthat::skip_if(length(dir) == 0, "R/ source not available")
  files <- list.files(dir[1], pattern = "^7-.*[.]R$", full.names = TRUE)
  testthat::expect_gt(length(files), 5L)

  mixed <- character()
  n_fun <- 0L
  for (f in files) {
    per_fun <- .theme_calls_by_function(f)
    n_fun <- n_fun + length(per_fun)
    for (fn in names(per_fun)) {
      looks <- per_fun[[fn]]
      if (length(looks) > 1L) {
        mixed <- c(mixed, sprintf("%s::%s uses %s", basename(f), fn,
                                  paste(looks, collapse = " and ")))
      }
    }
  }
  # Enumerated, not spot-checked: every plotter that sets a theme at all.
  testthat::expect_gt(n_fun, 20L)
  testthat::expect_equal(mixed, character())
})


testthat::test_that("the plotters use only the three package looks", {
  # Pins the convention the block above deliberately does not enforce. If a
  # fourth look appears, or .rceattle_theme() stops being the classic one, that
  # is a visual change to every figure and should arrive on purpose.
  dir <- c("R", testthat::test_path("..", "..", "R"))
  dir <- dir[dir.exists(dir)]
  testthat::skip_if(length(dir) == 0, "R/ source not available")
  files <- list.files(dir[1], pattern = "^7-.*[.]R$", full.names = TRUE)

  all_looks <- sort(unique(unlist(lapply(files, function(f) {
    unlist(.theme_calls_by_function(f), use.names = FALSE)
  }))))
  # The fourth entry is .rceattle_theme()'s own body passing its argument
  # through, not a fourth look.
  testthat::expect_equal(
    all_looks,
    c(".rceattle_theme(default)", "theme_bw(10)", "theme_classic(base_size)",
      "theme_classic(default)"))

  # .rceattle_theme() is the classic one, so the trajectory figures carry no
  # gridlines and the theme_bw() diagnostics do.
  looks <- .theme_calls_by_function(
    file.path(dir[1], "7-plot_helpers.R"))[[".rceattle_theme"]]
  testthat::expect_equal(looks, "theme_classic(base_size)")
  testthat::expect_true(inherits(.rceattle_theme()$panel.grid, "element_blank"))
  testthat::expect_true(inherits(ggplot2::theme_bw()$panel.grid,
                                 "element_line"))
  testthat::expect_equal(.rceattle_theme(base_size = 9)$text$size, 9)
})
