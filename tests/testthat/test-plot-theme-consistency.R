# A figure's theme sets whether it has gridlines, whether faceted panels are
# boxed, and at what font size it is legible. Two plotters each used two of
# them INSIDE ONE FUNCTION, so a single call drew panels a reader is meant to
# compare at different sizes, or with and without gridlines:
#
#   plot_comp         theme_bw(base_size = 10) on the Pearson residual bubbles,
#                     theme_bw(base_size = 9) on the composition fit panels
#   plot_diet_comp2   theme_bw() on the line and bar fits (cases 1-3),
#                     theme_classic() on the bubble panels (case 4)
#
# Both are now one look per function, each resolved on the arm the rest of the
# package already agreed with. plot_comp takes base_size 10, the size
# plot.rceattle_osa() draws the same residuals at; its own roxygen asks for "a
# consistent look with [plot.rceattle_osa()]". plot_diet_comp2 takes
# theme_classic, because its case-4 panels are the same observed / estimated /
# Pearson bubble triptych plot_diet_comp() draws and that function is
# theme_classic throughout.
#
# The package uses three looks, which is why these blocks assert consistency
# per function and per file rather than one theme everywhere:
#
#   .rceattle_theme()        theme_classic plus a panel border and bold,
#                            unboxed facet strips -- the trajectory figures,
#                            17 of the 24 functions that set a theme
#   theme_bw(base_size = 10) gridlines and a boxed panel -- the OSA panels in
#                            7-plot_osa.R and plot_comp
#   theme_classic()          the diet family, plot_diet_comp[2]
#
# Three limits of that, stated so none of it reads as tidier than it is.
#
# It is NOT "residual diagnostics get gridlines": plot_indexresidual() and
# plot_profile() draw residuals and profiles at .rceattle_theme(), with none.
#
# Pearson residual bubbles are drawn at theme_bw(10) in plot_comp and
# .osa_bubble_plot() but at theme_classic() in the diet family, which a
# per-function rule cannot see.
#
# And what these blocks assert is one look per set of theme CALLS, not per
# figure. plot_indexresidual(residual_type = "osa") and plot_comp() both
# delegate to plot.rceattle_osa(), so one function can still render two looks
# by argument; a look reached through a delegation is invisible here.
#
# theme_classic() also leaves panel.border blank, so the faceted diet panels
# (cases 1-3) are delimited by their strip boxes -- black at linewidth 1 -- and
# not by a box around each panel.
#
# .rceattle_theme() is NOT available as a third resolution here: it sets
# legend.title = element_blank(), which would delete the "Source", "95% CI",
# "Prop." and "Abs(Resid)" legend names the diet plotters set on every panel.
#
# Parsed from the call objects, not the source lines, so a reformatted or
# line-wrapped call cannot blind it.

.THEME_FNS <- c("theme_bw", "theme_classic", "theme_minimal", "theme_grey",
                "theme_gray", "theme_light", "theme_dark", "theme_linedraw",
                "theme_void", ".rceattle_theme")

# `ggplot2::theme_bw` and a bare `theme_bw` are the same call to a reader. The
# operator has to be matched as a NAME: as.character() on a call flattens it
# operator-first, so `(grDevices::colorRampPalette(x))(n)` -- a real call in
# 7-plot_helpers.R -- reads as "::" one level too deep and then has no third
# element.
.theme_call_name <- function(head) {
  if (is.name(head)) return(as.character(head))
  if (is.call(head) && length(head) == 3L && is.name(head[[1]]) &&
      as.character(head[[1]]) %in% c("::", ":::")) {
    return(as.character(head[[3]]))
  }
  NA_character_
}

# Collect the looks a single expression sets, as `theme(size)` strings.
#
# Two spellings count as setting a base size. A theme_*() call takes it as its
# first formal -- true of all nine ggplot2 themes and of .rceattle_theme(). And
# `theme(text = element_text(size = 9))` overrides it afterwards, which is the
# most natural way to reintroduce the defect above while leaving the theme_*()
# call alone; without it a plotter could draw two sizes with this guard green.
.theme_looks <- function(e) {
  # 10 and 10L are the same font size, so compare the VALUE where the argument
  # is a literal; anything else (a variable, an expression) keeps its text.
  as_label <- function(x) {
    v <- tryCatch(eval(x, baseenv()), error = function(e) NULL)
    if (is.numeric(v) && length(v) == 1L) format(v) else
      paste(deparse(x), collapse = "")
  }

  base_size_of <- function(call) {
    if (length(call) < 2L) return("default")
    args <- as.list(call)[-1L]
    nms <- names(args)
    if (!is.null(nms) && "base_size" %in% nms) {
      return(as_label(args[["base_size"]]))
    }
    # A positional base_size can follow a named argument, since R matches names
    # first: theme_bw(base_family = "x", 10) is base_size 10.
    unnamed <- if (is.null(nms)) seq_along(args) else which(!nzchar(nms))
    named <- if (is.null(nms)) character() else nms[nzchar(nms)]
    if (length(unnamed) && !("base_size" %in% named)) {
      return(as_label(args[[unnamed[1]]]))
    }
    "default"
  }

  # theme(text = element_text(size = <n>)) -- a base-size override.
  text_size_of <- function(call) {
    args <- as.list(call)[-1L]
    nms <- names(args)
    if (is.null(nms) || !("text" %in% nms)) return(NULL)
    el <- args[["text"]]
    if (!is.call(el) ||
        !identical(.theme_call_name(el[[1]]), "element_text")) return(NULL)
    el_args <- as.list(el)[-1L]
    if (is.null(names(el_args)) || !("size" %in% names(el_args))) return(NULL)
    paste0("theme(text ", as_label(el_args[["size"]]), ")")
  }

  found <- character()
  walk <- function(x) {
    if (!is.call(x)) return(invisible(NULL))
    nm <- .theme_call_name(x[[1]])
    if (!is.na(nm) && nm %in% .THEME_FNS) {
      found <<- c(found, paste0(nm, "(", base_size_of(x), ")"))
    } else if (identical(nm, "theme")) {
      found <<- c(found, text_size_of(x))
    }
    for (k in seq_along(x)) {
      if (!identical(x[[k]], quote(expr = ))) walk(x[[k]])
    }
    invisible(NULL)
  }
  walk(e)
  unique(found)
}

# Looks per top-level function. A plotter written any other way -- an alias, a
# `.ts_wrapper()` call, a theme held in a top-level constant -- has no function
# body here to attribute a look to, so the file-level block below is what
# covers it.
.theme_calls_by_function <- function(path) {
  exprs <- parse(path, keep.source = FALSE)
  out <- list()
  for (i in seq_along(exprs)) {
    e <- exprs[[i]]
    if (!(is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=") &&
          is.name(e[[2]]) && is.call(e[[3]]) &&
          as.character(e[[3]][[1]])[1] == "function")) next
    looks <- .theme_looks(e[[3]])
    if (length(looks)) out[[as.character(e[[2]])]] <- looks
  }
  out
}

.plot_source_dir <- function(candidates) {
  candidates[dir.exists(candidates)]
}


testthat::test_that("no plotter draws one figure family at two looks", {
  dir <- .plot_source_dir(c("R", testthat::test_path("..", "..", "R")))
  testthat::skip_if(length(dir) == 0, "R/ source not available")
  files <- list.files(dir[1], pattern = "^7-.*[.]R$", full.names = TRUE)
  testthat::expect_gt(length(files), 5L)

  mixed <- character()
  n_fun <- 0L
  for (f in files) {
    per_fun <- .theme_calls_by_function(f)
    n_fun <- n_fun + length(per_fun)
    for (fn in names(per_fun)) {
      if (length(per_fun[[fn]]) > 1L) {
        mixed <- c(mixed, sprintf("%s::%s uses %s", basename(f), fn,
                                  paste(per_fun[[fn]], collapse = " and ")))
      }
    }
  }
  # Every function in R/7-*.R written as `name <- function(...)` that sets a
  # theme, which is 24 of them. It is not every exported plotter: 13 exports
  # are aliases or `.ts_wrapper()` calls with no body of their own, and they
  # set no theme -- the file-level block is what holds them.
  testthat::expect_equal(n_fun, 24L)
  testthat::expect_equal(mixed, character())
})


testthat::test_that("each plot file keeps the look its figures are drawn at", {
  # A per-function rule alone cannot see a whole family changing look together:
  # switching plot_timeseries() from .rceattle_theme() to theme_bw() puts
  # gridlines on all seven exported trajectory plotters and leaves that
  # function internally consistent. This pins the assignment per FILE.
  #
  # Theme FUNCTION names only, deliberately. base_size is an argument with a
  # documented default, so .rceattle_theme(base_size = 14) for a poster figure
  # is a legitimate per-call choice and must not fail here; a size mismatch
  # WITHIN one function is the previous block's job.
  dir <- .plot_source_dir(c("R", testthat::test_path("..", "..", "R")))
  testthat::skip_if(length(dir) == 0, "R/ source not available")
  files <- list.files(dir[1], pattern = "^7-.*[.]R$", full.names = TRUE)

  by_file <- vapply(files, function(f) {
    looks <- unlist(.theme_calls_by_function(f), use.names = FALSE)
    paste(sort(unique(sub("[(].*$", "", looks))), collapse = " + ")
  }, character(1))
  names(by_file) <- basename(files)

  testthat::expect_equal(as.list(by_file), list(
    "7-plot_ceattle.R"       = ".rceattle_theme",
    "7-plot_comp.R"          = "theme_bw + theme_classic",
    "7-plot_data.R"          = ".rceattle_theme",
    "7-plot_diagnostics.R"   = ".rceattle_theme",
    "7-plot_helpers.R"       = "theme_classic",
    "7-plot_osa.R"           = "theme_bw",
    "7-plot_profile.R"       = ".rceattle_theme",
    "7-plot_sel_vs_mat.R"    = ".rceattle_theme",
    "7-plot_stock_recruit.R" = ".rceattle_theme"))

  # 7-plot_comp.R carries two because plot_comp is a composition diagnostic and
  # the diet plotters are not; each function in it is internally consistent.
  comp <- .theme_calls_by_function(file.path(dir[1], "7-plot_comp.R"))
  testthat::expect_equal(comp[["plot_comp"]], "theme_bw(10)")
  testthat::expect_equal(comp[["plot_diet_comp"]], "theme_classic(default)")
  testthat::expect_equal(comp[["plot_diet_comp2"]], "theme_classic(default)")
})


testthat::test_that(".rceattle_theme() is the classic one, with a panel border", {
  # The house theme's identity, so a change to the look of every trajectory
  # figure arrives on purpose. Read with [[ ]]: a theme list carries
  # panel.grid.major, panel.grid.minor and four more alongside `panel.grid`, so
  # a `$` read here is one ggplot2 release away from partial-matching.
  dir <- .plot_source_dir(c("R", testthat::test_path("..", "..", "R")))
  testthat::skip_if(length(dir) == 0, "R/ source not available")
  looks <- .theme_calls_by_function(
    file.path(dir[1], "7-plot_helpers.R"))[[".rceattle_theme"]]
  testthat::expect_equal(looks, "theme_classic(base_size)")

  th <- .rceattle_theme()
  testthat::expect_true(inherits(th[["panel.grid"]], "element_blank"))
  testthat::expect_true(inherits(th[["panel.border"]], "element_rect"))
  # Bold, unboxed facet strips: the panel border is what delimits a facet, so
  # the strip box is removed rather than drawn twice.
  testthat::expect_true(inherits(th[["strip.background"]], "element_blank"))
  testthat::expect_equal(.rceattle_theme(base_size = 9)[["text"]]$size, 9)

  # The other two looks differ from it in exactly the way the header says.
  testthat::expect_true(
    inherits(ggplot2::theme_bw()[["panel.grid"]], "element_line"))
  testthat::expect_true(
    inherits(ggplot2::theme_classic()[["panel.border"]], "element_blank"))
})
