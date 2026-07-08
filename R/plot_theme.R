# plot_theme.R
# ---------------------------------------------------------------------------
# Shared figure style for the whole framework (see FIGURE_STYLE.md). Every
# rendered figure — R/ plot helpers, vignettes, simulation — uses ggplot2 with
# this theme and the Okabe-Ito colourblind-safe palette, so the manuscript's
# figures read as one visual family. No base-R graphics anywhere.
# ---------------------------------------------------------------------------

#' Okabe-Ito colourblind-safe palette. Pick series colours by name; never
#' introduce ad-hoc colours. Reference/guide lines are grey (see below).
#' @export
bc_palette <- c(grey = "#999999", orange = "#E69F00", skyblue = "#56B4E9",
                green = "#009E73", blue = "#0072B2", vermillion = "#D55E00",
                purple = "#CC79A7", yellow = "#F0E442")

#' The house theme: theme_classic() with a top legend and bold facet strips.
#' @param base_size base font size.
#' @return a ggplot2 theme.
#' @export
theme_bayesval <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "top",
      strip.text.x    = ggplot2::element_text(face = "bold"),
      strip.text.y    = ggplot2::element_text(face = "bold"))
}

#' pretty()-based axis ticks. Use as `scale_*_continuous(breaks = number_ticks(6))`.
#' @param n target number of ticks.
#' @export
number_ticks <- function(n) function(limits) pretty(limits, n)
