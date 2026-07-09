# plot_theme.R
# ---------------------------------------------------------------------------
# Shared figure style for the whole framework (see FIGURE_STYLE.md). Every
# rendered figure---R/ plot helpers, vignettes, simulation---uses ggplot2 with
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

#' Write a figure at house resolution through a device that emits plain sRGB.
#'
#' R's default bitmap device on macOS (`type = "quartz"`) bakes an Apple
#' "Generic RGB" ICC profile into the PNG *and* converts the pixel values into
#' it. LaTeX ignores the embedded profile and renders the raw values, so an
#' Okabe-Ito orange (`#E69F00`) prints as `#DE8E08` and the blue (`#0072B2`) as
#' `#0B5DA2`---darker, hue-shifted, and with the palette's contrast flattened.
#' Pin an sRGB device so the written pixels are the colours we asked for.
#'
#' @param filename output path.
#' @param plot a ggplot object.
#' @param width,height size in inches.
#' @param dpi resolution (house style: 600).
#' @param ... passed to [ggplot2::ggsave()].
#' @export
bc_ggsave <- function(filename, plot, width, height, dpi = 600, ...) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(filename, plot, device = ragg::agg_png,
                    width = width, height = height, dpi = dpi, ...)
  } else if (isTRUE(capabilities("cairo"))) {
    ggplot2::ggsave(filename, plot, device = "png", type = "cairo",
                    width = width, height = height, dpi = dpi, ...)
  } else {
    stop("Need `ragg` or a cairo-capable R to write sRGB figures; the quartz ",
         "device would embed an Apple Generic RGB profile.", call. = FALSE)
  }
}
