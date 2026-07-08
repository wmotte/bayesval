# posterior_summary.R
# ---------------------------------------------------------------------------
# Unified reporting for the framework: turn posterior_metric objects into the
# manuscript's headline table (posterior median [95% CrI]) and provide the
# calibration / decision-curve plots as ggplot2 objects.
# ---------------------------------------------------------------------------

#' Format a posterior_metric (or several) as a median `[95% CrI]` table.
#'
#' @param ...   named posterior_metric objects, or a single named list of them.
#' @param probs credible-interval probabilities.
#' @param digits rounding for the formatted `"median [lo, hi]"` column.
#' @return tibble with columns metric, median, lower, upper, and a
#'         pre-formatted `summary` string.
#' @export
posterior_table <- function(..., probs = c(0.025, 0.5, 0.975), digits = 3) {
  args <- list(...)
  if (length(args) == 1L && is.list(args[[1]]) &&
      !inherits(args[[1]], "posterior_metric")) args <- args[[1]]
  fmt <- paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]")
  purrr::map(args, function(pm) summary(pm, probs = probs)) |>
    dplyr::bind_rows() |>
    dplyr::mutate(summary = sprintf(fmt, median, lower, upper))
}

#' Side-by-side Bayesian vs frequentist headline table (Part D).
#' @param bayes  a posterior_table() (or a posterior_metric) with the same
#'               metric names as `freq`.
#' @param freq   data.frame with columns metric, estimate, lower, upper.
#' @param digits rounding.
#' @export
posterior_vs_frequentist <- function(bayes, freq, digits = 3) {
  if (inherits(bayes, "posterior_metric")) bayes <- posterior_table(bayes)
  fmt <- paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]")
  dplyr::select(bayes, metric, median, lower, upper) |>
    dplyr::inner_join(tibble::as_tibble(freq), by = "metric",
                      suffix = c("_bayes", "_freq")) |>
    dplyr::transmute(
      metric,
      bayesian    = sprintf(fmt, median, lower_bayes, upper_bayes),
      frequentist = sprintf(fmt, estimate, lower_freq, upper_freq))
}

#' Plot a posterior calibration band from ppc_calibration() (ggplot2).
#'
#' House style (FIGURE_STYLE.md): Okabe-Ito palette, grey dashed diagonal,
#' credible band as a filled ribbon under the median line.
#' @param ppc output of ppc_calibration().
#' @param main plot title; default summarises the horizon and posterior
#'        predictive p-value.
#' @return a ggplot object.
#' @export
plot_calibration_band <- function(ppc, main = NULL) {
  b <- ppc$band
  if (is.null(main)) main <- sprintf("Bayesian calibration (t = %s, ppp = %.2f)",
                                     signif(ppc$horizon, 3), ppc$ppp)
  ggplot2::ggplot(b, ggplot2::aes(x = pred)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = bc_palette[["orange"]], alpha = 0.25) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2,
                         colour = "gray50") +
    ggplot2::geom_line(ggplot2::aes(y = median), colour = bc_palette[["orange"]],
                       linewidth = 1) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(breaks = number_ticks(6)) +
    ggplot2::scale_y_continuous(breaks = number_ticks(6)) +
    ggplot2::labs(x = "Predicted risk", y = "Observed risk", title = main) +
    theme_bayesval()
}

#' Plot posterior net-benefit bands from posterior_net_benefit() (ggplot2).
#'
#' House style (FIGURE_STYLE.md): model curve in the primary palette colour with
#' its credible ribbon, treat-all in grey, treat-none the grey zero line.
#' @param pnb output of posterior_net_benefit().
#' @param main plot title.
#' @return a ggplot object.
#' @export
plot_net_benefit <- function(pnb, main = "Posterior decision curve") {
  m <- pnb$model
  lines_df <- dplyr::bind_rows(
    dplyr::mutate(dplyr::select(m, threshold, nb = median), strategy = "Model"),
    dplyr::mutate(dplyr::select(pnb$all, threshold, nb = median),
                  strategy = "Treat all"))
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = m, ggplot2::aes(x = threshold, ymin = lower, ymax = upper),
      fill = bc_palette[["orange"]], alpha = 0.2) +
    ggplot2::geom_hline(yintercept = 0, linetype = 3, colour = "gray40") +
    ggplot2::geom_line(
      data = lines_df,
      ggplot2::aes(x = threshold, y = nb, colour = strategy,
                   linetype = strategy), linewidth = 1) +
    ggplot2::scale_colour_manual(
      values = c(Model = bc_palette[["orange"]], "Treat all" = "gray40")) +
    ggplot2::scale_linetype_manual(
      values = c(Model = 1, "Treat all" = 2)) +
    ggplot2::scale_x_continuous(breaks = number_ticks(6)) +
    ggplot2::scale_y_continuous(breaks = number_ticks(6)) +
    ggplot2::labs(x = "Threshold probability", y = "Net benefit",
                  title = main, colour = NULL, linetype = NULL,
                  caption = "Dotted line: treat none (net benefit 0)") +
    theme_bayesval()
}
