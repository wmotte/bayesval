# worked_example_composite.R
# ===========================================================================
# COMPOSITE FIGURES: the Rotterdam-to-GBSG worked example, drawn once per
# backend and faceted.
#
# `worked_example.R` (Cox-Laplace) and `worked_example_stan.R` (full-Bayes
# piecewise-exponential MCMC) validate the SAME model on the SAME cohort at the
# SAME 5-year horizon, and each already draws its own calibration band and its
# own decision curve. Shown alone, either one illustrates the method. Shown side
# by side they also carry the appendix's claim -- the fast approximate backend
# and the genuine full-Bayes fit agree to within Monte-Carlo error -- as
# something the reader can see instead of a sentence they must take on trust.
# Two panels per figure, at the page width one panel used to take.
#
# Prerequisites -- each source script writes the results file consumed here:
#
#   Rscript vignettes/worked_example.R        # worked_example_results.rds
#   Rscript vignettes/worked_example_stan.R   # worked_example_stan_results.rds
#   Rscript vignettes/worked_example_composite.R
#
# Companion paper: Otte WM, Kampman SL, [...] A Bayesian framework for
# validating survival prediction models: posterior credible intervals for
# discrimination, calibration, and clinical utility. in prep.
# ===========================================================================

source(file.path("..", "R", "load_all.R"), chdir = TRUE)

figdir <- "figures"

read_results <- function(file) {
  path <- file.path(figdir, file)
  if (!file.exists(path))
    stop("Missing `", path, "`. Run the two worked-example scripts first; see ",
         "the header of this file for the order.", call. = FALSE)
  readRDS(path)
}

lap  <- read_results("worked_example_results.rds")
bayes <- read_results("worked_example_stan_results.rds")
if (is.null(bayes$pnb))
  stop("`worked_example_stan_results.rds` predates the `pnb` field. Re-run ",
       "vignettes/worked_example_stan.R.", call. = FALSE)

## Backend order is fixed: the approximation first, the reference it approximates
## second, so the eye reads left to right as "does the shortcut hold up".
bk <- c("Cox-Laplace approximation", "Full Bayes (MCMC)")
by_backend <- function(x_lap, x_bayes) {
  dplyr::bind_rows(dplyr::mutate(tibble::as_tibble(x_lap),   backend = bk[1]),
                   dplyr::mutate(tibble::as_tibble(x_bayes), backend = bk[2])) |>
    dplyr::mutate(backend = factor(backend, levels = bk))
}

## ---- Posterior predictive calibration (B4) --------------------------------
## House style: see FIGURE_STYLE.md (Okabe-Ito palette + theme_bayesval()).
## Fixed scales across panels, so the two bands can be compared by position and
## not only by shape; the diagonal is perfect calibration in both.
band <- by_backend(lap$ppc$band, bayes$ppc$band)

p_cal <- ggplot2::ggplot(band, ggplot2::aes(x = pred)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                       fill = bc_palette[["orange"]], alpha = 0.25) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2,
                       colour = "gray50") +
  ggplot2::geom_line(ggplot2::aes(y = median), colour = bc_palette[["orange"]],
                     linewidth = 1) +
  ggplot2::facet_wrap(~ backend) +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  ggplot2::scale_x_continuous(breaks = number_ticks(6)) +
  ggplot2::scale_y_continuous(breaks = number_ticks(6)) +
  # No panel title: the two posterior predictive p-values belong in the caption.
  ggplot2::labs(x = "Predicted risk", y = "Observed risk") +
  theme_bayesval() +
  # expand = FALSE puts the "1.0" tick label exactly on the panel edge, so the
  # gutter has to hold two half-labels, not the usual whitespace.
  ggplot2::theme(panel.spacing.x = grid::unit(2.4, "lines"))
bc_ggsave(file.path(figdir, "calibration_ppc_backends.png"), p_cal,
          width = 6.5, height = 3.6)

## ---- Posterior decision curve (B6) ----------------------------------------
## Treat-none is the horizontal zero line, drawn dotted; treat-all is the grey
## dashed curve. Only the model curve carries a credible band, because treat-all
## is a strategy and not an estimate of this model's performance.
strat <- c("Model", "Treat all")
nb_lines <- dplyr::bind_rows(
  dplyr::mutate(by_backend(lap$pnb$model, bayes$pnb$model), strategy = strat[1]),
  dplyr::mutate(by_backend(lap$pnb$all,   bayes$pnb$all),   strategy = strat[2])) |>
  dplyr::mutate(strategy = factor(strategy, levels = strat))
nb_band <- by_backend(lap$pnb$model, bayes$pnb$model)

p_dca <- ggplot2::ggplot(mapping = ggplot2::aes(x = threshold)) +
  ggplot2::geom_ribbon(data = nb_band,
                       ggplot2::aes(ymin = lower, ymax = upper),
                       fill = bc_palette[["orange"]], alpha = 0.2) +
  ggplot2::geom_hline(yintercept = 0, linetype = 3, colour = "gray40") +
  ggplot2::geom_line(data = nb_lines,
                     ggplot2::aes(y = median, colour = strategy,
                                  linetype = strategy), linewidth = 1) +
  ggplot2::facet_wrap(~ backend) +
  ggplot2::scale_colour_manual(
    values = stats::setNames(c(bc_palette[["orange"]], "gray40"), strat)) +
  ggplot2::scale_linetype_manual(values = stats::setNames(c(1, 2), strat)) +
  ggplot2::scale_x_continuous(breaks = number_ticks(6)) +
  ggplot2::scale_y_continuous(breaks = number_ticks(6)) +
  # Treat-none is the dotted zero line; that is caption text, not a panel legend.
  ggplot2::labs(x = "Threshold probability", y = "Net benefit",
                colour = NULL, linetype = NULL) +
  theme_bayesval()
bc_ggsave(file.path(figdir, "decision_curve_backends.png"), p_dca,
          width = 6.5, height = 3.6)

## ---- Numbers the captions quote -------------------------------------------
cat("=====================================================================\n")
cat("WORKED-EXAMPLE BACKEND COMPOSITES (Rotterdam -> GBSG, 5-year horizon)\n")
cat(sprintf("  posterior predictive p-value: %s %.3f | %s %.3f\n",
            bk[1], lap$ppc$ppp, bk[2], bayes$ppc$ppp))
nb_at <- function(pnb, t) {
  i <- which.min(abs(pnb$model$threshold - t))
  sprintf("%.3f [%.3f, %.3f]", pnb$model$median[i], pnb$model$lower[i],
          pnb$model$upper[i])
}
cat(sprintf("  net benefit @20%%:            %s %s | %s %s\n",
            bk[1], nb_at(lap$pnb, 0.20), bk[2], nb_at(bayes$pnb, 0.20)))
cat(sprintf("  figures -> vignettes/%s/{calibration_ppc,decision_curve}_backends.png\n",
            figdir))
cat("=====================================================================\n")
