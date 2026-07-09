# showcase_joint_composite.R
# ===========================================================================
# COMPOSITE FIGURE: the joint gate posterior for both showcase cohorts.
#
# The simulated small cohort (main text) and the Mayo PBC external validation
# (Supplemental Digital Content) run the SAME three-criterion deployment gate on
# the SAME two axes. They belong in one faceted panel, not two figures: side by
# side the reader sees a gate that almost always fails (simulated, blue fraction
# 0.05) next to one that is genuinely undecided (PBC, 0.27), at no extra cost in
# page space. Each panel also carries its own C bar, so the figure shows that the
# gate is a pre-specified rule and not a fixed formula.
#
# Prerequisites -- each source script writes the small draws file consumed here:
#
#   Rscript vignettes/showcase_small_cohort.R --full   # showcase_joint_draws.rds
#   Rscript vignettes/showcase_pbc.R                   # showcase_pbc_joint_draws.rds
#   Rscript vignettes/showcase_joint_composite.R       # -> showcase_joint.png
#
# Companion paper: Otte WM, Kampman SL, [...] A Bayesian framework for
# validating survival prediction models: posterior credible intervals for
# discrimination, calibration, and clinical utility. in prep.
# ===========================================================================

source(file.path("..", "R", "load_all.R"), chdir = TRUE)

figdir <- "figures"

read_gate_draws <- function(file) {
  path <- file.path(figdir, file)
  if (!file.exists(path))
    stop("Missing `", path, "`. Run the two showcase scripts first; see the ",
         "header of this file for the order.", call. = FALSE)
  readRDS(path)
}

sim <- read_gate_draws("showcase_joint_draws.rds")
pbc <- read_gate_draws("showcase_pbc_joint_draws.rds")

## ===========================================================================
## The gate rule lives HERE, once, and is applied to both cohorts' raw draws.
## Each source script stores only its draws and its own thresholds, so the two
## panels cannot disagree on what "passes" means. Recomputing the pass flag from
## the stored thresholds must reproduce the P_gate each script reported; a
## mismatch means the draws and the reported probability came from different
## runs, which is exactly the drift this check exists to catch.
## ===========================================================================
gate_lab <- c("gate passes (also dNB > 0)", "gate fails")

flag_gate <- function(x) {
  g <- x$gate
  out <- tibble::as_tibble(x$draws) |>
    dplyr::mutate(
      cohort = x$label,
      gate = factor(dplyr::if_else(
        C_harrell >= g$C_min & cal_slope >= g$slope_lo &
          cal_slope <= g$slope_hi & dNB > 0,
        gate_lab[1], gate_lab[2]), levels = gate_lab))
  stopifnot(abs(mean(out$gate == gate_lab[1]) - x$P_gate) < 1e-12)
  out
}

lev <- c(sim$label, pbc$label)
joint_df <- dplyr::bind_rows(flag_gate(sim), flag_gate(pbc)) |>
  dplyr::mutate(cohort = factor(cohort, levels = lev)) |>
  # the passing draws are a minority inside the failing cloud; draw them last so
  # they are not painted over by later failing draws.
  dplyr::arrange(gate == gate_lab[1])

## Thresholds per panel, kept data-driven: the slope band happens to be the same
## in both cohorts but the C bar is not (0.65 vs 0.70), and a future change to
## one gate must not silently mislabel the other panel's reference lines.
thr <- dplyr::bind_rows(
  tibble::tibble(cohort = sim$label, C_min = sim$gate$C_min,
                 slope_lo = sim$gate$slope_lo, slope_hi = sim$gate$slope_hi),
  tibble::tibble(cohort = pbc$label, C_min = pbc$gate$C_min,
                 slope_lo = pbc$gate$slope_lo, slope_hi = pbc$gate$slope_hi)) |>
  dplyr::mutate(cohort = factor(cohort, levels = lev))
thr_v <- thr |>
  tidyr::pivot_longer(c(slope_lo, slope_hi), values_to = "x", names_to = NULL)

## ---- Figure --------------------------------------------------------------
## House style: see FIGURE_STYLE.md (Okabe-Ito palette + theme_bayesval()).
## Fixed (shared) scales, the facet_wrap default. The two cohorts sit in
## different regions of the plane -- PBC discriminates better, so its cloud rides
## higher -- and that vertical offset is itself a finding. Free scales would
## re-centre each panel and hide it, and the two C bars would then sit at
## visually identical heights while standing for different numbers.
p_joint <- ggplot2::ggplot(
    joint_df, ggplot2::aes(cal_slope, C_harrell, fill = gate)) +
  ggplot2::geom_vline(data = thr_v, ggplot2::aes(xintercept = x),
                      linetype = 2, colour = "gray40") +
  ggplot2::geom_hline(data = thr, ggplot2::aes(yintercept = C_min),
                      linetype = 2, colour = "gray40") +
  ggplot2::geom_point(shape = 21, colour = "gray40", stroke = 0.3,
                      size = 1.4, alpha = 0.7) +
  ggplot2::facet_wrap(~ cohort) +
  ggplot2::scale_fill_manual(
    values = stats::setNames(bc_palette[c("blue", "orange")], gate_lab)) +
  ggplot2::scale_x_continuous(breaks = number_ticks(5)) +
  ggplot2::scale_y_continuous(breaks = number_ticks(5)) +
  # No panel title: the gate probabilities belong in the caption. The strips name
  # the cohort; sample sizes go in the caption too.
  ggplot2::labs(x = "Calibration slope", y = "Harrell C", fill = NULL) +
  theme_bayesval()
# Authored at the printed width (\linewidth = 6.5in at 1in margins) so the 12pt
# figure text is not shrunk relative to the 12pt body text by \includegraphics.
bc_ggsave(file.path(figdir, "showcase_joint.png"), p_joint,
          width = 6.5, height = 3.8)

cat("=====================================================================\n")
cat("COMPOSITE JOINT-GATE FIGURE\n")
for (x in list(sim, pbc))
  cat(sprintf("  %-24s n = %3d, %3d events, C bar %.2f, P(gate) = %.2f\n",
              x$label, x$n, x$events, x$gate$C_min, x$P_gate))
cat(sprintf("  figure -> vignettes/%s/showcase_joint.png\n", figdir))
cat("=====================================================================\n")
