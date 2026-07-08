# Plot a posterior calibration band from ppc_calibration() (ggplot2).

House style (FIGURE_STYLE.md): Okabe-Ito palette, grey dashed diagonal,
credible band as a filled ribbon under the median line.

## Usage

``` r
plot_calibration_band(ppc, main = NULL)
```

## Arguments

- ppc:

  output of ppc_calibration().

- main:

  plot title; default summarises the horizon and posterior predictive
  p-value.

## Value

a ggplot object.
