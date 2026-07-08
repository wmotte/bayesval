# Plot posterior net-benefit bands from posterior_net_benefit() (ggplot2).

House style (FIGURE_STYLE.md): model curve in the primary palette colour
with its credible ribbon, treat-all in grey, treat-none the grey zero
line.

## Usage

``` r
plot_net_benefit(pnb, main = "Posterior decision curve")
```

## Arguments

- pnb:

  output of posterior_net_benefit().

- main:

  plot title.

## Value

a ggplot object.
