# Side-by-side Bayesian vs frequentist headline table (Part D).

Side-by-side Bayesian vs frequentist headline table (Part D).

## Usage

``` r
posterior_vs_frequentist(bayes, freq, digits = 3)
```

## Arguments

- bayes:

  a posterior_table() (or a posterior_metric) with the same metric names
  as `freq`.

- freq:

  data.frame with columns metric, estimate, lower, upper.

- digits:

  rounding.
