# Format a posterior_metric (or several) as a median `[95% CrI]` table.

Format a posterior_metric (or several) as a median `[95% CrI]` table.

## Usage

``` r
posterior_table(..., probs = c(0.025, 0.5, 0.975), digits = 3)
```

## Arguments

- ...:

  named posterior_metric objects, or a single named list of them.

- probs:

  credible-interval probabilities.

- digits:

  rounding for the formatted `"median [lo, hi]"` column.

## Value

tibble with columns metric, median, lower, upper, and a pre-formatted
`summary` string.
