# Antolini's time-dependent concordance (per draw) — A1.4. Uses the model's time-specific survival: a usable pair (i event at Ytilde_i, j later) is concordant when S_i(Ytilde_i) \< S_j(Ytilde_i), i.e. the subject who failed had the lower predicted survival AT THAT TIME. For a PH model this reduces to Harrell's C; for flexible / non-PH models whose curves may cross it is the natural single summary. Survival at the (arbitrary) event times is obtained by interpolating each subject's cumulative hazard across the horizon knots.

Antolini's time-dependent concordance (per draw) — A1.4. Uses the
model's time-specific survival: a usable pair (i event at Ytilde_i, j
later) is concordant when S_i(Ytilde_i) \< S_j(Ytilde_i), i.e. the
subject who failed had the lower predicted survival AT THAT TIME. For a
PH model this reduces to Harrell's C; for flexible / non-PH models whose
curves may cross it is the natural single summary. Survival at the
(arbitrary) event times is obtained by interpolating each subject's
cumulative hazard across the horizon knots.

## Usage

``` r
metric_antolini_c(draw, w = NULL, max_n = 2000L)
```

## Arguments

- draw:

  a
  [`psurv_draw()`](https://wmotte.github.io/bayesval/reference/psurv_draw.md)
  list for a single posterior draw.

- w:

  length-n case-weight vector (sums to 1); `NULL` uses uniform weights.

- max_n:

  guard: subsample the pair set above this n (kept tractable in the
  posterior loop; noted, not silent — see the returned NA path).

## Value

c(C_antolini = value)
