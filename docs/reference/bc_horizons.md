# Resolve evaluation horizons from data (or return the configured horizons).

Resolve evaluation horizons from data (or return the configured
horizons).

## Usage

``` r
bc_horizons(time, status, cfg = bc_config())
```

## Arguments

- time, status:

  observed follow-up and event indicator.

- cfg:

  a bc_config() list.

## Value

numeric horizons (event-time quantiles when cfg\$horizons is NULL).
