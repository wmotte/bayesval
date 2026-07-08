# Convergence gate: TRUE only if Rhat \< max_rhat, ESS \>= min_ess, no divergences.

Metrics computed on a fit that fails this gate must not be reported
(plan "MCMC diagnostics gate"). For a non-MCMC backend the gate passes
vacuously with a note.

## Usage

``` r
check_mcmc_gate(fit, max_rhat = 1.01, min_ess = 400, error = FALSE)

mcmc_gate(fit, max_rhat = 1.01, min_ess = 400, error = FALSE)
```

## Arguments

- fit:

  as in mcmc_diagnostics().

- max_rhat:

  maximum acceptable split-Rhat (default 1.01).

- min_ess:

  minimum acceptable bulk & tail ESS (default 400).

- error:

  if TRUE, stop() on failure instead of returning FALSE.

## Value

logical (invisibly), with attribute "diagnostics".
