# Pre-rendered articles

These take a long time to run (several minutes each to generate MCMC samples, requiring an rstan toolchain) and so are **not** built by `pkgdown::build_site()` with the rest of the site, and **not** shipped as package vignettes.

## Workflow

1. Edit the relevant `.Rmd` file or files here. All contain all necessary code, including `library(bayesval)`.
2. Temporarily move to `vignettes/`, build, and copy the resulting .html file to `articles/_rendered` (to inherit the `pkgdown` style, top bar, etc.):
   ```sh
   cp articles/small_sample_stan_vs_laplace.Rmd vignettes/
   Rscript -e "pkgdown::build_article('small_sample_stan_vs_laplace')"
   cp docs/articles/small_sample_stan_vs_laplace.html articles/_rendered/
   rm vignettes/small_sample_stan_vs_laplace.Rmd
   ```
3. Once `pkgdown::build_site()` has (re)built the remaining pages under `docs/`, copy the relevant .html files into `docs/articles/`, alongside the live-built vignettes:
   ```sh
   cp articles/_rendered/small_sample_stan_vs_laplace.html docs/articles/
   ```
   `pkgdown::build_site()` wipes `docs/` before rebuilding, so this copy must
   happen **after** every rebuild, not before.
4. The top bar links directly to `articles/small_sample_stan_vs_laplace.html` (see `_pkgdown.yml`); it does not appear under the auto-generated list of Articles.
