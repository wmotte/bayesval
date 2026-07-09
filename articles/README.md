# Pre-rendered articles

Content here is **not** built by `pkgdown::build_site()` and **not** shipped as
a package vignette - it needs an rstan toolchain and several minutes of MCMC,
which is too heavy to run on every site rebuild.

## Workflow

1. Edit the `.Rmd` source in this folder (uses `library(bayesval)`, self-contained).
2. To render it with the correct `pkgdown` layout (so it gets the navbar and styling), temporarily move it to `vignettes/`, build it, and copy the HTML to `articles/_rendered`:
   ```sh
   cp articles/small_sample_stan_vs_laplace.Rmd vignettes/
   Rscript -e "pkgdown::build_article('small_sample_stan_vs_laplace')"
   cp docs/articles/small_sample_stan_vs_laplace.html articles/_rendered/
   rm vignettes/small_sample_stan_vs_laplace.Rmd
   ```
3. After `pkgdown::build_site()` has (re)built `docs/`, copy the rendered HTML
   into `docs/articles/` so it is served alongside the live-built vignettes:
   ```sh
   cp articles/_rendered/small_sample_stan_vs_laplace.html docs/articles/
   ```
   `pkgdown::build_site()` wipes `docs/` before rebuilding, so this copy must
   happen **after** every rebuild, not before.
4. The site navbar links to `articles/small_sample_stan_vs_laplace.html`
   directly (see `_pkgdown.yml`); it is not part of the auto-generated
   "Articles" index.
