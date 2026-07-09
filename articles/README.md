# Pre-rendered articles

Content here is **not** built by `pkgdown::build_site()` and **not** shipped as
a package vignette---it needs an rstan toolchain and several minutes of MCMC,
which is too heavy to run on every site rebuild.

## Workflow

1. Edit the `.Rmd` source in this folder (uses `library(bayesval)`, self-contained).
2. Render it once locally:
   ```r
   rmarkdown::render("articles/small_sample_stan_vs_laplace.Rmd",
                     output_dir = "articles/_rendered")
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
