# worked_example_stan.R
# ---------------------------------------------------------------------------
# Part D with the PRIMARY, genuine full-Bayes model. Because CRAN rstanarm ships
# no stan_surv, we use the framework's rstan piecewise-exponential PH model
# (stan/pwexp_ph.stan) as the Bayesian analogue of Cox. Train on Rotterdam,
# externally validate on GBSG at 5 years, and compare the MCMC posterior to the
# Cox-Laplace approximation (from worked_example.R) and the frequentist estimate.
#
#   Rscript vignettes/worked_example_stan.R
#
# Companion paper: Otte WM, Kampman SL, [...] A Bayesian framework for
# validating survival prediction models: posterior credible intervals for
# discrimination, calibration, and clinical utility. in prep.
# ---------------------------------------------------------------------------
suppressMessages({library(survival); library(rstan)})
source(file.path("..", "R", "load_all.R"), chdir = TRUE)
rstan::rstan_options(auto_write = TRUE)
YEAR <- 365.25; figdir <- "figures"; dir.create(figdir, showWarnings = FALSE)

prep <- function(df, time, status) {
  size_cat <- if (is.factor(df$size)) df$size else
    cut(df$size, breaks = c(-Inf, 20, 50, Inf), labels = c("<=20","20-50",">50"))
  tibble::tibble(time = df[[time]] / YEAR, status = df[[status]], age = df$age,
                 nodes = df$nodes, grade3 = as.integer(df$grade >= 3),
                 size = factor(size_cat, levels = c("<=20","20-50",">50")),
                 logpgr = log(df$pgr + 1), hormon = df$hormon, meno = df$meno)
}
rot <- dplyr::mutate(rotterdam,
                     rfs = as.integer(recur == 1 | death == 1),
                     rfstime = ifelse(recur == 1, rtime, dtime))
train <- prep(rot, "rfstime", "rfs")  |> tidyr::drop_na()
test  <- prep(gbsg, "rfstime", "status") |> tidyr::drop_na()
horizons <- c(2, 5); h5 <- 2L; tau <- 5
form <- Surv(time, status) ~ age + nodes + grade3 + size + logpgr + hormon + meno

## ---- Genuine Bayesian fit (MCMC) -------------------------------------------
cat("Fitting full-Bayes piecewise-exponential PH (rstan)...\n")
t0 <- Sys.time()
obj <- fit_stan_pwexp(form, train, horizons = horizons, n_pieces = 8L,
                      beta_scale = 2.5, chains = 4L, iter = 1000L, seed = 2026)
cat(sprintf("fitted in %.0f s\n", difftime(Sys.time(), t0, units = "secs")))
print(obj)
stopifnot(obj$diagnostics$max_rhat < 1.05)          # MCMC diagnostics gate

## ---- External validation posterior on GBSG ---------------------------------
ps <- stan_pwexp_predict(obj, test, time = "time", status = "status")
print(ps)

suite <- list(
  harrell = metric_harrell_c,
  uno     = function(d, w) metric_uno_c(d, w, tau = tau),
  gonen   = metric_gonen_k,
  slope   = metric_cal_slope,
  ici     = function(d, w) metric_ici(d, w, horizon = h5, method = "auto"),
  brier   = function(d, w) metric_brier(d, w)[c(paste0("Brier_t", h5), "IBS")])
post <- propagate_many(ps, suite, weights = "bootstrap", seed = 1)
tab <- posterior_table(post) |>
  dplyr::mutate(metric = sub("^[^.]+\\.", "", metric))
cat("\n== Full-Bayes (stan_pwexp) posterior median [95% CrI], GBSG @5y ==\n")
print(as.data.frame(dplyr::select(tab, metric, summary)), row.names = FALSE)

## ---- Compare MCMC vs Cox-Laplace vs frequentist ----------------------------
lp <- colMeans(ps$eta); dfv <- tibble::tibble(time = test$time, status = test$status)
ch <- concordance(Surv(time, status) ~ lp, data = dfv)
cu <- concordance(Surv(time, status) ~ lp, data = dfv, timewt = "n/G2", ymax = tau)
cs <- coxph(Surv(time, status) ~ lp, data = dfv)
fmt <- function(m, l, u) sprintf("%.3f [%.3f, %.3f]", m, l, u)
pick <- function(nm) dplyr::filter(tab, metric == nm)
lap <- tryCatch(readRDS(file.path(figdir, "worked_example_results.rds"))$headline,
                error = function(e) NULL)
cmp <- tibble::tibble(
  metric = c("C_harrell", "C_uno", "cal_slope"),
  full_bayes = c(fmt(pick("C_harrell")$median, pick("C_harrell")$lower, pick("C_harrell")$upper),
                 fmt(pick("C_uno")$median, pick("C_uno")$lower, pick("C_uno")$upper),
                 fmt(pick("cal_slope")$median, pick("cal_slope")$lower, pick("cal_slope")$upper)),
  cox_laplace = if (!is.null(lap)) lap$bayesian[match(c("C_harrell","C_uno","cal_slope"), lap$metric)] else NA,
  frequentist = c(fmt(1 - ch$concordance, 1 - ch$concordance - 1.96*sqrt(ch$var), 1 - ch$concordance + 1.96*sqrt(ch$var)),
                  fmt(1 - cu$concordance, 1 - cu$concordance - 1.96*sqrt(cu$var), 1 - cu$concordance + 1.96*sqrt(cu$var)),
                  fmt(coef(cs)[[1]], coef(cs)[[1]] - 1.96*sqrt(vcov(cs)[1,1]), coef(cs)[[1]] + 1.96*sqrt(vcov(cs)[1,1]))))
cat("\n== Full-Bayes vs Cox-Laplace vs frequentist ==\n")
print(as.data.frame(cmp), row.names = FALSE)

## ---- PPC + decision curve (full-Bayes) -------------------------------------
## Drawn by `worked_example_composite.R`, which faces these curves with the
## Cox-Laplace ones from `worked_example.R`. Storing the objects here means the
## comparison figure never needs a Stan toolchain, nor a second MCMC fit.
ppc <- ppc_calibration(ps, horizon = h5, method = "auto", draws = 200L, seed = 3)
pnb <- posterior_net_benefit(ps, thresholds = seq(0.1, 0.7, 0.05), horizon = h5,
                             weights = "bootstrap", draws = 300L, seed = 4)
saveRDS(list(tab = tab, cmp = cmp, diagnostics = obj$diagnostics, ppc = ppc,
             pnb = pnb),
        file.path(figdir, "worked_example_stan_results.rds"))
cat(sprintf("\nMCMC diagnostics: max Rhat %.3f, min ESS %.0f, divergences %d\n",
            obj$diagnostics$max_rhat, obj$diagnostics$min_ess, obj$diagnostics$divergences))
cat("Done. Figures + results in vignettes/figures/\n")
