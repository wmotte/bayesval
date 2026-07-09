# Showcase: real-world deployment decision on the Mayo PBC cohort

This vignette mirrors
[`vignette("showcase_small_cohort")`](https://wmotte.github.io/bayesval/articles/showcase_small_cohort.md)
but on real clinical data. We use the historic Mayo Clinic Primary
Biliary Cholangitis (PBC) dataset to develop a standard Cox model on the
randomized trial participants (n = 312, 125 events) and externally
validate it on the independent non-randomized eligible patients (n =
106, 36 events).

With only 36 events in the validation cohort, classical frequentist
confidence intervals for separate performance measures are wide, often
leading to ambiguous local-use decisions. The posterior workflow yields
a coherent deployment gate probability instead.

``` r

suppressMessages(library(survival))
library(bayesval)
set.seed(2028)
```

## Prepare data and fit the published model

``` r

data(pbc)
pbc$status <- ifelse(pbc$status == 2, 1, 0)   # overall survival: status==2 (dead) is the event
vars_pred <- c("time", "status", "age", "edema", "bili", "protime", "albumin")
pbc_sub <- pbc[complete.cases(pbc[, vars_pred]), ]

train <- pbc_sub[!is.na(pbc_sub$trt), ]   # development: randomized patients
val   <- pbc_sub[is.na(pbc_sub$trt), ]    # validation: non-randomized patients

cat(sprintf("Development cohort: n = %d, events = %d\n", nrow(train), sum(train$status)))
```

    ## Development cohort: n = 312, events = 125

``` r

cat(sprintf("Validation cohort:  n = %d, events = %d\n", nrow(val), sum(val$status)))
```

    ## Validation cohort:  n = 104, events = 35

``` r

# Published Mayo model (coxph-Laplace backend for speed)
published <- coxph(Surv(time, status) ~ age + edema + log(bili) + log(protime) + log(albumin),
                   data = train, x = TRUE, y = TRUE)
```

## Horizons, deployment gate, and validation suite

``` r

HZ <- c(1096, 1826); H5 <- 2L; TAU <- 1826; PT <- 0.20   # days: 3y and 5y horizons
C_MIN <- 0.70; SLOPE_LO <- 0.80; SLOPE_HI <- 1.25

metric_dnb <- function(pt, h) function(draw, w) {
  nb <- net_benefit_draw(draw, thresholds = pt, horizon = h, w = w)
  c(NB_model = nb$model, dNB = nb$model - max(nb$all, 0))
}
suite <- list(
  harrell = metric_harrell_c,
  uno     = function(d, w) metric_uno_c(d, w, tau = TAU),
  slope   = metric_cal_slope,
  ici     = function(d, w) metric_ici(d, w, horizon = H5, method = "auto")["ICI"],
  ibs     = function(d, w) metric_brier(d, w)["IBS"],
  nb      = metric_dnb(PT, H5))
```

## Frequentist read-out on the validation cohort

``` r

Terms <- stats::delete.response(stats::terms(published))
mf <- stats::model.frame(Terms, val, xlev = published$xlevels)
X  <- stats::model.matrix(Terms, mf, contrasts.arg = published$contrasts)
X  <- X[, names(stats::coef(published)), drop = FALSE]
lp <- as.numeric(X %*% stats::coef(published))

dfv <- tibble::tibble(time = val$time, status = val$status, lp = lp)
wald <- function(est, se) c(est, est - 1.96 * se, est + 1.96 * se)
ch <- concordance(Surv(time, status) ~ lp, data = dfv)
fr_C <- wald(1 - ch$concordance, sqrt(ch$var))
cs <- coxph(Surv(time, status) ~ lp, data = dfv)
fr_slope <- wald(unname(coef(cs)), sqrt(vcov(cs)[1, 1]))

ps <- psurv_external(published, val, time = "time", status = "status",
                     horizons = HZ, S = 1000L, seed = 2028)

mean_draw <- list(surv2d = matrix(colMeans(ps$surv), nrow = ps$n), eta = lp,
                  time = val$time, status = val$status, horizons = HZ,
                  n = ps$n, H = length(HZ))
boot_ci <- function(fn, B = 400L, seed = 11) {
  set.seed(seed); n <- ps$n
  est <- as.numeric(fn(mean_draw, rep(1 / n, n)))
  v <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    db <- list(surv2d = mean_draw$surv2d[idx, , drop = FALSE], eta = lp[idx],
               time = val$time[idx], status = val$status[idx],
               horizons = HZ, n = n, H = length(HZ))
    as.numeric(fn(db, rep(1 / n, n)))
  })
  c(est, stats::quantile(v, c(0.025, 0.975), names = FALSE))
}
fr_ici <- boot_ci(function(d, w) metric_ici(d, w, horizon = H5, method = "auto")["ICI"])
fr_dnb <- boot_ci(function(d, w) metric_dnb(PT, H5)(d, w)["dNB"])

frtab <- rbind(`Harrell C` = fr_C, `Calibration slope` = fr_slope,
              `ICI (5y)` = fr_ici, `dNB @20% (vs default)` = fr_dnb)
colnames(frtab) <- c("est", "lo", "hi")
knitr::kable(round(frtab, 3), caption = "Frequentist point estimate [95% CI]")
```

|                       |   est |    lo |    hi |
|:----------------------|------:|------:|------:|
| Harrell C             | 0.786 | 0.714 | 0.857 |
| Calibration slope     | 0.857 | 0.554 | 1.159 |
| ICI (5y)              | 0.117 | 0.092 | 0.450 |
| dNB @20% (vs default) | 0.028 | 0.011 | 0.053 |

Frequentist point estimate \[95% CI\] {.table}

``` r

verdict <- function(lo, hi, need_lo = NA, need_hi = NA) {
  ok_lo <- if (is.na(need_lo)) TRUE else lo >= need_lo
  ok_hi <- if (is.na(need_hi)) TRUE else hi <= need_hi
  bad_lo <- if (is.na(need_lo)) FALSE else hi < need_lo
  bad_hi <- if (is.na(need_hi)) FALSE else lo > need_hi
  if (ok_lo && ok_hi) "PASS" else if (bad_lo || bad_hi) "FAIL" else "AMBIGUOUS"
}
cat(sprintf("discrimination C>=%.2f: %s\n", C_MIN, verdict(fr_C[2], fr_C[3], need_lo = C_MIN)))
```

    ## discrimination C>=0.70: PASS

``` r

cat(sprintf("calibration slope in [%.2f,%.2f]: %s\n", SLOPE_LO, SLOPE_HI,
            verdict(fr_slope[2], fr_slope[3], need_lo = SLOPE_LO, need_hi = SLOPE_HI)))
```

    ## calibration slope in [0.80,1.25]: AMBIGUOUS

``` r

cat(sprintf("net benefit dNB>0: %s\n", verdict(fr_dnb[2], fr_dnb[3], need_lo = 0)))
```

    ## net benefit dNB>0: PASS

## Bayesian read-out: posteriors and deployment probabilities

``` r

J <- propagate_many(ps, suite, weights = "bootstrap", seed = 1)
Jm <- unclass(J); colnames(Jm) <- sub("^[^.]+\\.", "", colnames(Jm))

qi <- function(x) stats::quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
show <- c(C_harrell = "Harrell C", C_uno = "Uno C", cal_slope = "Calibration slope",
          ICI = "ICI (5y)", IBS = "IBS", NB_model = "Net benefit @20%",
          dNB = "dNB @20% (vs default)")
post_tab <- purrr::map_dfr(names(show), function(col) {
  q <- qi(Jm[, col])
  tibble::tibble(metric = show[[col]], median = q[2], lo = q[1], hi = q[3])
})
knitr::kable(post_tab, digits = 3, caption = "Posterior median [95% CrI]")
```

| metric                | median |     lo |    hi |
|:----------------------|-------:|-------:|------:|
| Harrell C             |  0.786 |  0.696 | 0.850 |
| Uno C                 |  0.786 |  0.698 | 0.854 |
| Calibration slope     |  0.846 |  0.541 | 1.184 |
| ICI (5y)              |  0.274 |  0.036 | 0.691 |
| IBS                   |  0.216 |  0.112 | 0.736 |
| Net benefit @20%      |  0.131 |  0.000 | 0.292 |
| dNB @20% (vs default) |  0.000 | -0.222 | 0.092 |

Posterior median \[95% CrI\] {.table}

``` r

P_disc <- mean(Jm[, "C_harrell"] >= C_MIN)
P_cal  <- mean(Jm[, "cal_slope"] >= SLOPE_LO & Jm[, "cal_slope"] <= SLOPE_HI)
P_ici  <- mean(Jm[, "ICI"] <= 0.05)
P_bene <- mean(Jm[, "dNB"] > 0)
P_harm <- mean(Jm[, "dNB"] < 0)
P_gate <- mean(Jm[, "C_harrell"] >= C_MIN &
               Jm[, "cal_slope"] >= SLOPE_LO & Jm[, "cal_slope"] <= SLOPE_HI &
               Jm[, "dNB"] > 0)
cat(sprintf("P(discrimination adequate) = %.2f\n", P_disc))
```

    ## P(discrimination adequate) = 0.97

``` r

cat(sprintf("P(calibration acceptable)  = %.2f\n", P_cal))
```

    ## P(calibration acceptable)  = 0.60

``` r

cat(sprintf("P(ICI <= 0.05)              = %.2f\n", P_ici))
```

    ## P(ICI <= 0.05)              = 0.05

``` r

cat(sprintf("P(net clinical benefit)    = %.2f\n", P_bene))
```

    ## P(net clinical benefit)    = 0.43

``` r

cat(sprintf("P(net harm)                = %.2f\n", P_harm))
```

    ## P(net harm)                = 0.54

``` r

cat(sprintf("P(JOINT deployment gate)   = %.2f\n", P_gate))
```

    ## P(JOINT deployment gate)   = 0.27

Even in this small, real validation cohort, the coherent joint posterior
turns several individually ambiguous frequentist intervals into a single
deployment probability.

### Citation

`bayesval` is the software companion to:

> Otte WM, Kampman SL, \[…\] *A Bayesian framework for validating
> survival prediction models: posterior credible intervals for
> discrimination, calibration, and clinical utility.* in prep.
