# showcase_pbc.R
# ===========================================================================
# REAL-WORLD SHOWCASE: The Mayo Clinic Primary Biliary Cholangitis (PBC) cohort
#
# Setting: This script mirrors `showcase_small_cohort.R` but uses real clinical
# data. We use the historic Mayo PBC dataset to develop a standard Cox model
# on the randomized trial participants (n = 312, 125 events) and externally
# validate it on the independent non-randomized eligible patients (n = 104,
# 35 events).
#
# With only 35 events in the validation cohort, the classical frequentist
# confidence intervals for validation metrics are wide, often leading to
# ambiguous deployment decisions. The Bayesian framework yields a coherent
# deployment gate probability.
#
#   Rscript vignettes/showcase_pbc.R
#
# Companion paper: Otte WM, Kampman SL, [...] A Bayesian framework for
# validating survival prediction models: posterior credible intervals for
# discrimination, calibration, and clinical utility. in prep.
# ===========================================================================

suppressMessages(library(survival))
source(file.path("..", "R", "load_all.R"), chdir = TRUE)

figdir <- "figures"
dir.create(figdir, showWarnings = FALSE)
set.seed(2028)

## ===========================================================================
## 1. Prepare data and fit published model
## ===========================================================================
# PBC dataset from the survival package
data(pbc)
# Status: 0=censored, 1=transplant, 2=dead. We analyze overall survival (status==2 is event)
pbc$status <- ifelse(pbc$status == 2, 1, 0)
# Use the classical predictors
vars_pred <- c("time", "status", "age", "edema", "bili", "protime", "albumin")
pbc_sub <- pbc[complete.cases(pbc[, vars_pred]), ]

# Development cohort: randomized patients (trt is not NA)
train <- pbc_sub[!is.na(pbc_sub$trt), ]
# Validation cohort: non-randomized patients
val <- pbc_sub[is.na(pbc_sub$trt), ]

cat("=====================================================================\n")
cat("REAL-WORLD SHOWCASE: Mayo PBC Dataset\n")
cat(sprintf("Development cohort: n = %d, events = %d\n", nrow(train), sum(train$status)))
cat(sprintf("Validation cohort:  n = %d, events = %d\n", nrow(val), sum(val$status)))
cat("=====================================================================\n")

# Fit the published Mayo model
# (Note: we use the coxph Laplace backend for speed)
published <- coxph(Surv(time, status) ~ age + edema + log(bili) + log(protime) + log(albumin),
                   data = train, x = TRUE, y = TRUE)

## ===========================================================================
## 2. Define horizons, gate, and metrics
## ===========================================================================
# Horizons in days (3 years and 5 years)
HZ <- c(1096, 1826)
H5 <- 2L
TAU <- 1826
PT <- 0.20 # 20% risk threshold for decision

C_MIN <- 0.70
SLOPE_LO <- 0.80
SLOPE_HI <- 1.25

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
  nb      = metric_dnb(PT, H5)
)

## ===========================================================================
## 3. Frequentist read-out
## ===========================================================================
# Apply the model to the validation cohort
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

# Bootstrap for ICI and dNB
# Need psurv object first to get the mean draw baseline hazard
ps <- psurv_external(published, val, time = "time", status = "status",
                     horizons = HZ, S = 1000L, seed = 2028)

mean_draw <- list(surv2d = matrix(colMeans(ps$surv), nrow = ps$n),
                  eta = lp, time = val$time, status = val$status,
                  horizons = HZ, n = ps$n, H = length(HZ))
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

cat("\n-- Frequentist read-out (point estimate [95% CI]) -------------------\n")
frtab <- rbind(
  `Harrell C`            = fr_C,
  `Calibration slope`    = fr_slope,
  `ICI (5y)`             = fr_ici,
  `dNB @20% (vs default)`= fr_dnb)
colnames(frtab) <- c("est", "lo", "hi")
print(round(frtab, 3))

verdict <- function(lo, hi, need_lo = NA, need_hi = NA) {
  ok_lo <- if (is.na(need_lo)) TRUE else lo >= need_lo
  ok_hi <- if (is.na(need_hi)) TRUE else hi <= need_hi
  bad_lo <- if (is.na(need_lo)) FALSE else hi < need_lo
  bad_hi <- if (is.na(need_hi)) FALSE else lo > need_hi
  if (ok_lo && ok_hi) "PASS" else if (bad_lo || bad_hi) "FAIL" else "AMBIGUOUS"
}
cat(sprintf("  discrimination C>=%.2f : %s\n", C_MIN,
            verdict(fr_C[2], fr_C[3], need_lo = C_MIN)))
cat(sprintf("  calibration slope in [%.2f,%.2f] : %s\n", SLOPE_LO, SLOPE_HI,
            verdict(fr_slope[2], fr_slope[3], need_lo = SLOPE_LO, need_hi = SLOPE_HI)))
cat(sprintf("  net benefit dNB>0 : %s\n",
            verdict(fr_dnb[2], fr_dnb[3], need_lo = 0)))

## ===========================================================================
## 4. Bayesian read-out
## ===========================================================================
J <- propagate_many(ps, suite, weights = "bootstrap", seed = 1)
Jm <- unclass(J); colnames(Jm) <- sub("^[^.]+\\.", "", colnames(Jm))
post_summary <- summary(J) |>
  dplyr::mutate(metric = sub("^[^.]+\\.", "", metric))

qi <- function(x) stats::quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
fmt <- function(v, d = 3) sprintf(paste0("%.", d, "f"), v)

cat("\n-- Bayesian read-out (posterior median [95% CrI]) -------------------\n")
show <- c(C_harrell = "Harrell C", C_uno = "Uno C", cal_slope = "Calibration slope",
          ICI = "ICI (5y)", IBS = "IBS", NB_model = "Net benefit @20%",
          dNB = "dNB @20% (vs default)")
for (col in names(show)) {
  q <- qi(Jm[, col])
  cat(sprintf("  %-22s %s [%s, %s]\n", show[[col]],
              fmt(q[2]), fmt(q[1]), fmt(q[3])))
}

P_disc  <- mean(Jm[, "C_harrell"] >= C_MIN)
P_cal   <- mean(Jm[, "cal_slope"] >= SLOPE_LO & Jm[, "cal_slope"] <= SLOPE_HI)
P_ici   <- mean(Jm[, "ICI"] <= 0.05)
P_bene  <- mean(Jm[, "dNB"] > 0)
P_harm  <- mean(Jm[, "dNB"] < 0)
P_gate  <- mean(Jm[, "C_harrell"] >= C_MIN &
                Jm[, "cal_slope"] >= SLOPE_LO & Jm[, "cal_slope"] <= SLOPE_HI &
                Jm[, "dNB"] > 0)
# What an analyst would get by multiplying the three marginal probabilities, as
# if the criteria were independent. It misses the joint gate in EITHER direction
# depending on the sign of the dependence, which is the point: the product is not
# a conservative approximation of the joint event, it is the wrong quantity.
P_gate_indep <- P_disc * P_cal * P_bene

cat("\n-- Deployment probabilities -----------------------------------------\n")
cat(sprintf("  P(discrimination adequate, C>=%.2f)         = %.2f\n", C_MIN, P_disc))
cat(sprintf("  P(calibration acceptable, slope in[%.2f,%.2f]) = %.2f\n",
            SLOPE_LO, SLOPE_HI, P_cal))
cat(sprintf("  P(ICI <= 0.05)                             = %.2f\n", P_ici))
cat(sprintf("  P(net clinical benefit, dNB>0)             = %.2f\n", P_bene))
cat(sprintf("  P(net HARM, dNB<0)                         = %.2f\n", P_harm))
cat(sprintf("  --------------------------------------------------------\n"))
cat(sprintf("  P(JOINT deployment gate passes)            = %.2f\n", P_gate))
cat(sprintf("  naive product of the three marginals        = %.2f\n", P_gate_indep))

## ---- Joint posterior of the two gate-critical metrics --------------------
## Consumed by `showcase_joint_composite.R`, which faces this cohort with the
## simulated showcase. Store raw draws and thresholds only; the gate rule is
## applied once, there, to both cohorts.
saveRDS(list(
  draws  = tibble::tibble(cal_slope = Jm[, "cal_slope"],
                          C_harrell = Jm[, "C_harrell"],
                          dNB       = Jm[, "dNB"]),
  gate   = list(C_min = C_MIN, slope_lo = SLOPE_LO, slope_hi = SLOPE_HI),
  P_gate = P_gate,
  label  = "PBC external validation",
  n      = nrow(val),
  events = sum(val$status)),
  file.path(figdir, "showcase_pbc_joint_draws.rds"))

saveRDS(list(freq = frtab, post = post_summary,
             P = c(disc = P_disc, cal = P_cal, ici = P_ici, bene = P_bene,
                   harm = P_harm, gate = P_gate, gate_indep = P_gate_indep)),
        file.path(figdir, "showcase_pbc_results.rds"))

cat("\n=====================================================================\n")
cat("PBC SHOWCASE COMPLETE\n")
cat("=====================================================================\n")
