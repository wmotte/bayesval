// pwexp_ph.stan
// Full-Bayes piecewise-exponential proportional-hazards survival model — a
// genuine Bayesian analogue of Cox (B1). The baseline hazard is piecewise
// constant on K pieces with log-hazards log_h0; a first-order random-walk prior
// smooths them (the Bayesian counterpart of a spline/M-spline baseline). No
// person-time row expansion: the likelihood uses each subject's exposure in
// each piece (matrix E) and the piece containing their event (kev).
data {
  int<lower=1> N;                       // subjects
  int<lower=1> P;                       // covariates
  int<lower=1> K;                       // baseline-hazard pieces
  matrix[N, P] X;                       // covariates (no intercept)
  matrix[N, K] E;                       // E[i,k] = time subject i spent in piece k
  array[N] int<lower=0, upper=1> d;     // event indicator
  array[N] int<lower=1, upper=K> kev;   // piece containing the event time (events only)
  real<lower=0> beta_scale;             // prior SD on log-HRs (e.g. 2.5)
}
parameters {
  vector[P] beta;                       // log hazard ratios
  vector[K] log_h0;                     // log baseline hazard per piece
  real<lower=0> sigma_rw;               // random-walk smoothing SD
}
model {
  beta ~ normal(0, beta_scale);         // weakly-informative HR prior
  log_h0[1] ~ normal(-2, 3);            // weakly-informative level
  for (k in 2:K) log_h0[k] ~ normal(log_h0[k - 1], sigma_rw);  // RW1 smoothness
  sigma_rw ~ normal(0, 1);              // half-normal (sigma_rw >= 0)

  vector[N] eta = X * beta;
  vector[N] cumhaz = (E * exp(log_h0)) .* exp(eta);   // H_i(t_i) = H0(t_i) e^{eta_i}
  target += -cumhaz;                                   // survival to t_i
  for (n in 1:N)
    if (d[n] == 1) target += log_h0[kev[n]] + eta[n];  // hazard at the event
}
