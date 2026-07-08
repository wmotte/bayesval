# globals.R
# ---------------------------------------------------------------------------
# Declares NSE column names used unquoted inside dplyr/tidyr verbs, so R CMD
# check's "no visible binding for global variable" NOTE does not fire. These
# are tidy-eval data-mask symbols, not real globals.
# ---------------------------------------------------------------------------
utils::globalVariables(c(
  "beta_scale", "estimate", "lower", "lower_bayes", "lower_freq", "median",
  "metric", "nb", "pred", "prior", "strategy", "threshold", "upper",
  "upper_bayes", "upper_freq", "variable"
))
