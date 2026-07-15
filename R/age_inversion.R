# =============================================================================
# Stage 2: two-acid Bayesian age inversion for a single sample
# =============================================================================

# Log-likelihood for two-acid age inversion.
# Converts a sample's observed Asp and Glu D/L into two correlated point
# estimates of log(age) via the Stage 1 calibration curves (predict_age_asp,
# predict_age_glu), then evaluates how consistent a candidate age is with
# both, using the same residual covariance (sigma_Asp, sigma_Glu, rho)
# estimated during calibration. This is the Bayesian analogue of a GLS
# combination of Asp and Glu ages: instead of algebraically combining two
# point estimates and variances (with an ad hoc floor to keep GLS weights
# positive), the two acids simply contribute correlated evidence to one joint
# likelihood, which combines coherently for any value of rho.
# calib_params: list with elements a_Asp, e_Asp, sigma_Asp, a_Glu, e_Glu,
# sigma_Glu, rho (see summarize_calibration()).
log_lik_two_acid_age <- function(log_age, calib_params,
                                 DL_Asp_obs, DL_Glu_obs) {
  log_t_Asp_hat <- log(predict_age_asp(DL_Asp_obs, calib_params$a_Asp,
                                       calib_params$e_Asp))
  log_t_Glu_hat <- log(predict_age_glu(DL_Glu_obs, calib_params$a_Glu,
                                       calib_params$e_Glu))
  if (!is.finite(log_t_Asp_hat) || !is.finite(log_t_Glu_hat)) return(-Inf)

  resid_Asp <- log_t_Asp_hat - log_age
  resid_Glu <- log_t_Glu_hat - log_age

  .dbvnorm_log(resid_Asp, resid_Glu,
              calib_params$sigma_Asp, calib_params$sigma_Glu,
              calib_params$rho)
}


# Log-prior for sample age: a weak, bounded-uniform prior on log(age), wide
# enough to cover the full calibration range (default c(0.5, 3000) ka)
# without favouring any age within it. Because this is an ordinary Bayesian
# prior, it can be tightened using independent information (e.g. a
# stratigraphic age range, or a known minimum/maximum age) simply by
# narrowing range_ka or swapping in an informative density.
log_prior_age <- function(log_age, range_ka = c(0.5, 3000)) {
  if (log_age < log(range_ka[1]) || log_age > log(range_ka[2])) return(-Inf)
  0
}


# Log-posterior for two-acid age inversion (Stage 2).
# params: named numeric vector with element log_age.
log_posterior_two_acid_age <- function(params, calib_params,
                                       DL_Asp_obs, DL_Glu_obs,
                                       range_ka = c(0.5, 3000)) {
  log_age <- params[["log_age"]]
  lp <- log_prior_age(log_age, range_ka = range_ka)
  if (!is.finite(lp)) return(-Inf)

  lp + log_lik_two_acid_age(log_age, calib_params, DL_Asp_obs, DL_Glu_obs)
}


# Run the two-acid Bayesian age inversion for one sample.
# Convenience wrapper around run_mcmc() for log_posterior_two_acid_age().
# Initialises the chain at the average of the two single-acid point-age
# estimates.
invert_two_acid_age <- function(DL_Asp_obs, DL_Glu_obs, calib_params,
                                n_iter = 20000, proposal_sd = 0.3,
                                range_ka = c(0.5, 3000)) {
  t_Asp <- predict_age_asp(DL_Asp_obs, calib_params$a_Asp, calib_params$e_Asp)
  t_Glu <- predict_age_glu(DL_Glu_obs, calib_params$a_Glu, calib_params$e_Glu)
  init  <- c(log_age = mean(log(c(t_Asp, t_Glu))))

  run_mcmc(log_posterior_two_acid_age, init = init, n_iter = n_iter,
           proposal_sd = c(log_age = proposal_sd),
           calib_params = calib_params,
           DL_Asp_obs = DL_Asp_obs, DL_Glu_obs = DL_Glu_obs,
           range_ka = range_ka)
}


# Posterior summary for a two-acid age inversion.
# mcmc_out: list returned by invert_two_acid_age().
# Returns a one-row data frame with columns median_ka, lo50_ka, hi50_ka,
# lo66_ka, hi66_ka, lo95_ka, hi95_ka.
summarize_age_posterior <- function(mcmc_out, burnin = 1000) {
  age_ka <- exp(utils::tail(mcmc_out$samples[, "log_age"],
                            nrow(mcmc_out$samples) - burnin))
  data.frame(
    median_ka = stats::median(age_ka),
    lo50_ka   = stats::quantile(age_ka, 0.25),
    hi50_ka   = stats::quantile(age_ka, 0.75),
    lo66_ka   = stats::quantile(age_ka, 0.17),
    hi66_ka   = stats::quantile(age_ka, 0.83),
    lo95_ka   = stats::quantile(age_ka, 0.025),
    hi95_ka   = stats::quantile(age_ka, 0.975),
    row.names = NULL
  )
}


# Plot the posterior age distribution for a two-acid age inversion.
# true_age_ka: optional reference vertical line (e.g. an independently known
# age, for validation).
# compare_ages_ka: optional named numeric vector of point age estimates to
# compare against (e.g. single-acid or GLS-combined ages), drawn as dashed
# vertical lines.
plot_age_posterior <- function(mcmc_out, burnin = 1000, true_age_ka = NULL,
                               compare_ages_ka = NULL) {
  age_ka <- exp(utils::tail(mcmc_out$samples[, "log_age"],
                            nrow(mcmc_out$samples) - burnin))
  df <- data.frame(age_ka = age_ka)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$age_ka)) +
    ggplot2::geom_density(fill = "steelblue", alpha = 0.4) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Age (ka BP)", y = "Posterior density",
                  title = "Two-acid Bayesian age posterior") +
    ggplot2::theme_bw()

  if (!is.null(true_age_ka)) {
    p <- p + ggplot2::geom_vline(xintercept = true_age_ka, colour = "red",
                                 linetype = "dashed", linewidth = 0.8)
  }
  if (!is.null(compare_ages_ka)) {
    ref <- data.frame(age_ka = as.numeric(compare_ages_ka),
                      label  = names(compare_ages_ka))
    p <- p + ggplot2::geom_vline(data = ref,
                                 ggplot2::aes(xintercept = .data$age_ka),
                                 colour = "black", linetype = "dotted",
                                 linewidth = 0.6)
  }
  p
}
