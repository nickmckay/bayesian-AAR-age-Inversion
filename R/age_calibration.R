# =============================================================================
# Two-acid AAR age calibration: forward transforms + Stage 1 Bayesian fit
# =============================================================================

# Predictor transform for the time-dependent kinetics (TDK) model
# f(D/L) = atanh(D/L), used for aspartic acid (Asp).
f_tdk <- function(DL) atanh(DL)

# Predictor transform for the simple power-law kinetics (SPK) model
# f(D/L) = D/L, used for glutamic acid (Glu).
f_spk <- function(DL) DL

# Predict age from D/L using the TDK (Asp) calibration curve
# t = a * atanh(D/L)^e
predict_age_asp <- function(DL, a, e) a * f_tdk(DL)^e

# Predict age from D/L using the SPK (Glu) calibration curve
# t = a * (D/L)^e
predict_age_glu <- function(DL, a, e) a * f_spk(DL)^e

# Predict D/L from age using the TDK (Asp) calibration curve (inverse of predict_age_asp)
predict_DL_asp <- function(age_ka, a, e) tanh((age_ka / a)^(1 / e))

# Predict D/L from age using the SPK (Glu) calibration curve (inverse of predict_age_glu)
predict_DL_glu <- function(age_ka, a, e) (age_ka / a)^(1 / e)


# Bivariate normal log-density with zero means (avoids an mvtnorm dependency)
.dbvnorm_log <- function(x1, x2, sigma1, sigma2, rho) {
  if (sigma1 <= 0 || sigma2 <= 0 || abs(rho) >= 1) return(-Inf)
  z1  <- x1 / sigma1
  z2  <- x2 / sigma2
  q   <- (z1^2 - 2 * rho * z1 * z2 + z2^2) / (1 - rho^2)
  -log(2 * pi * sigma1 * sigma2 * sqrt(1 - rho^2)) - 0.5 * q
}


# =============================================================================
# Stage 1: Bayesian calibration fit (Asp + Glu jointly, correlated residuals)
# =============================================================================

# Log-likelihood for the joint two-acid age calibration.
# Models log(t) for each acid as a linear function of log(f(D/L)) (TDK for
# Asp, SPK for Glu), with residuals for the two acids drawn from a correlated
# bivariate normal at each calibration sample (both acids are measured on the
# same dated sample, so their departures from the mean trend are expected to
# covary).
# calibration_data: data frame with columns age_ka, DL_Asp, DL_Glu.
log_lik_age_calibration <- function(log_a_Asp, e_Asp, log_sigma_Asp,
                                    log_a_Glu, e_Glu, log_sigma_Glu,
                                    z_rho, calibration_data) {
  if (e_Asp <= 0 || e_Glu <= 0) return(-Inf)

  a_Asp <- exp(log_a_Asp); sigma_Asp <- exp(log_sigma_Asp)
  a_Glu <- exp(log_a_Glu); sigma_Glu <- exp(log_sigma_Glu)
  rho   <- tanh(z_rho)

  mu_Asp <- log_a_Asp + e_Asp * log(f_tdk(calibration_data$DL_Asp))
  mu_Glu <- log_a_Glu + e_Glu * log(f_spk(calibration_data$DL_Glu))

  resid_Asp <- log(calibration_data$age_ka) - mu_Asp
  resid_Glu <- log(calibration_data$age_ka) - mu_Glu

  if (any(!is.finite(c(resid_Asp, resid_Glu)))) return(-Inf)

  sum(.dbvnorm_log(resid_Asp, resid_Glu, sigma_Asp, sigma_Glu, rho))
}


# Log-prior for the joint two-acid age calibration. Weakly informative
# priors centred on literature-scale values, wide enough that a calibration
# dataset of ~100+ samples dominates the posterior.
log_prior_age_calibration <- function(log_a_Asp, e_Asp, log_sigma_Asp,
                                      log_a_Glu, e_Glu, log_sigma_Glu,
                                      z_rho) {
  if (e_Asp <= 0 || e_Glu <= 0) return(-Inf)
  stats::dnorm(log_a_Asp,    log(3000), 3,   log = TRUE) +
  stats::dnorm(e_Asp,        2.5,       2,   log = TRUE) +
  stats::dnorm(log_sigma_Asp, log(0.4), 1.5, log = TRUE) +
  stats::dnorm(log_a_Glu,    log(9000), 3,   log = TRUE) +
  stats::dnorm(e_Glu,        2.2,       2,   log = TRUE) +
  stats::dnorm(log_sigma_Glu, log(0.4), 1.5, log = TRUE) +
  stats::dnorm(z_rho,        0,         1,   log = TRUE)
}


# Log-posterior for the joint two-acid age calibration (Stage 1).
# params: named numeric vector log_a_Asp, e_Asp, log_sigma_Asp, log_a_Glu,
# e_Glu, log_sigma_Glu, z_rho.
log_posterior_age_calibration <- function(params, calibration_data) {
  lp <- log_prior_age_calibration(
    log_a_Asp     = params[["log_a_Asp"]],
    e_Asp         = params[["e_Asp"]],
    log_sigma_Asp = params[["log_sigma_Asp"]],
    log_a_Glu     = params[["log_a_Glu"]],
    e_Glu         = params[["e_Glu"]],
    log_sigma_Glu = params[["log_sigma_Glu"]],
    z_rho         = params[["z_rho"]]
  )
  if (!is.finite(lp)) return(-Inf)

  lp + log_lik_age_calibration(
    log_a_Asp     = params[["log_a_Asp"]],
    e_Asp         = params[["e_Asp"]],
    log_sigma_Asp = params[["log_sigma_Asp"]],
    log_a_Glu     = params[["log_a_Glu"]],
    e_Glu         = params[["e_Glu"]],
    log_sigma_Glu = params[["log_sigma_Glu"]],
    z_rho         = params[["z_rho"]],
    calibration_data = calibration_data
  )
}


# Extract posterior-mean calibration parameters from a Stage 1 MCMC run.
# Summarises a run_mcmc() calibration chain into a single named list of
# calibration constants (posterior means after burn-in), on the natural
# parameter scale, ready to pass to log_posterior_two_acid_age() for Stage 2
# inversion.
summarize_calibration <- function(mcmc_out, burnin = 1000) {
  post <- utils::tail(mcmc_out$samples, nrow(mcmc_out$samples) - burnin)
  m    <- colMeans(post)
  list(
    a_Asp     = exp(m[["log_a_Asp"]]),
    e_Asp     = m[["e_Asp"]],
    sigma_Asp = exp(m[["log_sigma_Asp"]]),
    a_Glu     = exp(m[["log_a_Glu"]]),
    e_Glu     = m[["e_Glu"]],
    sigma_Glu = exp(m[["log_sigma_Glu"]]),
    rho       = tanh(m[["z_rho"]])
  )
}


# Plot fitted Asp and Glu calibration curves against calibration data.
# Posterior-mean TDK (Asp) and SPK (Glu) curves with 66% prediction bands
# (using posterior residual SDs), overlaid on the calibration D/L-age data.
plot_age_calibration_fit <- function(mcmc_out, calibration_data,
                                     burnin = 1000) {
  cp <- summarize_calibration(mcmc_out, burnin = burnin)

  dl_seq <- seq(0.01, 0.99, length.out = 200)

  curve_asp <- data.frame(
    acid   = "Aspartic acid (TDK)",
    DL     = dl_seq,
    age_ka = predict_age_asp(dl_seq, cp$a_Asp, cp$e_Asp)
  )
  curve_glu <- data.frame(
    acid   = "Glutamic acid (SPK)",
    DL     = dl_seq,
    age_ka = predict_age_glu(dl_seq, cp$a_Glu, cp$e_Glu)
  )
  curves <- rbind(curve_asp, curve_glu)
  # 66% prediction band in log-age space (t_crit for the "likely" band)
  t_crit <- 0.9578
  curves$lo66 <- curves$age_ka * exp(-t_crit * ifelse(
    curves$acid == "Aspartic acid (TDK)", cp$sigma_Asp, cp$sigma_Glu))
  curves$hi66 <- curves$age_ka * exp(+t_crit * ifelse(
    curves$acid == "Aspartic acid (TDK)", cp$sigma_Asp, cp$sigma_Glu))

  obs <- rbind(
    data.frame(acid = "Aspartic acid (TDK)",
               DL = calibration_data$DL_Asp, age_ka = calibration_data$age_ka),
    data.frame(acid = "Glutamic acid (SPK)",
               DL = calibration_data$DL_Glu, age_ka = calibration_data$age_ka)
  )

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = curves,
                         ggplot2::aes(x = .data$DL, ymin = .data$lo66,
                                      ymax = .data$hi66),
                         fill = "steelblue", alpha = 0.25) +
    ggplot2::geom_line(data = curves,
                       ggplot2::aes(x = .data$DL, y = .data$age_ka),
                       colour = "steelblue", linewidth = 1) +
    ggplot2::geom_point(data = obs,
                        ggplot2::aes(x = .data$DL, y = .data$age_ka),
                        size = 1.2, alpha = 0.6) +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~acid, scales = "free_x") +
    ggplot2::labs(x = "D/L", y = "Age (ka BP)",
                  title = "Np AAR age calibration (posterior mean fit)",
                  subtitle = "Shading: 66% prediction interval") +
    ggplot2::theme_bw()
}
