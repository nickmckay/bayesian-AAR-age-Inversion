# =============================================================================
# Two-acid AAR age calibration: forward/inverse transforms
# =============================================================================
#
# Calibration constants (a, e per acid, residual SD per acid, and their
# residual correlation rho) are taken directly from the manuscript's
# published OLS fit (see R/gls_calculator.R), not re-estimated here. Bayesian
# machinery in this package (R/mcmc.R, R/age_inversion.R) is used only to
# invert those known, fixed calibration curves for a sample's age -- it is
# not used to refit the curves themselves.

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
