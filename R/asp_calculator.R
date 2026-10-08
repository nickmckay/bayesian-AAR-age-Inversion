# =============================================================================
# Asp-only age calculator (reference spreadsheet method)
# =============================================================================
#
# A faithful R implementation of the revised Np AAR age calculator described
# in the manuscript (AAR_Age_Calculator_NP.xlsx, "Age Calculator" sheet): the
# age is computed from aspartic acid alone through the TDK0 equation, with a
# Student-t prediction interval that combines calibration scatter, leverage,
# and measurement error propagated through the model derivative. Glutamic
# acid is not combined into the age; its SPK0 age is reported only as a
# quality check. Included so the Bayesian two-acid inversion built elsewhere
# in this repository (age_calibration.R, age_inversion.R) can be compared
# directly against the manuscript's published method and constants.
#
# The earlier blended Asp + Glu GLS/softplus calculator, which the manuscript
# no longer uses, is preserved in release v0.1-presubmission of this
# repository.

# Calibration constants, copied from the "Constants" sheet of the revised
# spreadsheet calculator (n = 140 Neogloboquadrina pachyderma samples, with
# the 22 hr hydrolysis samples from core MD99-2289 normalized to 6 hr
# equivalents). The spreadsheet does not carry xbar_Glu or Sxx_Glu, because
# it never computes a Glu prediction interval; those two values are computed
# here from data/np_calibration.csv so that the Bayesian inversion can give
# Glu the same per-sample leverage treatment as Asp.
asp_calibration_constants <- list(
  a_Asp      = 2580.171720155672,
  e_Asp      = 2.537515540180582,
  MSE_Asp    = 0.1077082153969766,
  n          = 140,
  xbar_Asp   = -1.431371812733824,
  Sxx_Asp    = 47.89650969613734,
  a_Glu      = 7592.070431548645,
  e_Glu      = 2.011676937361452,
  MSE_Glu    = 0.1706291661730284,
  xbar_Glu   = -2.3425025051,  # computed from data/np_calibration.csv
  Sxx_Glu    = 74.1709261240,  # computed from data/np_calibration.csv
  rho        = 0.7944985872107606,
  sigma_disc = 0.2508450876806995,
  delta_Asp  = 0.05058,  # 22 hr -> 6 hr hydrolysis offset in atanh(D/L)
  delta_Glu  = 0.02201,
  t50        = 0.676272,
  t66        = 0.957479,
  t90        = 1.65597
)


# Normalize a D/L value measured after 22 hr hydrolysis to its 6 hr
# equivalent: DL_6 = tanh(atanh(DL_22) - delta). Undefined (NA) for
# DL_22 <= tanh(delta).
normalize_22hr <- function(DL, delta) {
  out <- tanh(atanh(DL) - delta)
  out[out <= 0] <- NA_real_
  out
}


# Per-sample log-age prediction variance for one acid: calibration MSE
# inflated for the sample's leverage on the fitted line, plus the standard
# error of the mean D/L propagated through d ln(t) / d(D/L).
.prediction_variance <- function(x, slope, SE, MSE, n, xbar, Sxx) {
  MSE * (1 + 1 / n + (x - xbar)^2 / Sxx) + (slope * SE)^2
}


# Asp-only age with prediction intervals, plus the Glu quality check.
# DL_Asp, DL_Glu: observed D/L ratios (sample means, already 6 hr equivalent).
# SD_Asp, SD_Glu: standard deviation across subsamples.
# N: number of subsamples used for the mean.
# Returns a data frame with columns t_Asp, V_Asp, sigma_Asp, lo50, hi50,
# lo66, hi66, lo90, hi90 (the spreadsheet outputs), t_Glu, ratio and flag
# (the spreadsheet quality check), and V_Glu (the analogous Glu prediction
# variance, not reported by the spreadsheet but used by the Bayesian
# inversion).
predict_age_asp_only <- function(DL_Asp, SD_Asp, DL_Glu, SD_Glu, N,
                                 constants = asp_calibration_constants) {
  cst <- constants

  x_Asp     <- log(atanh(DL_Asp))
  slope_Asp <- cst$e_Asp / (atanh(DL_Asp) * (1 - DL_Asp^2))
  V_Asp     <- .prediction_variance(x_Asp, slope_Asp, SD_Asp / sqrt(N),
                                    cst$MSE_Asp, cst$n, cst$xbar_Asp, cst$Sxx_Asp)
  sigma_Asp <- sqrt(V_Asp)
  t_Asp     <- cst$a_Asp * atanh(DL_Asp)^cst$e_Asp

  x_Glu     <- log(DL_Glu)
  slope_Glu <- cst$e_Glu / DL_Glu
  V_Glu     <- .prediction_variance(x_Glu, slope_Glu, SD_Glu / sqrt(N),
                                    cst$MSE_Glu, cst$n, cst$xbar_Glu, cst$Sxx_Glu)
  t_Glu     <- cst$a_Glu * DL_Glu^cst$e_Glu

  disc <- abs(log(t_Asp / t_Glu))
  flag <- ifelse(disc > 2.576 * cst$sigma_disc, "FAIL (>99%)",
                 ifelse(disc > 1.96 * cst$sigma_disc, "check (>95%)", "ok"))

  data.frame(
    t_Asp = t_Asp, V_Asp = V_Asp, sigma_Asp = sigma_Asp,
    lo50 = t_Asp * exp(-cst$t50 * sigma_Asp),
    hi50 = t_Asp * exp(+cst$t50 * sigma_Asp),
    lo66 = t_Asp * exp(-cst$t66 * sigma_Asp),
    hi66 = t_Asp * exp(+cst$t66 * sigma_Asp),
    lo90 = t_Asp * exp(-cst$t90 * sigma_Asp),
    hi90 = t_Asp * exp(+cst$t90 * sigma_Asp),
    t_Glu = t_Glu, V_Glu = V_Glu, ratio = exp(disc), flag = flag
  )
}
