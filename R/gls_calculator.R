# =============================================================================
# GLS/softplus two-acid age calculator (reference spreadsheet method)
# =============================================================================
#
# A faithful R implementation of the spreadsheet age calculator described in
# the Np AAR calibration study: generalized least squares (GLS) combination
# of independent TDK (Asp) and SPK (Glu) age estimates, using a
# softplus-smoothed weight floor to keep both weights positive even when the
# raw GLS weight would go negative. Included so the Bayesian inversion built
# elsewhere in this repository (age_calibration.R, age_inversion.R) can be
# compared directly against the original point-estimate-plus-prediction-
# interval method, using its own published formulas and calibration
# constants.

# Calibration constants for the GLS age calculator.
# Table 2 / calibration constants from the Np AAR calibration study
# (n = 140 Neogloboquadrina pachyderma samples, the manuscript's current
# dataset), refit here by OLS from data/np_calibration.csv following the
# same TDK(Asp)/SPK(Glu) methodology as the spreadsheet's Constants sheet.
# Cross-checked against the manuscript text: MSE_Asp = 0.166 vs. published
# 0.166, MSE_Glu = 0.214 vs. published 0.214, rho = 0.839 vs. published
# 0.840 -- confirms this refit reproduces the manuscript's actual Table 2.
# The prior version of this file (n = 128) is archived at
# data/np_calibration_128_archive.csv; it predates several Nordic Seas and
# Southern Ocean samples added to the manuscript's calibration dataset.
gls_calibration_constants <- list(
  a_Asp    = 2755.1216120445,
  e_Asp    = 2.7086094456,
  MSE_Asp  = 0.1660775086,
  n_Asp    = 140,
  xbar_Asp = -1.3647882802,
  Sxx_Asp  = 40.9855516223,
  a_Glu    = 8928.6968578976,
  e_Glu    = 2.1528457896,
  MSE_Glu  = 0.2142008076,
  n_Glu    = 140,
  xbar_Glu = -2.2632774360,
  Sxx_Glu  = 63.4451634208,
  rho      = 0.8389317275,
  SE_rho   = 0.0463285022,
  t50      = 0.6762717115,
  t66      = 0.9574788976,
  t90      = 1.6559703824
)


# GLS-combined two-acid age with softplus weight smoothing.
# Reference implementation of the spreadsheet age calculator: independent
# TDK (Asp) and SPK (Glu) age predictions are combined by generalized least
# squares, with the (possibly negative) raw GLS weights passed through a
# softplus floor before combination.
# DL_Asp, DL_Glu: observed D/L ratios (sample means).
# SD_Asp, SD_Glu: standard deviation across subsamples.
# N_Asp, N_Glu: number of subsamples used for the mean.
# Returns a one-row data frame with columns t_Asp, t_Glu, V_Asp, V_Glu, C,
# w_Asp, w_Glu, D, V_int, V_excess, V_total, sigma_comb, t_comb, lo50, hi50,
# lo66, hi66, lo90, hi90.
predict_age_gls <- function(DL_Asp, SD_Asp, N_Asp, DL_Glu, SD_Glu, N_Glu,
                            constants = gls_calibration_constants) {
  cst <- constants
  SE_Asp <- SD_Asp / sqrt(N_Asp)
  SE_Glu <- SD_Glu / sqrt(N_Glu)

  t_Asp <- cst$a_Asp * atanh(DL_Asp)^cst$e_Asp
  t_Glu <- cst$a_Glu * DL_Glu^cst$e_Glu

  # Predictor-derivative terms propagate D/L measurement uncertainty
  # through the log(f(D/L)) transform via the chain rule.
  slope_Asp <- cst$e_Asp / (atanh(DL_Asp) * (1 - DL_Asp^2))
  slope_Glu <- cst$e_Glu / DL_Glu

  V_Asp <- cst$MSE_Asp * (1 + 1 / cst$n_Asp +
                           (log(atanh(DL_Asp)) - cst$xbar_Asp)^2 / cst$Sxx_Asp) +
    (slope_Asp * SE_Asp)^2
  V_Glu <- cst$MSE_Glu * (1 + 1 / cst$n_Glu +
                           (log(DL_Glu) - cst$xbar_Glu)^2 / cst$Sxx_Glu) +
    (slope_Glu * SE_Glu)^2

  C <- cst$rho * sqrt(V_Asp * V_Glu)
  s <- sqrt(V_Asp * V_Glu) * cst$SE_rho

  .softplus_floor <- function(raw, s) {
    pmax(raw, 0) + s * log(1 + exp(-abs(raw / s)))
  }
  w_Asp <- .softplus_floor(V_Glu - C, s)
  w_Glu <- .softplus_floor(V_Asp - C, s)
  D <- w_Asp + w_Glu

  V_int    <- (V_Asp * V_Glu - C^2) / D
  V_excess <- w_Asp * w_Glu * (log(t_Asp) - log(t_Glu))^2 / D^2
  V_total  <- V_int + V_excess
  sigma_comb <- sqrt(V_total)

  t_comb <- t_Asp^(w_Asp / D) * t_Glu^(w_Glu / D)

  data.frame(
    t_Asp = t_Asp, t_Glu = t_Glu, V_Asp = V_Asp, V_Glu = V_Glu, C = C,
    w_Asp = w_Asp, w_Glu = w_Glu, D = D, V_int = V_int, V_excess = V_excess,
    V_total = V_total, sigma_comb = sigma_comb, t_comb = t_comb,
    lo50 = t_comb * exp(-cst$t50 * sigma_comb),
    hi50 = t_comb * exp(+cst$t50 * sigma_comb),
    lo66 = t_comb * exp(-cst$t66 * sigma_comb),
    hi66 = t_comb * exp(+cst$t66 * sigma_comb),
    lo90 = t_comb * exp(-cst$t90 * sigma_comb),
    hi90 = t_comb * exp(+cst$t90 * sigma_comb)
  )
}
