# bayesian-AAR-age-Inversion

Independent Bayesian validation of the Np (*Neogloboquadrina pachyderma*)
two-acid AAR age calculator, in support of section 4.1 of the accompanying
manuscript. Given the manuscript's own fitted TDK (Asp) / SPK (Glu)
calibration curves, a Bayesian MCMC inversion (bivariate-normal correlated
residual likelihood) is compared against the manuscript's GLS/softplus
spreadsheet calculator across the full 128-sample Np calibration set.

The analysis is a single self-contained R Markdown notebook,
[`np_bayesian_age_validation.Rmd`](np_bayesian_age_validation.Rmd), backed
by plain R scripts in [`R/`](R/) and the calibration dataset in
[`data/np_calibration.csv`](data/np_calibration.csv). No package
installation is required.

## Rendering locally

```r
install.packages(c("rmarkdown", "knitr", "ggplot2", "dplyr"))
rmarkdown::render("np_bayesian_age_validation.Rmd")
```

## Continuous integration

A GitHub Actions workflow (`.github/workflows/render.yml`) renders the
notebook on every push/PR to `main` and publishes the rendered HTML to
GitHub Pages on pushes to `main`.
