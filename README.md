# bayesian-AAR-age-Inversion

Independent Bayesian validation of the Np (*Neogloboquadrina pachyderma*)
AAR age calculator, in support of the Validation section of the
accompanying manuscript. The manuscript computes ages from aspartic acid
alone (TDK0), keeping glutamic acid (SPK0) only as a quality check. Here a
Bayesian MCMC inversion that uses both amino acids (bivariate-normal
correlated residual likelihood) is compared against the Asp-only
spreadsheet calculator across the full 140-sample Np calibration set,
first with the manuscript's published constants (Test 1) and then with an
independent joint Bayesian refit of the calibration (Test 2).

The analysis is a single self-contained R Markdown notebook,
[`np_bayesian_age_validation.Rmd`](np_bayesian_age_validation.Rmd), backed
by plain R scripts in [`R/`](R/) and the calibration dataset in
[`data/np_calibration.csv`](data/np_calibration.csv) (the Np subset of the
manuscript's Table S1, with 22 hr hydrolysis samples normalized to 6 hr
equivalents). No package installation is required.

Release v0.1-presubmission validated an earlier, blended Asp + Glu
calculator that the manuscript no longer uses; its dataset is archived at
`data/np_calibration_140_presubmission_archive.csv`.

## Rendering locally

```r
install.packages(c("rmarkdown", "knitr", "ggplot2", "dplyr"))
rmarkdown::render("np_bayesian_age_validation.Rmd")
```

## Continuous integration

A GitHub Actions workflow (`.github/workflows/render.yml`) renders the
notebook on every push/PR to `main` and publishes the rendered HTML to
GitHub Pages on pushes to `main`.
