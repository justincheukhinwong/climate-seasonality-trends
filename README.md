# Seasonality and temperature trends in climate station data

**Separating seasonal cycles, gradual temperature change, and daily variability with linear models in R.**

This independent SDS 439: Linear Statistical Models project fits station-specific harmonic regressions to daily average temperatures from a course-provided NOAA climate-station dataset. It demonstrates model specification, coefficient interpretation, uncertainty assessment, residual diagnostics, and careful limits on inference.

**Start with the [rendered analysis](analysis.md).** It includes the statistical explanation, code, results, and figures; you can read it on GitHub without installing R. The [R Markdown source](analysis.Rmd) reproduces the report.

![Fitted seasonal cycles at four selected stations](figures/selected-cycles-1.png)

## Main findings

| Finding | Result in the supplied snapshot |
| --- | --- |
| Most station trends are positive | 130 of 131 contiguous-U.S. stations |
| Typical fitted change | Mean station slope: 0.0483 degrees C/year; median: 0.0492 |
| Residuals remain dependent | Median consecutive-day residual correlation: 0.71 |
| Uncertainty depends on error assumptions | 122 classical 95% intervals are above zero; 95 are above zero with the 30-day HAC sensitivity calculation |
| Seasonal variation differs substantially | Fitted peak-to-trough ranges run from 4.14 to 33.54 degrees C |

The main sample contains **941,130 usable daily observations**, collectively spanning **June 4, 2001–April 18, 2025**. Station windows differ. The original file contains 132 stations; the Ontario station is excluded from the contiguous-U.S. analysis.

The mean slope is an unweighted descriptive average, not an area-weighted national climate estimate. Intervals are pointwise, with no multiple-comparison adjustment. HAC intervals are an approximate sensitivity analysis; they do not resolve spatial dependence or mean-model misspecification.

## Statistical approach

Each station gets six coefficients: an intercept, a linear time trend, and annual and semiannual sine/cosine pairs. The nonlinear time functions are fixed predictors, so estimation remains ordinary linear regression. Separate fits allow each station to have its own trend, seasonal shape, and residual scale.

The report explains the intercept and annual trend, separates seasonal range from within-year drift, and distinguishes residual standard deviation from coefficient uncertainty. It compares classical t intervals with calendar-aware Bartlett HAC intervals using 7-, 30-, and 90-day lags. A custom calendar-aware calculation preserves gaps in the daily record and is checked against an established R implementation on complete data.

## Run it in RStudio

1. Install [R](https://cran.r-project.org/) and, if you want a graphical editor, [RStudio](https://docs.posit.co/ide/user/).
2. Download or clone this repository. Open `climate-seasonality-trends.Rproj` in RStudio.
3. Run this in the **R Console** once to install packages:

   ```r
   source("setup.R")
   ```

4. Put your supplied `uscrn_daily.RData` file in `data/raw/`. Keep the filename unchanged. No path editing or unzipping is needed.
5. Run:

   ```r
   source("run_analysis.R")
   ```

The script rebuilds `analysis.md`, seven figures, and the CSV results. Open `analysis.md` on GitHub or in a Markdown preview to read it. You can also open `analysis.Rmd` in RStudio to inspect and run its chunks. Run the chunks from the beginning so their objects exist in the correct order.

If you prefer RStudio's **Knit** button, install the optional `rmarkdown` package first with `install.packages("rmarkdown")`; the supplied runner only needs `knitr` and `maps`.

## Run it in Terminal

From this project folder, after installing R:

```bash
Rscript setup.R
Rscript run_analysis.R
```

Run `setup.R` once, then repeat only `Rscript run_analysis.R` after changing the analysis. Run the numerical HAC checks when changing the covariance calculation:

```bash
Rscript tests/validate_hac.R
```

Tested with **R 4.3.3**, **knitr 1.45**, **maps 3.4.2**, and **sandwich 3.1.0** for the numerical checks. `setup.R` installs available package versions rather than locking them to this environment. The recorded run environment is in [results/session_info.txt](results/session_info.txt).

## Repository map

| File or folder | Purpose |
| --- | --- |
| `analysis.Rmd` / `analysis.md` | Reproducible source and rendered statistical report |
| `R/analysis_helpers.R` | Data checks, station models, summaries, and plots |
| `R/calendar_hac.R` | Covariance sensitivity calculation using calendar-day lags |
| `run_analysis.R` / `setup.R` | One-command analysis and package setup |
| `data/` | Input instructions; raw data are excluded from Git |
| `figures/` | Seven generated figures displayed in the report |
| `results/` | Data audit, station results, coefficients, seasonal curves, and provenance |
| `tests/validate_hac.R` | Numerical checks against independent covariance calculations |
| `VALIDATION.md` | Changes from the coursework version and validation evidence |

The [data instructions](data/README.md) distinguish the course-prepared snapshot from NOAA's live public archive. The exact RData snapshot is not bundled; readers need that file to reproduce the saved results.

## Project context

Original analysis: Justin Wong, SDS 439: Linear Statistical Models, Spring 2025. The course supplied the question, dataset, and two-harmonic model specification. This repository revises the implementation and presentation; the HAC sensitivity analysis and expanded diagnostics are later additions. Assignment prompts, the student ID, and the local Mac path have been removed from the public presentation.
