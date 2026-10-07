# Validation and changes from the coursework version

## What was retained

The original scientific question and station-specific OLS mean model are retained: an intercept, a linear trend, and annual and semiannual sine/cosine pairs. Scaling time from days to years changes coefficient units, not fitted values. The calendar origin remains January 1, 2020.

## Corrections and additions

| Original issue | Repository treatment |
| --- | --- |
| Early seasonal and station plots used uncleaned temperature sentinels | Clean the response before every plot and fit; export the sample audit |
| Description called all 132 stations lower-48 sites | Exclude the Ontario station explicitly; export all 132 stations' inclusion status |
| A user-specific absolute file path | Load the input from `data/raw/` using the project folder |
| Repeated fitting of identical models | Fit each station once and reuse named model objects |
| Slope reported in degrees per day | Report degrees per year and explain the seasonal adjustment |
| Intercept described as lacking interpretation | Explain the seasonally adjusted baseline and the nonzero cosine terms at the reference date |
| Only approximate independent-error intervals | Use exact classical t critical values under the teaching assumptions; add separate HAC sensitivity intervals |
| Autocorrelation mentioned without measurement | Add calendar-lag diagnostics and 7/30/90-day HAC sensitivity |
| Seasonal range included within-year linear drift | Hold trend at its reference baseline and vary harmonic phase |
| “Amplitude” meant maximum minus minimum | Label the quantity seasonal peak-to-trough range |
| `sd(residuals)` used denominator n-1 | Use `sigma(fit)`, based on residual degrees of freedom n-6 |
| Selected stations described by guessed locations | Identify stations by name and recompute their extrema |
| Regional explanations stated too strongly | Report observed station patterns and identify untested mechanisms as hypotheses |
| Assignment questions, due date, and student ID appeared in the document | Present a standalone research narrative with course attribution |

HAC inference is a substantive later extension. It changes the covariance estimate, not the OLS coefficients. It should not be presented as part of the original submitted coursework.

## Completed checks

- Executed the full R Markdown analysis in a fresh Rscript process with the supplied RData file.
- Verified 959,126 input rows, 132 station names, 10,433 missing-temperature sentinels, no duplicate station-date keys, valid calendar dates, and stable coordinates within every station.
- Verified 131 fitted contiguous-U.S. stations and 941,130 usable observations; every fitted design has rank six.
- Independently recomputed OLS coefficients, classical intervals, residual standard deviations, seasonal ranges, and consecutive-day correlations using Python/NumPy and SciPy. The R results agree within numerical tolerance.
- Checked the calendar-aware HAC estimator against `sandwich::NeweyWest(prewhite = FALSE)` on complete daily data and an independent full pair-weight calculation on irregular dates. The checks cover multiple lags, observed-row finite-sample correction, gap handling, row-order invariance, and the lag-zero HC1 identity.
- Inspected all seven generated figures for labels, scales, and rendering problems.

The checked environment was R 4.3.3, knitr 1.45, maps 3.4.2, and sandwich 3.1.0 for covariance checks. The analysis does not require Python; Python was used only for an independent validation during repository preparation.

## Reconciliation with the original results

The original 132-station specification gives an unweighted average slope of 0.04827643 degrees C/year, matching the original report's approximately 0.048. Excluding Ontario gives 0.04827048. The point-estimate conclusion is therefore essentially unchanged by the geographic correction.

Uncertainty interpretation does change: among the 131 included stations, 122 classical 95% intervals lie entirely above zero, compared with 95 under the HAC30 calculation. Median HAC30 slope standard error is 2.13 times the classical standard error. At lags 7 and 90 days, the counts are 105 and 88 respectively. These are pointwise sensitivity results, not a simultaneous or spatially adjusted national test.

Corrected seasonal ranges run from 4.14389 degrees C at CA_Bodega_6_WSW to 33.54088 at MN_Goodridge_12_NNW. The idealized-year grid places the earliest winter minimum at OR_Riley_10_WSW (day 361) and the latest at CA_Yosemite_Village_12_W (day 37). These corrected seasonal-phase summaries are intentionally not exact copies of the original trend-contaminated predictions.

## Re-run the numerical covariance checks

From the repository root after `setup.R`:

```bash
Rscript tests/validate_hac.R
```

The input data are not needed for these synthetic-data tests. The tests verify the estimator's implementation; they do not prove that any fixed lag choice captures all dependence in real temperature data.
