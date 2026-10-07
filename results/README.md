# Generated results

These files are rebuilt from the input snapshot by `run_analysis.R`.

| File | Contents |
| --- | --- |
| `data_audit.csv` | Raw counts, missing-response counts, geographic exclusions, and modeled sample |
| `station_coverage.csv` | All 132 input stations, observation windows, valid counts, gaps, coordinates, and inclusion status |
| `station_results.csv` | 131 fitted U.S. stations: trends, uncertainty, residual diagnostics, seasonal features, and coverage |
| `model_coefficients.csv` | Six OLS coefficients for each station; trend is per year |
| `seasonal_curves.csv` | 365 seasonal phase evaluations per station at the reference baseline |
| `selected_stations.csv` | Stations selected for the seasonal-curve figure and their features |
| `hac_sensitivity.csv` | Median uncertainty ratios and above-zero interval counts at three calendar lag choices |
| `summary.csv` | Readable headline results |
| `input_manifest.csv` | Input filename, byte count, and MD5 fingerprint |
| `session_info.txt` | R and loaded-package versions used for the saved run |

In `station_results.csv`, `iid_*` fields use the classical independent-error OLS covariance and t intervals. `hac7_se`, `hac30_se`, and `hac90_se` are slope standard errors from calendar-aware Bartlett HAC sensitivity calculations; `hac30_lower` and `hac30_upper` use a normal critical value. All trend coefficients, standard errors, and limits are in degrees Celsius per year.

`seasonal_range_c` is peak-to-trough range at a fixed trend baseline. `coldest_day` and `warmest_day` label 365 equally spaced phases over an idealized 365.25-day cycle. `winter_order` runs from July to June. `residual_sd_c` uses residual degrees of freedom. `residual_cor_day1` uses only pairs separated by exactly one calendar day.
