# Data source and exact reproduction

The project uses `uscrn_daily.RData`, a course-provided snapshot of daily climate-station observations derived from NOAA's USCRN/USRCRN daily product. It contains a data frame named `dat` with 959,126 rows, 12 columns, and 132 station names. The file supplied for this repository includes `ON_Egbert_1_W`; the main contiguous-U.S. analysis excludes that station explicitly.

## Use the supplied file

Copy `uscrn_daily.RData` into `data/raw/`, then run the project normally. Do not rename it, unpack it, or edit an absolute path in the code. The file is intentionally excluded from the GitHub upload.

To reproduce the exact saved results, obtain this same snapshot from the author or the course materials. An exact public download URL for the course-prepared RData was not supplied. This repository does not claim that downloading the current NOAA archive recreates that snapshot.

| Fingerprint | Value |
| --- | --- |
| Filename | `uscrn_daily.RData` |
| Size | 11,851,208 bytes |
| MD5 | `f513a11063c7bcadce096d219d0c98e4` |

The MD5 is a reproducibility identifier, not a security guarantee. Each execution records the actual input's fingerprint in `results/input_manifest.csv`.

## NOAA documentation

- [Network overview](https://www.ncei.noaa.gov/access/crn/)
- [Daily-product archive](https://www.ncei.noaa.gov/pub/data/uscrn/products/daily01/)
- [Field definitions and missing-value codes](https://www.ncei.noaa.gov/pub/data/uscrn/products/daily01/readme.txt)

NOAA's archive provides station-year text files and snapshots. Those are public source data; using them here would require a separate import step that constructs the expected `dat` object, matches the station collection, and fixes the historical snapshot date. Such an importer is not included.

## Fields used

| Field | Use |
| --- | --- |
| `site` | Station identifier; its state/province prefix defines the geographic sample |
| `LST_DATE` | Local-standard-time calendar date, stored as YYYYMMDD |
| `T_DAILY_AVG` | Daily average temperature, in degrees Celsius |
| `LONGITUDE`, `LATITUDE` | Observed station coordinates for maps |

The response's `-9999` sentinel becomes `NA` before exploration or modeling. Missing responses are omitted, with counts exported. We do not claim to clean unused precipitation, humidity, or other fields.

Raw dates span November 14, 2000–April 19, 2025; usable temperature dates span June 4, 2001–April 18, 2025. The contiguous-U.S. subset has 951,542 raw rows, of which 10,412 lack a usable temperature. Its 941,130 remaining observations come from 131 stations with unequal coverage. See `results/station_coverage.csv` for station-specific dates and counts.
