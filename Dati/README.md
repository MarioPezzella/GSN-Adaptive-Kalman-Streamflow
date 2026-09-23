# Daily discharge data

The 37 files in this directory are the only inputs to the code. Everything
else in the repository is derived from them.

All files are read by `Main_GSN_ADPT_Kalman_Camastra.m` and are never
modified. Do not re-save them from a spreadsheet application: the values use
`.` as the decimal mark and `NaN` for missing observations, and a round trip
through Excel with an Italian locale will silently change the decimal
separator and may replace `NaN` with text or empty cells.

## Format

| Property | Value |
|---|---|
| Files | `q_1984.txt` … `q_2020.txt`, one per calendar year |
| Rows | one per day, from 1 January to 31 December (365 or 366 rows) |
| Columns | one, no header |
| Quantity | daily inflow discharge of the Camastra reservoir |
| Units | m³/s |
| Missing values | `NaN` |
| Encoding | ASCII, Windows line endings (CRLF) |

The files are read with `readmatrix` and concatenated in chronological
order into a single record of N = 13515 days. The script checks the number
of rows of each file against the calendar, 366 for leap years and 365
otherwise, and issues a warning if they do not match: an off-by-one in a
single year would shift the rest of the record against the calendar.

## Missing observations

Of the 13515 daily values, 1327 (9.82%) are missing. They are concentrated
in a few years:

| Year | Missing days |
|---|---|
| 1984 | 214 |
| 1985 | 153 |
| 1986 | 365 (whole year) |
| 1988 | 316 |
| 1994 | 5 |
| 2002 | 18 |
| 2003 | 72 |
| 2017 | 184 |

All other years are complete. Missing days are not interpolated: the filter
handles them through its prediction-only branch, propagating the state and
its covariance until the next available observation, and they are excluded
from every metric.

## Summary of the valid observations

| Statistic | Value [m³/s] |
|---|---|
| Valid days | 12188 |
| Mean | 3.758 |
| Median | 0.926 |
| Standard deviation | 7.157 |
| Minimum | 0.001 |
| Maximum | 130.185 |

The strong skewness, with a median about four times smaller than the mean,
reflects the alternation of long dry-season recessions with short
wet-season flood events.

## Provenance

The daily inflow record of the Camastra reservoir (Basilicata, southern
Italy) was provided by the Agency for Development of Irrigation and Land
Transformation in Puglia, Lucania and Irpinia. The historical inflow
series of the same reservoir was used for the GSN calibration of
Cimorelli et al.
(Journal of Water Resources Planning and Management, 2021,
[10.1061/(ASCE)WR.1943-5452.0001307](https://doi.org/10.1061/(ASCE)WR.1943-5452.0001307)).

If you use these data, cite the manuscript and acknowledge the data
provider, not only this repository. See the root `README.md` for the
citation block.
