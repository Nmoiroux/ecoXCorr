# ecoXCorr 0.2.0

## New features

- `aggregate_lagged_intervals()` now supports multiple independent time series via a new `id_col` argument.
  When specified, lagged aggregations are performed separately within each group, preventing unintended mixing
  of values across distinct time series (e.g. sites, stations, individuals).

## Improvements

- Improved robustness of `fit_models_by_lag()` when computing R² using `performance::r2()`.
  Warnings are now handled more gracefully, and R² is re-evaluated when possible instead of returning `NA`.
- More informative warning messages are collected and returned to the user when issues occur during R² computation.

## Bug fixes

- Fixed aggregation interval: The function aggregate_lagged_intervals() has been corrected to properly include the last day of each interval in the aggregation. Previously, intervals were defined as start <= date < end, which systematically excluded the final day. This has now been changed to start <= date <= end, ensuring that all dates within the intended interval are considered when computing summary statistics. This fix guarantees that the aggregated results accurately reflect the intervals specified by the user.

- Fixed an error occurring in `fit_models_by_lag()` during result aggregation:
`
Error in match.names(clabs, names(xi)) : names do not match previous names
`
This error occurred when `performance::r2()` failed (error or warning) for some models,
leading to inconsistent column structures across elements of the results list.
The function now ensures consistent output structure across all lag windows by:
- standardizing R² extraction,
- returning properly formatted values even when warnings occur,
- preventing failures in `do.call(rbind, results)`.

## Notes

- This update improves stability when fitting large numbers of models or using complex model families
(e.g. zero-inflated or truncated distributions via `glmmTMB`).