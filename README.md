# ecoXCorr

**ecoXCorr** is an R package designed to explore **lagged associations between environmental time series and ecological or epidemiological responses**.  

It provides a coherent workflow to:

1. Aggregate environmental variables over multiple lagged time windows  
2. Fit regression models across all lag structures  
3. Visualise effect strength and direction using **cross-correlation maps** (CCM)

The package is particularly suited for studying **delayed environmental effects**, such as the influence of meteorological conditions on insect abundance or disease dynamics.

---

## Installation

Install the development version from GitHub:

```r
devtools::install_github("Nmoiroux/ecoXCorr")
```

## Overview of the workflow

The typical workflow in **ecoXCorr** follows three steps:

1. **Lagged aggregation of environmental data**  
   (`aggregate_lagged_intervals()`)

2. **Model fitting across lag windows**  
   (`fit_models_by_lag()`)

3. **Visualisation using cross-correlation maps**  
   (`plotCCM()`)

The package ships with two example datasets to illustrate this workflow:

- `meteoMPL2023`: daily meteorological data (Montpellier, France, 2023)
- `albopictusMPL2023`: mosquito sampling data (*Aedes albopictus*, Montpellier, 2023)

## Example

Load required packages and example data
```r
library(ecoXCorr)

# Meteorological daily time series
?meteoMPL2023
data(meteoMPL2023)

# Response data: tiger mosquito collections
?albopictusMPL2023
data(albopictusMPL2023)

```

### Aggregate meteorological variables over lagged intervals

Meteorological variables are aggregated over all possible lag windows defined by:

- a reference date (`d`),  
- a base interval length (`i`),
- a maximum lag (`m`)

For each reference date `d`, all intervals $[d - k \times i,\; d - (l-1) \times i)]$ are generated,
where `i` is the interval length (in days) and `k, l` range from 1 to `m` with `k >= l`.


In the example below `i` is set to 7 days, indicting that the unit for time intervals and lags will be a week.
`m` is set to 8, indicating that the maximum lag (between response and predictor variables) considered will be
8 weeks (or $m \times i$ days)

```r
?aggregate_lagged_intervals

# get sampling days from our entomological dataset
sampling_dates <- unique(albopictusMPL2023$date)

# perform data aggregation for cumulated rainfalls and mean daily temperatures
met_agg <- aggregate_lagged_intervals(
  data       = meteoMPL2023,
  date_col   = "date",
  value_cols = c("rain_sum", "temp_mean"),
  d          = sampling_dates,
  i          = 7,               
  m          = 8
)

head(met_agg)
```

### Join environmental and response data

Each reference date is associated with multiple lag windows, resulting in a many-to-many join:

```r
data <- full_join(
  met_agg,
  albopictusMPL2023,
  by = "date",
  relationship = "many-to-many"
)

```
### Fit models across lag windows
Simple GLM example

Here, a Poisson GLM is fitted independently for each lag window:

```r
?fit_models_by_lag

res_glm <- fit_models_by_lag(
  data       = data,
  response   = "individualCount",
  predictors = "temp_mean_mean",
  model      = "GLM",
  family     = poisson(link = "log")
)
```

### Visualise results as a cross-correlation map

```r
?plotCCM

plotCCM(res_glm, threshold_p = 0.2)
```
Each tile represents a lag window, with colour indicating the signed R²
(effect strength × direction). Non-significant associations (p>0.2) are masked.


### Mixed-effects model example

Because mosquito counts in the `albopictusMPL2023` dataset are zero-truncated and observations are not independent
(repeated measurements, spatial structure), a mixed-effects model may be more appropriate.
This example illustrates how ecoXCorr supports such models:

```r
?fit_models_by_lag

res_glmm <- fit_models_by_lag(
  data       = data,
  response   = "individualCount",
  predictors = "temp_mean_mean",
  model      = "GLM",
  random     = "(1|area/trap)",
  family     = truncated_nbinom2(link = "log")
)

plotCCM(res_glmm, threshold_p = 0.2)
```

The modelling function used depends on the combination of `model` and `random` arguments:

- `model = "LM"`  and `random = ""`:  `stats::lm()`
- `model = "GLM"` and `random = ""`:  `stats::glm()`
- `model = "LM"`  and `random != ""`: `lme4::lmer()`
- `model = "GLM"` and `random != ""`: `glmmTMB::glmmTMB()`


## When should I use ecoXCorr?

ecoXCorr is useful when:

- environmental drivers are expected to have delayed effects
- the relevant time scale of these effects is unknown
- you want a global view of lagged associations rather than testing a single lag

Typical applications include:

- vector ecology
- disease ecology
- environmental epidemiology
- climate–biology interactions


### References 
This package builds upon : 

- Curriero FC, Shone SM, Glass GE. (2005) *Cross correlation maps: a tool for visualizing and modeling time lagged associations.* [Vector Borne Zoonotic Dis.](https://doi.org/10.1089/vbz.2005.5.267)


### License 
This package is released under the [GPL-3 License](https://www.gnu.org/licenses/gpl-3.0-standalone.html). 
