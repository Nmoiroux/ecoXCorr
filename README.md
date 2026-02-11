[![DOI](https://zenodo.org/badge/1147854443.svg)](https://doi.org/10.5281/zenodo.18600567)

# ecoXCorr

**ecoXCorr** (pronounce "*Eco-Cross-Corr*") is an R package designed to explore **lagged associations between environmental time series and ecological or epidemiological responses** based on the method proposed by Curriero *et al.* (2005)[^1].  

It provides a coherent workflow to:

1. Aggregate environmental variables over multiple lagged time windows  
2. Fit regression models across all lag structures  
3. Visualise effect strength and direction using **cross-correlation maps** (CCM)

The package is particularly suited for studying **delayed environmental effects**, such as the influence of meteorological conditions on insect abundance or disease dynamics.

`ecoXCorr` has similar features than amazing [`climwin`](https://cran.r-project.org/web/packages/climwin/index.html) package but use [`glmmTMB`](https://cran.r-project.org/web/packages/glmmTMB/index.html) allowing to fit (generalized) linear (mixed-)models using a large variety of error distribution (including negative-binomial, zero-inflated, zero-truncated... see [`?glmmTMB::family_glmmTMB`](https://glmmtmb.github.io/glmmTMB/reference/nbinom2.html)) and covariance structures (see [`vignette(glmmTMB::covstruct)`](https://glmmtmb.github.io/glmmTMB/articles/covstruct.html)). `ecoXCorr` is also more flexible for interval lengths allowing to specify interval in number of days, not restricted to standard time periods (e.g. "week" or "month") as in `climwin`. 

Below is an exemple of figure computed using `ecoXCorr`.
![plot](/man/figures/Rplot.jpg)

*Fig. 1: Cross correlation maps showing the lagged effect of rainfall on Ae. albopictus abundance. Time lags are expressed in days prior to sampling. The signed R² reflects the variance explained by the explanatory variable, multiplied by the sign of the estimated effect. Pink-bordered square highlight the time lag with the highest R². Grey squares represent correlations with adjusted (for multiple testing) p-values > 0.05.* 

---

## Installation

Install the development version from GitHub:

```r
devtools::install_github("Nmoiroux/ecoXCorr")
```

## Overview of the workflow

The typical workflow in **ecoXCorr** follows three steps with respective functions:

1. **Lagged aggregation of environmental data**  
   `aggregate_lagged_intervals()`

2. **Model fitting across lag windows**  
   `fit_models_by_lag()`

3. **Visualisation using cross-correlation maps**  
   `plotCCM()`

Steps 1 and 2 can be merged into a single step using the `ecoXCorr()` wrapper function.

The package ships with two example datasets to illustrate this workflow:

- `meteoMPL2023`: daily meteorological data (Montpellier, France, 2023)
- `albopictusMPL2023`: mosquito sampling data (*Aedes albopictus*, Montpellier, 2023)

The package includes a user-friendly Shiny application: [`ecoXCorrApp`](https://nicolas-moiroux.shinyapps.io/ecoXCorrApp/).

## Example

Load required packages and example data
```r
library(ecoXCorr)

# Meteorological daily time series
?meteoMPL2023
head(meteoMPL2023)

# Response data: tiger mosquito collections
?albopictusMPL2023
head(albopictusMPL2023)

```

### Aggregate meteorological variables over lagged intervals

Meteorological variables are aggregated over all possible lag windows defined by:

- a reference date (`d`),  
- a base interval length (`i`),
- a maximum lag (`m`)

For each reference date `d`, all intervals $[d - k \times i + 1; d - (l-1) \times i)]$ are generated,
where `i` is the interval length (in days) and `k, l` range from 1 to `m` with `k >= l`.


In the example below `i` (`interval`) is set to 7 days, indicting that the unit for time intervals and lags will be a week.
`m` (`max_lag`) is set to 8, indicating that the maximum lag (between response and predictor variables) considered will be
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
  ref_date   = sampling_dates,
  interval   = 7,               
  max_lag    = 8
)

head(met_agg)
```

### Join environmental and response data

Each reference date is associated with multiple lag windows, resulting in a many-to-many join:

```r
data <- merge(met_agg,
              albopictusMPL2023, by = "date", all = TRUE)

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
  family     = "poisson")

```

### One-step approach using the `ecoXCorr()` wrapper function:

This approach is allow to plot CCM in one step giving one dataset of meteorological time series and one dataset of response variable. Aggregation is performed for only one variable in the meteo dataset (argument `value_cols`) according to one function (`agg_fun`). The resulting variable is used as the predictor in the modeling process. 


```r
?ecoXCorr()

res_glm <- ecoXCorr(
  meteo_data    = meteoMPL2023,
  response_data = albopictusMPL2023,
  date_col_meteo = "date",
  date_col_resp = "date",
  value_cols    = "rain_sum",
  agg_fun       = "sum",
  response      = "individualCount",
  interval      = 7,
  max_lag       = 8,
  family        = "poisson"
)
```
 
### Visualise results as cross-correlation maps

```r
?plotCCM

plotCCM(res_glm, model_outcome ="R2sign", threshold_p = 0.2)
```
Each tile represents a lag window, with colour indicating the signed R²
(% of variance explained × direction). Non-significant associations (p>0.2) are masked.

Other outcomes can be plotted (R², delta-AIC [2], Akaike weight [2], beta parameters of the linear predictor):

```r
plotCCM(res_glm, model_outcome = "R2")
plotCCM(res_glm, model_outcome = "d_aic")
plotCCM(res_glm, model_outcome = "betas")
```

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
  random     = "(1|area/trap)",
  family     = "truncated_nbinom2")

plotCCM(res_glmm, model_outcome ="R2sign", threshold_p = 0.2)
```

The modelling function used depends on the `random` and `family` arguments:

- `random` is empty:  [`stats::glm()`](https://rdrr.io/r/stats/glm.html)
- `random` is specified OR `family` is a valid glmmTMB family: [`glmmTMB::glmmTMB()`](https://glmmtmb.github.io/glmmTMB/reference/glmmTMB.html)


## When should I use ecoXCorr?

ecoXCorr is useful when:

- environmental drivers are expected to have delayed effects
- the relevant time scale of these effects is unknown
- you want a global view of lagged associations rather than testing a single lag

Typical applications include:

- variable and feature selection in modelling
- vector ecology
- disease ecology
- environmental epidemiology
- climate–biology interactions


### References (Methods)

[^1].  Curriero FC, Shone SM, Glass GE. (2005) Cross correlation maps: a tool for visualizing and modeling time lagged associations. *Vector Borne Zoonotic Dis.* [doi:10.1089/vbz.2005.5.267](https://doi.org/10.1089/vbz.2005.5.267)
2.  van de Pol M, Bailey LD, McLean N, et al. (2016) Identifying the best climatic predictors in ecology and evolution. *Methods in Ecology and Evolution.* [doi:10.1111/2041-210X.12590](https://doi.org/10.1111/2041-210X.12590)

### References (Use examples)
3.  Bartholomée C, Taconet P, Mercat M, Grail C, Bouhsira E, Fournet F, et al. Investigating the role of urban vegetation alongside other environmental variables in shaping Aedes albopictus presence and abundance in Montpellier, France. PLOS ONE. 2025;20: e0335793. [doi:10.1371/journal.pone.0335793](https://doi.org/10.1371/journal.pone.0335793)
4.  Taconet P, Porciani A, Soma DD, Mouline K, Simard F, Koffi AA, et al. Data-driven and interpretable machine-learning modeling to explore the fine-scale environmental determinants of malaria vectors biting rates in rural Burkina Faso. Parasites & Vectors. 2021;14: 345. [doi:10.1186/s13071-021-04851-x](https://doi.org/10.1186/s13071-021-04851-x)


### License 
This package is released under the [GPL-3 License](https://www.gnu.org/licenses/gpl-3.0-standalone.html). 
