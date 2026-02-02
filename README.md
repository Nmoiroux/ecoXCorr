# ecoXCorr

**ecoXCorr** is an R package designed to explore **lagged associations between environmental time series and ecological or epidemiological responses**.  

It provides a coherent workflow to:

1. Aggregate environmental variables over multiple lagged time windows  
2. Fit regression models across all lag structures  
3. Visualise effect strength and direction using **cross-correlation maps**

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
- a time unit in days (`lag_unit`)

For each reference date \code{d}, all intervals \eqn{[d - k \times i \times u,\; d - (l-1) \times i \times u)} are generated,
where \code{i} is the interval length (in units of \code{lag_unit}), \code{u} is the time unit in days, and \code{k, l} range from 1 to \code{m} with
\code{k >= l}.

```r
sampling_dates <- unique(albopictusMPL2023$date)

met_agg <- aggregate_lagged_intervals(
  data       = meteoMPL2023,
  date_col   = "date",
  value_cols = c("rain_sum", "temp_mean"),
  d          = sampling_dates,
  i          = 1,
  m          = 8,
  lag_unit   = 7   # weekly intervals
)

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

Because mosquito counts are zero-truncated and observations are not independent
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
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
