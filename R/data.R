#' Entomological records of female *Aedes albopictus* in Montpellier (2023)
#'
#' This dataset contains entomological sampling records of female
#' *Aedes albopictus* collected in Montpellier (France) during 2023.
#' The data originate from GBIF and correspond to a zero-truncated subset
#' (i.e. only positive captures) of adult female mosquitoes collected using
#' fixed traps across different areas of the city.
#'
#' The dataset was derived from a GBIF annotated archive
#' (DOI: 10.15468/4qafbu) and processed to extract sampling dates, trap
#' identifiers, and spatial grouping variables.
#'
#' @format A data.frame with the following variables:
#' \describe{
#'   \item{species}{Scientific name of the species (*Aedes albopictus*).}
#'   \item{individualCount}{Number of female individuals captured during the sampling event.}
#'   \item{eventDate}{Original event date as provided by GBIF.}
#'   \item{trap}{Identifier of the trap used for sampling.}
#'   \item{area}{Identifier of the area where the trap is located.}
#'   \item{date}{Sampling date (class \code{Date}).}
#'   \item{ID}{Unique sampling event identifier.}
#' }
#'
#' @details
#' Only records corresponding to adult females of *Aedes albopictus*
#' were retained. Observations with zero counts were removed, resulting
#' in a zero-truncated dataset suitable for abundance modelling.
#'
#' The dataset can be used to illustrate lagged associations between
#' environmental variables and mosquito abundance, for example in
#' conjunction with the functions \code{aggregate_lagged_intervals()},
#' \code{fit_models_by_lag()}, and \code{plotCCM()}.
#'
#' @source
#' GBIF occurrence data:
#' \doi{10.15468/4qafbu}
#'
#' @docType data
#' @keywords datasets
"albopictusMPL2023"

#' Daily meteorological conditions in Montpellier (2023)
#'
#' This dataset contains daily meteorological variables for Montpellier
#' (Fréjorgues Airport, WMO station 07643) during the year 2023. The data
#' were derived from 3-hourly SYNOP observations provided by Météo-France
#' and aggregated to daily summaries.
#'
#' The dataset includes temperature, humidity, wind speed, atmospheric
#' pressure and precipitation, and is intended for use in studies of
#' lagged associations between environmental conditions and ecological or
#' epidemiological time series.
#'
#' @format A data.frame with one row per day and the following variables:
#' \describe{
#'   \item{date}{Date of observation (class \code{Date}).}
#'   \item{wind_mean}{Daily mean wind speed (m·s⁻¹).}
#'   \item{wind_min}{Daily minimum wind speed (m·s⁻¹).}
#'   \item{wind_max}{Daily maximum wind speed (m·s⁻¹).}
#'   \item{temp_mean}{Daily mean air temperature (°C).}
#'   \item{temp_min}{Daily minimum air temperature (°C).}
#'   \item{temp_max}{Daily maximum air temperature (°C).}
#'   \item{dew.p_mean}{Daily mean dew point temperature (°C).}
#'   \item{dew.p_min}{Daily minimum dew point temperature (°C).}
#'   \item{dew.p_max}{Daily maximum dew point temperature (°C).}
#'   \item{rh_mean}{Daily mean relative humidity (\%).}
#'   \item{rh_min}{Daily minimum relative humidity (\%).}
#'   \item{rh_max}{Daily maximum relative humidity (\%).}
#'   \item{pres_mean}{Daily mean atmospheric pressure (Pa).}
#'   \item{pres_min}{Daily minimum atmospheric pressure (Pa).}
#'   \item{pres_max}{Daily maximum atmospheric pressure (Pa).}
#'   \item{rain_mean}{Daily mean precipitation (mm).}
#'   \item{rain_min}{Daily minimum precipitation (mm).}
#'   \item{rain_max}{Daily maximum precipitation (mm).}
#'   \item{rain_sum}{Daily total precipitation (mm).}
#' }
#'
#' @details
#' Original 3-hourly SYNOP observations were filtered to retain data from
#' Montpellier–Fréjorgues Airport (WIGOS station 0-20000-0-07643). Air temperature and
#' dew point temperature were converted from Kelvin to degrees Celsius.
#' Negative precipitation values were set to zero prior to aggregation.
#'
#' Daily statistics (mean, minimum and maximum) were computed for all
#' variables. For precipitation, daily totals are also
#' provided.
#'
#' This dataset is designed to be used in combination with
#' \code{aggregate_lagged_intervals()} to generate lagged environmental
#' predictors for ecological or epidemiological modelling.
#'
#' @source
#' Météo-France SYNOP data, distributed via data.gouv.fr:
#' \url{https://meteo.data.gouv.fr/datasets/686f8595b351c06a3a790867}
#'
#' @docType data
#' @keywords datasets
"meteoMPL2023"
