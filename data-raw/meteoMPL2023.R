## code to prepare `meteoMPL2023` dataset

# meteoMPL2023----
library(tidyverse)

# load meteo data from the SYNOP network (downloaded at https://meteo.data.gouv.fr/datasets/686f8595b351c06a3a790867) and unziped
# 3-hr data for several station in France

meteoMPL2023_3h <- read_delim('inst/extdata/synop_2023.csv', delim = ";",na = c("", "NA", "mq")) %>%
  filter(geo_id_wmo == "07643") %>% # Montpellier (Frejorgues) Airport
  select(validity_time, ff, t, td, u, pres, rr3) %>% # useful variables (ff: wind speed, t: temperature, td: dew point, u: relative humidity, pres: atmo pressure and rr3: precipitation)
  mutate(wind = as.numeric(ff),
         temp = as.numeric(t) -273.15 , # convert K to celsius
         dew.p = as.numeric(td) -273.15,
         rh = as.integer(u),
         pres = as.integer(pres),
         rain = as.numeric(rr3)) %>%
  select(-c(ff,t,td,u,rr3)) %>%
  mutate(rain = ifelse(rain<0,0,rain)) # some data are -0.1, converted to 0


# compute daily statistics (mean, min and max + sum for rainfall)
meteoMPL2023 <- meteoMPL2023_3h %>%
  mutate(date = date(validity_time)) %>%
  group_by(date) %>%
  summarise(across(!(validity_time),list(mean=mean, min = min, max=max)), rain_sum = sum(rain, na.rm = TRUE)) %>%
  ungroup()

# save
usethis::use_data(meteoMPL2023, overwrite = TRUE)
