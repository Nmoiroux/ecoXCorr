## code to prepare `albopictusMPL2023` dataset

# entomological sampling in Montpellier
# (from GBIF: https://doi.org/10.15468/4qafbu) downloaded as GBIF anotated archive in .csv format, then unzipped
# female Aedes albopictus dataset

library(tidyverse)

albopictusMPL2023 <- read_tsv('inst/extdata/0004826-260129131611470.csv') %>%
  select(species,individualCount,eventDate,occurrenceID) %>%
  mutate(sex = str_sub(occurrenceID,nchar(occurrenceID),nchar(occurrenceID))) %>%
  mutate(trap = str_split_i(occurrenceID,"_",1)) %>% # extract trap ID
  mutate(area = str_split_i(occurrenceID,"_",2)) %>% # extract area ID where trap is located
  mutate(area = if_else(area=="ZPLSou","ZPSou",area)) %>% # corrects ID
  filter(species == "Aedes albopictus", eventDate != "2023", sex == "F") %>% # filter the dataset for female Ae. aedes albopictus
  select(-c(sex, occurrenceID)) %>%
  mutate(date = str_split_i(eventDate,"/",1) %>% ymd_hm() %>% date()) %>% # extract sampling date
  filter(individualCount>0) # zero-truncated dataset


usethis::use_data(albopictusMPL2023, overwrite = TRUE)
