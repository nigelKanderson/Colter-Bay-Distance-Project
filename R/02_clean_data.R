library(tidyverse)
library(lubridate)
library(purrr)

#First I make a table of the sites and coordinates

sites <- tibble(
  site = c("GRTE01", "GRTE02", "GRTE03", "GRTE05","GRTE06", "GRTE08", "GRTE09", "GRTE10", "GRTE11", "GRTE13", "GRTE14", "GRTE15", "GRTE16", "GRTE17"), #, "GRTF02"),
  lon = c(-110.64153, -110.64346, -110.64275, -110.63919, -110.63607, -110.62608, -110.6221, -110.620952, -110.62, -110.60068, -110.60471, -110.60903, -110.6405, -110.645134), #, -110.64346),
  lat = c(43.90109, 43.89614, 43.89061, 43.90586, 43.90257, 43.90206, 43.90013, 43.897718, 43.90569, 43.84588, 43.84504, 43.84558, 43.90519, 43.904879), #, 43.89614)
  
)

#This function ensures that the species name is clean.

clean_species <- function(x) {
  x %>%
    stringr::str_split("/") %>%
    purrr::map_chr(1) %>%
    stringr::str_trim()
}

#This function filters out rows with NA's and ensures the variables are the right class.

clean_data <- function(data) {
  
  data_clean <- data_raw %>%
    
    mutate(
      species = clean_species(species)
    ) %>%
    
    filter(
      !is.na(detections),
      !is.na(site),
      !is.na(date),
    ) %>%
    
    left_join(sites, by = 'site') %>%
    
    mutate(
      detections = as.numeric(detections),
      weighted_detections = as.numeric(weighted_detections),
      year = as.integer(year),
      jd = as.integer(jd)
    ) %>%
    
    filter(
      !is.na(jd),
      !is.na(lon),
      !is.na(lat)
    ) %>%
    
    arrange(site, date, species)
  
  return(data_clean)
}

