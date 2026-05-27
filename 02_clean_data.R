library(tidyverse)
library(lubridate)
library(purrr)

sites <- tibble(
  site = c("CORA1", "DALO1", "CRNR1", "NORO1", "AMLA1","SHLA1", "GSLA1"),
  lon = c(-110.6354, -110.64717, -110.63882, -110.63902, -110.64412, -110.64197, -110.64065),
  lat = c(43.90172, 43.90791, 43.90151, 43.90642, 43.90422, 43.90194, 43.90496)
  
)

clean_species <- function(x) {
  x %>%
    stringr::str_split("/") %>%
    purrr::map_chr(1) %>%
    stringr::str_trim()
}

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

#saveRDS(data_clean, "data_clean.rds")