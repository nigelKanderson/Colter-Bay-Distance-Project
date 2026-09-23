library(tidyverse)
library(lubridate)
library(purrr)

#First I make a table of the sites and coordinates

sites <- tibble(
  site = c("GRTE01", "GRTE02", "GRTE03", "GRTE05","GRTE06", "GRTE08", "GRTE09", "GRTE10", "GRTE11", "GRTE13", "GRTE14", "GRTE15", "GRTE16", "GRTE17"), #, "GRTF02"),
  lon = c(-110.64153, -110.64346, -110.64275, -110.63919, -110.63607, -110.62608, -110.6221, -110.620952, -110.62, -110.60068, -110.60471, -110.60903, -110.6405, -110.645134), #, -110.64346),
  lat = c(43.90109, 43.89614, 43.89061, 43.90586, 43.90257, 43.90206, 43.90013, 43.897718, 43.90569, 43.84588, 43.84504, 43.84558, 43.90519, 43.904879), #, 43.89614)

)

# sanity: the hard-coded site reference table must be well-formed and in-region
stopifnot(
  "sites: duplicate site codes"       = !any(duplicated(sites$site)),
  "sites: coordinates missing"        = all(!is.na(sites$lon)) && all(!is.na(sites$lat)),
  "sites: longitudes outside GRTE"    = all(sites$lon > -110.70 & sites$lon < -110.55),
  "sites: latitudes outside GRTE"     = all(sites$lat >  43.80  & sites$lat <  43.95)
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

  # Species labels (single or compound) are finalised in 01_import_sonobat.R,
  # so we do NOT collapse compounds here. Operate on the passed `data`.
  data_clean <- data %>%

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

  # surface silently-dropped rows from unknown site codes (not in the sites table)
  .unknown <- setdiff(unique(data$site), sites$site)
  if (length(.unknown) > 0)
    warning("clean_data: dropped rows for site codes not in the sites table: ",
            paste(.unknown, collapse = ", "))

  # sanity: the cleaned output must be usable by the downstream pipeline
  stopifnot(
    "clean_data: output has no rows"                 = nrow(data_clean) > 0,
    "clean_data: missing required columns"           =
      all(c("site","date","species","detections","lon","lat","jd","year") %in% names(data_clean)),
    "clean_data: detections must be non-negative"    = all(data_clean$detections >= 0, na.rm = TRUE),
    "clean_data: key fields must be complete"        =
      all(!is.na(data_clean$detections)) && all(!is.na(data_clean$date)) &&
      all(!is.na(data_clean$jd)) && all(!is.na(data_clean$lon)) && all(!is.na(data_clean$lat)),
    "clean_data: coordinates fell outside GRTE"      =
      all(data_clean$lon > -110.70 & data_clean$lon < -110.55) &&
      all(data_clean$lat >  43.80  & data_clean$lat <  43.95),
    "clean_data: Julian day out of range"            = all(data_clean$jd >= 1 & data_clean$jd <= 366, na.rm = TRUE)
  )

  return(data_clean)
}

