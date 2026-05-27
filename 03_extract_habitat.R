library(sf)
library(terra)
library(exactextractr)
library(tidyverse)

nlcd_forest <- terra::rast('Annual_NLCD_LndCov_2021_CU_C1V1.tif')



add_habitat <- function(data, forest_raster, buffer = 500) {
  
  sites <- data %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(
      lon = first(lon),
      lat = first(lat),
      .groups = "drop"
    )
  
  
  pts <- st_as_sf(
    sites,
    coords = c("lon", "lat"),
    crs = 4326
  )
  
  pts <- st_transform(pts, terra::crs(forest_raster))
  
  pts_buf <- st_buffer(pts, dist = 500)
  
  openness <- exact_extract(forest_raster, pts_buf, function(values, coverage_fraction) {
    
    tibble(
      pct_forest = mean(values %in% c(41,42,43), na.rm = TRUE),
      pct_nonforest = 1 - pct_forest
      #pct_open = mean(values %in% c(31,41,42,43,52,71,81,82), na.rm = TRUE),
      #pct_developed = mean(values %in% c(21:24), na.rm = TRUE),
    ) 
    
  })
  
  openness_df <- dplyr::bind_rows(openness)
  
  sites_env <- dplyr::bind_cols(
    sites, openness_df
  )
  
  data_out <- data %>%
    left_join(sites_env, by = "site")
  
  saveRDS(data_out, "data_out.rds")
  
  return(data_out)
}

add_moonlight <- function(data,
                          timezone = "America/Denver") {
  library(dplyr)
  library(moonlit)
  library(purrr)
  
  site_night <- data_env %>%
    mutate(date = as.Date(datetime)) %>%
    distinct(site, date, lat.y, lon.y)
  
  moon_env <- site_night %>%
    rowwise () %>%
    mutate(
      moon_stats = list(
        calculateMoonlightStatistics(
          lat = lat.y,
          lon = lon.y,
          e = 0.16,
          date = as.POSIXct(date, tz = timezone),
          timezone = timezone,
          t = "15 mins"
        )
      )
    ) %>%
    
    ungroup()

  print(names(site_night))
  
  moon_env <- moon_env %>%
    mutate(
      mean_moonlight = map_dbl(moon_stats, "meanMoonlightIntensity"),
      max_moonlight = map_dbl(moon_stats, "maxMoonlightIntensity"),
      mean_phase = map_dbl(moon_stats, "meanMoonPhase")
    ) %>%
    
    select(
      site,
      date,
      mean_moonlight,
      max_moonlight,
      mean_phase
    )
  
  data_out <- data %>%
    mutate(date = as.Date(datetime)) %>%
    left_join(moon_env, by = c("site", "date"))
  
  return(data_out)
    
}
               