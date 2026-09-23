library(sf)
library(terra)
library(exactextractr)
library(tidyverse)

nlcd_forest <- terra::rast('data/Annual_NLCD_LndCov_2021_CU_C1V1.tif')
stopifnot("NLCD raster failed to load" = inherits(nlcd_forest, "SpatRaster"))

#This function takes in a raster file and buffer value and then extracts percent forest
# and non-forest values. 


add_habitat <- function(data, forest_raster, buffer = 50, cache_path = "data/data_out.rds") {

  stopifnot(
    "add_habitat: forest_raster must be a SpatRaster" = inherits(forest_raster, "SpatRaster"),
    "add_habitat: input has no rows"                  = nrow(data) > 0,
    "add_habitat: input needs site/lon/lat columns"   = all(c("site","lon","lat") %in% names(data))
  )

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
  
  pts_buf <- st_buffer(pts, dist = buffer)
  
  openness <- exact_extract(forest_raster, pts_buf, function(values, coverage_fraction) {
    
    tibble(
      pct_forest = mean(values %in% c(41,42,43), na.rm = TRUE),
      pct_developed = mean(values %in% c(21, 22, 23, 24), na.rm = TRUE),
      pct_water = mean(values %in% 11, na.rm = TRUE),
      pct_wetland = mean(values %in% c(90,95), na.rm = TRUE),
      pct_shrub = mean(values %in% c(51,52,71), na.rm = TRUE),
      pct_nonforest = 1 - pct_forest - pct_shrub
    )

  })

  openness_df <- dplyr::bind_rows(openness)

  # sanity: every habitat proportion must be a valid fraction in [0, 1]
  .pcts <- c("pct_forest","pct_developed","pct_water","pct_wetland","pct_shrub","pct_nonforest")
  stopifnot(
    "add_habitat: habitat proportions must lie in [0,1]" =
      all(vapply(.pcts, function(cc) all(openness_df[[cc]] >= 0 & openness_df[[cc]] <= 1, na.rm = TRUE), logical(1)))
  )

  sites_env <- dplyr::bind_cols(
    sites, openness_df
  )

  data_out <- data %>%
    left_join(sites_env, by = "site")

  stopifnot(
    "add_habitat: join changed the row count"      = nrow(data_out) == nrow(data),
    "add_habitat: some sites got no habitat values" = all(!is.na(data_out$pct_nonforest))
  )

  saveRDS(data_out, cache_path)

  return(data_out)
}

#This funciton uses the date and time to extract the moon characteristics at the point of
#interest.

add_moonlight <- function(data,
                          timezone = "America/Denver") {
  library(dplyr)
  library(moonlit)
  library(purrr)
  
  site_night <- data %>%
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

  # sanity: the join must not drop/duplicate rows, and moon values must be sane
  stopifnot(
    "add_moonlight: join changed the row count"       = nrow(data_out) == nrow(data),
    "add_moonlight: moon phase must be a fraction in [0,1]" =
      all(data_out$mean_phase >= 0 & data_out$mean_phase <= 1, na.rm = TRUE),
    "add_moonlight: moonlight intensity must be non-negative" =
      all(data_out$mean_moonlight >= 0, na.rm = TRUE) && all(data_out$max_moonlight >= 0, na.rm = TRUE)
  )

  return(data_out)

}
               