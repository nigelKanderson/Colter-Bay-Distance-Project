# =============================================================================
# 31_alt_distance.R
# Alternative "nearest-source" distance gradient (to try alongside the existing
# single-reference dist_km).
#
# Three sites are treated as light sources / parking lot (distance = 0):
#   GRTE17, GRTE16, GRTE01
# Every other site's distance is measured to an ASSIGNED source site (GRTE01 or
# GRTE16), per the project owner's mapping below.
#
# Sourcing this file:
#   - builds `site_dist_alt` (site-level lookup) and caches it to
#     data/site_dist_alt.csv
#   - defines add_alt_distance(df): attaches a `dist_km_alt` column to any data
#     frame that has a `site` column (GRTE01 … GRTE17 style).
#
# Usage later, e.g.:
#   source("R/31_alt_distance.R")
#   data_env <- add_alt_distance(data_env)   # then model with dist_km_alt
# =============================================================================

suppressMessages({ library(dplyr) })

# site -> assigned reference site; NA = the site is itself a source (dist 0)
.alt_ref_map <- c(
  GRTE17 = NA_character_, GRTE16 = NA_character_, GRTE01 = NA_character_,
  GRTE05 = "GRTE16", GRTE06 = "GRTE16", GRTE08 = "GRTE16", GRTE11 = "GRTE16",
  GRTE02 = "GRTE01", GRTE03 = "GRTE01", GRTE09 = "GRTE01", GRTE10 = "GRTE01",
  GRTE15 = "GRTE01", GRTE14 = "GRTE01", GRTE13 = "GRTE01"
)

# Measured distances (meters) to the assigned source, provided by project owner.
# Sources (GRTE17/16/01) = 0.
.alt_dist_m <- c(
  GRTE17 = 0,        GRTE16 = 0,        GRTE01 = 0,
  GRTE05 = 115.872,  GRTE06 = 398.01,   GRTE08 = 1188.111, GRTE11 = 1646.61,
  GRTE02 = 469.13,   GRTE03 = 1085.588, GRTE09 = 1412.75,  GRTE10 = 1587.811,
  GRTE15 = 6685.86,  GRTE14 = 6886.48,  GRTE13 = 6952.187
)

build_alt_distance_lookup <- function() {
  sites <- names(.alt_dist_m)
  tibble(site = sites,
         assigned_ref = ifelse(is.na(.alt_ref_map[sites]), "(source)",
                               .alt_ref_map[sites]),
         dist_m_alt  = unname(.alt_dist_m[sites]),
         dist_km_alt = unname(.alt_dist_m[sites]) / 1000) %>%
    arrange(dist_km_alt)
}

site_dist_alt <- build_alt_distance_lookup()
dir.create("data", showWarnings = FALSE)
readr::write_csv(site_dist_alt, "data/site_dist_alt.csv")

# Attach dist_km_alt to any df with a `site` column (robust to factor/character)
add_alt_distance <- function(df) {
  lk <- setNames(site_dist_alt$dist_km_alt, site_dist_alt$site)
  df$dist_km_alt <- unname(lk[as.character(df$site)])
  df
}

message("Alt-distance lookup built for ", nrow(site_dist_alt),
        " sites -> data/site_dist_alt.csv")
