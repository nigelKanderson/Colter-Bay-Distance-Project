# 13_import_insects.R
# Import + parse moth PIC_ID metadata (site_date_treatment, e.g. "GRTE01_2aug22_w100D")
# and aggregate specimen records to total- and family/subfamily-level detection counts.
#
# PIC_ID parsing rules (per project owner):
#  - trailing d/D is a processing marker, not a treatment level -> dropped
#  - light-meter intensity readings that deviate from the bat project's
#    10/30/50/70/100 schedule are rounded to the nearest of those levels
#  - repeat-photo suffixes like " (2)", " (3)", "_b" mark duplicate photos of
#    the same site/date/treatment sampling event -> merged by grouping on
#    site+date+color+intensity (the parsed key drops those suffixes entirely)
#
# Input:  path to grte_distance_insectID.xlsx (identifications on `sheet`)
# Output: list(specimens, total, family, subfamily)

import_insects <- function(path, sheet = "Sheet2") {
  library(readxl)
  library(dplyr)
  library(stringr)

  raw <- read_excel(path, sheet = sheet)

  intensity_levels <- c(10, 30, 50, 70, 100)
  round_intensity  <- function(x) intensity_levels[which.min(abs(x - intensity_levels))]

  # Fix inconsistent case, e.g. "noctuidae" / "Noctuidae"
  fix_case <- function(x) {
    ifelse(is.na(x), NA_character_,
           paste0(toupper(str_sub(x, 1, 1)), tolower(str_sub(x, 2))))
  }

  parts <- str_split_fixed(raw$PIC_ID, "_", 3)

  specimens <- raw %>%
    mutate(
      site          = parts[, 1],
      date          = as.Date(parts[, 2], format = "%d%b%y"),
      treat_raw     = parts[, 3],
      color         = toupper(str_sub(treat_raw, 1, 1)),
      intensity_raw = as.numeric(str_extract(treat_raw, "(?<=^[a-zA-Z])[0-9]+")),
      intensity     = vapply(intensity_raw, round_intensity, numeric(1)),
      Family        = fix_case(Family),
      Subfamily     = fix_case(Subfamily),
      # "Noctuoidea" is a superfamily, not a family -- mislabeled in a few rows
      Family        = ifelse(Family == "Noctuoidea", NA_character_, Family)
    ) %>%
    filter(color %in% c("R", "W"), !is.na(date), !is.na(intensity))

  n_dropped <- nrow(raw) - nrow(specimens)
  if (n_dropped > 0) cat("Dropped", n_dropped, "rows with unparseable PIC_ID\n")
  cat("Parsed", nrow(specimens), "specimen records from",
      n_distinct(raw$PIC_ID), "photos (", n_distinct(specimens$site, specimens$date),
      "site-nights )\n")

  total <- specimens %>%
    group_by(site, date, color, intensity) %>%
    summarise(detections = n(), .groups = "drop")

  family <- specimens %>%
    filter(!is.na(Family)) %>%
    group_by(site, date, color, intensity, Family) %>%
    summarise(detections = n(), .groups = "drop")

  subfamily <- specimens %>%
    filter(!is.na(Subfamily)) %>%
    group_by(site, date, color, intensity, Subfamily) %>%
    summarise(detections = n(), .groups = "drop")

  list(specimens = specimens, total = total, family = family, subfamily = subfamily)
}
