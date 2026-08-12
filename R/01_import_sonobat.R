#| warning: false
#| message: false

library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)
library(purrr)
library(stringr)

# -----------------------------------------------------------------------------
# SonoBat v30.1 detections for the distance project.
# The complete v30.1 set is staged locally (matches the prior analysis cache and
# avoids the flaky Google Drive mount). To re-pull from Drive instead, point
# `sonobat_dir` at the Drive .../grandteton_distanceproject v30.1 folder.
# -----------------------------------------------------------------------------
sonobat_dir <- "data/sonobat_compare/v30"

files <- list.files(sonobat_dir, pattern = "\\.xlsx$", full.names = TRUE,
                    recursive = TRUE)
files <- files[grepl("v30\\.1", basename(files), ignore.case = TRUE)]

if (length(files) == 0)
  stop("No v30.1 SonoBatch .xlsx files found in '", sonobat_dir,
       "'. Is the folder present / the Google Drive mounted?")
message("Importing ", length(files), " v30.1 SonoBatch file(s) from ", sonobat_dir)

# This function reads one file, applies the species/probability rule, extracts
# metadata from the filename, and aggregates to detection level.
#
# Probability rule (v30.1): SppAccp can be a single species ("Mylu") or a
# compound of candidates ("Lano/Epfu") with matching vote-share probabilities
# ("0.32/0.68") that sum to ~1 within a call. Read Prob as text and parse it
# numerically BEFORE thresholding (a numeric read + the old text comparison had
# let ambiguous calls slip through). Keep a recording when the COMBINED
# probability mass is >= 0.9 -- for a single-species call that is just its own
# probability; for a compound it is essentially always true -- and keep the
# recording's full (possibly compound) species label. Blank / NoID accepts are
# dropped.
process_file <- function(file) {
  df <- tryCatch(
    readxl::read_excel(file, col_types = "text"),
    error = function(e) return(NULL)
  )
  if (is.null(df)) return(NULL)

  required_cols <- c("Prob", "SppAccp", "Filename")
  if (!all(required_cols %in% names(df))) return(NULL)

  # --- parse compound species/probabilities ---------------------------------
  spp <- stringr::str_split(df$SppAccp, "/")
  prb <- lapply(stringr::str_split(df$Prob, "/"),
                function(x) suppressWarnings(as.numeric(x)))

  df$combined_prob <- vapply(prb, function(p) sum(p, na.rm = TRUE), numeric(1))
  # Keep the full (possibly compound) label, dropping only empty tokens such as
  # the trailing slash in "Laci/". Single-candidate calls collapse to one name.
  df$label <- vapply(spp, function(s) {
    s <- stringr::str_trim(s); keep <- !is.na(s) & s != ""
    if (!any(keep)) return(NA_character_)
    paste(s[keep], collapse = "/")
  }, character(1))

  df <- df %>%
    filter(!is.na(label), label != "NoID", combined_prob >= 0.9) %>%
    mutate(SppAccp = label, Prob = combined_prob)

  if (nrow(df) == 0) return(NULL)

  df %>%
    mutate(

      Prob = as.numeric(Prob),

      file_name_clean =
        tools::file_path_sans_ext(Filename),

      site =
        stringr::word(file_name_clean, 1, sep = "_"),

      site = recode(site, "GRTF02" = "GRTE02"),

      date_string =
        stringr::str_extract(file_name_clean, "\\d{8}"),

      time_string =
        stringr::str_extract(file_name_clean, "\\d{6}$"),

      time = lubridate::hms(
        paste0(
          substr(time_string, 1, 2), ":",
          substr(time_string, 3, 4), ":",
          substr(time_string, 5, 6)
        )
      ),

      time_clean = str_replace_all(time, "H|M|S", " "),
      time_clean = str_squish(time_clean),
      hours = as.numeric(str_extract(time, "\\d+(?=H)")),
      minutes = as.numeric(str_extract(time, "\\d+(?=M)")),
      seconds = as.numeric(str_extract(time, "\\d+(?=S)")),

      hours = ifelse(is.na(hours), 0, hours),
      minutes = ifelse(is.na(minutes), 0, minutes),
      seconds = ifelse(is.na(seconds), 0, seconds),

      time_sec = hours * 3600 + minutes * 60 + seconds,

      date =
        lubridate::ymd(date_string),

      datetime =
        as.POSIXct(
          paste(date, time),
          tz = "UTC"
        ),

      year =
        lubridate::year(date),

      month =
        lubridate::month(date),

      day =
        lubridate::day(date),

      jd = lubridate::yday(date),
      jd2 = jd^2

    ) %>%

    group_by(site, year, date, jd, jd2, SppAccp, datetime, time_sec) %>%

    summarise(

      detections = n(),

      weighted_detections =
        sum(Prob, na.rm = TRUE),

      mean_confidence =
        mean(Prob, na.rm = TRUE),

      .groups = 'drop'

    ) %>%

    rename(species = SppAccp) %>%

    mutate(
      source_file = basename(file)
    )


}


data_raw <- map_dfr(files, process_file)

if (nrow(data_raw) == 0 || !"species" %in% names(data_raw))
  stop("Import produced no usable rows (no 'species' column). ",
       "Check that the input files are SonoBatch outputs with Prob/SppAccp/Filename.")

message("Imported ", nrow(data_raw), " detection rows across ",
        dplyr::n_distinct(data_raw$species), " species from ",
        dplyr::n_distinct(data_raw$site), " sites.")
