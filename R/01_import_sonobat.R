#| warning: false
#| message: false

library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)
library(purrr)
library(stringr)

files_all <- list.files(
  path = "/Users/nanderson/Library/CloudStorage/GoogleDrive-nigel_anderson@brown.edu/.shortcut-targets-by-id/1sSdpOAdUOgAVbJGpTKgB3-CJAfnjsvKJ/grandteton_distanceproject",
  pattern = "\\.xlsx$",
  full.names = TRUE,
  recursive = TRUE
)

# I want to exclude all files that are version 4.4.5 because we want to ensure
#that all files are the same version for sonobat for consistencey/accuracy

excluded <- files_all[grepl('v4.4.5', basename(files_all), ignore.case = TRUE)]

files <- files_all[!grepl('v4.4.5', basename(files_all), ignore.case = TRUE)]

cat(
  "Found ", length(files_all), " Excel files.\n",
  "Excluded ", length(excluded), " v4.4.5 files. \n",
  "Importing ", length(files), " files. \n"
)
  

#This function simply reads in a file and ensures it has the right columns. It also
#filters based on probability of 90% or higher and extracts all of the metadata
#from the file name. 

process_file <- function(file) {
  df <- tryCatch(
    readxl::read_excel(file),
    error = function(e) return(NULL)
  )
  
  if (is.null(df)) return(NULL)
  
  required_cols <- c("Prob","SppAccp","Filename")
  
  if(!all(required_cols %in% names(df))) return(NULL)
  
  df <- df %>%
    filter(Prob >= 0.9)
  
  df <- df %>%
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

write_csv(
  all_counts,
  "data/glmm_ready_data.csv"
)

