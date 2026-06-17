library(readr)
library(readxl)
library(dplyr)
library(lubridate)

lighting <- read_excel("~/Library/CloudStorage/GoogleDrive-nigel_anderson@brown.edu/.shortcut-targets-by-id/1sSdpOAdUOgAVbJGpTKgB3-CJAfnjsvKJ/grandteton_distanceproject/grte_distance_lightingschedule_2022.xlsx")

names(lighting)

lighting <- lighting %>%
  mutate(date = as.Date(date))

saveRDS(lighting, "data_full.rds")
