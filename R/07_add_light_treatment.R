library(readr)
library(readxl)
library(dplyr)
library(lubridate)

lighting <- read_excel("/Users/nanderson/Desktop/grte_distance_lightingschedule_2022.xlsx")

names(lighting)

lighting <- lighting %>%
  mutate(date = as.Date(date))

saveRDS(lighting, "data/data_full.rds")
