library(readr)
library(dplyr)

sky <- read_csv("data/sky_brightness_measurements.csv", show_col_types = FALSE) %>%
  mutate(site = toupper(site))

data_env <- data_env %>%
  left_join(sky %>% select(site, brightness_dark, brightness_white100), by = "site")

n_matched <- sum(!is.na(data_env$brightness_dark))
message("Sky brightness joined: ", n_matched, "/", nrow(data_env), " rows matched")
message("  brightness_dark range:     ", round(min(data_env$brightness_dark,     na.rm=TRUE), 3),
        " – ", round(max(data_env$brightness_dark,     na.rm=TRUE), 3))
message("  brightness_white100 range: ", round(min(data_env$brightness_white100, na.rm=TRUE), 3),
        " – ", round(max(data_env$brightness_white100, na.rm=TRUE), 3))
