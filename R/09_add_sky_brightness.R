library(readr)
library(dplyr)

sky <- read_csv("data/sky_brightness_measurements.csv", show_col_types = FALSE) %>%
  mutate(site = toupper(site))

sky_cols <- intersect(
  c("brightness_dark", "brightness_white100",
    "brightness_dark_median", "brightness_white100_median"),
  names(sky)
)

data_env <- data_env %>%
  left_join(sky %>% select(site, all_of(sky_cols)), by = "site")

n_matched <- sum(!is.na(data_env$brightness_dark))
message("Sky brightness joined: ", n_matched, "/", nrow(data_env), " rows matched")
for (col in sky_cols) {
  vals <- data_env[[col]]
  message("  ", col, " range: ",
          round(min(vals, na.rm=TRUE), 4), " – ", round(max(vals, na.rm=TRUE), 4))
}
