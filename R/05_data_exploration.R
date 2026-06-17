library(tidyverse)
library(janitor)
library(ggplot2)

data <- readRDS('data_out.rds')

glimpse(data)
summary(data)

data %>%
  summarise(
    n = n(),
    n_sites = n_distinct(site),
    n_species = n_distinct(species)
  )

data %>%
  count(site, year) %>%
  arrange(site)

data %>%
  group_by(site, date) %>%
  summarise(detections = sum(detections), .groups = "drop") %>%
  ggplot(aes(date, detections, color = site)) +
  geom_line()

data %>%
  ggplot(aes(jd, weighted_detections)) +
  geom_point(alpha = 0.3) +
  geom_smooth()

data %>%
  ggplot(aes(weighted_detections)) +
  geom_histogram(bins = 50)

data %>%
  group_by(site) %>%
  summarise(
    mean_det = mean(detections),
    mean_wdet = mean(weighted_detections),
    #pct_forest = first(openness$pct_forest)
  ) %>%
  ggplot(aes(site, mean_wdet)) +
  geom_col()

data %>%
  select(weighted_detections, jd, jd2) %>%
  cor(use = 'complete.obs')

data %>%
  summarise(
    
    n = n(),
    n_sites = n_distinct(site),
    n_species = n_distinct(species),
    min_date = min(date),
    max_date = max(date),
    missing_lons = sum(is.na(lon)),
    missing_lats = sum(is.na(lat))
  )


