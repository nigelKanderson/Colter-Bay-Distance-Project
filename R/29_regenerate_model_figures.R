# =============================================================================
# 29_regenerate_model_figures.R
# Regenerate the model-based figures used in the panels with the unified
# theme (theme_classic + Arial + larger axis fonts), WITHOUT the Google Drive
# import step. Reconstructs data_env / insect_total_env from the cached
# intermediates (data/data_out.rds = cleaned + habitat; data/data_full.rds =
# lighting schedule; data/insect_data_out.rds), re-derives moonlight/distance/
# sky/traits exactly as the qmd does, refits the two general models, and
# re-runs 08_figures.R (bat), 15 (moth), 17 (bat vs moth).
#
# 08_figures.R defines its own bat_theme internally (already edited to the
# unified style); 15 and 17 receive bat_theme as an argument, defined below.
# =============================================================================

suppressMessages({
  library(tidyverse); library(readxl); library(sf); library(glmmTMB)
})

# ---- unified theme + palette (matches the edited bat_theme in 08_figures.R) --
bat_theme <- theme_classic(base_size = 15, base_family = "Arial") +
  theme(panel.grid.major.y = element_line(color = "#EBEBEB", linewidth = 0.4),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(color = "#444444", linewidth = 0.5),
        axis.text  = element_text(color = "#333333", size = 16),
        axis.title = element_text(color = "#222222", size = 18),
        plot.title = element_text(color = "#111111", size = 14, face = "bold"),
        plot.subtitle = element_text(color = "#666666", size = 10),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold"),
        strip.text = element_text(face = "bold", size = 13))
col_red <- "#7A1C2E"; col_white <- "#C9A84C"; col_navy <- "#1B2A4A"
col_purple <- "#4A1259"; col_ivory <- "#F5F0E3"
intensity_pal <- c("10" = col_ivory, "30" = col_white, "50" = col_purple,
                   "70" = col_red, "100" = col_navy)

source("R/03_extract_habitat.R")   # add_moonlight()
lighting <- readRDS("data/data_full.rds")
sky <- read_csv("data/sky_brightness_measurements.csv", show_col_types = FALSE) %>%
  mutate(site = toupper(site))
sky_cols <- intersect(c("brightness_dark", "brightness_white100",
                        "brightness_dark_median", "brightness_white100_median"),
                      names(sky))

ref_point <- st_as_sf(data.frame(lon = -110.645134, lat = 43.904879),
                      coords = c("lon", "lat"), crs = 4326) %>% st_transform(32612)
site_distance <- function(df) {
  s <- df %>% distinct(site, lon = lon.x, lat = lat.x)
  s_sf <- st_as_sf(s, coords = c("lon", "lat"), crs = 4326) %>% st_transform(32612)
  s$dist_m <- as.numeric(st_distance(s_sf, ref_point))
  s %>% select(site, dist_m)
}

# ---- BAT data_env ------------------------------------------------------------
message("Reconstructing bat data_env ...")
data <- readRDS("data/data_out.rds")
data <- add_moonlight(data)                                   # adds mean_phase
data_env <- data %>%
  left_join(lighting %>% select(date, color, intensity), by = "date") %>%
  left_join(site_distance(data), by = "site") %>%
  mutate(dist_km = dist_m / 1000,
         dist_cat = case_when(dist_km < 0.52 ~ "Close", dist_km < 1 ~ "Medium",
                              dist_km < 2.1 ~ "Further", TRUE ~ "Far"),
         jd_c = jd - mean(jd, na.rm = TRUE)) %>%
  filter(!is.na(dist_km))
focal_spp <- c("Laci", "Lano", "Mylu", "Epfu", "Myev", "Myvo")
data_env <- data_env %>% filter(species %in% focal_spp) %>%
  left_join(sky %>% select(site, all_of(sky_cols)), by = "site")
bat_traits <- read_excel("data/bat_traits.xlsx") %>% rename(species = Species)
data_env <- data_env %>% left_join(bat_traits, by = "species") %>%
  mutate(color = factor(color), site = factor(site), intensity = factor(intensity))

message("Fitting bat model ...")
simple_model1 <- glmmTMB(
  detections ~ color * intensity + intensity * dist_km + mean_phase +
    pct_nonforest + brightness_dark + jd_c + I(jd_c^2) + (1 | site),
  data = data_env, ziformula = ~1, family = nbinom2())

message("Regenerating bat figures (08_figures.R) ...")
source("R/08_figures.R")   # uses data_env + simple_model1; saves figNN PNGs

# ---- MOTH insect_total_env ---------------------------------------------------
message("Reconstructing moth insect_total_env ...")
insect <- readRDS("data/insect_data_out.rds")
insect$datetime <- as.POSIXct(paste(insect$date, "12:00:00"), tz = "America/Denver")
insect <- add_moonlight(insect)
insect_total_env <- insect %>%
  left_join(sky %>% select(site, brightness_dark), by = "site") %>%
  mutate(jd_c = jd - mean(jd, na.rm = TRUE),
         color = factor(color), site = factor(site), intensity = factor(intensity))

message("Fitting moth model + figures (15) ...")
source("R/15_insect_models.R")
insect_general <- run_insect_general_model(insect_total_env, bat_theme, intensity_pal)

message("Regenerating bat-vs-moth relative-effect figures (17) ...")
source("R/17_bat_moth_effects.R")
bm_effects <- compare_bat_moth_effects(simple_model1, insect_general$model,
                                       data_env, insect_total_env,
                                       bat_theme, col_navy, col_purple)

message("Done regenerating model figures.")
