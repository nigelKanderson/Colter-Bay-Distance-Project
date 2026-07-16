# =============================================================================
# panel_subfigs.R
# Shared builder for the "generated" (non-model) panel sub-figures, so the
# raster montage (28_figure_panels.R) and the vector patchwork panels
# (30_figure_panels_vector.R) draw the SAME figures from one definition.
#
# build_panel_subfigs() returns a named list of ggplot objects: g1..g7.
# =============================================================================

suppressMessages({ library(tidyverse) })

build_panel_subfigs <- function() {
  col_red <- "#7A1C2E"; col_white <- "#C9A84C"; col_navy <- "#1B2A4A"
  col_purple <- "#4A1259"; col_ivory <- "#F5F0E3"
  intensity_pal <- c("10" = col_ivory, "30" = col_white, "50" = col_purple,
                     "70" = col_red, "100" = col_navy)
  panel_theme <- theme_classic(base_size = 16, base_family = "Arial") +
    theme(panel.grid.major.y = element_line(color = "#EBEBEB", linewidth = 0.4),
          plot.title    = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(color = "#666666", size = 10),
          axis.title    = element_text(size = 18, color = "#222222"),
          axis.text     = element_text(size = 16, color = "#333333"),
          legend.title  = element_text(size = 13),
          legend.text   = element_text(size = 13),
          legend.position = "top")

  moth <- readRDS("data/insect_data_out.rds") %>%
    mutate(intensity = factor(intensity, levels = c(10, 30, 50, 70, 100)))
  site_dist <- moth %>% distinct(site, dist_km)

  ## G1  moth detections by intensity
  g1 <- moth %>%
    group_by(intensity) %>%
    summarise(mean_det = mean(detections), se = sd(detections)/sqrt(n()), .groups = "drop") %>%
    ggplot(aes(intensity, mean_det, fill = intensity)) +
    geom_col(color = "grey30") +
    geom_errorbar(aes(ymin = mean_det - se, ymax = mean_det + se), width = 0.2) +
    scale_fill_manual(values = intensity_pal, guide = "none") +
    labs(title = "Moth detections by intensity", x = "Light intensity (%)",
         y = "Mean detections per site-night") + panel_theme

  ## G2  moth detections by distance
  g2 <- moth %>%
    group_by(site, dist_km) %>%
    summarise(mean_det = mean(detections), .groups = "drop") %>%
    ggplot(aes(dist_km, mean_det)) +
    geom_smooth(method = "lm", se = TRUE, color = col_purple, fill = "#E7DDEE") +
    geom_point(size = 3, color = "#333333") +
    labs(title = "Moth detections by distance", x = "Distance from Colter Bay parking lot (km)",
         y = "Mean detections per site-night") + panel_theme

  ## Moth family data (for panel 6)
  source("R/13_import_insects.R")
  fam <- import_insects("data/grte_distance_insectID.xlsx")$family %>%
    mutate(intensity = factor(intensity, levels = c(10, 30, 50, 70, 100))) %>%
    left_join(site_dist, by = "site")
  top_fam <- fam %>% count(Family, wt = detections, sort = TRUE) %>% slice_head(n = 6) %>% pull(Family)
  famt <- fam %>% filter(Family %in% top_fam)
  fam_pal <- setNames(
    c(col_red, col_navy, col_purple, col_white, "#2E86AB", "#6B8E23")[seq_along(top_fam)],
    top_fam)

  ## G3  family detections by intensity
  g3 <- famt %>%
    group_by(Family, intensity) %>%
    summarise(mean_det = mean(detections), .groups = "drop") %>%
    ggplot(aes(intensity, mean_det, color = Family, group = Family)) +
    geom_line(linewidth = 1) + geom_point(size = 2.5) +
    scale_color_manual(values = fam_pal, name = NULL) +
    labs(title = "Moth family detections by intensity", x = "Light intensity (%)",
         y = "Mean detections per event") + panel_theme

  ## G4  family detections by distance
  g4 <- famt %>%
    group_by(Family, site, dist_km) %>%
    summarise(mean_det = mean(detections), .groups = "drop") %>%
    ggplot(aes(dist_km, mean_det, color = Family)) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    scale_color_manual(values = fam_pal, name = NULL) +
    labs(title = "Moth family detections by distance", x = "Distance from Colter Bay parking lot (km)",
         y = "Mean detections per event") + panel_theme

  ## G4b  family activity heatmap: z-scored within family, sites ordered by distance
  site_km  <- site_dist %>% arrange(dist_km)
  fam_heat <- famt %>%
    group_by(Family, site) %>%
    summarise(mean_det = mean(detections), .groups = "drop") %>%
    group_by(Family) %>%
    mutate(z = as.numeric(scale(mean_det))) %>%
    ungroup() %>%
    mutate(site   = factor(site, levels = site_km$site),
           Family = fct_reorder(Family, mean_det, .fun = sum))
  g4b <- ggplot(fam_heat, aes(site, Family, fill = z)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient2(low = col_red, mid = col_ivory, high = col_navy,
                         midpoint = 0, name = "Relative\nactivity (z)") +
    scale_x_discrete(labels = setNames(paste0(round(site_km$dist_km, 2), " km"),
                                       site_km$site)) +
    labs(title = "Moth family activity by site",
         subtitle = "Z-scored within family; sites ordered by increasing distance from the Colter Bay parking lot",
         x = "Distance from Colter Bay parking lot", y = NULL) +
    panel_theme +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

  ## G5  family detections by moon phase (illuminated fraction via suncalc)
  g5 <- famt %>%
    mutate(moon = suncalc::getMoonIllumination(date)$fraction) %>%
    ggplot(aes(moon, detections, color = Family)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
    scale_color_manual(values = fam_pal, name = NULL) +
    labs(title = "Moth family detections by moon phase",
         x = "Moon illuminated fraction", y = "Detections per event") + panel_theme

  ## G6  single-trait morphology figure (forearm length)
  g6 <- read_csv("output/posthoc_trait_dist_pred.csv", show_col_types = FALSE) %>%
    filter(trait == "forearm_length") %>%
    mutate(trait_level = factor(trait_level, levels = c("Low", "Mean", "High"))) %>%
    ggplot(aes(dist_km, response, color = trait_level, fill = trait_level)) +
    geom_ribbon(aes(ymin = lower_CL, ymax = upper_CL), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = c(Low = col_navy, Mean = col_red, High = col_purple),
                       name = "Forearm length") +
    scale_fill_manual(values = c(Low = col_navy, Mean = col_red, High = col_purple),
                      guide = "none") +
    labs(title = "Bat activity by distance × forearm length",
         subtitle = "Predicted detections ± 95% CI at low/mean/high (±1 SD)",
         x = "Distance from Colter Bay parking lot (km)", y = "Predicted detections per night") +
    panel_theme

  ## G7  human-survey component scores by light condition
  g7 <- read_csv("output/people_pca_component_scores.csv", show_col_types = FALSE) %>%
    filter(StreetlightCondition %in% c(1, 2)) %>%
    mutate(light = factor(StreetlightCondition, labels = c("Red", "White"))) %>%
    pivot_longer(c(PC1_score, PC2_score), names_to = "component", values_to = "score") %>%
    mutate(component = recode(component,
                              PC1_score = "PC1: wildlife &\nexperience benefit",
                              PC2_score = "PC2: discomfort\nin the dark")) %>%
    group_by(light, component) %>%
    summarise(mean = mean(score), se = sd(score)/sqrt(n()), .groups = "drop") %>%
    ggplot(aes(component, mean, fill = light)) +
    geom_col(position = position_dodge(0.8), width = 0.7, color = "grey30") +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(0.8), width = 0.2) +
    scale_fill_manual(values = c(Red = col_red, White = col_white), name = "Light seen") +
    labs(title = "Survey component scores by light condition",
         x = NULL, y = "Mean component score (1-5)") + panel_theme

  list(g1 = g1, g2 = g2, g3 = g3, g4 = g4, g4b = g4b, g5 = g5, g6 = g6, g7 = g7)
}
