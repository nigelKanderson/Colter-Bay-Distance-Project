library(ggplot2)
library(dplyr)
library(patchwork)
library(forcats)

# ==============================================================================
# ── Shared theme (matches ColterBay_Practice style) ───────────────────────────
# ==============================================================================

bat_theme <- theme_classic(base_size = 12) +
  theme(
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "#EBEBEB", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    axis.line          = element_line(color = "#444444", linewidth = 0.5),
    axis.ticks         = element_line(color = "#444444"),
    axis.text          = element_text(color = "#333333", size = 10),
    axis.title         = element_text(color = "#222222", size = 11),
    plot.title         = element_text(color = "#111111", size = 13,
                                      face = "bold", margin = margin(b = 3)),
    plot.subtitle      = element_text(color = "#666666", size = 9,
                                      margin = margin(b = 8)),
    legend.background  = element_rect(fill = "white", color = NA),
    legend.key         = element_rect(fill = "white", color = NA),
    legend.text        = element_text(color = "#333333", size = 10),
    legend.title       = element_text(color = "#444444", size = 9,
                                      face = "bold"),
    strip.background   = element_rect(fill = "#F2F2F2", color = NA),
    strip.text         = element_text(color = "#222222", face = "bold",
                                      size = 10),
    plot.caption       = element_text(color = "#999999", size = 8,
                                      hjust = 0, margin = margin(t = 8)),
    plot.margin        = margin(12, 14, 10, 12)
  )

# Catholic-inspired palette
col_red    <- "#7A1C2E"   # cardinal crimson
col_white  <- "#C9A84C"   # gilded gold
col_navy   <- "#1B2A4A"   # nave navy
col_purple <- "#4A1259"   # papal purple
col_ivory  <- "#F5F0E3"   # illuminated manuscript ivory

intensity_pal <- c("10"  = col_ivory,
                   "30"  = col_white,
                   "50"  = col_purple,
                   "70"  = col_red,
                   "100" = col_navy)

# Distance category order
dist_order <- c("Close", "Medium", "Further", "Far")

# Species to include (matches ColterBay project)
focal_spp <- c("Laci", "Lano", "Mylu", "Epfu", "Myev", "Myvo", "Myyu", "Myci")

# Ensure jd_c exists (centered Julian day)
if (!"jd_c" %in% names(data_env)) {
  data_env$jd_c <- data_env$jd - mean(data_env$jd, na.rm = TRUE)
}

# ==============================================================================
# ── Fig 1: Red vs White — headline color effect (p < 0.001) ──────────────────
# ==============================================================================

color_summary <- data_env %>%
  group_by(color) %>%
  summarise(
    mean_det = mean(detections, na.rm = TRUE),
    se       = sd(detections,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

fig1 <- ggplot(color_summary,
               aes(x = color, y = mean_det, fill = color)) +
  geom_col(width = 0.5, color = NA) +
  geom_errorbar(aes(ymin = mean_det - se, ymax = mean_det + se),
                width = 0.12, linewidth = 0.9, color = "#333333") +
  scale_fill_manual(values = c(R = col_red, W = col_white), guide = "none") +
  scale_x_discrete(labels = c(R = "Red", W = "White")) +
  labs(
    title    = "Effect of Light Color on Bat Detections",
    subtitle = "Mean ± SE  |  White light p < 0.001",
    x        = "Light color",
    y        = "Mean detections per night"
  ) +
  bat_theme

# ==============================================================================
# ── Fig 2: Distance from light source ─────────────────────────────────────────
# ==============================================================================

dist_summary <- data_env %>%
  group_by(site, dist_km) %>%
  summarise(
    mean_det = mean(detections, na.rm = TRUE),
    se       = sd(detections, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  arrange(dist_km) %>%
  mutate(site = fct_reorder(site, dist_km))

fig2 <- ggplot(dist_summary,
               aes(x = site, y = mean_det)) +
  geom_col(fill = col_ivory, color = col_navy,
           linewidth = 0.6, width = 0.7) +
  geom_errorbar(aes(ymin = mean_det - se, ymax = mean_det + se),
                width = 0.2, linewidth = 0.8, color = col_navy) +
  geom_text(aes(label = round(dist_km, 2), y = -0.4),
            size = 2.8, color = "#666666", vjust = 1) +
  labs(
    title    = "Bat Detections by Distance from Light Source",
    subtitle = "Site-level means ± SE, ordered by increasing distance (km shown below bars)",
    x        = "Site",
    y        = "Mean detections per night"
  ) +
  bat_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ==============================================================================
# ── Fig 3: Color × Intensity — marginal interaction (p = 0.09) ───────────────
# Mean ± SE at each discrete intensity level by color
# ==============================================================================

intensity_summary <- data_env %>%
  group_by(intensity, color) %>%
  summarise(
    mean_det = mean(detections, na.rm = TRUE),
    se       = sd(detections,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(intensity = factor(intensity))

fig3 <- ggplot(intensity_summary,
               aes(x = intensity, y = mean_det,
                   color = color, group = color)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = mean_det - se, ymax = mean_det + se),
                width = 0.15, linewidth = 0.8,
                position = position_dodge(0.1)) +
  scale_color_manual(values = c(R = col_red, W = col_white),
                     labels = c(R = "Red", W = "White"),
                     name   = "Light color") +
  labs(
    title    = "Light Intensity Effect by Color",
    subtitle = "Mean ± SE at each intensity level",
    x        = "Light intensity (%)",
    y        = "Mean detections per night"
  ) +
  bat_theme +
  theme(legend.position = c(0.88, 0.88))

# ==============================================================================
# ── Fig 4: Seasonal patterns by light color ───────────────────────────────────
# ==============================================================================

fig4 <- ggplot(data_env,
               aes(x = jd, y = detections, color = color)) +
  geom_smooth(method = "loess", span = 0.4, se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = c(R = col_red, W = col_white),
                     labels = c(R = "Red", W = "White"),
                     name   = "Light color") +
  labs(
    title    = "Seasonal Bat Activity by Light Color",
    subtitle = "LOESS smooth across Julian day",
    x        = "Julian day",
    y        = "Detections per night"
  ) +
  bat_theme +
  theme(legend.position = c(0.12, 0.88))

# ==============================================================================
# ── Fig 5: Species community heatmap by treatment ─────────────────────────────
# ==============================================================================

heat_data <- data_env %>%
  filter(species %in% focal_spp) %>%
  group_by(species, color) %>%
  summarise(
    mean_det = mean(detections, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(species = fct_reorder(species, mean_det, .fun = sum))

fig5 <- ggplot(heat_data,
               aes(x = color, y = species, fill = mean_det)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient(low = col_ivory, high = col_navy,
                      name = "Mean\ndetections") +
  scale_x_discrete(labels = c(R = "Red", W = "White")) +
  labs(
    title    = "Bat Community Response by Light Color",
    subtitle = "Mean detections per species",
    x        = "Light color",
    y        = NULL
  ) +
  bat_theme +
  theme(
    axis.line  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

# ==============================================================================
# ── Fig 6: Moon phase effect by color (mean_phase p = 0.008) ─────────────────
# ==============================================================================

fig6 <- ggplot(data_env %>% filter(!is.na(mean_phase)),
               aes(x = mean_phase, y = detections, color = color)) +
  geom_smooth(method = "loess", span = 0.8, se = FALSE, linewidth = 1.2) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("New", "1st Q", "Full", "3rd Q", "New")) +
  scale_color_manual(values = c(R = col_red, W = col_white),
                     labels = c(R = "Red", W = "White"),
                     name   = "Light color") +
  labs(
    title    = "Effect of Moon Phase on Bat Detections",
    subtitle = "LOESS smooth by light color (p = 0.008)",
    x        = "Moon phase",
    y        = "Detections per night"
  ) +
  bat_theme +
  theme(legend.position = c(0.12, 0.88))

# ==============================================================================
# ── Fig 7: Detections vs continuous distance (site-level means) ───────────────
# ==============================================================================

library(ggrepel)

site_means <- data_env %>%
  group_by(site, dist_km) %>%
  summarise(mean_det = mean(detections, na.rm = TRUE), .groups = "drop")

fig7 <- ggplot(site_means, aes(x = dist_km, y = mean_det)) +
  geom_smooth(method = "lm", se = TRUE,
              color = col_navy, fill = col_ivory, linewidth = 1) +
  geom_point(size = 3.5, color = col_red) +
  geom_text_repel(aes(label = site), size = 3, color = "#555555",
                  box.padding = 0.4) +
  labs(
    title    = "Bat Detections vs Distance from Light Source",
    subtitle = "Site-level means with linear trend (Spearman correlation)",
    x        = "Distance from light (km)",
    y        = "Mean detections per night"
  ) +
  bat_theme

# Print correlation in console for reference
message("Distance-detection Spearman correlation:")
print(cor.test(site_means$dist_km, site_means$mean_det, method = "spearman"))

# ==============================================================================
# ── Fig 8: Species × Distance — community heatmap ────────────────────────────
# ==============================================================================

spp_dist <- data_env %>%
  filter(species %in% focal_spp) %>%
  mutate(dist_cat = factor(dist_cat, levels = dist_order)) %>%
  group_by(species, dist_cat) %>%
  summarise(mean_det = mean(detections, na.rm = TRUE), .groups = "drop") %>%
  mutate(species = fct_reorder(species, mean_det, .fun = sum))

fig8 <- ggplot(spp_dist,
               aes(x = dist_cat, y = species, fill = mean_det)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient(low = col_ivory, high = col_navy,
                      name = "Mean\ndetections") +
  labs(
    title    = "Species Detections by Distance from Light Source",
    subtitle = "Mean detections per species × distance category",
    x        = "Distance category",
    y        = NULL
  ) +
  bat_theme +
  theme(
    axis.line  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

# ==============================================================================
# ── Fig 9: Species × Site heatmap (z-scored within species) ──────────────────
# ==============================================================================

# Site order by increasing distance
site_order <- data_env %>%
  distinct(site, dist_km) %>%
  arrange(dist_km) %>%
  pull(site)

spp_site <- data_env %>%
  filter(species %in% focal_spp) %>%
  group_by(species, site) %>%
  summarise(mean_det = mean(detections, na.rm = TRUE), .groups = "drop") %>%
  group_by(species) %>%
  mutate(mean_det_z = as.numeric(scale(mean_det))) %>%
  ungroup() %>%
  mutate(
    site    = factor(site,    levels = site_order),
    species = fct_reorder(species, mean_det, .fun = sum)
  )

dist_labels <- data_env %>%
  distinct(site, dist_km) %>%
  arrange(dist_km) %>%
  mutate(site  = factor(site, levels = site_order),
         label = paste0(round(dist_km, 2), " km"))

fig9 <- ggplot(spp_site, aes(x = site, y = species, fill = mean_det_z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low      = col_red,
    mid      = col_ivory,
    high     = col_navy,
    midpoint = 0,
    name     = "Relative\nactivity (z)"
  ) +
  scale_x_discrete(labels = setNames(dist_labels$label, dist_labels$site)) +
  labs(
    title    = "Species Activity by Site",
    subtitle = "Z-scored within species · sites ordered by increasing distance from light source",
    x        = "Distance from light source",
    y        = NULL
  ) +
  bat_theme +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid   = element_blank(),
    axis.line    = element_blank(),
    axis.ticks   = element_blank(),
    legend.position = "right"
  )

# ==============================================================================
# ── Fig 10: Mean detections by intensity level ────────────────────────────────
# ==============================================================================

intensity_bar <- data_env %>%
  group_by(intensity) %>%
  summarise(
    mean_det = mean(detections, na.rm = TRUE),
    se       = sd(detections, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(intensity = factor(intensity, levels = c("10","30","50","70","100")))

fig10 <- ggplot(intensity_bar, aes(x = intensity, y = mean_det,
                                    fill = intensity)) +
  geom_col(width = 0.6, color = NA) +
  geom_errorbar(aes(ymin = mean_det - se, ymax = mean_det + se),
                width = 0.18, linewidth = 0.8, color = "#333333") +
  scale_fill_manual(values = intensity_pal, guide = "none") +
  labs(
    title    = "Bat Detections by Light Intensity",
    subtitle = "Mean ± SE across all sites and colours",
    x        = "Light intensity (%)",
    y        = "Mean detections per night"
  ) +
  bat_theme

# ==============================================================================
# ── Fig 11: Mean detections by intensity × color (grouped bar) ───────────────
# ==============================================================================

int_color_bar <- data_env %>%
  mutate(intensity = factor(intensity, levels = c("10","30","50","70","100"))) %>%
  group_by(intensity, color, site, dist_km) %>%
  summarise(
    mean_det = mean(detections, na.rm = TRUE),
    se       = sd(detections, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(site_label = paste0(site, "\n(", round(dist_km, 2), " km)"),
         site_label = fct_reorder(site_label, dist_km))

fig11 <- ggplot(int_color_bar,
                aes(x = intensity, y = mean_det, fill = color)) +
  geom_col(position = position_dodge(0.7), width = 0.6, color = NA) +
  geom_errorbar(aes(ymin = mean_det - se, ymax = mean_det + se),
                position = position_dodge(0.7),
                width = 0.18, linewidth = 0.7, color = "#333333") +
  scale_fill_manual(values = c(R = col_red, W = col_white),
                    labels  = c(R = "Red", W = "White"),
                    name    = "Light color") +
  facet_wrap(~ site_label, nrow = 2) +
  labs(
    title    = "Bat Detections by Intensity and Light Color at Each Site",
    subtitle = "Mean ± SE · sites ordered by increasing distance from light source",
    x        = "Light intensity (%)",
    y        = "Mean detections per night"
  ) +
  bat_theme +
  theme(legend.position = "top",
        strip.text      = element_text(size = 8))

# ==============================================================================
# ── Fig 12: Detections by color × intensity, faceted by species ───────────────
# ==============================================================================

spp_int_color <- data_env %>%
  filter(species %in% focal_spp) %>%
  mutate(intensity = factor(intensity, levels = c("10","30","50","70","100"))) %>%
  group_by(species, intensity, color) %>%
  summarise(
    mean_det = mean(detections, na.rm = TRUE),
    se       = sd(detections, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

fig12 <- ggplot(spp_int_color,
                aes(x = intensity, y = mean_det, fill = color)) +
  geom_col(position = position_dodge(0.7), width = 0.6, color = NA) +
  geom_errorbar(aes(ymin = mean_det - se, ymax = mean_det + se),
                position = position_dodge(0.7),
                width = 0.18, linewidth = 0.6, color = "#333333") +
  scale_fill_manual(values = c(R = col_red, W = col_white),
                    labels  = c(R = "Red", W = "White"),
                    name    = "Light color") +
  facet_wrap(~ species, nrow = 2, scales = "free_y") +
  labs(
    title    = "Bat Detections by Light Color and Intensity per Species",
    subtitle = "Mean ± SE · y-axes free within species",
    x        = "Light intensity (%)",
    y        = "Mean detections per night"
  ) +
  bat_theme +
  theme(legend.position = "top",
        strip.text      = element_text(size = 9, face = "bold"))

# ==============================================================================
# ── Fig 13: Detections vs pct_nonforest (site-level means) ───────────────────
# ==============================================================================

nonforest_means <- data_env %>%
  group_by(site, pct_nonforest) %>%
  summarise(mean_det = mean(detections, na.rm = TRUE), .groups = "drop")

fig13 <- ggplot(nonforest_means,
                aes(x = pct_nonforest, y = mean_det)) +
  geom_smooth(method = "lm", se = TRUE,
              color = col_navy, fill = col_ivory, linewidth = 1) +
  geom_point(size = 3.5, color = col_red) +
  geom_text_repel(aes(label = site), size = 3, color = "#555555",
                  box.padding = 0.4) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title    = "Bat Detections vs. Non-forest Cover",
    subtitle = "Site-level means with linear trend",
    x        = "Non-forest cover within 50 m buffer (%)",
    y        = "Mean detections per night"
  ) +
  bat_theme

# ==============================================================================
# ── Fig 14: Predicted detections ~ distance by intensity (ggpredict) ──────────
# ==============================================================================

library(emmeans)

dist_grid <- seq(0, 2, by = 0.05)

pred_dist <- emmeans(
  simple_model1,
  ~ intensity | dist_km,
  at   = list(dist_km = dist_grid),
  type = "response"
) %>%
  as.data.frame() %>%
  rename(conf.low  = asymp.LCL,
         conf.high = asymp.UCL) %>%
  mutate(intensity = factor(intensity, levels = c("10","30","50","70","100")))

fig14 <- ggplot(pred_dist, aes(x = dist_km, y = response,
                                color = intensity, fill = intensity)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = intensity_pal, name = "Intensity (%)") +
  scale_fill_manual(values  = intensity_pal, guide = "none") +
  scale_x_continuous(breaks = c(0, 0.3, 0.7, 1.0, 1.5, 2.0),
                     name   = "Distance from light (km)") +
  geom_vline(xintercept = c(0.3, 0.7, 2.0),
             linetype = "dashed", color = "#BBBBBB", linewidth = 0.4) +
  labs(
    title    = "Bat Activity vs. Distance from Light by Intensity",
    subtitle = "Marginal predictions ± 95% CI · dashed lines = distances tested in contrasts",
    x        = "Distance from light (km)",
    y        = "Predicted detections per night"
  ) +
  bat_theme +
  theme(legend.position = "right")


# ==============================================================================
# ── Fig 15: Species detections by color across the distance gradient ──────────
# ==============================================================================

fig15_data <- data_env %>%
  filter(species %in% focal_spp) %>%
  group_by(species, color, site, dist_km) %>%
  summarise(mean_det = mean(detections, na.rm = TRUE),
            se       = sd(detections, na.rm = TRUE) / sqrt(n()),
            .groups  = "drop")

fig15 <- ggplot(fig15_data,
                aes(x = dist_km, y = mean_det,
                    color = color, fill = color, group = color)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.0, alpha = 0.15) +
  scale_color_manual(values = c(R = col_red, W = col_white),
                     labels  = c(R = "Red", W = "White"),
                     name    = "Light color") +
  scale_fill_manual(values  = c(R = col_red, W = col_white),
                    guide   = "none") +
  facet_wrap(~ species, scales = "free_y", nrow = 2) +
  labs(
    title    = "Bat Detections by Light Color Across the Distance Gradient",
    subtitle = "Linear trend ± 95% CI · one observation per site",
    x        = "Distance from light source (km)",
    y        = "Mean detections per night"
  ) +
  bat_theme +
  theme(
    legend.position = "top",
    strip.text      = element_text(size = 9, face = "bold")
  )

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

ggsave("output/figures/fig1_color_effect.png",       fig1, width = 5,  height = 5,   dpi = 300)
ggsave("output/figures/fig2_distance_effect.png",    fig2, width = 6,  height = 5,   dpi = 300)
ggsave("output/figures/fig3_color_intensity.png",    fig3, width = 7,  height = 5,   dpi = 300)
ggsave("output/figures/fig4_seasonal.png",           fig4, width = 9,  height = 5,   dpi = 300)
ggsave("output/figures/fig5_community_heatmap.png",  fig5, width = 5,  height = 5.5, dpi = 300)
ggsave("output/figures/fig6_moon_phase.png",         fig6, width = 7,  height = 5,   dpi = 300)
ggsave("output/figures/fig7_distance_gradient.png",  fig7, width = 7,  height = 5.5, dpi = 300)
ggsave("output/figures/fig8_species_distance.png",   fig8, width = 9,  height = 6,   dpi = 300)
ggsave("output/figures/fig9_species_site.png",       fig9,  width = 11, height = 6,   dpi = 300)
ggsave("output/figures/fig10_intensity_bar.png",     fig10, width = 6,  height = 5,   dpi = 300)
ggsave("output/figures/fig11_intensity_color_bar.png", fig11, width = 16, height = 8,  dpi = 300)
ggsave("output/figures/fig12_species_intensity_color.png", fig12, width = 14, height = 7, dpi = 300)
ggsave("output/figures/fig13_nonforest_detections.png",   fig13, width = 7,  height = 5,  dpi = 300)
ggsave("output/figures/fig14_dist_intensity_pred.png",    fig14, width = 10, height = 6,  dpi = 300, bg = "white")
# ==============================================================================
# ── Fig 16: Species detections by intensity across the distance gradient ───────
# ==============================================================================

fig16_data <- data_env %>%
  filter(species %in% focal_spp) %>%
  mutate(intensity = factor(intensity, levels = c("10","30","50","70","100"))) %>%
  group_by(species, intensity, site, dist_km) %>%
  summarise(mean_det = mean(detections, na.rm = TRUE), .groups = "drop")

fig16 <- ggplot(fig16_data,
                aes(x = dist_km, y = mean_det,
                    color = intensity, fill = intensity, group = intensity)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.0, alpha = 0.12) +
  scale_color_manual(values = intensity_pal, name = "Intensity (%)") +
  scale_fill_manual(values  = intensity_pal, guide = "none") +
  facet_wrap(~ species, scales = "free_y", nrow = 2) +
  labs(
    title    = "Bat Detections by Light Intensity Across the Distance Gradient",
    subtitle = "Linear trend ± 95% CI · one observation per site",
    x        = "Distance from light source (km)",
    y        = "Mean detections per night"
  ) +
  bat_theme +
  theme(
    legend.position = "top",
    strip.text      = element_text(size = 9, face = "bold")
  )

ggsave("output/figures/fig15_species_color_by_site.png",  fig15, width = 14, height = 8,  dpi = 300, bg = "white")
ggsave("output/figures/fig16_species_intensity_dist.png", fig16, width = 14, height = 8,  dpi = 300, bg = "white")

# ==============================================================================
# ── Fig 17: Total detections by species ───────────────────────────────────────
# ==============================================================================

fig17_data <- data_env %>%
  filter(species %in% focal_spp) %>%
  group_by(species) %>%
  summarise(
    mean_det = mean(detections, na.rm = TRUE),
    se       = sd(detections, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(species = fct_reorder(species, mean_det))

fig17 <- ggplot(fig17_data,
                aes(x = mean_det, y = species, fill = species)) +
  geom_col(width = 0.65, color = NA) +
  geom_errorbar(aes(xmin = mean_det - se, xmax = mean_det + se),
                height = 0.25, linewidth = 0.7, color = "#333333") +
  scale_fill_manual(values = colorRampPalette(c(col_ivory, col_navy))(length(focal_spp)),
                    guide  = "none") +
  labs(
    title    = "Mean Bat Detections by Species",
    subtitle = "Mean ± SE across all sites, nights, and treatments",
    x        = "Mean detections per night",
    y        = NULL
  ) +
  bat_theme

ggsave("output/figures/fig17_species_detections.png",
       fig17, width = 7, height = 5, dpi = 300, bg = "white")

# ==============================================================================
# ── Fig 18: Species bar charts — light color × detections ─────────────────────
# ==============================================================================

fig18_data <- data_env %>%
  filter(species %in% focal_spp) %>%
  group_by(species, color) %>%
  summarise(
    mean_det = mean(detections, na.rm = TRUE),
    se       = sd(detections,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

fig18 <- ggplot(fig18_data,
                aes(x = color, y = mean_det, fill = color)) +
  geom_col(width = 0.6, color = NA) +
  geom_errorbar(aes(ymin = mean_det - se, ymax = mean_det + se),
                width = 0.18, linewidth = 0.7, color = "#333333") +
  scale_fill_manual(values = c(R = col_red, W = col_white),
                    labels  = c(R = "Red", W = "White"),
                    name    = "Light color") +
  scale_x_discrete(labels = c(R = "Red", W = "White")) +
  facet_wrap(~ species, nrow = 2, scales = "free_y") +
  labs(
    title    = "Bat Detections by Light Color per Species",
    subtitle = "Mean ± SE · y-axes free within species",
    x        = "Light color",
    y        = "Mean detections per night"
  ) +
  bat_theme +
  theme(legend.position = "none",
        strip.text      = element_text(size = 9, face = "bold"))

ggsave("output/figures/fig18_species_color_bar.png",
       fig18, width = 12, height = 7, dpi = 300, bg = "white")

# ==============================================================================
# ── Fig 19: Species bar charts — intensity × detections ───────────────────────
# ==============================================================================

fig19_data <- data_env %>%
  filter(species %in% focal_spp) %>%
  mutate(intensity = factor(intensity, levels = c("10","30","50","70","100"))) %>%
  group_by(species, intensity) %>%
  summarise(
    mean_det = mean(detections, na.rm = TRUE),
    se       = sd(detections,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

fig19 <- ggplot(fig19_data,
                aes(x = intensity, y = mean_det, fill = intensity)) +
  geom_col(width = 0.65, color = NA) +
  geom_errorbar(aes(ymin = mean_det - se, ymax = mean_det + se),
                width = 0.18, linewidth = 0.7, color = "#333333") +
  scale_fill_manual(values = intensity_pal, guide = "none") +
  facet_wrap(~ species, nrow = 2, scales = "free_y") +
  labs(
    title    = "Bat Detections by Light Intensity per Species",
    subtitle = "Mean ± SE · y-axes free within species",
    x        = "Light intensity (%)",
    y        = "Mean detections per night"
  ) +
  bat_theme +
  theme(legend.position = "none",
        strip.text      = element_text(size = 9, face = "bold"))

ggsave("output/figures/fig19_species_intensity_bar.png",
       fig19, width = 12, height = 7, dpi = 300, bg = "white")

message("All figures saved to output/figures/")
