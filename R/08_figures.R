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

# Color palette
col_red   <- "#C0392B"
col_white <- "#2E86AB"   # blue proxy for white light (visible on plots)

# Distance category order
dist_order <- c("Close", "Medium", "Further", "Far")

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
  mutate(dist_cat = factor(dist_cat, levels = dist_order)) %>%
  group_by(dist_cat) %>%
  summarise(
    mean_det = mean(detections, na.rm = TRUE),
    se       = sd(detections, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

fig2 <- ggplot(dist_summary,
               aes(x = dist_cat, y = mean_det)) +
  geom_col(fill = "#DDEAF5", color = "#2E86AB",
           linewidth = 0.6, width = 0.55) +
  geom_errorbar(aes(ymin = mean_det - se, ymax = mean_det + se),
                width = 0.18, linewidth = 0.9, color = "#1A5276") +
  labs(
    title    = "Bat Detections by Distance from Light Source",
    subtitle = "Mean ± SE across all treatments",
    x        = "Distance category",
    y        = "Mean detections per night"
  ) +
  bat_theme

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
  group_by(species, color) %>%
  summarise(
    mean_det = mean(detections, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(species = fct_reorder(species, mean_det, .fun = sum))

fig5 <- ggplot(heat_data,
               aes(x = color, y = species, fill = mean_det)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C",
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
              color = "#2E86AB", fill = "#DDEAF5", linewidth = 1) +
  geom_point(size = 3.5, color = "#333333") +
  geom_text_repel(aes(label = site), size = 3, color = "#666666",
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
  mutate(dist_cat = factor(dist_cat, levels = dist_order)) %>%
  group_by(species, dist_cat) %>%
  summarise(mean_det = mean(detections, na.rm = TRUE), .groups = "drop") %>%
  mutate(species = fct_reorder(species, mean_det, .fun = sum))

fig8 <- ggplot(spp_dist,
               aes(x = dist_cat, y = species, fill = mean_det)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C",
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
# ── Save figures ──────────────────────────────────────────────────────────────
# ==============================================================================

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

ggsave("output/figures/fig1_color_effect.png",       fig1, width = 5,  height = 5,   dpi = 300)
ggsave("output/figures/fig2_distance_effect.png",    fig2, width = 6,  height = 5,   dpi = 300)
ggsave("output/figures/fig3_color_intensity.png",    fig3, width = 7,  height = 5,   dpi = 300)
ggsave("output/figures/fig4_seasonal.png",           fig4, width = 9,  height = 5,   dpi = 300)
ggsave("output/figures/fig5_community_heatmap.png",  fig5, width = 5,  height = 5.5, dpi = 300)
ggsave("output/figures/fig6_moon_phase.png",         fig6, width = 7,  height = 5,   dpi = 300)
ggsave("output/figures/fig7_distance_gradient.png",  fig7, width = 7,  height = 5.5, dpi = 300)
ggsave("output/figures/fig8_species_distance.png",   fig8, width = 9,  height = 6,   dpi = 300)

message("All figures saved to output/figures/")
