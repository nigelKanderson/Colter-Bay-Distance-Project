library(ggplot2)
library(dplyr)
library(ggeffects)
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
col_red    <- "#C0392B"
col_white  <- "#2E86AB"   # blue proxy for white light (visible on plots)
col_pts    <- "#555555"

# Distance category order
dist_order <- c("Close", "Medium", "Further", "Far")

# Helper: ggpredict with type="fixed" (conditional mean, bypasses ZI dominance),
# returned as plain data frame
gp_df <- function(model, terms, ...) {
  as.data.frame(ggpredict(model, terms = terms, type = "fixed", ...))
}

# ==============================================================================
# ── Fig 1: Color × Intensity interaction (headline result) ────────────────────
# Raw data means + SE — most reliable for factor × factor comparisons
# ==============================================================================

ci_summary <- data_env %>%
  group_by(intensity, color) %>%
  summarise(
    mean_det = mean(detections, na.rm = TRUE),
    se       = sd(detections, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(intensity = factor(intensity))

fig1 <- ggplot(ci_summary,
               aes(x = intensity, y = mean_det,
                   color = color, group = color)) +
  geom_errorbar(aes(ymin = mean_det - se, ymax = mean_det + se),
                width = 0.15, linewidth = 0.8,
                position = position_dodge(0.35)) +
  geom_point(size = 3.5, position = position_dodge(0.35)) +
  geom_line(position = position_dodge(0.35), linewidth = 0.9) +
  scale_color_manual(values = c(R = col_red, W = col_white),
                     labels = c(R = "Red", W = "White"),
                     name   = "Light color") +
  labs(
    title    = "Light Color × Intensity Effect on Bat Detections",
    subtitle = "Mean ± SE",
    x        = "Light intensity",
    y        = "Mean detections per night"
  ) +
  bat_theme +
  theme(legend.position = c(0.88, 0.88))

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
# ── Fig 3: Distance × Treatment interaction ───────────────────────────────────
# ==============================================================================

dist_treat <- data_env %>%
  mutate(dist_cat = factor(dist_cat, levels = dist_order)) %>%
  group_by(dist_cat, color, intensity) %>%
  summarise(
    mean_det = mean(detections, na.rm = TRUE),
    se       = sd(detections, na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

fig3 <- ggplot(dist_treat,
               aes(x = factor(intensity), y = mean_det,
                   color = color, group = color)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_det - se, ymax = mean_det + se),
                width = 0.15, linewidth = 0.7) +
  facet_wrap(~ dist_cat, nrow = 1) +
  scale_color_manual(values = c(R = col_red, W = col_white),
                     labels = c(R = "Red", W = "White"),
                     name   = "Light color") +
  labs(
    title    = "Detections by Distance, Color, and Intensity",
    subtitle = "Mean ± SE",
    x        = "Light intensity",
    y        = "Mean detections per night"
  ) +
  bat_theme +
  theme(legend.position = "top")

# ==============================================================================
# ── Fig 4: Seasonal patterns by light color ───────────────────────────────────
# ==============================================================================

pred_jd  <- gp_df(simple_model1, terms = c("jd [all]", "color"))

fig4 <- ggplot(pred_jd,
               aes(x = x, y = predicted, color = group)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c(R = col_red, W = col_white),
                     labels = c(R = "Red", W = "White"),
                     name   = "Light color") +
  labs(
    title    = "Seasonal Bat Activity by Light Color",
    subtitle = "Predicted detections (conditional mean) across Julian day",
    x        = "Julian day",
    y        = "Predicted detections"
  ) +
  bat_theme +
  theme(legend.position = c(0.12, 0.88))

# ==============================================================================
# ── Fig 5: Species community heatmap by treatment ─────────────────────────────
# ==============================================================================

heat_data <- data_env %>%
  group_by(species, color, intensity) %>%
  summarise(
    mean_det = mean(detections, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(
    treatment = paste(color, intensity, sep = "_"),
    species   = fct_reorder(species, mean_det, .fun = sum)
  )

fig5 <- ggplot(heat_data,
               aes(x = treatment, y = species, fill = mean_det)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C",
                      name = "Mean\ndetections") +
  labs(
    title    = "Bat Community Response to Light Treatment",
    subtitle = "Mean detections per species × treatment combination",
    x        = "Treatment (color_intensity)",
    y        = NULL
  ) +
  bat_theme +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1),
    axis.line    = element_blank(),
    axis.ticks   = element_blank(),
    panel.grid   = element_blank()
  )

# ==============================================================================
# ── Fig 6: Moon phase effect ──────────────────────────────────────────────────
# ==============================================================================

pred_moon <- gp_df(simple_model1, terms = c("mean_phase [all]", "intensity"))

fig6 <- ggplot(pred_moon,
               aes(x = x, y = predicted, color = group)) +
  geom_line(linewidth = 1.1) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25),
                     labels = c("New", "1st Q", "Full", "3rd Q", "New")) +
  scale_color_brewer(palette = "Dark2", name = "Intensity") +
  labs(
    title    = "Effect of Moon Phase on Bat Detections",
    subtitle = "Predicted detections (conditional mean) by intensity level",
    x        = "Moon phase",
    y        = "Predicted detections"
  ) +
  bat_theme +
  theme(legend.position = "right")

# ==============================================================================
# ── Save figures ──────────────────────────────────────────────────────────────
# ==============================================================================

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

ggsave("output/figures/fig1_color_intensity.png",    fig1, width = 7,  height = 5,   dpi = 300)
ggsave("output/figures/fig2_distance_effect.png",    fig2, width = 6,  height = 5,   dpi = 300)
ggsave("output/figures/fig3_distance_treatment.png", fig3, width = 10, height = 4.5, dpi = 300)
ggsave("output/figures/fig4_seasonal.png",           fig4, width = 9,  height = 5,   dpi = 300)
ggsave("output/figures/fig5_community_heatmap.png",  fig5, width = 8,  height = 5.5, dpi = 300)
ggsave("output/figures/fig6_moon_phase.png",         fig6, width = 7,  height = 5,   dpi = 300)

message("All figures saved to output/figures/")
