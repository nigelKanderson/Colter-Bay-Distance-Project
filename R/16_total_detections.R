# 16_total_detections.R
# Total detection bar charts comparing bats and moths — by site, light color,
# and light intensity.
# Input:  data_env (bat), insect_total_env (moth), bat_theme, col_navy, col_purple,
#         col_red, col_white, intensity_pal
# Output: list(fig_site, fig_color, fig_intensity, data frames behind each)

plot_total_detections <- function(data_env, insect_total_env, bat_theme,
                                   col_navy, col_purple, col_red, col_white,
                                   intensity_pal) {
  library(dplyr)
  library(ggplot2)
  library(forcats)

  taxon_pal <- c(Bat = col_navy, Moth = col_purple)

  # ── By site ────────────────────────────────────────────────────────────────
  bat_site  <- data_env %>%
    group_by(site) %>%
    summarise(total = sum(detections, na.rm = TRUE),
              dist_km = first(dist_km), .groups = "drop") %>%
    mutate(taxon = "Bat")

  moth_site <- insect_total_env %>%
    group_by(site) %>%
    summarise(total = sum(detections, na.rm = TRUE),
              dist_km = first(dist_km), .groups = "drop") %>%
    mutate(taxon = "Moth")

  site_totals <- bind_rows(bat_site, moth_site) %>%
    mutate(
      site_label = paste0(site, "\n(", round(dist_km, 2), " km)"),
      site_label = fct_reorder(site_label, dist_km)
    )

  fig_site <- ggplot(site_totals, aes(x = site_label, y = total, fill = taxon)) +
    geom_col(width = 0.65, color = NA) +
    scale_fill_manual(values = taxon_pal, guide = "none") +
    facet_wrap(~ taxon, scales = "free_y", ncol = 1) +
    labs(
      title    = "Total Detections by Site — Bat vs. Moth",
      subtitle = "Sites ordered by increasing distance from light source",
      x        = "Site",
      y        = "Total detections"
    ) +
    bat_theme +
    theme(axis.text.x = element_text(size = 8))

  # ── By light color ─────────────────────────────────────────────────────────
  bat_color  <- data_env %>%
    group_by(color) %>%
    summarise(total = sum(detections, na.rm = TRUE), .groups = "drop") %>%
    mutate(taxon = "Bat")

  moth_color <- insect_total_env %>%
    group_by(color) %>%
    summarise(total = sum(detections, na.rm = TRUE), .groups = "drop") %>%
    mutate(taxon = "Moth")

  color_totals <- bind_rows(bat_color, moth_color) %>%
    mutate(color = factor(color, levels = c("R", "W")))

  fig_color <- ggplot(color_totals, aes(x = color, y = total, fill = color)) +
    geom_col(width = 0.55, color = NA) +
    geom_text(aes(label = total), vjust = -0.4, size = 3.5, color = "#333333") +
    scale_x_discrete(labels = c(R = "Red", W = "White")) +
    scale_fill_manual(values = c(R = col_red, W = col_white), guide = "none") +
    facet_wrap(~ taxon, scales = "free_y") +
    labs(
      title    = "Total Detections by Light Color — Bat vs. Moth",
      x        = "Light color",
      y        = "Total detections"
    ) +
    bat_theme

  # ── By light intensity ─────────────────────────────────────────────────────
  bat_intensity  <- data_env %>%
    group_by(intensity) %>%
    summarise(total = sum(detections, na.rm = TRUE), .groups = "drop") %>%
    mutate(taxon = "Bat")

  moth_intensity <- insect_total_env %>%
    group_by(intensity) %>%
    summarise(total = sum(detections, na.rm = TRUE), .groups = "drop") %>%
    mutate(taxon = "Moth")

  intensity_totals <- bind_rows(bat_intensity, moth_intensity) %>%
    mutate(intensity = factor(intensity, levels = c("10","30","50","70","100")))

  fig_intensity <- ggplot(intensity_totals, aes(x = intensity, y = total, fill = intensity)) +
    geom_col(width = 0.65, color = NA) +
    geom_text(aes(label = total), vjust = -0.4, size = 3.5, color = "#333333") +
    scale_fill_manual(values = intensity_pal, guide = "none") +
    facet_wrap(~ taxon, scales = "free_y") +
    labs(
      title    = "Total Detections by Light Intensity — Bat vs. Moth",
      x        = "Light intensity (%)",
      y        = "Total detections"
    ) +
    bat_theme

  dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
  ggsave("output/figures/fig_bat_vs_moth_site_totals.png",
         fig_site, width = 9, height = 8, dpi = 300, bg = "white")
  ggsave("output/figures/fig_bat_vs_moth_color_totals.png",
         fig_color, width = 8, height = 5, dpi = 300, bg = "white")
  ggsave("output/figures/fig_bat_vs_moth_intensity_totals.png",
         fig_intensity, width = 9, height = 5, dpi = 300, bg = "white")

  list(
    site_totals      = site_totals,
    color_totals      = color_totals,
    intensity_totals  = intensity_totals,
    fig_site          = fig_site,
    fig_color         = fig_color,
    fig_intensity     = fig_intensity
  )
}
