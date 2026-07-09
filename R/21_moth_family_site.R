# 21_moth_family_site.R
# Moth counterpart to fig9 in 08_figures.R: family x site heatmap
# (z-scored within family), sites ordered by increasing distance from the
# light source.
#
# Input:  insect_family_env (site, dist_km, species [= Family], detections),
#         bat_theme, col_red, col_navy, col_ivory
# Output: fig_family_site (ggplot object; also saved to output/figures/)

plot_moth_family_site <- function(insect_family_env, bat_theme, col_red, col_navy, col_ivory) {
  library(dplyr)
  library(forcats)
  library(ggplot2)

  site_order <- insect_family_env %>%
    distinct(site, dist_km) %>%
    arrange(dist_km) %>%
    pull(site)

  fam_site <- insect_family_env %>%
    group_by(species, site) %>%
    summarise(mean_det = mean(detections, na.rm = TRUE), .groups = "drop") %>%
    group_by(species) %>%
    mutate(mean_det_z = as.numeric(scale(mean_det))) %>%
    ungroup() %>%
    mutate(
      site    = factor(site,    levels = site_order),
      species = fct_reorder(species, mean_det, .fun = sum)
    )

  dist_labels <- insect_family_env %>%
    distinct(site, dist_km) %>%
    arrange(dist_km) %>%
    mutate(site  = factor(site, levels = site_order),
           label = paste0(round(dist_km, 2), " km"))

  fig_family_site <- ggplot(fam_site, aes(x = site, y = species, fill = mean_det_z)) +
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
      title    = "Moth Family Activity by Site",
      subtitle = "Z-scored within family · sites ordered by increasing distance from light source",
      x        = "Distance from light source",
      y        = NULL
    ) +
    bat_theme +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
      panel.grid      = element_blank(),
      axis.line       = element_blank(),
      axis.ticks      = element_blank(),
      legend.position = "right"
    )

  dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
  ggsave("output/figures/moth_fig_family_site.png", fig_family_site,
         width = 11, height = 6, dpi = 300, bg = "white")

  fig_family_site
}
