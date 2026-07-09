# 22_moth_moon_phase.R
# Moth counterpart to fig6 in 08_figures.R: detections vs. moon phase by
# light color.
#
# Input:  insect_total_env (mean_phase, color, detections), bat_theme,
#         col_red, col_white
# Output: fig_moon_phase (ggplot object; also saved to output/figures/)

plot_moth_moon_phase <- function(insect_total_env, bat_theme, col_red, col_white) {
  library(dplyr)
  library(ggplot2)

  fig_moon_phase <- ggplot(insect_total_env %>% filter(!is.na(mean_phase)),
                           aes(x = mean_phase, y = detections, color = color)) +
    geom_smooth(method = "loess", span = 0.8, se = FALSE, linewidth = 1.2) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("New", "1st Q", "Full", "3rd Q", "New")) +
    scale_color_manual(values = c(R = col_red, W = col_white),
                       labels = c(R = "Red", W = "White"),
                       name   = "Light color") +
    labs(
      title    = "Effect of Moon Phase on Moth Detections",
      subtitle = "LOESS smooth by light color (moon phase p = 0.059, marginal)",
      x        = "Moon phase",
      y        = "Detections per night"
    ) +
    bat_theme +
    theme(legend.position = c(0.12, 0.88))

  dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
  ggsave("output/figures/moth_fig_moon_phase.png", fig_moon_phase,
         width = 8, height = 5, dpi = 300, bg = "white")

  fig_moon_phase
}
