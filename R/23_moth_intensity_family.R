# 23_moth_intensity_family.R
# Diagnostic figure for the red-100% moth suppression: (A) mean moth
# detections per night across intensity, split by color, to show the decline
# is specific to red at maximum intensity rather than red generally being
# lower; (B) family composition at 100% intensity, red vs. white, to show
# the effect is broad-based across families rather than driven by one taxon.
#
# Input:  insect_total_env, specimens (from import_insects()$specimens),
#         bat_theme, col_red, col_white
# Output: list(fig_intensity, fig_family_100, data_intensity, data_family_100)

plot_moth_100red_diagnostic <- function(insect_total_env, specimens, bat_theme, col_red, col_white) {
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(patchwork)

  # ── Panel A: mean detections per night by color x intensity ──────────────
  intensity_summary <- insect_total_env %>%
    group_by(color, intensity) %>%
    summarise(
      n_nights = n(),
      mean_det = mean(detections, na.rm = TRUE),
      se       = sd(detections, na.rm = TRUE) / sqrt(n()),
      .groups  = "drop"
    ) %>%
    mutate(intensity = factor(intensity, levels = c("10", "30", "50", "70", "100")))

  fig_intensity <- ggplot(intensity_summary,
                          aes(x = intensity, y = mean_det, color = color, group = color)) +
    geom_line(linewidth = 0.9, position = position_dodge(0.2)) +
    geom_point(size = 3, position = position_dodge(0.2)) +
    geom_errorbar(aes(ymin = mean_det - se, ymax = mean_det + se),
                  width = 0.15, linewidth = 0.6, position = position_dodge(0.2)) +
    scale_color_manual(values = c(R = col_red, W = col_white),
                       labels = c(R = "Red", W = "White"), name = "Light color") +
    labs(
      title    = "Mean Moth Detections by Intensity",
      subtitle = "Mean ± SE per night · red peaks at 30%, bottoms out at 100%",
      x        = "Light intensity (%)",
      y        = "Mean detections per night"
    ) +
    bat_theme +
    theme(legend.position = "top")

  # ── Panel B: family composition at 100% intensity, red vs. white ─────────
  fam_100 <- specimens %>%
    filter(intensity == 100, !is.na(Family)) %>%
    count(color, Family, name = "n") %>%
    group_by(Family) %>%
    mutate(total = sum(n)) %>%
    ungroup() %>%
    mutate(Family = fct_reorder(Family, total))

  fig_family_100 <- ggplot(fam_100, aes(x = n, y = Family, fill = color)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    scale_fill_manual(values = c(R = col_red, W = col_white),
                      labels = c(R = "Red", W = "White"), name = "Light color") +
    labs(
      title    = "Family Composition at 100% Intensity",
      subtitle = "Specimens identified per family · red vs. white at max intensity",
      x        = "Specimens identified",
      y        = NULL
    ) +
    bat_theme +
    theme(legend.position = "top")

  fig_combined <- (fig_intensity | fig_family_100) +
    plot_annotation(
      title = "What's Happening at 100% Red: Moths",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )

  dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
  ggsave("output/figures/moth_fig_100red_diagnostic.png", fig_combined,
         width = 12, height = 6, dpi = 300, bg = "white")

  list(
    fig_intensity    = fig_intensity,
    fig_family_100   = fig_family_100,
    fig_combined     = fig_combined,
    data_intensity   = intensity_summary,
    data_family_100  = fam_100
  )
}
