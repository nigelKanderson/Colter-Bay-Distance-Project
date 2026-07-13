# 15_insect_models.R
# General (taxon-agnostic) moth activity model + core figures — Distance Project
# Mirrors the bat `simple_model1` GLMM so effect estimates are directly comparable.
# Input:  insect_env — one row per site/date/color/intensity sampling event, with
#         total detections and the same covariates as the bat pipeline
#         (dist_km, mean_phase, pct_nonforest, jd_c, site)
# Output: list(model, em_color_intensity, em_dist_trend, fig_color_intensity, fig_dist_intensity)

run_insect_general_model <- function(insect_env, bat_theme, intensity_pal) {

  library(glmmTMB)
  library(emmeans)
  library(tidyverse)

  if (!"jd_c" %in% names(insect_env)) {
    insect_env$jd_c <- insect_env$jd - mean(insect_env$jd, na.rm = TRUE)
  }

  model <- glmmTMB(
    detections ~
      color * intensity +
      intensity * dist_km +
      mean_phase +
      pct_nonforest +
      jd_c + I(jd_c^2) +
      (1 | site),
    data      = insect_env,
    ziformula = ~1,
    family    = nbinom2()
  )

  # ── Color x intensity marginal means ──────────────────────────────────────
  em_color_intensity <- emmeans(model, ~ color | intensity, type = "response")
  ci_df <- as.data.frame(em_color_intensity) %>%
    rename(lower_CL = asymp.LCL, upper_CL = asymp.UCL) %>%
    mutate(intensity = factor(intensity, levels = c("10","30","50","70","100")))

  fig_color_intensity <- ggplot(ci_df,
                                aes(x = intensity, y = response,
                                    color = color, group = color)) +
    geom_line(linewidth = 0.8, position = position_dodge(0.3)) +
    geom_point(size = 3.5, position = position_dodge(0.3)) +
    geom_errorbar(aes(ymin = lower_CL, ymax = upper_CL),
                  width = 0.15, linewidth = 0.7,
                  position = position_dodge(0.3)) +
    scale_color_manual(values = c(R = "#7A1C2E", W = "#C9A84C"),
                       labels = c(R = "Red", W = "White"),
                       name   = "Light color") +
    labs(
      title    = "Light Color × Intensity Effect on Moth Detections",
      subtitle = "Marginal means ± 95% CI from negative binomial GLMM",
      x        = "Light intensity (%)",
      y        = "Predicted detections per night"
    ) +
    bat_theme +
    theme(legend.position = "top")

  # ── Distance gradient by intensity ────────────────────────────────────────
  em_dist_trend <- emtrends(model, ~ intensity, var = "dist_km", type = "response")

  dist_grid <- seq(min(insect_env$dist_km, na.rm = TRUE),
                    max(insect_env$dist_km, na.rm = TRUE), length.out = 60)

  dist_df <- emmeans(model, ~ intensity | dist_km,
                      at = list(dist_km = dist_grid), type = "response") %>%
    as.data.frame() %>%
    rename(lower_CL = asymp.LCL, upper_CL = asymp.UCL)

  fig_dist_intensity <- ggplot(dist_df,
                               aes(x = dist_km, y = response,
                                   color = factor(intensity), fill = factor(intensity))) +
    geom_ribbon(aes(ymin = lower_CL, ymax = upper_CL), alpha = 0.12, color = NA) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = intensity_pal, name = "Intensity (%)") +
    scale_fill_manual(values  = intensity_pal, name = "Intensity (%)") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    labs(
      title    = "Predicted Moth Detections by Distance and Light Intensity",
      subtitle = "Marginal means ± 95% CI from negative binomial GLMM",
      x        = "Distance from Colter Bay parking lot (km)",
      y        = "Predicted detections per night",
      caption  = "Other covariates held at their means"
    ) +
    bat_theme +
    theme(legend.position = "right")

  dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
  ggsave("output/figures/moth_fig_color_intensity.png",
         fig_color_intensity, width = 8, height = 5, dpi = 300, bg = "white")
  ggsave("output/figures/moth_fig_dist_intensity.png",
         fig_dist_intensity, width = 8, height = 5, dpi = 300, bg = "white")

  list(
    model               = model,
    em_color_intensity  = em_color_intensity,
    em_dist_trend       = em_dist_trend,
    fig_color_intensity = fig_color_intensity,
    fig_dist_intensity  = fig_dist_intensity
  )
}
