# 17_bat_moth_effects.R
# Bat vs. moth comparison figures on shared, comparable scales.
# Raw totals differ hugely between taxa in magnitude and effort, so these
# figures instead show *relative* change (% change from each taxon's own
# baseline) for color, intensity, and distance, plus a combined rate-ratio
# forest plot with proper CIs for a rigorous side-by-side effect comparison.
#
# Input:  bat_model, moth_model — glmmTMB objects with matching
#         color * intensity + intensity * dist_km structure
#         data_env, insect_total_env — for each taxon's observed dist_km range
# Output: list(fig_color, fig_intensity, fig_distance, fig_forest, forest_df)

compare_bat_moth_effects <- function(bat_model, moth_model, data_env, insect_total_env,
                                      bat_theme, col_navy, col_purple) {
  library(emmeans)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(stringr)
  library(purrr)

  taxon_pal        <- c(Bat = col_navy, Moth = col_purple)
  intensity_levels <- c("10", "30", "50", "70", "100")

  # ── A. Color effect (% diff White vs Red) across intensity levels ────────
  # Use revpairwise (White / Red) to match the forest plot's contrast direction;
  # plain pairs() defaults to Red / White since factor levels are alphabetical.
  color_pct <- function(model, taxon) {
    em <- emmeans(model, ~ color | intensity, type = "response")
    as.data.frame(contrast(em, method = "revpairwise", adjust = "none")) %>%
      transmute(
        intensity = factor(intensity, levels = intensity_levels),
        pct_diff  = (ratio - 1) * 100,
        taxon     = taxon
      )
  }
  fig_color_df <- bind_rows(color_pct(bat_model, "Bat"), color_pct(moth_model, "Moth"))

  fig_color <- ggplot(fig_color_df, aes(x = intensity, y = pct_diff,
                                        color = taxon, group = taxon)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#999999") +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = taxon_pal, name = NULL) +
    labs(
      title    = "Color Effect Relative to Red — Bat vs. Moth",
      subtitle = "% change in detections under white vs. red light, by intensity",
      x        = "Light intensity (%)",
      y        = "% change (White vs. Red)"
    ) +
    bat_theme +
    theme(legend.position = "top")

  # ── B. Intensity effect (% change from 10% baseline), marginal over color ─
  intensity_pct <- function(model, taxon) {
    em <- as.data.frame(emmeans(model, ~ intensity, type = "response")) %>%
      mutate(intensity = factor(intensity, levels = intensity_levels)) %>%
      arrange(intensity)
    base <- em$response[em$intensity == "10"]
    em %>% transmute(intensity = intensity,
                      pct_change = (response / base - 1) * 100,
                      taxon = taxon)
  }
  fig_intensity_df <- bind_rows(intensity_pct(bat_model, "Bat"), intensity_pct(moth_model, "Moth"))

  fig_intensity <- ggplot(fig_intensity_df, aes(x = intensity, y = pct_change,
                                                color = taxon, group = taxon)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#999999") +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = taxon_pal, name = NULL) +
    labs(
      title    = "Intensity Effect Relative to 10% Baseline — Bat vs. Moth",
      subtitle = "% change in detections vs. the dimmest (10%) treatment, marginal over color",
      x        = "Light intensity (%)",
      y        = "% change from 10% baseline"
    ) +
    bat_theme +
    theme(legend.position = "top")

  # ── C. Distance gradient (% change from the nearest shared distance) ─────
  bat_range    <- range(data_env$dist_km, na.rm = TRUE)
  moth_range   <- range(insect_total_env$dist_km, na.rm = TRUE)
  shared_range <- c(max(bat_range[1], moth_range[1]), min(bat_range[2], moth_range[2]))
  dist_grid    <- seq(shared_range[1], shared_range[2], length.out = 60)

  dist_pct <- function(model, taxon) {
    em   <- as.data.frame(emmeans(model, ~ dist_km, at = list(dist_km = dist_grid), type = "response"))
    base <- em$response[1]
    em %>% transmute(dist_km = dist_km,
                      pct_change = (response / base - 1) * 100,
                      taxon = taxon)
  }
  fig_distance_df <- bind_rows(dist_pct(bat_model, "Bat"), dist_pct(moth_model, "Moth"))

  fig_distance <- ggplot(fig_distance_df, aes(x = dist_km, y = pct_change,
                                              color = taxon, group = taxon)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#999999") +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = taxon_pal, name = NULL) +
    labs(
      title    = "Distance Gradient Relative to the Nearest Shared Site — Bat vs. Moth",
      subtitle = paste0("% change in detections vs. ", round(shared_range[1], 2),
                        " km · marginal over color and intensity"),
      x        = "Distance from light source (km)",
      y        = paste0("% change from ", round(shared_range[1], 2), " km")
    ) +
    bat_theme +
    theme(legend.position = "top")

  # ── D. Combined effect-size forest plot (rate ratios, log scale) ─────────
  ratio_row <- function(term, estimate, se, taxon) {
    tibble(term = term, ratio = exp(estimate),
           lower_CL = exp(estimate - 1.96 * se),
           upper_CL = exp(estimate + 1.96 * se),
           taxon = taxon)
  }

  effect_rows <- function(model, taxon) {
    color_c   <- as.data.frame(contrast(emmeans(model, ~ color), method = "revpairwise", adjust = "none"))
    color_row <- ratio_row("Color (White vs Red)", color_c$estimate, color_c$SE, taxon)

    int_c <- as.data.frame(contrast(emmeans(model, ~ intensity), method = "trt.vs.ctrl1", adjust = "none"))
    int_rows <- pmap_dfr(int_c, function(contrast, estimate, SE, ...) {
      lvl <- str_extract(contrast, "\\d+")
      ratio_row(paste0("Intensity ", lvl, "% vs 10%"), estimate, SE, taxon)
    })

    dist_c   <- as.data.frame(emtrends(model, ~ 1, var = "dist_km"))
    dist_row <- ratio_row("Distance (per km)", dist_c$dist_km.trend, dist_c$SE, taxon)

    bind_rows(color_row, int_rows, dist_row)
  }

  forest_df <- bind_rows(effect_rows(bat_model, "Bat"), effect_rows(moth_model, "Moth")) %>%
    mutate(term = factor(term, levels = c(
      "Color (White vs Red)",
      "Intensity 30% vs 10%", "Intensity 50% vs 10%",
      "Intensity 70% vs 10%", "Intensity 100% vs 10%",
      "Distance (per km)"
    )))

  fig_forest <- ggplot(forest_df, aes(x = ratio, y = fct_rev(term), color = taxon)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#999999") +
    geom_errorbar(aes(xmin = lower_CL, xmax = upper_CL), orientation = "y",
                  width = 0.2, position = position_dodge(0.5), linewidth = 0.7) +
    geom_point(size = 3, position = position_dodge(0.5)) +
    scale_x_log10() +
    scale_color_manual(values = taxon_pal, name = NULL) +
    labs(
      title    = "Effect Sizes: Bat vs. Moth",
      subtitle = "Rate ratios ± 95% CI on a log scale · dashed line = no effect",
      x        = "Rate ratio (log scale)",
      y        = NULL
    ) +
    bat_theme +
    theme(legend.position = "top")

  dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
  ggsave("output/figures/fig_bat_vs_moth_color_relative.png",     fig_color,     width = 7, height = 5, dpi = 300, bg = "white")
  ggsave("output/figures/fig_bat_vs_moth_intensity_relative.png", fig_intensity, width = 7, height = 5, dpi = 300, bg = "white")
  ggsave("output/figures/fig_bat_vs_moth_distance_relative.png",  fig_distance,  width = 8, height = 5, dpi = 300, bg = "white")
  ggsave("output/figures/fig_bat_vs_moth_forest.png",             fig_forest,    width = 8, height = 6, dpi = 300, bg = "white")

  list(
    fig_color     = fig_color,
    fig_intensity = fig_intensity,
    fig_distance  = fig_distance,
    fig_forest    = fig_forest,
    forest_df     = forest_df
  )
}
