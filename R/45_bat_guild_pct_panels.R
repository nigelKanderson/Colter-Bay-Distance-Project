# 45_bat_guild_pct_panels.R
# Bat call-frequency guilds: % change in detections, 6-panel plate for the storyboard.
#   Rows: distance | distance x intensity | distance x color
#   Cols: full data | far sites (GRTE13/14/15, ~7 km) removed
# Each guild (low frequency = Epfu/Lano/Laci, Myev, high frequency = Mylu/Myvo/Myci/Myca/Myyu)
# gets its own zero-inflated negative binomial GLMM on the zero-filled guild-night data
# (`cluster_night`), with the same structure as the bat vs. moth plate:
#   detections ~ color x intensity + intensity x distance + color x distance
#                + moon phase + non-forest + sky brightness + season + (1 | site), zi ~ 1
# Panels, CIs and asterisks as in 42_bat_moth_pct_panels.R (shared pct_panel_plate()), but on a
# log y-axis labeled in % change: guild CIs are very asymmetric (e.g. Myev up to +3000%).
# Guild colors match the species-composition figure (Okabe-Ito blue / pink + darkened gold,
# color-blind checked); linetypes differ on the distance curves.
# Output: output/figures/fig_bat_guild_pct_panels.png

source("R/42_bat_moth_pct_panels.R")

make_bat_guild_pct_panels <- function(cluster_night, bat_theme,
                                      drop_sets = default_drop_sets, col_labels = default_col_labels,
                                      title = "Bat call-frequency guilds: % change in detections",
                                      file = "output/figures/fig_bat_guild_pct_panels.png") {
  library(glmmTMB)
  library(tidyverse)

  guilds <- c("Low frequency" = "LowFreq", "Myev" = "Myev", "High frequency" = "HighFreq")
  f <- detections ~ color * intensity + intensity * dist_km + color * dist_km +
    mean_phase + pct_nonforest + brightness_dark + jd_c + I(jd_c^2) + (1 | site)

  fit_guild <- function(d, g) {
    glmmTMB(f, data = droplevels(filter(d, cluster == g)), ziformula = ~1, family = nbinom2())
  }
  sets <- map(drop_sets, ~ droplevels(filter(cluster_night, !site %in% .x)))
  fits <- map(sets, function(d) map(guilds, ~ fit_guild(d, .x)))
  conv <- map(fits, ~ map_lgl(.x, ~ isTRUE(.x$sdr$pdHess)))
  if (!all(unlist(conv))) warning("guild model(s) without a positive-definite Hessian: ",
                                  paste(names(which(!unlist(conv))), collapse = ", "))

  dist_range <- map(sets, ~ range(.x$dist_km))

  out <- pct_panel_plate(
    fits, dist_range,
    pal = c("Low frequency" = "#0072B2", "Myev" = "#CC79A7", "High frequency" = "#9A7400"),
    lty = c("Low frequency" = "solid", "Myev" = "13", "High frequency" = "42"),
    bat_theme = bat_theme, col_labels = col_labels,
    title    = title,
    subtitle = paste0("Each guild's zero-inflated negative binomial GLMM with color × intensity, intensity × distance ",
                      "and color × distance terms; 95% CIs. Rows share a log-scaled y-axis.\n",
                      "Low frequency = Epfu, Lano, Laci; High frequency = Mylu, Myvo, Myci, Myca, Myyu. ",
                      "Asterisks: estimate differs from 0 (* p < 0.05, ** p < 0.01, *** p < 0.001); C–F Holm-adjusted within each guild and panel."),
    file = file, log_y = TRUE)
  out$converged <- conv
  out
}
