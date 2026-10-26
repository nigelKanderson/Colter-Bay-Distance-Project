# 47_bat_moth_white_vs_red_panel.R
# % change in detections under white vs. red light at each intensity, two panels:
#   A  full data (14 sites)
#   B  far sites removed (11 sites)
# make_bat_moth_white_vs_red_panel(): bats vs. moths, from the storyboard bat and moth GLMMs
#   (color x intensity + intensity x distance + color x distance + moon, non-forest, sky
#   brightness, season, site RE). Okabe-Ito blue / vermilion.
# make_guild_white_vs_red_panel(): the three call-frequency guilds, from the per-guild
#   zero-inflated NB GLMMs fitted in 45_bat_guild_pct_panels.R (same structure); guild
#   colors as in the guild figure; log y-axis labeled in % change (wide guild CIs).
# Contrasts at mean distance. Points = % change white vs red, 95% CI; asterisks = differs
# from 0 (* p<0.05, ** p<0.01, *** p<0.001; Wald z, Holm-adjusted across the 5 intensities within
# each taxon / guild and panel). Panels share a y-axis.
# Output: output/figures/fig_bat_moth_white_vs_red.png, fig_bat_guild_white_vs_red.png

# shared engine: `fits` = list(full = list(<series> = model, ...), nofar = list(...))
white_vs_red_plate <- function(fits, pal, bat_theme, title, subtitle, file,
                               log_y = FALSE, width = 11, height = 5.6,
                               col_labels = c(full = "Full data (14 sites)",
                                              nofar = "Far sites removed (11 sites)")) {
  library(emmeans)
  library(tidyverse)
  library(patchwork)

  int_levels <- c("10", "30", "50", "70", "100")
  stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")

  # r = rate ratio white / red; plotted as % change on a linear or log scale
  tabs <- map(fits, ~ imap_dfr(.x, function(m, sr) {
    as.data.frame(contrast(emmeans(m, ~ color | intensity), "revpairwise", adjust = "none")) %>%
      transmute(series = factor(sr, names(pal)), intensity = factor(intensity, int_levels),
                r = exp(estimate), lo = exp(estimate - 1.96 * SE), hi = exp(estimate + 1.96 * SE),
                p_raw = p.value, p = p.adjust(p.value, "holm"))   # Holm within this series' intensities
  }))
  if (!log_y) tabs <- map(tabs, ~ mutate(.x, across(c(r, lo, hi), ~ (.x - 1) * 100)))
  base <- if (log_y) 1 else 0

  lims <- range(c(base, unlist(map(tabs, ~ c(.x$lo, .x$hi)))))
  lims <- if (log_y) lims * c(1, (lims[2] / lims[1])^0.08) else lims + c(0, 0.08 * diff(lims))
  y_scale <- if (log_y) {
    scale_y_log10(breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 10, 20),
                  labels = function(r) { v <- round((r - 1) * 100); paste0(ifelse(v > 0, "+", ""), v, "%") })
  } else scale_y_continuous()

  dodge <- position_dodge(0.25 * length(pal) - 0.05 * (length(pal) - 2))
  panel <- function(key) {
    ggplot(tabs[[key]], aes(intensity, r, color = series, group = series)) +
      geom_hline(yintercept = base, color = "#999999", linewidth = 0.4) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.6, position = dodge) +
      geom_point(size = 2.8, position = dodge) +
      geom_text(aes(y = hi, label = stars(p)), position = dodge, vjust = -0.15, size = 4.2,
                color = "#222222", show.legend = FALSE) +
      scale_color_manual(values = pal, name = NULL) +
      y_scale + coord_cartesian(ylim = lims) +
      labs(title = col_labels[[key]], x = "Light intensity (%)", y = "% change, white vs red") +
      bat_theme +
      theme(plot.title = element_text(size = 11), axis.text = element_text(size = 10),
            axis.title = element_text(size = 11))
  }

  plate <- (panel(names(fits)[1]) | panel(names(fits)[2])) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A", title = title, subtitle = subtitle,
                    theme = theme(plot.title = element_text(face = "bold", size = 15),
                                  plot.subtitle = element_text(size = 10, color = "#555555"))) &
    theme(legend.position = "top", plot.tag = element_text(face = "bold", size = 14))

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    ggsave(file, plate, width = width, height = height, dpi = 300, bg = "white")
  }

  list(plate = plate, panels = map(set_names(names(fits)), panel), tables = tabs, models = fits,
       converged = map(fits, ~ map_lgl(.x, ~ isTRUE(.x$sdr$pdHess))))
}

make_bat_moth_white_vs_red_panel <- function(bat_model, moth_model, site_night, insect_total_env,
                                             bat_theme, far_sites = c("GRTE13", "GRTE14", "GRTE15"),
                                             col_bat = "#0072B2", col_moth = "#D55E00",
                                             file = "output/figures/fig_bat_moth_white_vs_red.png") {
  library(tidyverse)
  add_dc <- . ~ . + color:dist_km
  drop <- function(df) droplevels(filter(df, !site %in% far_sites))
  fits <- list(
    full  = list(Bat  = update(bat_model,  add_dc, data = site_night),
                 Moth = update(moth_model, add_dc, data = insect_total_env)),
    nofar = list(Bat  = update(bat_model,  add_dc, data = drop(site_night)),
                 Moth = update(moth_model, add_dc, data = drop(insect_total_env)))
  )
  white_vs_red_plate(
    fits, pal = c(Bat = col_bat, Moth = col_moth), bat_theme = bat_theme,
    title = "Bat vs. moth: white vs. red light at each intensity",
    subtitle = paste0("Storyboard bat and moth GLMMs, contrasts at the mean distance; 95% CIs. ",
                      "Asterisks: differs from 0 (* p < 0.05, ** p < 0.01, *** p < 0.001; Holm-adjusted across intensities within each taxon)."),
    file = file)
}

# `guild_fits`: the `models` element returned by make_bat_guild_pct_panels() (45_*.R)
make_guild_white_vs_red_panel <- function(guild_fits, bat_theme,
                                          file = "output/figures/fig_bat_guild_white_vs_red.png") {
  # file = NULL: build only (e.g. for the combined guild + species figure in 48_*.R)
  white_vs_red_plate(
    guild_fits,
    pal = c("Low frequency" = "#0072B2", "Myev" = "#CC79A7", "High frequency" = "#9A7400"),
    bat_theme = bat_theme, log_y = TRUE, width = 12, height = 5.8,
    title = "Bat guilds: white vs. red light at each intensity",
    subtitle = paste0("Each guild's zero-inflated NB GLMM (as in the guild % change figure), contrasts at the mean distance; ",
                      "95% CIs; log-scaled y-axis shared across panels.\n",
                      "Low frequency = Epfu, Lano, Laci; High frequency = Mylu, Myvo, Myci, Myca, Myyu. ",
                      "Asterisks: differs from 0 (* p < 0.05, ** p < 0.01, *** p < 0.001; Holm-adjusted across intensities within each guild)."),
    file = file)
}
