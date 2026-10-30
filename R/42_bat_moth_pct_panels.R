# 42_bat_moth_pct_panels.R
# Bat vs. moth % change comparison, 6-panel plate for the storyboard.
#   Rows: distance | distance x intensity | distance x color
#   Cols: all 14 sites | GRTE15 excluded (13) | far sites GRTE13/14/15 removed (11) (R/storyboard_sites.R)
# Each taxon's general model (bat `simple_model1`, moth `insect_general$model`) is
# refit with an added color x distance term so every panel, including the
# distance x color row, comes from one model per taxon and dataset.
#   distance            : % change in detections vs the nearest site, along the gradient
#   distance x intensity: % change per km, at each intensity (intensity_by_fill = TRUE puts the
#                         series on the x-axis with a purple intensity ramp; used for bat vs moth and the guilds)
#   distance x color   : % change per km, under red and under white (color_by_light = TRUE puts
#                        the series on the x-axis with red / white points; used for bat vs moth and the guilds)
# All bands/bars are 95% CIs from the model's log-scale estimate and SE; asterisks mark
# estimates that differ from 0 (* p<0.05, ** p<0.01, *** p<0.001; Wald z). Per-intensity and per-color
# slopes (C-F) are Holm-adjusted within each series and panel; the single gradient test (A-B) is not.
# Taxon colors are the Okabe-Ito blue/orange pair (color-blind safe); both taxa use
# circles, with solid vs dashed lines as the redundant cue on the distance curves.
# Output: output/figures/fig_bat_moth_pct_panels.png (+ list of plots and tables)
#
# pct_panel_plate() is the shared engine (also used by 45_bat_guild_pct_panels.R): it takes
# fitted models per series for the full and far-removed data and draws the 6-panel plate.

source("R/storyboard_sites.R")   # SB_FULL_DROP, SB_FAR_DROP, SB_LABELS
source("R/pct_arrows.R")         # direction arrows next to the 0 line

pct_panel_plate <- function(fits, dist_range, pal, lty, bat_theme, title, subtitle, file,
                            height = 12.5, log_y = FALSE, extra_panels = NULL,
                            arrow_up = "Higher farther\nfrom the lot", arrow_down = "Higher near\nthe lot",
                            col_labels = SB_LABELS, color_by_light = FALSE,
                            intensity_by_fill = FALSE) {
  library(emmeans)
  library(tidyverse)
  library(patchwork)

  shape      <- setNames(rep(16, length(pal)), names(pal))
  int_levels <- c("10", "30", "50", "70", "100")
  pct    <- function(x) (exp(x) - 1) * 100
  wald_p <- function(est, se) 2 * pnorm(-abs(est / se))
  stars  <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
  dodge  <- position_dodge(0.25 * length(pal) - 0.05 * (length(pal) - 2))

  # ── per-dataset tables ──────────────────────────────────────────────────────
  tabs <- imap(fits, function(ms, ds) {
    d0   <- dist_range[[ds]][1]
    grid <- seq(dist_range[[ds]][1], dist_range[[ds]][2], length.out = 80)

    gradient <- imap_dfr(ms, function(m, sr) {
      tr <- as.data.frame(emtrends(m, ~ 1, var = "dist_km"))
      tibble(series = sr, dist_km = grid, dd = grid - d0,
             pct = pct(tr$dist_km.trend * dd),
             lo  = pct((tr$dist_km.trend - 1.96 * tr$SE) * dd),
             hi  = pct((tr$dist_km.trend + 1.96 * tr$SE) * dd),
             p   = wald_p(tr$dist_km.trend, tr$SE)) %>%
        mutate(lo2 = pmin(lo, hi), hi2 = pmax(lo, hi), lo = lo2, hi = hi2) %>%
        select(-lo2, -hi2)
    })

    slope_by <- function(by) imap_dfr(ms, function(m, sr) {
      as.data.frame(emtrends(m, as.formula(paste("~", by)), var = "dist_km")) %>%
        transmute(series = sr, level = as.character(.data[[by]]),
                  pct = pct(dist_km.trend),
                  lo  = pct(dist_km.trend - 1.96 * SE),
                  hi  = pct(dist_km.trend + 1.96 * SE),
                  p   = wald_p(dist_km.trend, SE)) %>%
        mutate(p_raw = p, p = p.adjust(p, "holm"))   # Holm within this series' levels
    })

    list(gradient  = gradient,
         dist_int  = slope_by("intensity") %>% mutate(level = factor(level, int_levels)),
         dist_col  = slope_by("color") %>%
                       mutate(level = factor(recode(level, R = "Red", W = "White"),
                                             c("Red", "White"))))
  })
  tabs <- map(tabs, ~ map(.x, function(t) mutate(t, series = factor(series, names(pal)))))

  # y-axis: % change on a linear scale, or (log_y) the underlying ratio on a log scale
  # labeled as % change, so that e.g. -50% and +100% sit equally far from 0
  if (log_y) {
    tabs <- map(tabs, ~ map(.x, function(t) mutate(t, across(c(pct, lo, hi), ~ 1 + .x / 100))))
    base <- 1
    y_scale <- scale_y_log10(
      breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 10, 20, 40),
      labels = function(r) { v <- round((r - 1) * 100); paste0(ifelse(v > 0, "+", ""), v, "%") })
  } else {
    base <- 0
    y_scale <- scale_y_continuous()
  }

  # ── plotting helpers ────────────────────────────────────────────────────────
  series_scales <- list(
    scale_color_manual(values = pal, name = NULL),
    scale_fill_manual(values = pal, name = NULL),
    scale_shape_manual(values = shape, name = NULL),
    scale_linetype_manual(values = lty, name = NULL)
  )
  zero <- geom_hline(yintercept = base, color = "#999999", linewidth = 0.4)
  row_lims <- function(key, cols = c("lo", "hi")) {
    v <- unlist(map(tabs, ~ unlist(.x[[key]][cols])))
    r <- range(c(v, base), na.rm = TRUE)
    if (log_y) r * c(1, (r[2] / r[1])^0.08)   # headroom for the asterisks (log scale)
    else r + c(0, 0.08 * diff(r))
  }

  p_gradient <- function(ds, ttl) {
    g    <- tabs[[ds]]$gradient
    ends <- g %>% group_by(series) %>% slice_tail(n = 1) %>% ungroup()
    ggplot(g, aes(dist_km, pct, color = series, fill = series)) +
      zero +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, color = NA) +
      geom_line(aes(linetype = series), linewidth = 1) +
      geom_point(data = ends, aes(shape = series), size = 4) +
      ggrepel::geom_text_repel(data = ends, aes(label = paste0(series, " ", stars(p))),
                               size = 3.6, color = "#222222", show.legend = FALSE,
                               direction = "y", hjust = 0, nudge_x = 0.04 * diff(range(g$dist_km)),
                               xlim = c(max(g$dist_km), NA), segment.color = "#BBBBBB",
                               min.segment.length = 0, seed = 1) +
      scale_x_continuous(expand = expansion(mult = c(0.02, 0.32)),       # room for the labels,
                         breaks = function(l) { b <- pretty(c(l[1], max(g$dist_km)))  # but no ticks
                                                b[b <= max(g$dist_km) + 1e-9] }) +     # past the data
      series_scales + y_scale +
      guides(color = "none", fill = "none", linetype = "none", shape = "none") +  # direct-labeled
      coord_cartesian(ylim = row_lims("gradient"), clip = "off") +
      pct_arrows_y(arrow_up, arrow_down, zero = base) +
      labs(title = ttl, x = SB_DIST_LABEL,
           y = "% change vs nearest site") +
      bat_theme
  }
  p_points <- function(ds, key, ttl, xlab, ylab) {
    ggplot(tabs[[ds]][[key]], aes(level, pct, color = series, shape = series, group = series)) +
      zero +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.6, position = dodge) +
      geom_point(size = 4, position = dodge) +
      geom_text(aes(y = hi, label = stars(p)), position = dodge,
                vjust = -0.15, size = 4.2, color = "#222222", show.legend = FALSE) +
      series_scales + y_scale +
      coord_cartesian(ylim = row_lims(key), clip = "off") +
      pct_arrows_y(arrow_up, arrow_down, zero = base) +
      (if (key == "dist_col") red_tick_box(levels(tabs[[ds]][[key]]$level), bottom = if (log_y) 0 else -Inf)) +
      labs(title = ttl, x = xlab, y = ylab) +
      bat_theme
  }

  # distance x color with the series on the x-axis and light color as the points
  # (red = red point, white = white point with a black outline)
  p_by_light <- function(ds, ttl, ylab) {
    light_fill <- c(Red = "#B2182B", White = "white")
    light_line <- c(Red = "#B2182B", White = "black")
    ld <- position_dodge(0.45)
    ggplot(tabs[[ds]]$dist_col, aes(series, pct, group = level)) +
      zero +
      geom_errorbar(aes(ymin = lo, ymax = hi, color = level), width = 0.18, linewidth = 0.6, position = ld) +
      geom_point(aes(color = level, fill = level), shape = 21, size = 4.4, stroke = 1, position = ld) +
      geom_text(aes(y = hi, label = stars(p)), position = ld,
                vjust = -0.15, size = 4.2, color = "#222222", show.legend = FALSE) +
      scale_color_manual(values = light_line, name = "Light color") +
      scale_fill_manual(values = light_fill, name = "Light color") +
      y_scale +
      coord_cartesian(ylim = row_lims("dist_col"), clip = "off") +
      pct_arrows_y(arrow_up, arrow_down, zero = base) +
      labs(title = ttl, x = NULL, y = ylab) +
      bat_theme
  }

  # distance x intensity with the series on the x-axis and intensity as a light-to-dark purple
  # ramp (color-blind safe: lightness rises steadily with intensity); dark outlines and error bars
  # keep the palest level visible
  p_by_intensity <- function(ds, ttl, ylab) {
    idodge <- position_dodge(0.75)
    ggplot(tabs[[ds]]$dist_int, aes(series, pct, group = level)) +
      zero +
      geom_errorbar(aes(ymin = lo, ymax = hi), color = "#444444", width = 0.12, linewidth = 0.5,
                    position = idodge) +
      geom_point(aes(fill = level), shape = 21, color = "#222222", size = 4.4, stroke = 0.8,
                 position = idodge) +
      geom_text(aes(y = hi, label = stars(p)), position = idodge,
                vjust = -0.15, size = 4.2, color = "#222222", show.legend = FALSE) +
      scale_fill_manual(values = SB_INT_FILL, name = "Light intensity (%)") +
      y_scale +
      coord_cartesian(ylim = row_lims("dist_int"), clip = "off") +
      pct_arrows_y(arrow_up, arrow_down, zero = base) +
      labs(title = ttl, x = NULL, y = ylab) +
      bat_theme
  }

  # one column per dataset in `fits` (in order); rows = distance | x intensity | x color
  ks <- names(fits)
  panels <- c(
    lapply(ks, function(k) p_gradient(k, paste("Distance ·", col_labels[[k]]))),
    lapply(ks, function(k) {
      ttl <- paste("Distance × intensity ·", col_labels[[k]])
      if (intensity_by_fill) p_by_intensity(k, ttl, "% change per km")
      else p_points(k, "dist_int", ttl, "Light intensity (%)", "% change per km")
    }),
    lapply(ks, function(k) {
      ttl <- paste("Distance × color ·", col_labels[[k]])
      if (color_by_light) p_by_light(k, ttl, "% change per km")
      else p_points(k, "dist_col", ttl, "Light color", "% change per km")
    }))
  panels <- c(panels, extra_panels)   # optional extra rows appended below (tagged in order)
  panels <- map(panels, ~ .x + theme(plot.title = element_text(size = 11),
                                     axis.text  = element_text(size = 10),
                                     axis.title = element_text(size = 11)) +
                  pct_arrow_margin_y())

  plate <- wrap_plots(panels, ncol = length(ks), guides = "collect") +
    plot_annotation(
      tag_levels = "A", title = title, subtitle = subtitle,
      theme = theme(plot.title = element_text(face = "bold", size = 15),
                    plot.subtitle = element_text(size = 10, color = "#555555"))) &
    theme(legend.position = "top", plot.tag = element_text(face = "bold", size = 14))

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plate, width = max(7.3 * length(ks), 12), height = height, dpi = 300, bg = "white")

  list(plate = plate, panels = panels, tables = tabs, models = fits)
}

# `drop_sets`: named site vectors to drop, one per column; `col_labels` the matching column
# titles. Default = the storyboard site sets (R/storyboard_sites.R).
default_drop_sets  <- SB_SETS
default_col_labels <- SB_LABELS

make_bat_moth_pct_panels <- function(bat_model, moth_model, site_night, insect_total_env,
                                     bat_theme, col_bat = "#0072B2", col_moth = "#D55E00",
                                     drop_sets = default_drop_sets, col_labels = default_col_labels,
                                     title = "Bat vs. moth: % change in detections",
                                     file = "output/figures/fig_bat_moth_pct_panels.png") {
  library(tidyverse)
  drop <- function(df, sites) droplevels(filter(df, !site %in% sites))
  add_dc <- . ~ . + color:dist_km

  fits <- map(drop_sets, function(sites)
    list(Bat  = update(bat_model,  add_dc, data = drop(site_night, sites)),
         Moth = update(moth_model, add_dc, data = drop(insect_total_env, sites))))
  dist_range <- map(drop_sets, ~ range(drop(site_night, .x)$dist_km))

  pct_panel_plate(
    fits, dist_range,
    pal = c(Bat = col_bat, Moth = col_moth), lty = c(Bat = "solid", Moth = "22"),
    bat_theme = bat_theme, col_labels = col_labels, color_by_light = TRUE, intensity_by_fill = TRUE,
    title    = title,
    subtitle = paste0("Each taxon's negative binomial GLMM with color × intensity, intensity × distance ",
                      "and color × distance terms; 95% CIs. Rows share a y-axis.\n",
                      "Asterisks: estimate differs from 0 (* p < 0.05, ** p < 0.01, *** p < 0.001); per-intensity and per-color rows Holm-adjusted within each taxon and panel."),
    file = file)
}
