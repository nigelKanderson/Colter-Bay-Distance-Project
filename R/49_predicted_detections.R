# 49_predicted_detections.R
# White vs. red light and predicted detections by light color and intensity, from the
# storyboard bat and moth GLMMs (the same models as the bat vs. moth % change plate in
# 42_bat_moth_pct_panels.R):
#   Rows: white vs red % change at each intensity, bats vs moths (A-B; 47_*.R engine) |
#         predicted bats, detections per night (C-D) | predicted moths, per sampling event (E-F)
#   Cols: full data (14 sites) | far sites (GRTE13/14/15) removed (11 sites)
# Points = marginal means on the response scale (other covariates at their means, mean
# distance), 95% CIs; asterisks = white vs red differs at that intensity (* p<0.05,
# ** p<0.01, *** p<0.001; Holm-adjusted across the 5 intensities). y-axis shared within each row.
# Colors: storyboard crimson (red light) and darkened gold (white light; color-blind checked).
# Output: output/figures/fig_predicted_detections.png

make_predicted_detections <- function(bat_model, moth_model, site_night, insect_total_env, bat_theme,
                                      far_sites = c("GRTE13", "GRTE14", "GRTE15"),
                                      col_red = "#7A1C2E", col_white = "#A8842A",
                                      file = "output/figures/fig_predicted_detections.png") {
  library(emmeans)
  library(tidyverse)
  library(patchwork)

  int_levels <- c("10", "30", "50", "70", "100")
  stars  <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
  add_dc <- . ~ . + color:dist_km
  drop   <- function(df) droplevels(filter(df, !site %in% far_sites))
  fits <- list(
    full  = list(Bat  = update(bat_model,  add_dc, data = site_night),
                 Moth = update(moth_model, add_dc, data = insect_total_env)),
    nofar = list(Bat  = update(bat_model,  add_dc, data = drop(site_night)),
                 Moth = update(moth_model, add_dc, data = drop(insect_total_env)))
  )

  tabs <- map(fits, ~ imap(.x, function(m, tx) {
    em  <- emmeans(m, ~ color | intensity, type = "response")
    est <- as.data.frame(em) %>%
      transmute(color = recode(color, R = "Red", W = "White"), intensity = factor(intensity, int_levels),
                y = response, lo = asymp.LCL, hi = asymp.UCL)
    ct  <- as.data.frame(contrast(em, "revpairwise", adjust = "none")) %>%
      transmute(intensity = factor(intensity, int_levels), p_raw = p.value,
                p = p.adjust(p.value, "holm"))   # Holm across the 5 intensities
    list(est = est, ct = ct)
  }))

  row_max <- function(tx) max(unlist(map(tabs, ~ .x[[tx]]$est$hi)))
  dodge <- position_dodge(0.45)
  ylab  <- c(Bat = "Predicted bat detections per night", Moth = "Predicted moths per sampling event")
  col_labels <- c(full = "Full data (14 sites)", nofar = "Far sites removed (11 sites)")

  panel <- function(ds, tx) {
    r   <- tabs[[ds]][[tx]]
    top <- r$est %>% group_by(intensity) %>% summarise(top = max(hi), .groups = "drop") %>%
      left_join(r$ct, by = "intensity")
    ggplot(r$est, aes(intensity, y, color = color, group = color)) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.6, position = dodge) +
      geom_line(linewidth = 0.5, alpha = 0.5, position = dodge) +
      geom_point(size = 2.8, position = dodge) +
      geom_text(data = top, aes(intensity, top, label = stars(p)), inherit.aes = FALSE,
                vjust = -0.3, size = 4.2, color = "#222222") +
      scale_color_manual(values = c(Red = col_red, White = col_white), name = "Light color") +
      scale_y_continuous(limits = c(0, row_max(tx) * 1.1), expand = expansion(mult = c(0, 0.02))) +
      labs(title = paste(ifelse(tx == "Bat", "Bats", "Moths"), "·", col_labels[[ds]]),
           x = "Light intensity (%)", y = ylab[[tx]]) +
      bat_theme +
      theme(plot.title = element_text(size = 11), axis.text = element_text(size = 10),
            axis.title = element_text(size = 11))
  }

  # top row: white vs red % change, bats vs moths, from the same models
  source("R/47_bat_moth_white_vs_red_panel.R")
  wr <- white_vs_red_plate(fits, pal = c(Bat = "#0072B2", Moth = "#D55E00"), bat_theme = bat_theme,
                           title = NULL, subtitle = NULL, file = NULL, col_labels = col_labels)$panels
  wr <- imap(wr, ~ .x + labs(title = paste("White vs. red, bats vs. moths ·", col_labels[[.y]])) +
               scale_color_manual(values = c(Bat = "#0072B2", Moth = "#D55E00"), name = "Taxon") +
               theme(plot.title = element_text(size = 11), axis.text = element_text(size = 10),
                     axis.title = element_text(size = 11)))

  plate <- (wr$full | wr$nofar) /
           (panel("full", "Bat") | panel("nofar", "Bat")) /
           (panel("full", "Moth") | panel("nofar", "Moth")) +
    plot_layout(guides = "collect") +
    plot_annotation(
      tag_levels = "A",
      title = "Bat vs. moth: white vs. red light and predicted detections",
      subtitle = paste0("Same GLMMs as the % change figure. A–B: % change in detections under white vs. red at each intensity. ",
                        "C–F: predicted detections (marginal means, other covariates\nat their averages and the mean distance). ",
                        "95% CIs; y-axis shared within each row. ",
                        "Asterisks: white vs. red differs at that intensity (* p < 0.05, ** p < 0.01, *** p < 0.001; Holm-adjusted across intensities)."),
      theme = theme(plot.title = element_text(face = "bold", size = 15),
                    plot.subtitle = element_text(size = 10, color = "#555555"))) &
    theme(legend.position = "top", plot.tag = element_text(face = "bold", size = 14))

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plate, width = 12, height = 13.5, dpi = 300, bg = "white")

  list(plate = plate, tables = tabs, models = fits,
       converged = map(fits, ~ map_lgl(.x, ~ isTRUE(.x$sdr$pdHess))))
}

# Bat guilds: predicted detections per guild-night by light color and intensity, from the
# per-guild zero-inflated NB GLMMs of 45_bat_guild_pct_panels.R (`guild_fits` =
# make_bat_guild_pct_panels()$models). Expected count = conditional mean x (1 - zero-inflation
# probability) (zi is intercept-only, so this rescales the conditional mean and its CI).
#   Rows: low frequency | Myev | high frequency;  Cols: full data | far sites removed
# Output: output/figures/fig_guild_predicted_detections.png
make_guild_predicted_detections <- function(guild_fits, bat_theme,
                                            col_red = "#7A1C2E", col_white = "#A8842A",
                                            file = "output/figures/fig_guild_predicted_detections.png") {
  library(glmmTMB)
  library(emmeans)
  library(tidyverse)
  library(patchwork)

  int_levels <- c("10", "30", "50", "70", "100")
  stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
  col_labels <- c(full = "Full data (14 sites)", nofar = "Far sites removed (11 sites)")

  tabs <- map(guild_fits, ~ imap(.x, function(m, g) {
    keep_p <- if (length(fixef(m)$zi)) 1 - plogis(fixef(m)$zi[["(Intercept)"]]) else 1
    em  <- emmeans(m, ~ color | intensity, type = "response")
    est <- as.data.frame(em) %>%
      transmute(color = recode(color, R = "Red", W = "White"), intensity = factor(intensity, int_levels),
                y = response * keep_p, lo = asymp.LCL * keep_p, hi = asymp.UCL * keep_p)
    ct  <- as.data.frame(contrast(em, "revpairwise", adjust = "none")) %>%
      transmute(intensity = factor(intensity, int_levels), p_raw = p.value,
                p = p.adjust(p.value, "holm"))   # Holm across the 5 intensities
    list(est = est, ct = ct)
  }))

  guilds  <- names(guild_fits[[1]])
  row_max <- function(g) max(unlist(map(tabs, ~ .x[[g]]$est$hi)))
  dodge   <- position_dodge(0.45)

  panel <- function(ds, g) {
    r   <- tabs[[ds]][[g]]
    top <- r$est %>% group_by(intensity) %>% summarise(top = max(hi), .groups = "drop") %>%
      left_join(r$ct, by = "intensity")
    ggplot(r$est, aes(intensity, y, color = color, group = color)) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.6, position = dodge) +
      geom_line(linewidth = 0.5, alpha = 0.5, position = dodge) +
      geom_point(size = 2.8, position = dodge) +
      geom_text(data = top, aes(intensity, top, label = stars(p)), inherit.aes = FALSE,
                vjust = -0.3, size = 4.2, color = "#222222") +
      scale_color_manual(values = c(Red = col_red, White = col_white), name = "Light color") +
      scale_y_continuous(limits = c(0, row_max(g) * 1.1), expand = expansion(mult = c(0, 0.02))) +
      labs(title = paste(g, "·", col_labels[[ds]]), x = "Light intensity (%)",
           y = "Predicted detections per night") +
      bat_theme +
      theme(plot.title = element_text(size = 11), axis.text = element_text(size = 10),
            axis.title = element_text(size = 11))
  }

  rows  <- map(guilds, ~ panel("full", .x) | panel("nofar", .x))
  plate <- wrap_plots(rows, ncol = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(
      tag_levels = "A",
      title = "Bat guilds: predicted detections by light color and intensity",
      subtitle = paste0("Each guild's zero-inflated NB GLMM (as in the guild % change figure); expected detections per night ",
                        "with other covariates at their averages\nand the mean distance; 95% CIs; y-axis shared within each row. ",
                        "Low frequency = Epfu, Lano, Laci; High frequency = Mylu, Myvo, Myci, Myca, Myyu.\n",
                        "Asterisks: white vs. red differs at that intensity (* p < 0.05, ** p < 0.01, *** p < 0.001; Holm-adjusted across intensities)."),
      theme = theme(plot.title = element_text(face = "bold", size = 15),
                    plot.subtitle = element_text(size = 10, color = "#555555"))) &
    theme(legend.position = "top", plot.tag = element_text(face = "bold", size = 14))

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plate, width = 12, height = 13, dpi = 300, bg = "white")
  list(plate = plate, tables = tabs)
}
