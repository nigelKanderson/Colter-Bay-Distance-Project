# 43_moth_color_intensity_distance.R
# Moth distance trend (% change per km) under red vs white light at each intensity,
# far sites (GRTE13/14/15, ~7 km) removed. Fits the storyboard moth model
# (color x intensity + intensity x distance + color x distance + covariates) plus a
# color x intensity x distance term, and reports a likelihood-ratio test of that term.
# Points = per-km % change with 95% CI; asterisks = differs from 0
# (* p<0.05, ** p<0.01, *** p<0.001; Wald z, unadjusted).
# Colors: storyboard crimson (red light) and a darkened gold (white light; the
# storyboard's #C9A84C is too light on white for thin marks).
# Output: output/figures/fig_moth_color_intensity_distance.png (+ plot, table, LRT)

make_moth_color_intensity_panel <- function(moth_model, insect_total_env, bat_theme,
                                            far_sites = c("GRTE13", "GRTE14", "GRTE15"),
                                            col_red = "#7A1C2E", col_white = "#A8842A",
                                            file = "output/figures/fig_moth_color_intensity_distance.png") {
  library(emmeans)
  library(tidyverse)

  d  <- droplevels(filter(insect_total_env, !site %in% far_sites))
  m2 <- update(moth_model, . ~ . + color:dist_km, data = d)
  m3 <- update(m2, . ~ . + color:intensity:dist_km, data = d)
  lrt <- anova(m2, m3)

  pct   <- function(x) (exp(x) - 1) * 100
  stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")

  tab <- as.data.frame(emtrends(m3, ~ color | intensity, var = "dist_km")) %>%
    transmute(color = factor(recode(color, R = "Red", W = "White"), c("Red", "White")),
              intensity = factor(intensity, c("10", "30", "50", "70", "100")),
              pct = pct(dist_km.trend),
              lo  = pct(dist_km.trend - 1.96 * SE),
              hi  = pct(dist_km.trend + 1.96 * SE),
              p   = 2 * pnorm(-abs(dist_km.trend / SE)))

  lrt_lab <- sprintf("Color × intensity × distance term: LRT χ²(%d) = %.1f, p = %.3f; ΔAIC = %s",
                     lrt$Df[2] - lrt$Df[1], lrt$Chisq[2], lrt$`Pr(>Chisq)`[2],
                     sub("-", "\u2212", sprintf("%.1f", AIC(m3) - AIC(m2))))

  p <- ggplot(tab, aes(intensity, pct, color = color, group = color)) +
    geom_hline(yintercept = 0, color = "#999999", linewidth = 0.4) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.6,
                  position = position_dodge(0.5)) +
    geom_point(size = 2.8, position = position_dodge(0.5)) +
    geom_text(aes(y = hi, label = stars(p)), position = position_dodge(0.5),
              vjust = -0.15, size = 4.2, color = "#222222", show.legend = FALSE) +
    scale_color_manual(values = c(Red = col_red, White = col_white), name = "Light color") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    labs(title = "Moth distance trend by light color and intensity",
         subtitle = paste0("Far sites removed (11 sites, 0–2.1 km); % change in moth detections per km from the lot, 95% CI.\n",
                           lrt_lab, ". Asterisks: differs from 0."),
         x = "Light intensity (%)", y = "% change per km") +
    bat_theme +
    theme(legend.position = "top", legend.justification = "left",
          axis.text = element_text(size = 10), axis.title = element_text(size = 11),
          plot.subtitle = element_text(size = 9.5, color = "#555555"))

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, p, width = 8.5, height = 5.5, dpi = 300, bg = "white")

  list(plot = p, table = tab, lrt = lrt, models = list(two_way = m2, three_way = m3))
}
