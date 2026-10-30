# 54_nofar_distance_panels.R
# Figures for the far-sites-removed distance storyboard (distance_storyboard_nofar_distance.qmd),
# 11 sites within 0-2.1 km of the lot:
#   make_nofar_distance_figure(): all distance results on one plate
#     A  bat vs moth % change along the distance gradient   (R/42 panel, far sites removed)
#     B  moth taxa responding to distance (BH tree)          (PNG from the moth_sig_hierarchy chunk)
#     C  bat species responding to distance                  (R/53 tree, far sites removed)
#   make_nofar_int_color_figure(): the rest of the bat vs moth plate
#     A  distance x intensity   B  distance x color
# `panels` = make_bat_moth_pct_panels(...)$panels run with SB_KEEP <- "nofar" (gradient, x intensity,
# x color); `bat_hier` = fit_bat_species_hierarchy(...)["nofar"] or the saved fits.

source("R/53_bat_species_hierarchy_panel.R")

make_nofar_distance_figure <- function(panels, bat_hier,
                                       moth_png = "output/figures/moth_sig_hierarchy_distance_nofar_BH.png",
                                       file = "output/figures/nofar_distance/fig_distance_results.png") {
  library(tidyverse); library(patchwork)
  ttl <- function(p, t) p + labs(title = t) + theme(plot.title = element_text(face = "bold", size = 13))
  pA <- ttl(panels[[1]], "Bat vs. moth: % change along the distance gradient") + labs(tag = "A")
  moth <- as.raster(magick::image_read(moth_png))
  pB <- wrap_elements(full = grid::rasterGrob(moth, interpolate = TRUE)) + labs(tag = "B")
  sub <- paste0("each inset = model-predicted bat detections per night vs distance (0–max km, auto-scaled) ± 95% CI; ",
                "gray text = % change per km\nsolid colored box = p < 0.05 (all bats: single test; guilds and species: ",
                "Benjamini–Hochberg-adjusted across those 12 nodes); dashed gray box = not significant")
  hr <- bat_hier[[1]]; hr$sumry <- bat_hier_adjust(hr$sumry)
  pC <- bat_hier_plot(hr, "Bat species hierarchy (distance) — Far sites removed (11 sites)", sub) + labs(tag = "C")

  W <- 16
  hB <- W * nrow(moth) / ncol(moth)
  top <- wrap_plots(plot_spacer(), pA, plot_spacer(), widths = c(0.18, 0.64, 0.18))
  plate <- wrap_plots(top, pB, pC, ncol = 1, heights = c(5.2, hB, 5.8)) +
    plot_annotation() &
    theme(plot.tag = element_text(face = "bold", size = 18))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plate, width = W, height = 5.2 + hB + 5.8 + 0.6, dpi = 300, bg = "white", limitsize = FALSE)
  invisible(plate)
}

make_nofar_int_color_figure <- function(panels, file = "output/figures/nofar_distance/fig_bat_moth_int_color.png") {
  library(tidyverse); library(patchwork)
  plate <- wrap_plots(panels[[2]], panels[[3]], nrow = 1, guides = "collect") +
    plot_annotation(
      tag_levels = "A", title = "Bat vs. moth: distance by light intensity and color",
      subtitle = paste0("Each taxon's negative binomial GLMM (color × intensity, intensity × distance and color × distance terms); ",
                        "95% CIs.\nAsterisks: estimate differs from 0 (* p < 0.05, ** p < 0.01, *** p < 0.001); ",
                        "Holm-adjusted within each taxon and panel."),
      theme = theme(plot.title = element_text(face = "bold", size = 15),
                    plot.subtitle = element_text(size = 10, color = "#555555"))) &
    theme(legend.position = "top", plot.tag = element_text(face = "bold", size = 14))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plate, width = 15, height = 5.6, dpi = 300, bg = "white")
  invisible(plate)
}
