# =============================================================================
# 58_overall_v1_figures.R
# Figures of distance_storyboard_overall_v1.qmd, built from existing panels:
#   Figure 1  A  white-lights-on sky brightness heat map (R/57_white_light_heatmap.R)
#             B  bats along the distance gradient, red vs white (manager summary A)
#             C  moths along the distance gradient, red vs white (manager summary B)
#   Figure 2  A-D = manager summary C-F (distance slope and white vs red, by intensity)
#   Figure 3  bat species hierarchy, far sites removed, one overall line per node (R/53; saved fits)
#   Figure 4  moth significant-taxa hierarchy, far sites removed, BH-adjusted, one overall line per
#             node (copied from the `moth_sig_hierarchy` chunk of Distance_Project.qmd)
#   Figure 5  as Figure 1 but with the far stations (GRTE15 excluded, 13 stations): the map gains
#             an inset of GRTE13/14, and B-C come from make_manager_summary(drop_sites = "GRTE15")
# Figures 1-2 use the 11 sites within 2.1 km of the lot (make_manager_summary() default).
# Output: output/figures/overall_v1/fig1_overview.png ... fig5_overview_far.png
# =============================================================================

# heat map styled for a page-width figure (larger station / lot labels, short title)
overall_map <- function(heat) {
  for (i in seq_along(heat$layers)) {
    g <- class(heat$layers[[i]]$geom)[1]
    if (g %in% c("GeomLabelRepel", "GeomText") && !is.null(heat$layers[[i]]$aes_params$size)) {
      # layers are shared ggproto objects: enlarge a child copy, not the caller's plot
      ly <- rlang::env_clone(heat$layers[[i]]); class(ly) <- class(heat$layers[[i]])
      ly$aes_params$size <- ly$aes_params$size * 1.45
      heat$layers[[i]] <- ly
    }
  }
  heat + labs(title = "Sky brightness with the white lights on (100%)",
              subtitle = "Fisheye sky photos, interpolated between stations (DSU after each name)",
              caption = NULL) +
    theme(plot.title = element_text(face = "bold", size = 13.5, color = "#222222"),
          plot.subtitle = element_text(size = 10.5, color = "#666666"),
          legend.title = element_text(size = 11), legend.text = element_text(size = 10))
}

make_overall_v1_figures <- function(mgr, heat, dir = "output/figures/overall_v1") {
  library(ggplot2); library(patchwork)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tag_theme <- theme(plot.tag = element_text(face = "bold", size = 16))

  map <- overall_map(heat)

  fig1 <- (map | (mgr$panels$A / mgr$panels$B)) + plot_layout(widths = c(1.25, 1)) +
    plot_annotation(tag_levels = "A") & tag_theme
  ggsave(file.path(dir, "fig1_overview.png"), fig1, width = 15.5, height = 9.4, dpi = 300, bg = "white")

  p <- mgr$panels
  fig2 <- ((p$C | p$D) / (p$E | p$F)) + plot_annotation(tag_levels = "A") & tag_theme
  ggsave(file.path(dir, "fig2_intensity.png"), fig2, width = 15, height = 11, dpi = 300, bg = "white")

  invisible(list(fig1 = fig1, fig2 = fig2))
}

make_overall_v1_fig3_4 <- function(bat_fits = readRDS("output/bat_species_hierarchy_fits.rds"),
                                   moth_tree = "output/figures/moth_sig_hierarchy_distance_nofar_BH_overall.png",
                                   dir = "output/figures/overall_v1") {
  source("R/53_bat_species_hierarchy_panel.R")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  # single tree, so build it with bat_hier_plot() directly (make_bat_species_hierarchy() adds panel tags)
  res <- bat_fits[["nofar"]]; res$sumry <- bat_hier_adjust(res$sumry)
  sub <- paste0("each inset = model-predicted bat detections per night vs distance (0–max km, auto-scaled) ± 95% CI; ",
                "gray text = % change per km\nsolid colored box = p < 0.05 (all bats: single test; guilds and species: ",
                "Benjamini–Hochberg-adjusted across those 12 nodes); dashed gray box = not significant")
  ggsave(file.path(dir, "fig3_bat_hierarchy.png"),
         bat_hier_plot(res, paste("Bat species hierarchy (distance) —", SB_LABELS[["nofar"]]), sub, lines = "overall"),
         width = 16, height = 5.6, dpi = 300, bg = "white")
  stopifnot("moth tree not found; run the moth_sig_hierarchy chunk first" = file.exists(moth_tree))
  file.copy(moth_tree, file.path(dir, "fig4_moth_hierarchy.png"), overwrite = TRUE)
}

# mgr_full: make_manager_summary(..., drop_sites = SB_FULL_DROP); far_inset: p_far_inset (R/57)
make_overall_v1_fig5 <- function(mgr_full, heat, far_inset, dir = "output/figures/overall_v1") {
  library(ggplot2); library(patchwork); library(emmeans)
  pct <- function(x) (exp(x) - 1) * 100
  slope <- function(m) { e <- as.data.frame(test(emtrends(m, ~ color, var = "dist_km")))
                         setNames(Map(c, pct(e$dist_km.trend), e$p.value), e$color) }
  bt <- slope(mgr_full$bat); mt <- slope(mgr_full$moth)
  km <- scale_x_continuous("Distance from the parking lot (km)", breaks = 0:7,
                           expand = expansion(mult = c(0.02, 0.2)))
  ns <- function(p) if (p < 0.05) "" else " (n.s.)"

  pB <- mgr_full$panels$A + km +
    labs(title = if (max(bt$R[2], bt$W[2]) >= 0.05) "With the far stations, bats show no distance trend"
                 else "Bats are most active close to the lit lot, under red and white light",
         subtitle = sprintf("%+.0f%% per km under red, %+.0f%% under white%s, 0–7.5 km",
                            bt$R[1], bt$W[1], ns(max(bt$R[2], bt$W[2]))))
  pC <- mgr_full$panels$B + km +
    labs(title = "White light depletes moths near the lot; red light does not",
         subtitle = sprintf("Under white light, moths rise %.0f%% per km away from the lot; under red, %+.0f%%/km%s",
                            mt$W[1], mt$R[1], ns(mt$R[2])))
  map <- overall_map(heat) + inset_element(far_inset, left = 0.63, bottom = 0.01, right = 0.99, top = 0.34)

  fig5 <- (map | (pB / pC)) + plot_layout(widths = c(1.25, 1)) +
    plot_annotation(tag_levels = list(c("A", "", "B", "C"))) &
    theme(plot.tag = element_text(face = "bold", size = 16))
  ggsave(file.path(dir, "fig5_overview_far.png"), fig5, width = 15.5, height = 9.4, dpi = 300, bg = "white")
  invisible(fig5)
}
