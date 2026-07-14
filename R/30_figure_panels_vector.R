# =============================================================================
# 30_figure_panels_vector.R
# Canonical panel builder: composes the 8 panels from the LIVE ggplot objects
# with patchwork, and writes true-vector PDFs (cairo_pdf) plus matching
# high-res PNGs. Because it uses the ggplot objects (not rasterized images),
# the PDFs are vector — crisp at any zoom / print-ready.
#
# Object sources:
#   - model figures  : fig3, fig7, fig10, fig14, fig16 (08_figures.R),
#                      insect_general (15), bm_effects (17)
#   - generated figs : build_panel_subfigs() -> g1..g7  (panel_subfigs.R)
#   - people cluster : p  (25_people_factor_cluster.R)
#
# If those objects are already in the environment (e.g. when run at the end of
# Distance_Project.qmd) they are reused; otherwise the needed scripts are
# sourced. Standalone, this reconstructs data + models via 29 (cache-based, no
# Google Drive needed) and takes ~1-2 min.
# =============================================================================

suppressMessages({ library(tidyverse); library(patchwork) })

OUT <- "output/figures/panels"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- ensure all sub-figure objects exist ------------------------------------
model_objs <- c("fig3", "fig7", "fig10", "fig14", "fig16",
                "insect_general", "bm_effects")
if (!all(vapply(model_objs, exists, logical(1)))) {
  message("Model figure objects missing - reconstructing via 29 ...")
  source("R/29_regenerate_model_figures.R")
}
if (!exists("p")) source("R/25_people_factor_cluster.R")
people_fig <- p                                   # cluster plot from 25

source("R/panel_subfigs.R")
sub <- build_panel_subfigs()

# ---- compose + save ----------------------------------------------------------
save_panel <- function(name, plt, w, h, title_size = 14) {
  plt <- plt + plot_annotation(tag_levels = "A") &
    theme(plot.tag   = element_text(size = 20, face = "bold", family = "Arial"),
          plot.title = element_text(size = title_size))
  ggsave(file.path(OUT, paste0(name, ".pdf")), plt, width = w, height = h,
         device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), plt, width = w, height = h,
         dpi = 200, bg = "white", device = ragg::agg_png)
  cat("wrote", name, ".pdf + .png\n")
}

save_panel("panel1_bat_marginal",             fig10 | fig7,                     15, 6)
save_panel("panel2_bat_color_intensity_dist", fig3  | fig14,                    15, 6)
save_panel("panel3_bat_species_trait",        fig16 | sub$g6,                   15, 6)
save_panel("panel4_moth_marginal",            sub$g1 | sub$g2,                  15, 6)
save_panel("panel5_moth_color_intensity_dist",
           insect_general$fig_color_intensity | insect_general$fig_dist_intensity, 15, 6)
save_panel("panel6_moth_family",              sub$g3 | sub$g4,                 15, 6)
save_panel("panel7_bat_vs_moth_pct",
           bm_effects$fig_color | bm_effects$fig_intensity | bm_effects$fig_distance, 23, 6, title_size = 11)
save_panel("panel8_human_survey",             people_fig | sub$g7,             15, 6)

cat("\nVector panels (PDF + PNG) written to", OUT, "\n")
