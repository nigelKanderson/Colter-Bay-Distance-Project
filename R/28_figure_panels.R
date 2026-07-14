# =============================================================================
# 28_figure_panels.R
# Assemble multi-part figure panels (plates) from the project's figures.
#
# Most sub-figures already exist as saved PNGs in output/figures/. A few
# requested marginal figures don't, so they are generated here from the saved
# processed data (no model refitting):
#   - moth detections by intensity / by distance      (data/insect_data_out.rds)
#   - moth-family detections by intensity / dist / moon (import_insects + suncalc)
#   - single-trait morphology figure                  (output/posthoc_trait_dist_pred.csv)
#   - human-survey component scores by light           (output/people_pca_component_scores.csv)
#
# Panels are composed with magick (each sub-figure resized to a common height
# and tagged A/B/C), then written to output/figures/panels/.
# =============================================================================

suppressMessages({
  library(tidyverse)
  library(magick)
})

FIG <- "output/figures"
SRC <- file.path(FIG, "panel_src")     # generated sub-figures
OUT <- file.path(FIG, "panels")        # final panels
dir.create(SRC, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ggsave2 <- function(f, p, w = 7, h = 5.5)
  ggsave(file.path(SRC, f), p, width = w, height = h, dpi = 200, bg = "white",
         device = ragg::agg_png)

# ---- generate the missing sub-figures (shared builder) ----------------------
source("R/panel_subfigs.R")
sub <- build_panel_subfigs()
ggsave2("g1_moth_intensity.png", sub$g1)
ggsave2("g2_moth_distance.png",  sub$g2)
ggsave2("g3_family_intensity.png", sub$g3)
ggsave2("g4_family_distance.png",  sub$g4)
ggsave2("g5_family_moon.png",      sub$g5)
ggsave2("g6_trait_forearm.png",    sub$g6)
ggsave2("g7_people_scores.png",    sub$g7)

# ---- panel composition helpers ----------------------------------------------
H <- 1100L  # common sub-figure height (px)
tag <- function(path, letter) {
  image_read(path) |>
    image_resize(paste0("x", H)) |>
    image_border("white", "10x10") |>
    image_annotate(letter, size = 60, weight = 700, color = "black",
                   font = "Arial", location = "+18+10", boxcolor = "white")
}
panel <- function(name, specs) {
  imgs <- image_join(lapply(specs, function(s) tag(s[[1]], s[[2]])))
  out <- image_append(imgs)                       # side by side
  image_write(out, file.path(OUT, name))          # PNG montage (raster fallback)
  cat("wrote", file.path(OUT, name), "\n")
}
f  <- function(p) file.path(FIG, p)        # existing saved figure
s  <- function(p) file.path(SRC, p)        # generated sub-figure

# ---- the 8 panels ------------------------------------------------------------
panel("panel1_bat_marginal.png", list(
  list(f("fig10_intensity_bar.png"),  "A"),
  list(f("fig7_distance_gradient.png"), "B")))

panel("panel2_bat_color_intensity_dist.png", list(
  list(f("fig3_color_intensity.png"),     "A"),
  list(f("fig14_dist_intensity_pred.png"), "B")))

panel("panel3_bat_species_trait.png", list(
  list(f("fig16_species_intensity_dist.png"), "A"),
  list(s("g6_trait_forearm.png"),             "B")))

panel("panel4_moth_marginal.png", list(
  list(s("g1_moth_intensity.png"), "A"),
  list(s("g2_moth_distance.png"),  "B")))

panel("panel5_moth_color_intensity_dist.png", list(
  list(f("moth_fig_color_intensity.png"), "A"),
  list(f("moth_fig_dist_intensity.png"),  "B")))

panel("panel6_moth_family.png", list(
  list(s("g3_family_intensity.png"), "A"),
  list(s("g4_family_distance.png"),  "B")))

panel("panel7_bat_vs_moth_pct.png", list(
  list(f("fig_bat_vs_moth_color_relative.png"),    "A"),
  list(f("fig_bat_vs_moth_intensity_relative.png"), "B"),
  list(f("fig_bat_vs_moth_distance_relative.png"),  "C")))

panel("panel8_human_survey.png", list(
  list(f("people_clusters_factorspace.png"), "A"),
  list(s("g7_people_scores.png"),            "B")))

cat("\nAll panels written to", OUT, "\n")
