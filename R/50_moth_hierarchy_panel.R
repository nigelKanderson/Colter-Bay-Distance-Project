# 50_moth_hierarchy_panel.R
# Storyboard plate of the moth significant-taxa distance hierarchy: A = full data (14 sites),
# B = far sites removed (11 sites). Defaults to the Benjamini-Hochberg-corrected versions
# (taxa with BH-adjusted p < 0.05 across all taxa tested in that dataset). Composes the two figures written by the
# `moth_sig_hierarchy` chunk of Distance_Project.qmd (no refitting), so run that chunk first.
#   output/figures/moth_sig_hierarchy_distance_BH.png        (full data)
#   output/figures/moth_sig_hierarchy_distance_nofar_BH.png  (far sites removed)
# Output: output/figures/fig_moth_hierarchy_full_nofar.png

make_moth_hierarchy_panel <- function(full = "output/figures/moth_sig_hierarchy_distance_BH.png",
                                      nofar = "output/figures/moth_sig_hierarchy_distance_nofar_BH.png",
                                      file = "output/figures/fig_moth_hierarchy_full_nofar.png") {
  library(magick)
  stopifnot("hierarchy figures not found; run the moth_sig_hierarchy chunk first" =
              file.exists(full) && file.exists(nofar))
  imgs <- c(image_read(full), image_read(nofar))
  w <- max(image_info(imgs)$width)
  # keep both at their native scale (same text size); pad the narrower one with white
  tag <- function(img, lab) {
    img <- image_extent(img, paste0(w, "x", image_info(img)$height), gravity = "west", color = "white")
    img <- image_border(img, "white", "0x110")   # room above for the panel letter
    image_annotate(img, lab, size = 130, weight = 700, color = "black",
                   gravity = "northwest", location = "+30+20")
  }
  plate <- image_append(c(tag(imgs[1], "A"), tag(imgs[2], "B")), stack = TRUE)
  plate <- image_background(plate, "white")
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  image_write(plate, file, format = "png", density = 300)
  invisible(plate)
}
