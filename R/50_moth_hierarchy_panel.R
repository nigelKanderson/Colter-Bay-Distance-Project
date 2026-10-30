# 50_moth_hierarchy_panel.R
# Storyboard plate of the moth significant-taxa distance hierarchy, one tree per storyboard site
# set: A = all 14 sites, B = GRTE15 excluded (13 sites), C = far sites removed (11 sites).
# Benjamini-Hochberg-corrected versions (taxa with BH-adjusted p < 0.05 across all taxa tested in
# that dataset), composed from the trees written by the `moth_sig_hierarchy` chunk of
# Distance_Project.qmd (no refitting; run that chunk first):
#   output/figures/moth_sig_hierarchy_distance_all14_BH.png  (all 14 sites)
#   output/figures/moth_sig_hierarchy_distance_BH.png        (GRTE15 excluded)
#   output/figures/moth_sig_hierarchy_distance_nofar_BH.png  (far sites removed)
# Output: output/figures/fig_moth_hierarchy_full_nofar.png

make_moth_hierarchy_panel <- function(files = c(A = "output/figures/moth_sig_hierarchy_distance_all14_BH.png",
                                                B = "output/figures/moth_sig_hierarchy_distance_BH.png",
                                                C = "output/figures/moth_sig_hierarchy_distance_nofar_BH.png"),
                                      file = "output/figures/fig_moth_hierarchy_full_nofar.png") {
  library(magick)
  stopifnot("hierarchy figures not found; run the moth_sig_hierarchy chunk first" = all(file.exists(files)))
  imgs <- lapply(files, image_read)
  w <- max(vapply(imgs, function(i) image_info(i)$width, numeric(1)))
  # keep all at their native scale (same text size); pad narrower ones with white
  tag <- function(img, lab) {
    img <- image_extent(img, paste0(w, "x", image_info(img)$height), gravity = "west", color = "white")
    img <- image_border(img, "white", "0x110")   # room above for the panel letter
    image_annotate(img, lab, size = 130, weight = 700, color = "black",
                   gravity = "northwest", location = "+30+20")
  }
  plate <- image_append(do.call(c, Map(tag, imgs, names(files))), stack = TRUE)
  plate <- image_background(plate, "white")
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  image_write(plate, file, format = "png", density = 300)
  invisible(plate)
}
