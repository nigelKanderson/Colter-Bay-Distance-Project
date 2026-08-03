# =============================================================================
# 32_sonobat_version_compare.R
# Compare bat species IDs between two SonoBat classifier versions (v4.4.5 vs
# v30.1) on the same acoustic recordings -- DISTANCE-PROJECT subset.
#
# Input : data/sonobat_compare/v445/*.xlsx  and  data/sonobat_compare/v30/*.xlsx
# Output: console summary + output/sonobat_compare/*.csv + confusion heatmap
#
# The comparison logic lives in R/sonobat_compare_fns.R and is shared with
# 30_vs_445.R (which runs the same comparison on the grandteton_wide dataset).
# See that file for the method (timestamp matching + confident single-species
# accept at Prob >= PROB_CUT).
# =============================================================================

source("R/sonobat_compare_fns.R")

PROB_CUT <- 0.9

compare_sonobat_versions(
  files_a  = list.files("data/sonobat_compare/v445", "\\.xlsx$", full.names = TRUE),
  files_b  = list.files("data/sonobat_compare/v30",  "\\.xlsx$", full.names = TRUE),
  out_dir  = "output/sonobat_compare",
  label_a  = "v4.4.5",
  label_b  = "v30.1",
  prob_cut = PROB_CUT)
