# =============================================================================
# 30_vs_445.R
# Compare bat species IDs between SonoBat v4.4.5 and v30.1 on the FULL
# grandteton_wide dataset (many more nights than the distance-project subset
# handled by 32_sonobat_version_compare.R).
#
# Reuses the validated comparison in R/sonobat_compare_fns.R:
#   recordings matched by timestamp key SITE_YYYYMMDD_HHMMSS; each version's
#   confident single-species accept (lone SppAccp species, first Prob >= 0.9)
#   is compared; compound / low-prob calls -> NoID.
#
# Output: console summary + output/sonobat_compare_wide/*.csv + confusion heatmap
# =============================================================================

source("R/sonobat_compare_fns.R")

wide_base <- paste0(
  "/Users/nanderson/Library/CloudStorage/",
  "GoogleDrive-nigel_anderson@brown.edu/.shortcut-targets-by-id/",
  "1QnPhsCxHagLxbcgwGcPhBMpzEksptTdA/barberlab_lightingprojects/grandteton_wide")

files_445 <- list.files(file.path(wide_base, "grte_wide_sonobatch_v4.4.5"),
                        pattern = "\\.xlsx$", full.names = TRUE, recursive = TRUE)
files_30  <- list.files(file.path(wide_base, "grte_wide_sonobatch_v30.1"),
                        pattern = "\\.xlsx$", full.names = TRUE, recursive = TRUE)

cat(sprintf("Found %d v4.4.5 files and %d v30.1 files in grandteton_wide.\n",
            length(files_445), length(files_30)))

compare_sonobat_versions(
  files_a  = files_445,
  files_b  = files_30,
  out_dir  = "output/sonobat_compare_wide",
  label_a  = "v4.4.5",
  label_b  = "v30.1",
  prob_cut = 0.9)
