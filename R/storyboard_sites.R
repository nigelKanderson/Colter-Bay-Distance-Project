# storyboard_sites.R
# Site sets shared by the storyboard figure scripts (R/42-R/52). Each set is one column of
# every storyboard figure, in this order:
#   "all"   = all 14 sites
#   "full"  = every site except GRTE15 (13 sites)
#   "nofar" = far sites GRTE13/14/15 removed (11 sites, 0-2.1 km from the lot)
SB_FULL_DROP <- "GRTE15"
SB_FAR_DROP  <- c("GRTE13", "GRTE14", "GRTE15")
SB_SETS      <- list(all = character(0), full = SB_FULL_DROP, nofar = SB_FAR_DROP)   # sites dropped
SB_LABELS    <- c(all = "All sites (14 sites)", full = "GRTE15 excluded (13 sites)",
                  nofar = "Far sites removed (11 sites)")
SB_SHORT     <- c(all = "All 14 sites", full = "GRTE15 excl.", nofar = "Far removed")
# a run can keep only some sets (e.g. SB_KEEP <- "nofar" for a far-sites-removed-only storyboard)
if (exists("SB_KEEP")) { SB_SETS <- SB_SETS[SB_KEEP]; SB_LABELS <- SB_LABELS[SB_KEEP]; SB_SHORT <- SB_SHORT[SB_KEEP] }
# light-intensity fill ramp (10-100%): pale to deep purple, lightness rising steadily so it reads
# for color-blind viewers too; points drawn with a dark outline so 10% stays visible
SB_INT_FILL <- c("10" = "#DADAEB", "30" = "#BCBDDC", "50" = "#9E9AC8", "70" = "#756BB1", "100" = "#54278F")
# x-axis label for distance gradients; a run with another distance measure can set this first
if (!exists("SB_DIST_LABEL")) SB_DIST_LABEL <- "Distance from lot (km)"
