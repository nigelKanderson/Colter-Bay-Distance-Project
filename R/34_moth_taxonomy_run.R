# =============================================================================
# 34_moth_taxonomy_run.R
# Standalone runner for the moth taxonomic (subfamily/genus/species) community
# analyses + richness models. Saves:
#   - output/figures/panels/panel9_moth_richness.{pdf,png}
#   - output/moth_taxonomy/richness_color_intensity.csv
#   - output/moth_taxonomy/community_dbrda_<level>.csv
# Reconstructs moth event covariates from caches (no Google Drive needed).
# =============================================================================

suppressMessages({
  library(tidyverse); library(glmmTMB); library(emmeans); library(patchwork)
})
source("R/03_extract_habitat.R"); source("R/13_import_insects.R")
source("R/33_moth_taxonomy.R");   source("R/10_community_analysis.R")
dir.create("output/moth_taxonomy", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures/panels", showWarnings = FALSE, recursive = TRUE)

col_red <- "#7A1C2E"; col_white <- "#C9A84C"; col_navy <- "#1B2A4A"
col_purple <- "#4A1259"; col_ivory <- "#F5F0E3"
intensity_pal <- c("10" = col_ivory, "30" = col_white, "50" = col_purple,
                   "70" = col_red, "100" = col_navy)
panel_theme <- theme_classic(base_size = 16, base_family = "Arial") +
  theme(panel.grid.major.y = element_line(color = "#EBEBEB", linewidth = 0.4),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "#666666", size = 10),
        axis.title = element_text(size = 18, color = "#222222"),
        axis.text = element_text(size = 16, color = "#333333"),
        legend.title = element_text(size = 13), legend.text = element_text(size = 13),
        strip.text = element_text(size = 13, face = "bold"),
        legend.position = "top")

# ---- reconstruct moth event covariates --------------------------------------
sky <- read_csv("data/sky_brightness_measurements.csv", show_col_types = FALSE) %>%
  mutate(site = toupper(site))
ie <- readRDS("data/insect_data_out.rds")
ie$datetime <- as.POSIXct(paste(ie$date, "12:00:00"), tz = "America/Denver")
ie <- add_moonlight(ie) %>% left_join(sky %>% select(site, brightness_dark), by = "site")
event_covariates <- ie %>%
  select(site, date, color, intensity, dist_km, dist_cat, jd,
         mean_phase, pct_nonforest, brightness_dark)
ins <- import_insects("data/grte_distance_insectID.xlsx")

# ---- richness models + panel figure -----------------------------------------
message("Fitting richness models ...")
rich <- compute_richness(ins$specimens, event_covariates)
rr <- run_richness_models(rich, panel_theme, intensity_pal, col_red, col_white)
write_csv(rr$ci, "output/moth_taxonomy/richness_color_intensity.csv")
ggsave("output/figures/panels/panel9_moth_richness.pdf", rr$fig_richness,
       width = 13, height = 5, device = cairo_pdf)
ggsave("output/figures/panels/panel9_moth_richness.png", rr$fig_richness,
       width = 13, height = 5, dpi = 200, bg = "white", device = ragg::agg_png)
message("Richness panel written.")

# ---- community analysis at each taxonomic level -----------------------------
for (lv in c("Subfamily", "Genus", "Species_bin")) {
  message("Community analysis: ", lv, " ...")
  env <- build_taxon_env(ins$specimens, event_covariates, lv)
  cc  <- run_community_analysis(env, panel_theme, col_red, col_white,
                                taxon = "Moth", unit = lv, nBoot = 999)
  as.data.frame(cc$db_margin) %>% rownames_to_column("term") %>%
    write_csv(paste0("output/moth_taxonomy/community_dbrda_", lv, ".csv"))
}
message("DONE moth taxonomy run")
