# =============================================================================
# 33_moth_taxonomy.R
# Extend the moth analysis to subfamily / genus / species resolution:
#   - build_taxon_env(): a site-night x taxon detection frame at any level,
#     shaped like insect_family_env so run_community_analysis() can use it.
#   - compute_richness(): per-site-night taxon richness (species/genus/subfamily).
#   - run_richness_models(): GLMMs of richness ~ color*intensity +
#     intensity*dist_km + covariates, adjusted for catch size, + a figure.
#
# Species is the Genus + species binomial (rows with species == "sp." are
# dropped from the species level and from species richness).
# =============================================================================

suppressMessages({ library(tidyverse); library(glmmTMB); library(emmeans) })

# Add a clean binomial column to the specimen table
add_species_bin <- function(specimens) {
  specimens %>%
    mutate(Species_bin = if_else(!is.na(Genus) & !is.na(Species) & Species != "sp.",
                                 paste(Genus, Species), NA_character_))
}

# Build a site-night x taxon env frame at `level` (e.g. "Subfamily","Genus",
# "Species_bin"); `event_covariates` is one row per site/date/color/intensity.
build_taxon_env <- function(specimens, event_covariates, level) {
  specimens <- add_species_bin(specimens)
  specimens %>%
    filter(!is.na(.data[[level]])) %>%
    group_by(site, date, color, intensity, taxon = .data[[level]]) %>%
    summarise(detections = n(), .groups = "drop") %>%
    left_join(event_covariates, by = c("site", "date", "color", "intensity")) %>%
    filter(!is.na(dist_km)) %>%
    rename(species = taxon) %>%
    mutate(color = factor(color), site = factor(site),
           intensity = factor(intensity, levels = c("10","30","50","70","100")))
}

# Per-site-night richness at each level + covariates
compute_richness <- function(specimens, event_covariates) {
  specimens <- add_species_bin(specimens)
  specimens %>%
    group_by(site, date, color, intensity) %>%
    summarise(rich_species = n_distinct(Species_bin, na.rm = TRUE),
              rich_genus   = n_distinct(Genus,       na.rm = TRUE),
              rich_subfam  = n_distinct(Subfamily,   na.rm = TRUE),
              n_specimens  = n(), .groups = "drop") %>%
    left_join(event_covariates, by = c("site", "date", "color", "intensity")) %>%
    filter(!is.na(dist_km)) %>%
    mutate(color = factor(color), site = factor(site),
           intensity = factor(intensity, levels = c("10","30","50","70","100")),
           jd_c = jd - mean(jd, na.rm = TRUE),
           log_catch = log1p(n_specimens))
}

# Fit richness GLMMs (one per level) + build a summary figure.
# log_catch (log catch size) is included so treatment effects are assessed
# beyond what raw abundance explains.
run_richness_models <- function(rich, bat_theme, intensity_pal, col_red, col_white) {
  levels_map <- c(rich_species = "Species", rich_genus = "Genus", rich_subfam = "Subfamily")
  has_bd <- "brightness_dark" %in% names(rich)
  rhs <- paste("~ color * intensity + intensity * dist_km + mean_phase +",
               "pct_nonforest +", if (has_bd) "brightness_dark +" else "",
               "log_catch + jd_c + I(jd_c^2) + (1 | site)")

  models <- list(); ci_all <- list()
  for (resp in names(levels_map)) {
    form <- as.formula(paste(resp, rhs))
    m <- glmmTMB(form, data = rich, family = nbinom2)
    # richness counts are often ~Poisson; nbinom2 dispersion can hit a boundary
    # (NaN SEs) on the low-count levels -> fall back to Poisson.
    bad <- !isTRUE(m$sdr$pdHess) ||
      any(is.na(suppressWarnings(sqrt(diag(vcov(m)$cond)))))
    if (bad) {
      m <- glmmTMB(form, data = rich, family = poisson)
      cat("  [", levels_map[[resp]], "] nbinom2 unstable -> refit with Poisson\n")
    }
    models[[resp]] <- m
    cat("\n=== Richness model:", levels_map[[resp]], "===\n"); print(summary(m))
    cat("\n-- Red vs White at each intensity (", levels_map[[resp]], ") --\n")
    print(pairs(emmeans(m, ~ color | intensity, type = "response"), adjust = "holm"))
    ci <- as.data.frame(emmeans(m, ~ color | intensity, type = "response")) %>%
      rename(lower = asymp.LCL, upper = asymp.UCL) %>%
      mutate(level = levels_map[[resp]],
             intensity = factor(intensity, levels = c("10","30","50","70","100")))
    ci_all[[resp]] <- ci
  }

  fig <- bind_rows(ci_all) %>%
    mutate(level = factor(level, levels = c("Subfamily","Genus","Species"))) %>%
    ggplot(aes(intensity, response, color = color, group = color)) +
    geom_line(linewidth = 0.8, position = position_dodge(0.3)) +
    geom_point(size = 2.5, position = position_dodge(0.3)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15,
                  position = position_dodge(0.3)) +
    scale_color_manual(values = c(R = col_red, W = col_white),
                       labels = c(R = "Red", W = "White"), name = "Light color") +
    facet_wrap(~ level, scales = "free_y") +
    labs(title = "Moth taxonomic richness by light color and intensity",
         subtitle = "Predicted richness per site-night ± 95% CI · NB GLMM (catch-size adjusted)",
         x = "Light intensity (%)", y = "Predicted taxa per site-night") +
    bat_theme + theme(legend.position = "top")

  list(models = models, fig_richness = fig, ci = bind_rows(ci_all))
}
