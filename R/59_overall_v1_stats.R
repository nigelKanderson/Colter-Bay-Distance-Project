# =============================================================================
# 59_overall_v1_stats.R
# Statistics behind each figure of distance_storyboard_overall_v1.qmd, as tables for the
# companion document distance_storyboard_overall_v1_stats.qmd:
#   fig1  station sky brightness (lights off / white 100%); bat and moth distance slopes by color
#         (11 stations within 2.1 km; same refits as R/56, incl. the night random effect)
#   fig2  per-intensity distance slopes and white-vs-red contrasts (Holm across the 5 intensities)
#   fig3  bat hierarchy nodes (output/bat_species_hierarchy_fits.rds, far sites removed)
#   fig4  moth hierarchy: every taxon tested and the BH-significant nodes
#         (output/overall_v1_moth_hierarchy_stats.rds, from the `moth_sig_hierarchy` chunk)
#   fig5  as fig1 with the far stations, GRTE15 excluded (13 stations)
# Output: output/overall_v1_stats.rds
# =============================================================================

make_overall_v1_stats <- function(bat_model, moth_model, site_night, insect_total_env,
                                  bat_fits = readRDS("output/bat_species_hierarchy_fits.rds"),
                                  moth_hier = readRDS("output/overall_v1_moth_hierarchy_stats.rds"),
                                  file = "output/overall_v1_stats.rds") {
  library(tidyverse); library(glmmTMB); library(emmeans)
  source("R/storyboard_sites.R"); source("R/02_clean_data.R")
  pct <- function(x) (exp(x) - 1) * 100

  refit <- function(m, d, drop) {
    d <- droplevels(filter(d, !site %in% drop))
    update(m, . ~ . + color:dist_km + (1 | date), data = d)   # as make_manager_summary(night_re = TRUE)
  }
  info <- function(m, resp) tibble(model = resp, n_obs = nobs(m), n_sites = n_distinct(m$frame$site),
                                   formula = paste(deparse(formula(m), width.cutoff = 500), collapse = " "),
                                   family = family(m)$family, converged = isTRUE(m$sdr$pdHess))
  # distance slope under each color
  slopes <- function(m) as.data.frame(test(emtrends(m, ~ color, var = "dist_km"))) %>%
    transmute(color = recode(as.character(color), R = "Red", W = "White"),
              pct_per_km = pct(dist_km.trend), lo = pct(dist_km.trend - 1.96 * SE), hi = pct(dist_km.trend + 1.96 * SE),
              z = z.ratio, p = p.value)
  by_int_dist <- function(m) as.data.frame(test(emtrends(m, ~ intensity, var = "dist_km"))) %>%
    transmute(intensity = paste0(intensity, "%"), pct_per_km = pct(dist_km.trend),
              lo = pct(dist_km.trend - 1.96 * SE), hi = pct(dist_km.trend + 1.96 * SE),
              z = z.ratio, p = p.value, p_holm = p.adjust(p.value, "holm"))
  by_int_wvr <- function(m) as.data.frame(summary(contrast(emmeans(m, ~ color | intensity),
                                                           list(WvR = c(-1, 1))), infer = TRUE)) %>%
    transmute(intensity = paste0(intensity, "%"), pct_white_vs_red = pct(estimate),
              lo = pct(asymp.LCL), hi = pct(asymp.UCL), z = z.ratio, p = p.value, p_holm = p.adjust(p.value, "holm"))

  # ---- station sky brightness ---------------------------------------------------
  sky <- read_csv("data/sky_brightness_measurements.csv", show_col_types = FALSE) %>%
    transmute(site, lights_off = brightness_dark,
              white100 = brightness_white100, white100_filled = is.na(brightness_white100),
              white100_map = coalesce(brightness_white100, brightness_dark))
  ref <- sf::st_transform(sf::st_as_sf(data.frame(lon = -110.645134, lat = 43.904879), coords = c("lon", "lat"), crs = 4326), 32612)
  st  <- sf::st_transform(sf::st_as_sf(sites, coords = c("lon", "lat"), crs = 4326), 32612)
  sky <- sky %>% left_join(tibble(site = sites$site, dist_km = as.numeric(sf::st_distance(st, ref)) / 1000), by = "site") %>%
    mutate(lot = site %in% c("GRTE01", "GRTE16", "GRTE17"), far = site %in% SB_FAR_DROP) %>% arrange(dist_km)

  # ---- figures 1-2 (far removed) and 5 (GRTE15 excluded) ---------------------------
  out <- list()
  for (set in c("nofar", "full")) {
    drop <- SB_SETS[[set]]
    b <- refit(bat_model, site_night, drop); m <- refit(moth_model, insect_total_env, drop)
    out[[set]] <- list(models = bind_rows(info(b, "Bats: total detections per site-night"),
                                          info(m, "Moths: specimens per trap night")),
                       bat_slopes = slopes(b), moth_slopes = slopes(m))
    if (set == "nofar") out$fig2 <- list(bat_dist = by_int_dist(b), moth_dist = by_int_dist(m),
                                         bat_wvr = by_int_wvr(b), moth_wvr = by_int_wvr(m))
  }

  # ---- figure 3: bat hierarchy ----------------------------------------------------
  source("R/53_bat_species_hierarchy_panel.R")
  bh <- bat_hier_adjust(bat_fits[["nofar"]]$sumry) %>%
    transmute(tier = recode(lev, All = "All bats", Guild = "Guild", Species = "Species"), taxon,
              pct_per_km = pct(est), lo = pct(est - 1.96 * SE), hi = pct(est + 1.96 * SE),
              z = est / SE, p_raw, p_adj = p, adjustment = ifelse(lev == "All", "none (single test)", "BH (12 nodes)"),
              covariates)

  # ---- figure 4: moth hierarchy ---------------------------------------------------
  lev_lab <- c(Family = "Family", Subfamily = "Subfamily", Genus = "Genus", Species_bin = "Species")
  mt <- moth_hier$tested %>% transmute(level = lev_lab[lev], taxon, pct_per_km = pct(est), p_raw, p_BH,
                                       significant = p_BH < 0.05) %>%
    arrange(factor(level, lev_lab), p_BH)

  res <- list(sky = sky, fig1 = out$nofar, fig2 = out$fig2, fig3 = bh, fig4 = mt, fig5 = out$full)
  saveRDS(res, file)
  invisible(res)
}
