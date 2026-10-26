# 48_species_white_vs_red_100.R
# White vs. red light at a given intensity (100%, and 50% for the combined figure), by bat
# species (forest plot):
#   A  full data (14 sites)      B  far sites removed (11 sites)
# Each focal species gets its own GLMM on the zero-filled species-night data (`spp_night`,
# single-species detections only), same structure as the guild models:
#   detections ~ color x intensity + intensity x distance + color x distance
#                + moon phase + non-forest + sky brightness + season + (1 | site)
#   (NB; zero-inflated where the species has zero nights)
# Points = % change white vs red at that intensity (contrast at mean distance), 95% CI, log x-axis
# labeled in %; asterisks = differs from 0 (* p<0.05, ** p<0.01, *** p<0.001; Wald z,
# Holm-adjusted across the species in each panel). Colors = call-frequency guild (as in the guild figures). Hollow points =
# rare species (< 150 detections at that intensity in that dataset).
# Output: output/figures/fig_species_white_vs_red_100.png; combined guild + species figure:
#         output/figures/fig_bat_guild_species_white_vs_red.png

species_guild <- c(Epfu = "Low frequency", Lano = "Low frequency", Laci = "Low frequency",
                   Myev = "Myev",
                   Mylu = "High frequency", Myvo = "High frequency", Myci = "High frequency",
                   Myca = "High frequency", Myyu = "High frequency")
species_common <- c(Epfu = "big brown", Lano = "silver-haired", Laci = "hoary",
                    Myev = "long-eared myotis", Mylu = "little brown", Myvo = "long-legged myotis",
                    Myci = "small-footed myotis", Myca = "California myotis", Myyu = "Yuma myotis")

# one row per dataset x species x intensity (white vs red contrast at each of `intensities`)
fit_species_white_vs_red_100 <- function(spp_night, far_sites = c("GRTE13", "GRTE14", "GRTE15"),
                                         intensities = c("50", "100")) {
  library(glmmTMB)
  library(emmeans)
  library(tidyverse)
  f <- detections ~ color * intensity + intensity * dist_km + color * dist_km +
    mean_phase + pct_nonforest + brightness_dark + jd_c + I(jd_c^2) + (1 | site)
  sets <- list(full = spp_night, nofar = droplevels(filter(spp_night, !site %in% far_sites)))
  imap_dfr(sets, function(d0, ds) map_dfr(names(species_guild), function(sp) {
    d  <- droplevels(filter(d0, species == sp))
    zi <- any(d$detections == 0)
    m  <- tryCatch(glmmTMB(f, data = d, ziformula = if (zi) ~1 else ~0, family = nbinom2()),
                   error = function(e) NULL)
    ok <- !is.null(m) && isTRUE(m$sdr$pdHess) && all(is.finite(sqrt(diag(vcov(m)$cond))))
    out <- tibble(dataset = ds, species = sp, guild = species_guild[[sp]],
                  detections = sum(d$detections), intensity = intensities,
                  det_at = map_dbl(intensities, ~ sum(d$detections[d$intensity == .x])),
                  converged = ok)
    if (!ok) return(out)
    ct <- as.data.frame(contrast(emmeans(m, ~ color | intensity), "revpairwise", adjust = "none")) %>%
      mutate(intensity = as.character(intensity)) %>% filter(intensity %in% intensities)
    out %>% left_join(ct %>% transmute(intensity, r = exp(estimate), lo = exp(estimate - 1.96 * SE),
                                       hi = exp(estimate + 1.96 * SE), p = p.value), by = "intensity")
  }))
}

species_white_vs_red_panels <- function(res, bat_theme, level = "100", rare_n = 150) {
  library(tidyverse)
  library(patchwork)
  pal   <- c("Low frequency" = "#0072B2", "Myev" = "#CC79A7", "High frequency" = "#9A7400")
  stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")

  # species order: by guild (as in the legend), then by total detections (most at top)
  ord <- res %>% filter(dataset == "full") %>% distinct(species, guild, detections) %>%
    mutate(guild = factor(guild, names(pal))) %>% arrange(desc(guild), detections) %>% pull(species)
  d <- res %>% filter(converged, intensity == level) %>%
    group_by(dataset) %>% mutate(p_raw = p, p = p.adjust(p, "holm")) %>% ungroup() %>%   # Holm across species
    mutate(guild = factor(guild, names(pal)),
           label = factor(paste0(species, " (", species_common[species], ")"),
                          paste0(ord, " (", species_common[ord], ")")),
           rare  = det_at < rare_n,
           n_lab = paste0("n = ", format(det_at, big.mark = ",", trim = TRUE)))

  # x-range from every converged contrast (all intensities), so panels for different
  # intensities share one axis and can be compared directly
  all_ok <- filter(res, converged)
  xl <- range(c(all_ok$lo, all_ok$hi, 1))
  xl <- xl * c(1 / 1.05, (xl[2] / xl[1])^0.32)   # room on the right for n and asterisks
  brks <- c(0.25, 0.5, 1, 2, 4)
  lab_pct <- function(r) { v <- round((r - 1) * 100); paste0(ifelse(v > 0, "+", ""), v, "%") }

  panel <- function(ds, ttl) {
    dd <- filter(d, dataset == ds)
    ggplot(dd, aes(r, label, color = guild)) +
      geom_vline(xintercept = 1, color = "#999999", linewidth = 0.4) +
      geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, linewidth = 0.6, orientation = "y") +
      geom_point(aes(shape = rare), size = 2.8, stroke = 1, fill = "white") +
      geom_text(aes(x = hi, label = stars(p)), hjust = -0.25, vjust = 0.75, size = 4.2,
                color = "#222222", show.legend = FALSE) +
      geom_text(aes(x = xl[2], label = n_lab), hjust = 1, size = 3, color = "#666666",
                show.legend = FALSE) +
      scale_x_log10(breaks = brks, labels = lab_pct) +
      scale_color_manual(values = pal, name = NULL, drop = FALSE) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21),
                         breaks = "TRUE", labels = paste0("rare (< ", rare_n, " detections at that intensity)"),
                         name = NULL) +
      coord_cartesian(xlim = xl) +
      labs(title = ttl, x = paste0("% change at ", level, "% intensity, white vs red"), y = NULL) +
      bat_theme +
      theme(plot.title = element_text(size = 11), axis.text = element_text(size = 10),
            axis.title = element_text(size = 11),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(color = "#EBEBEB", linewidth = 0.4))
  }

  list(full = panel("full", "Full data (14 sites)"), nofar = panel("nofar", "Far sites removed (11 sites)"))
}

plot_species_white_vs_red_100 <- function(res, bat_theme, rare_n = 150,
                                          file = "output/figures/fig_species_white_vs_red_100.png") {
  library(patchwork)
  pn <- species_white_vs_red_panels(res, bat_theme, level = "100", rare_n = rare_n)
  plate <- (pn$full | pn$nofar) +
    plot_layout(guides = "collect") +
    plot_annotation(
      tag_levels = "A",
      title = "Bat species: white vs. red light at 100% intensity",
      subtitle = paste0("Each species' GLMM (as for the guilds; single-species detections only), contrast at the mean distance; ",
                        "95% CIs, log x-axis. n = detections at 100% intensity.\n",
                        "Asterisks: differs from 0 (* p < 0.05, ** p < 0.01, *** p < 0.001; Holm-adjusted across the 9 species)."),
      theme = theme(plot.title = element_text(face = "bold", size = 15),
                    plot.subtitle = element_text(size = 10, color = "#555555"))) &
    theme(legend.position = "top", plot.tag = element_text(face = "bold", size = 14))

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plate, width = 13, height = 5.8, dpi = 300, bg = "white")
  plate
}

# Combined figure: A-B guilds, white vs red at each intensity (from the per-guild models of
# 45_*.R via 47_*.R); C-D species at 100% and E-F species at 50% (fit_species_white_vs_red_100()).
# `guild_fits` = make_bat_guild_pct_panels()$models; `res` = fit_species_white_vs_red_100() output.
make_guild_species_white_vs_red <- function(guild_fits, res, bat_theme,
                                            file = "output/figures/fig_bat_guild_species_white_vs_red.png") {
  library(tidyverse)
  library(patchwork)
  source("R/47_bat_moth_white_vs_red_panel.R")
  gp <- make_guild_white_vs_red_panel(guild_fits, bat_theme, file = NULL)$panels
  no_col <- function(pn) map(pn, ~ .x + guides(color = "none"))   # guild colors already in the legend from A-B
  sp  <- no_col(species_white_vs_red_panels(res, bat_theme, level = "100"))
  sp50 <- no_col(species_white_vs_red_panels(res, bat_theme, level = "50"))
  retitle <- function(p, ttl) p + labs(title = ttl)

  plate <- (retitle(gp$full,  "Guilds, all intensities · Full data (14 sites)") |
            retitle(gp$nofar, "Guilds, all intensities · Far sites removed (11 sites)")) /
           (retitle(sp$full,  "Species at 100% intensity · Full data (14 sites)") |
            retitle(sp$nofar, "Species at 100% intensity · Far sites removed (11 sites)")) /
           (retitle(sp50$full,  "Species at 50% intensity · Full data (14 sites)") |
            retitle(sp50$nofar, "Species at 50% intensity · Far sites removed (11 sites)")) +
    plot_layout(guides = "collect", heights = c(1, 1.15, 1.15)) +
    plot_annotation(
      tag_levels = "A",
      title = "Bat guilds and species: white vs. red light",
      subtitle = paste0("A–B: each guild's zero-inflated NB GLMM, white vs. red at each intensity (log y-axis). ",
                        "C–F: each species' GLMM (single-species detections only), white vs. red at 100% (C–D)\n",
                        "and 50% (E–F) intensity (log x-axis; n = detections at that intensity; hollow = rare, < 150). ",
                        "Contrasts at the mean distance; 95% CIs. Asterisks: differs from 0 ",
                        "(* p < 0.05, ** p < 0.01, *** p < 0.001; Holm-adjusted across intensities within each guild (A–B)\nand across the 9 species in each panel (C–F))."),
      theme = theme(plot.title = element_text(face = "bold", size = 15),
                    plot.subtitle = element_text(size = 10, color = "#555555"))) &
    theme(legend.position = "top", plot.tag = element_text(face = "bold", size = 14))

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plate, width = 13, height = 16.5, dpi = 300, bg = "white")
  plate
}
