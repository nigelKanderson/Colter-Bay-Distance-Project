# 53_bat_species_hierarchy_panel.R
# Bat counterpart of the moth taxonomic distance hierarchy (R/50, Distance_Project.qmd
# `moth_sig_hierarchy`): one tree per storyboard site set (R/storyboard_sites.R), stacked A/B/C.
#   tiers: All bats -> call-frequency guild -> species (9 focal species)
#   each node: inset of model-predicted detections per night vs distance (auto-scaled) +/- 95% CI
# Unlike the moth tree (significant taxa only), every node is drawn, since there are only 13:
# a significant slope is a solid box colored by direction (teal = more active away from the lot,
# brown = toward the lot); a non-significant one is a dashed gray box.
# Models per tier, as in the moth hierarchy: zero-inflated NB GLMM
#   detections ~ taxon * dist_km + color + intensity + moon phase + non-forest + sky brightness
#                + season (jd_c + jd_c^2) + (1 | site)
# (the all-bats tier has a single taxon, so no taxon terms); covariates are dropped in order if
# a model doesn't converge. The all-bats node is a single test (the same question as the
# storyboard's bat vs moth figure), so its p is left unadjusted; the 3 guild and 9 species nodes
# are Benjamini-Hochberg-adjusted together within each dataset (bat_hier_adjust()).
# Output: output/figures/fig_bat_species_hierarchy.png

source("R/storyboard_sites.R")   # SB_SETS, SB_LABELS

# p = p_raw for the all-bats node; BH across the 12 guild + species nodes
bat_hier_adjust <- function(sumry) {
  sub <- sumry$lev != "All"
  sumry$p <- sumry$p_raw
  sumry$p[sub] <- p.adjust(sumry$p_raw[sub], "BH")
  sumry
}

bat_hier_fit <- function(dat, taxa = NULL) {
  library(glmmTMB); library(emmeans); library(tidyverse)
  stbl <- function(m) !is.null(m) && isTRUE(m$sdr$pdHess) &&
    !any(is.na(suppressWarnings(sqrt(diag(vcov(m)$cond)))))
  cand <- c("mean_phase + pct_nonforest + brightness_dark + jd_c + I(jd_c^2)",
            "mean_phase + brightness_dark + jd_c + I(jd_c^2)", "mean_phase + jd_c + I(jd_c^2)",
            "mean_phase + jd_c", "jd_c")
  lhs <- if (is.null(taxa)) "detections ~ dist_km" else "detections ~ taxon * dist_km"
  m <- NULL; used <- NA
  for (rhs in cand) {
    mi <- tryCatch(suppressWarnings(glmmTMB(
      as.formula(paste(lhs, "+ color + intensity +", rhs, "+ (1 | site)")),
      data = dat, ziformula = ~1, family = nbinom2())), error = function(e) NULL)
    if (is.null(m) && !is.null(mi)) { m <- mi; used <- rhs }
    if (stbl(mi)) { m <- mi; used <- rhs; break }
  }
  spec <- if (is.null(taxa)) ~ 1 else ~ taxon
  s  <- as.data.frame(emtrends(m, spec, var = "dist_km"))
  tx <- if (is.null(taxa)) "All bats" else as.character(s$taxon)
  dg <- seq(min(dat$dist_km), max(dat$dist_km), length.out = 40)
  cv <- as.data.frame(emmeans(m, if (is.null(taxa)) ~ dist_km else ~ taxon | dist_km,
                              at = list(dist_km = dg), type = "response"))
  list(sumry = tibble(taxon = tx, est = s$dist_km.trend, SE = s$SE,
                      p_raw = 2 * pnorm(-abs(s$dist_km.trend / s$SE)), covariates = used),
       curves = tibble(taxon = if (is.null(taxa)) "All bats" else as.character(cv$taxon),
                       dist_km = cv$dist_km, response = cv$response,
                       lower = cv$asymp.LCL, upper = cv$asymp.UCL))
}

fit_bat_species_hierarchy <- function(site_night, cluster_night, spp_night, sets = SB_SETS) {
  library(tidyverse)
  guild_lab <- c(LowFreq = "Low frequency", Myev = "Myev", HighFreq = "High frequency")
  drop <- function(d, s) droplevels(filter(d, !site %in% s))
  prep <- function(d) mutate(d, site = factor(site), color = factor(color),
                             intensity = factor(intensity, c("10", "30", "50", "70", "100")),
                             jd_c = if ("jd_c" %in% names(d)) jd_c else jd - mean(jd, na.rm = TRUE))
  map(sets, function(s) {
    all <- bat_hier_fit(prep(drop(site_night, s)) %>% rename(detections = total_det))
    gld <- bat_hier_fit(prep(drop(cluster_night, s)) %>%
                          mutate(taxon = factor(guild_lab[as.character(cluster)], guild_lab)), taxa = TRUE)
    spp <- bat_hier_fit(prep(drop(spp_night, s)) %>% mutate(taxon = factor(species)), taxa = TRUE)
    lv <- list(All = all, Guild = gld, Species = spp)
    sumry  <- bat_hier_adjust(imap_dfr(lv, ~ mutate(.x$sumry, lev = .y)))
    curves <- imap_dfr(lv, ~ mutate(.x$curves, lev = .y))
    list(sumry = sumry, curves = curves)
  })
}

bat_hier_plot <- function(res, title, subtitle) {
  library(tidyverse)
  spp_guild <- c(Epfu = "Low frequency", Lano = "Low frequency", Laci = "Low frequency",
                 Myev = "Myev", Mylu = "High frequency", Myvo = "High frequency",
                 Myci = "High frequency", Myca = "High frequency", Myyu = "High frequency")
  common <- c(Epfu = "big brown", Lano = "silver-haired", Laci = "hoary", Myev = "long-eared myotis",
              Mylu = "little brown", Myvo = "long-legged myotis", Myci = "small-footed myotis",
              Myca = "California myotis", Myyu = "Yuma myotis")
  stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
  SPX <- 1.6; BOXW <- 1.15; BOXH <- 0.56
  pal <- c("Away from lot" = "#018571", "Toward lot" = "#8C510A", "Not significant" = "#9A9A9A")   # teal / brown: not used elsewhere in the storyboard

  # layout: species leaves left to right by guild; guilds over their species; All bats over all
  sp  <- names(spp_guild); sx <- seq_along(sp) * SPX
  gx  <- tapply(sx, factor(spp_guild, unique(spp_guild)), mean)
  nd <- res$sumry %>%
    mutate(y = c(All = 3, Guild = 2, Species = 1)[lev],
           x = case_when(lev == "All" ~ mean(sx), lev == "Guild" ~ gx[taxon], TRUE ~ sx[match(taxon, sp)]),
           sig = p < 0.05,
           dir = ifelse(!sig, "Not significant", ifelse(est > 0, "Away from lot", "Toward lot")),
           key = paste(lev, taxon),
           parent = case_when(lev == "Guild" ~ "All All bats",
                              lev == "Species" ~ paste("Guild", spp_guild[taxon]), TRUE ~ NA_character_),
           lab = case_when(lev == "Species" ~ paste0(taxon, " (", common[taxon], ") ", stars(p)),
                           TRUE ~ paste0(taxon, " ", stars(p))),
           pct = sprintf("%+.0f%%/km", (exp(est) - 1) * 100))
  cv <- res$curves %>% inner_join(nd %>% select(lev, taxon, key, nx = x, ny = y, dir), by = c("lev", "taxon")) %>%
    group_by(key) %>%
    mutate(lo = min(lower), hi = max(upper), rng = ifelse(hi - lo <= 0, 1, hi - lo),
           xr = (dist_km - min(dist_km)) / (max(dist_km) - min(dist_km)),
           X = nx + (xr - 0.5) * BOXW, Yl = ny + ((response - lo) / rng - 0.5) * BOXH,
           Ylo = ny + ((lower - lo) / rng - 0.5) * BOXH, Yhi = ny + ((upper - lo) / rng - 0.5) * BOXH) %>%
    ungroup() %>% arrange(key, X)
  box <- nd %>% transmute(key, dir, sig, xmin = x - BOXW / 2, xmax = x + BOXW / 2,
                          ymin = y - BOXH / 2, ymax = y + BOXH / 2)
  seg <- nd %>% filter(!is.na(parent)) %>% left_join(nd %>% select(parent = key, px = x, py = y), by = "parent")
  tiers <- data.frame(y = 3:1)

  ggplot() +
    geom_hline(data = tiers, aes(yintercept = y), color = "#F2F2F2", linewidth = 10) +
    geom_segment(data = seg, aes(x = px, y = py - BOXH / 2 - 0.36, xend = x, yend = y + BOXH / 2),   # from below the parent's labels
                 color = "#C4C4C4", linewidth = 0.45) +
    geom_rect(data = box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = dir, linetype = sig),
              fill = "white", linewidth = 0.55, show.legend = TRUE) +   # TRUE: keys for every level
    scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"), guide = "none") +
    geom_ribbon(data = cv, aes(x = X, ymin = Ylo, ymax = Yhi, group = key, fill = dir), alpha = 0.22, show.legend = TRUE) +
    geom_line(data = cv, aes(x = X, y = Yl, group = key, color = dir), linewidth = 0.7, show.legend = TRUE) +
    geom_text(data = nd, aes(x = x, y = y - BOXH / 2 - 0.08, label = lab), size = 2.8, vjust = 1, color = "#222222") +
    geom_text(data = nd, aes(x = x, y = y - BOXH / 2 - 0.25, label = pct), size = 2.5, vjust = 1, color = "#666666") +
    scale_color_manual(values = pal, limits = names(pal), name = "More active:") +   # same legend in every panel
    scale_fill_manual(values = pal, limits = names(pal), name = "More active:") +
    scale_y_continuous(breaks = 3:1, labels = c("All bats", "Guild", "Species"), limits = c(0.45, 3.4)) +
    coord_cartesian(xlim = c(min(sx) - BOXW, max(sx) + BOXW / 2 + 0.2), clip = "off") +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_family = "Arial") +
    theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 9, color = "#555555"),
          axis.text.y = element_text(size = 10, face = "bold", color = "#333333", hjust = 1), legend.position = "top",
          legend.title = element_text(face = "bold", size = 9), legend.text = element_text(size = 9),
          plot.margin = margin(8, 14, 8, 8))
}

make_bat_species_hierarchy <- function(res, col_labels = SB_LABELS,
                                       file = "output/figures/fig_bat_species_hierarchy.png") {
  library(patchwork)
  res <- map(res, ~ { .x$sumry <- bat_hier_adjust(.x$sumry); .x })   # also updates fits saved earlier
  sub <- paste0("each inset = model-predicted bat detections per night vs distance (0–max km, auto-scaled) ± 95% CI; ",
                "gray text = % change per km\nsolid colored box = p < 0.05 (all bats: single test; guilds and species: ",
                "Benjamini–Hochberg-adjusted across those 12 nodes); dashed gray box = not significant")
  ps <- imap(res, ~ bat_hier_plot(.x, paste("Bat species hierarchy (distance) —", col_labels[[.y]]), sub))
  plate <- wrap_plots(ps, ncol = 1) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plate, width = 16, height = 5.6 * length(ps), dpi = 300, bg = "white")
  invisible(plate)
}
