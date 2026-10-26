# 52_community_panels.R
# Bat (9 focal species) and moth (9 families) community composition: what drives it and
# which taxa carry each effect, for the full data and with the far sites (GRTE13/14/15) removed.
# Uses run_community_analysis() (10_community_analysis.R): one sample = site x color x
# intensity (detections per night); predictors color, intensity, distance, moon phase,
# % non-forest, sky brightness.
#   A  Predictor importance: % of community variation explained by each predictor alone
#      (marginal dbRDA sum of squares / total, Bray-Curtis), asterisks from the marginal
#      permutation tests; bats vs moths, full data | far sites removed.
#   B  Taxon responses (mvabund NB multivariate GLM): % change in each taxon per 1 SD of each
#      continuous predictor (white vs red for color), asterisks from mvabund's per-taxon tests
#      (adjusted for multiple testing across taxa); intensity (5 levels) shown separately as
#      deviance explained.
# Output: output/figures/fig_community_panels.png; fits cached in output/community_fits.rds

fit_community_sets <- function(bat_data, moth_data, bat_theme, col_red, col_white,
                               far_sites = c("GRTE13", "GRTE14", "GRTE15"), nBoot = 9999,
                               cache = "output/community_fits.rds") {
  library(tidyverse)
  source("R/10_community_analysis.R")
  keep <- function(r) list(perm = r$perm, db_margin = r$db_margin, env = r$env,
                           coef = coef(r$mv_fit), mv_anova = r$mv_anova)
  run <- function(d, taxon, unit) keep(run_community_analysis(
    d, bat_theme, col_red, col_white, taxon = taxon, unit = unit, nBoot = nBoot))
  drop <- function(d) droplevels(filter(d, !site %in% far_sites))
  res <- list(
    full  = list(Bats  = run(bat_data, "Bat", "Species"),
                 Moths = run(moth_data, "Moth", "Family")),
    nofar = list(Bats  = run(drop(bat_data), "Bat", "Species"),
                 Moths = run(drop(moth_data), "Moth", "Family"))
  )
  dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, cache)
  res
}

make_community_panels <- function(res, bat_theme, file = "output/figures/fig_community_panels.png") {
  library(tidyverse)
  library(patchwork)

  pred_lab <- c(color = "Light color", intensity = "Light intensity", dist_km = "Distance",
                mean_phase = "Moon phase", pct_nonforest = "% non-forest",
                brightness_dark = "Sky brightness")
  ds_lab   <- c(full = "Full data (14 sites)", nofar = "Far sites removed (11 sites)")
  tx_pal   <- c(Bats = "#0072B2", Moths = "#D55E00")
  stars    <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
  bat_names  <- c(Epfu = "big brown", Lano = "silver-haired", Laci = "hoary", Myev = "long-eared myotis",
                  Mylu = "little brown", Myvo = "long-legged myotis", Myci = "small-footed myotis",
                  Myca = "California myotis", Myyu = "Yuma myotis")

  # ── A: variance explained by each predictor ─────────────────────────────────
  imp <- imap_dfr(res, function(r2, ds) imap_dfr(r2, function(r, tx) {
    tot <- sum(r$perm$SumOfSqs[rownames(r$perm) %in% c("Model", "Residual")])
    m <- as.data.frame(r$db_margin) %>% rownames_to_column("term") %>% filter(term %in% names(pred_lab))
    tibble(dataset = ds, taxon = tx, term = m$term, pct = 100 * m$SumOfSqs / tot, p = m$`Pr(>F)`,
           total_r2 = 100 * r$perm$R2[rownames(r$perm) == "Model"])
  })) %>%
    mutate(term = factor(pred_lab[term], rev(pred_lab)), dataset = factor(ds_lab[dataset], ds_lab),
           taxon = factor(taxon, names(tx_pal)))
  r2_lab <- imp %>% distinct(dataset, taxon, total_r2) %>%
    group_by(dataset) %>%
    summarise(lab = paste0("Total explained: bats ", round(total_r2[taxon == "Bats"]), "%, moths ",
                           round(total_r2[taxon == "Moths"]), "%"), .groups = "drop")
  dodge <- position_dodge(0.6)
  pA <- ggplot(imp, aes(pct, term, color = taxon)) +
    geom_linerange(aes(xmin = 0, xmax = pct), linewidth = 0.6, position = dodge, alpha = 0.5) +
    geom_point(size = 3, position = dodge) +
    # asterisks just past each dot; where the bat and moth dots are close, aligned past the longer one
    geom_text(data = imp %>% group_by(dataset, term) %>%
                mutate(xl = if (diff(range(pct)) < 1.5) max(pct) + 0.25 else pct + 0.25) %>% ungroup(),
              aes(x = xl, label = stars(p), group = taxon), position = dodge, hjust = 0, vjust = 0.75, size = 3.6,
              color = "#222222", show.legend = FALSE) +
    geom_text(data = r2_lab, aes(x = Inf, y = -Inf, label = lab), inherit.aes = FALSE,
              hjust = 1.02, vjust = -0.6, size = 3.3, color = "#555555") +
    facet_wrap(~ dataset) +
    scale_color_manual(values = tx_pal, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18)), labels = function(x) paste0(x, "%")) +
    labs(title = "What drives community composition",
         x = "Community variation explained by the predictor alone", y = NULL) +
    bat_theme +
    theme(plot.title = element_text(size = 12), axis.text = element_text(size = 10),
          axis.title = element_text(size = 11), strip.text = element_text(size = 10.5, face = "bold"),
          strip.background = element_rect(fill = "#F2F2F2", color = NA),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(color = "#EBEBEB", linewidth = 0.4),
          legend.position = "top")

  # ── B: per-taxon responses (mvabund) ────────────────────────────────────────
  uni <- function(r) {
    up <- r$mv_anova$uni.p; ud <- r$mv_anova$uni.test
    list(p = up, dev = ud)
  }
  resp <- imap_dfr(res, function(r2, ds) imap_dfr(r2, function(r, tx) {
    cf <- r$coef; u <- uni(r); env <- r$env
    sds <- c(dist_km = sd(env$dist_km), mean_phase = sd(env$mean_phase),
             pct_nonforest = sd(env$pct_nonforest), brightness_dark = sd(env$brightness_dark))
    signed <- bind_rows(
      tibble(term = "color", taxon_name = colnames(cf), est = cf["colorW", ]),
      map_dfr(names(sds), ~ tibble(term = .x, taxon_name = colnames(cf), est = cf[.x, ] * sds[[.x]])))
    p_long <- as.data.frame(u$p) %>% rownames_to_column("term") %>%
      pivot_longer(-term, names_to = "taxon_name", values_to = "p")
    d_long <- as.data.frame(u$dev) %>% rownames_to_column("term") %>%
      pivot_longer(-term, names_to = "taxon_name", values_to = "dev")
    full_join(signed, p_long, by = c("term", "taxon_name")) %>%
      left_join(d_long, by = c("term", "taxon_name")) %>%
      filter(term %in% names(pred_lab)) %>%
      mutate(dataset = ds, group = tx)
  })) %>%
    mutate(label = ifelse(group == "Bats", paste0(taxon_name, " (", bat_names[taxon_name], ")"), taxon_name),
           dataset = factor(ds_lab[dataset], ds_lab), group = factor(group, names(tx_pal)))
  # row order: by total abundance-independent alphabetical within group, bats first
  resp <- resp %>% mutate(label = factor(label, rev(sort(unique(label)))))

  lim <- log(4)   # colour scale capped at x4 / x0.25
  sgn <- resp %>% filter(term != "intensity") %>%
    mutate(term = factor(pred_lab[term], pred_lab[names(pred_lab) != "intensity"]),
           fillv = pmax(pmin(est, lim), -lim), sig = !is.na(p) & p < 0.05)
  pB1 <- ggplot(sgn, aes(term, label, fill = fillv)) +
    geom_tile(aes(alpha = sig), color = "white", linewidth = 0.8) +
    scale_alpha_manual(values = c(`FALSE` = 0.25, `TRUE` = 1), guide = "none") +
    geom_text(aes(label = stars(p)), size = 3.6, vjust = 0.75, color = "#111111") +
    facet_grid(group ~ dataset, scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#B35806", mid = "white", high = "#542788", midpoint = 0,
                         limits = c(-lim, lim), na.value = "grey90",
                         breaks = log(c(0.25, 0.5, 1, 2, 4)),
                         labels = c("-75%", "-50%", "0", "+100%", "+300%"),
                         name = "% change per SD\n(white vs red for color)") +
    labs(title = "Which taxa respond", x = NULL, y = NULL) +
    bat_theme +
    theme(plot.title = element_text(size = 12), axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 9.5, angle = 30, hjust = 1),
          strip.text = element_text(size = 10.5, face = "bold"),
          strip.background = element_rect(fill = "#F2F2F2", color = NA),
          axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          strip.text.y = element_blank(),
          legend.position = "right", legend.title = element_text(size = 9))

  int <- resp %>% filter(term == "intensity") %>%
    mutate(col = factor(ifelse(as.integer(dataset) == 1, "Full data", "Far removed"), c("Full data", "Far removed")),
           sig = !is.na(p) & p < 0.05)
  pB2 <- ggplot(int, aes(col, label, fill = dev)) +
    geom_tile(aes(alpha = sig), color = "white", linewidth = 0.8) +
    scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1), guide = "none") +
    geom_text(aes(label = stars(p)), size = 3.6, vjust = 0.75, color = "#111111") +
    facet_grid(group ~ "Light intensity", scales = "free_y", space = "free_y") +
    scale_fill_gradient(low = "#F2F2F2", high = "#4D4D4D", name = "Intensity:\ndeviance\nexplained") +
    labs(x = NULL, y = NULL) +
    bat_theme +
    theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 9.5, angle = 30, hjust = 1),
          strip.text = element_text(size = 10.5, face = "bold"),
          strip.background = element_rect(fill = "#F2F2F2", color = NA),
          axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          legend.position = "right", legend.title = element_text(size = 9))

  plate <- (pA / (pB1 + pB2 + plot_layout(widths = c(5, 1.1)))) +
    plot_layout(heights = c(1, 2.2)) +
    plot_annotation(
      tag_levels = list(c("A", "B", "")),
      title = "Bat and moth community composition",
      subtitle = paste0("Samples = site × color × intensity (detections per night); bats = 9 species, moths = 9 families. ",
                        "A: marginal dbRDA (Bray–Curtis), asterisks from permutation tests.\n",
                        "B: mvabund negative binomial multivariate GLM; % change per 1 SD of each predictor (white vs red for color), ",
                        "asterisks from per-taxon tests adjusted across taxa;\nnon-significant cells faded (rare moth families give large but unreliable estimates). ",
                        "* p < 0.05, ** p < 0.01, *** p < 0.001."),
      theme = theme(plot.title = element_text(face = "bold", size = 15),
                    plot.subtitle = element_text(size = 9.5, color = "#555555"))) &
    theme(plot.tag = element_text(face = "bold", size = 14))

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plate, width = 13, height = 14, dpi = 300, bg = "white")
  list(plate = plate, importance = imp, responses = resp)
}
