# 56_manager_summary_panel.R
# Four-panel summary plate of the main findings for park managers (far sites removed: the 11
# sites within 0-2.1 km of the Colter Bay lot, R/storyboard_sites.R). Plain-language panel
# titles state each finding; no subtitle or model terms on the page; a footnote gives the caveats.
#   A  bats per night along the distance gradient, under red vs white light (model prediction +/- 95% CI)
#   B  moths per trap night along the gradient, under red vs white light
#   C  bat distance effect (% change per km) at each light intensity
#   D  moth distance effect (% change per km) at each light intensity
#   E  bats under white vs red light (% difference) at each light intensity
#   F  moths under white vs red light (% difference) at each light intensity
# Models: each taxon's general model (bat `simple_model1`, moth `insect_general$model`) refit on
# the 11 sites with the color x distance term, exactly as in the storyboard bat vs moth figure
# (R/42). Filled points = differs from 0 (p < 0.05, Holm-adjusted across the 5 intensities);
# open points = not significant. Asterisks: * p < 0.05, ** p < 0.01, *** p < 0.001 (Holm in C-F; in A-B,
# on each line's end label, the test of that color's distance slope). The 100% level is flagged (shaded, †): its red nights fell on the July
# full moon and the two wettest nights of the study.
# Output: output/figures/fig_manager_summary.png; returns the plate, the six panels ($panels$A-F) and both refits.

source("R/storyboard_sites.R")   # SB_FAR_DROP, SB_INT_FILL

make_manager_summary <- function(bat_model, moth_model, site_night, insect_total_env,
                                 drop_sites = SB_FAR_DROP, night_re = TRUE, col_bat = "#0072B2",
                                 file = "output/figures/fig_manager_summary.png") {
  library(tidyverse); library(glmmTMB); library(emmeans); library(patchwork)
  drop <- function(df) droplevels(filter(df, !site %in% drop_sites))
  bat_dat <- drop(site_night); moth_dat <- drop(insect_total_env)
  # night random effect: light color/intensity were set per night, so the stations on one night share
  # the treatment (and moon, weather); (1 | date) makes nights, not site-nights, the treatment replicates
  add <- if (night_re) . ~ . + color:dist_km + (1 | date) else . ~ . + color:dist_km
  bat  <- update(bat_model,  add, data = bat_dat)
  moth <- update(moth_model, add, data = moth_dat)

  pct   <- function(x) (exp(x) - 1) * 100
  grid  <- seq(min(bat_dat$dist_km), max(bat_dat$dist_km), length.out = 60)
  sites <- distinct(bat_dat, site, dist_km)
  red   <- "#B2182B"; ink <- "#222222"; mute <- "#666666"

  theme_mgr <- theme_minimal(base_family = "Arial", base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold", size = 13.5, color = ink),
          plot.subtitle = element_text(size = 10.5, color = mute, margin = margin(b = 8)),
          axis.title = element_text(size = 11, color = "#333333"), axis.text = element_text(size = 10.5),
          legend.position = "none", plot.margin = margin(8, 14, 8, 8))
  rug <- geom_rug(data = sites, aes(x = dist_km), inherit.aes = FALSE, sides = "b",
                  color = "#888888", length = unit(2.2, "mm"))
  km_axis <- scale_x_continuous("Distance from the parking lot (km)", breaks = 0:2,
                                expand = expansion(mult = c(0.02, 0.2)))

  # ── shared helper for the distance panels (A, B): red vs white curves ──────
  stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
  lc <- c("Red light" = red, "White light" = "#333333")
  color_curves <- function(m, ytitle, title, subtitle) {
    d <- as.data.frame(emmeans(m, ~ color | dist_km, at = list(dist_km = grid), type = "response")) %>%
      mutate(light = factor(recode(as.character(color), R = "Red light", W = "White light"), names(lc)))
    sp <- as.data.frame(test(emtrends(m, ~ color, var = "dist_km")))
    ends <- d %>% group_by(light) %>% slice_tail(n = 1) %>% ungroup() %>%
      mutate(light_lab = paste0(light, " ", stars(sp$p.value[match(color, sp$color)])))
    # keep the two end labels at least 8% of the y-range apart (centered on their midpoint)
    gap <- 0.08 * max(d$asymp.UCL)
    if (abs(diff(ends$response)) < gap) {
      mid <- mean(ends$response); up <- which.max(ends$response)
      ends$label_y <- ifelse(seq_len(nrow(ends)) == up, mid + gap / 2, mid - gap / 2)
    } else ends$label_y <- ends$response
    ggplot(d, aes(dist_km, response, color = light, fill = light)) +
      geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), alpha = 0.10, color = NA) +
      geom_line(aes(linetype = light), linewidth = 1.3) +
      geom_text(data = ends, aes(y = label_y, label = light_lab), hjust = 0, nudge_x = 0.05, size = 4,
                fontface = "bold") +
      rug + km_axis +
      scale_color_manual(values = lc) + scale_fill_manual(values = lc) +
      scale_linetype_manual(values = c("Red light" = "solid", "White light" = "22")) +
      scale_y_continuous(ytitle, limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
      coord_cartesian(clip = "off") +                     # end labels may run past the panel edge
      labs(title = title, subtitle = subtitle) + theme_mgr + theme(plot.margin = margin(8, 30, 8, 8))
  }
  slopes <- function(m) { e <- as.data.frame(test(emtrends(m, ~ color, var = "dist_km"))); setNames(pct(e$dist_km.trend), e$color) }

  # ── A: bats along the gradient, red vs white ──────────────────────────────
  bt <- slopes(bat)
  pA <- color_curves(bat, "Bat detections per night",
                     "Bats are most active close to the lit lot, under red and white light",
                     if (round(bt[["R"]]) == round(bt[["W"]]))
                       sprintf("About %d%% fewer bat detections for every km away from the lot, under either color", round(-bt[["R"]]))
                     else sprintf("About %d%% fewer bat detections per km away from the lot under red, %d%% under white",
                                  round(-bt[["R"]]), round(-bt[["W"]])))

  # ── B: moths along the gradient, red vs white ─────────────────────────────
  mt <- slopes(moth)
  pB <- color_curves(moth, "Moths per trap night",
                     "White light depletes moths near the lot; red light does not",
                     sprintf("Under white light, moths rise %d%% per km away from the lot; under red, the gradient is flat",
                             round(mt[["W"]])))

  # ── shared helper for the intensity panels (C, D) ─────────────────────────
  int_panel <- function(d, ytitle, up, down, title, subtitle) {
    d <- d %>% mutate(intensity = factor(intensity, names(SB_INT_FILL)),
                      sig = p < 0.05, lab = paste0(intensity, "%", ifelse(intensity == "100", "†", "")))
    ggplot(d, aes(intensity, est)) +
      annotate("rect", xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf, fill = "#F0F0F0") +
      geom_hline(yintercept = 0, color = "#999999", linewidth = 0.5) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.14, linewidth = 0.6, color = "#444444") +
      geom_point(aes(fill = ifelse(sig, as.character(intensity), "ns")), shape = 21, size = 5,
                 stroke = 0.9, color = ink) +
      geom_text(aes(y = hi, label = paste0(sprintf("%+.0f%%", est), stars(p))), vjust = -0.7, size = 3.6, color = "#333333") +
      annotate("text", x = 0.45, y = Inf, label = up, hjust = 0, vjust = 1.3, size = 4.3, color = ink, fontface = "bold") +
      annotate("text", x = 0.45, y = -Inf, label = down, hjust = 0, vjust = -0.6, size = 4.3, color = ink, fontface = "bold") +
      scale_fill_manual(values = c(SB_INT_FILL, ns = "white")) +
      scale_x_discrete("Light intensity", labels = d$lab) +
      scale_y_continuous(ytitle, expand = expansion(mult = 0.16), labels = function(v) paste0(ifelse(v > 0, "+", ""), v, "%")) +
      labs(title = title, subtitle = subtitle) + theme_mgr
  }
  holm <- function(e) mutate(e, p = p.adjust(p.value, "holm"))

  dist_by_int <- function(m) as.data.frame(test(emtrends(m, ~ intensity, var = "dist_km"))) %>% holm() %>%
    mutate(est = pct(dist_km.trend), lo = pct(dist_km.trend - 1.96 * SE), hi = pct(dist_km.trend + 1.96 * SE))
  wvr_by_int <- function(m) as.data.frame(summary(contrast(emmeans(m, ~ color | intensity), list(WvR = c(-1, 1))), infer = TRUE)) %>%
    holm() %>% mutate(est = pct(estimate), lo = pct(asymp.LCL), hi = pct(asymp.UCL))

  # ── D: moth distance effect by intensity ──────────────────────────────────
  mc <- dist_by_int(moth)
  pD <- int_panel(mc, "Change in moths per km from the lot",
                  "↑ more moths farther from the lot", "↓ more moths near the lot",
                  "Dimming to 30% or less removed the moth gradient",
                  "Moths built up away from the lot only at 50% and 70% intensity")

  # ── C: bat distance effect by intensity ───────────────────────────────────
  bdist <- dist_by_int(bat)
  brng  <- range(round(bdist$est))
  pC <- int_panel(bdist, "Change in bats per km from the lot",
                  "↑ more bats farther from the lot", "↓ more bats near the lot",
                  "Dimming did not change how bats fall off with distance",
                  sprintf("%d to %d%% fewer bat detections per km away from the lot at every intensity", -brng[2], -brng[1]))

  # ── E: bats, white vs red, by intensity ───────────────────────────────────
  bd <- wvr_by_int(bat)
  mid <- filter(bd, intensity %in% c("30", "50", "70")); sigE <- filter(mid, p < 0.05)
  pE <- int_panel(bd, "Bat activity, white vs red light",
                  "↑ more bats under white", "↓ more bats under red",
                  if (nrow(sigE) == 3) "Bats were less active under white than red light"
                  else "Bats tended to be less active under white than red light",
                  sprintf("%.0f to %.0f%% fewer detections under white at 30–70%%%s", -max(mid$est), -min(mid$est),
                          if (nrow(sigE) %in% 1:2) paste0("; significant at ", paste0(sigE$intensity, "%", collapse = " and "))
                          else ""))

  # ── F: moths, white vs red, by intensity ──────────────────────────────────
  md <- wvr_by_int(moth)
  pF <- int_panel(md, "Moths per trap, white vs red light",
                  "↑ more moths under white", "↓ more moths under red",
                  "Moth catches did not differ between white and red light",
                  "No clear white vs. red difference at 10–70% intensity")

  note <- paste0(
    "Grand Teton National Park, Colter Bay, summer 2022; 11 acoustic detector and moth trap sites 0–2.1 km from the parking lot. ",
    "Lines and points are model predictions with 95% confidence intervals; ticks on the x-axis mark the sites.\n",
    "Filled points differ from zero; open points do not. * p < 0.05, ** p < 0.01, *** p < 0.001 (A, B: each line's change with distance; C–F: adjusted across the five intensities).\n",
    "Light color and intensity were set experimentally, night by night; ",
    "the distance pattern is observational, since sites farther from the lot also differ in habitat. ",
    "A night random effect accounts for conditions shared by all sites on a night (moon, weather).\n",
    "† 100% (shaded): the red-light nights at 100% included the July full moon and the two wettest nights of the study; ",
    "the moth white vs. red difference rests on four red-light trap nights.")

  plate <- (pA | pB) / (pC | pD) / (pE | pF) +
    plot_annotation(
      tag_levels = "A",
      title = "How the Colter Bay parking lot affects bats and moths at a distance",
      caption = note,
      theme = theme(plot.title = element_text(face = "bold", size = 20, family = "Arial"),
                    plot.caption = element_text(size = 9, color = mute, hjust = 0, lineheight = 1.15, family = "Arial"))) &
    theme(plot.tag = element_text(face = "bold", size = 16))

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plate, width = 15, height = 16.2, dpi = 300, bg = "white")
  invisible(list(plate = plate, panels = list(A = pA, B = pB, C = pC, D = pD, E = pE, F = pF), bat = bat, moth = moth))
}
