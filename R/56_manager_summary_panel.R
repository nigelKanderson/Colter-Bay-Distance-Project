# 56_manager_summary_panel.R
# Four-panel summary plate of the main findings for park managers (far sites removed: the 11
# sites within 0-2.1 km of the Colter Bay lot, R/storyboard_sites.R). Plain-language panel
# titles state each finding; no subtitle or model terms on the page; a footnote gives the caveats.
#   A  bats per night along the distance gradient (model prediction +/- 95% CI)
#   B  moths per trap night along the gradient, under red vs white light
#   C  moth distance effect (% change per km) at each light intensity
#   D  bats under white vs red light (% difference) at each light intensity
# Models: each taxon's general model (bat `simple_model1`, moth `insect_general$model`) refit on
# the 11 sites with the color x distance term, exactly as in the storyboard bat vs moth figure
# (R/42). Filled points = differs from 0 (p < 0.05, Holm-adjusted across the 5 intensities);
# open points = not significant. The 100% level is flagged (shaded, †): its red nights fell on the July
# full moon and the two wettest nights of the study.
# Output: output/figures/fig_manager_summary.png

source("R/storyboard_sites.R")   # SB_FAR_DROP, SB_INT_FILL

make_manager_summary <- function(bat_model, moth_model, site_night, insect_total_env,
                                 drop_sites = SB_FAR_DROP, col_bat = "#0072B2",
                                 file = "output/figures/fig_manager_summary.png") {
  library(tidyverse); library(glmmTMB); library(emmeans); library(patchwork)
  drop <- function(df) droplevels(filter(df, !site %in% drop_sites))
  bat_dat <- drop(site_night); moth_dat <- drop(insect_total_env)
  bat  <- update(bat_model,  . ~ . + color:dist_km, data = bat_dat)
  moth <- update(moth_model, . ~ . + color:dist_km, data = moth_dat)

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

  # ── A: bats along the gradient ────────────────────────────────────────────
  bA <- as.data.frame(emmeans(bat, ~ dist_km, at = list(dist_km = grid), type = "response"))
  bt <- as.data.frame(test(emtrends(bat, ~ 1, var = "dist_km")))
  pA <- ggplot(bA, aes(dist_km, response)) +
    geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), fill = col_bat, alpha = 0.15) +
    geom_line(color = col_bat, linewidth = 1.3) + rug + km_axis +
    scale_y_continuous("Bat detections per night", limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Bats are most active close to the lit lot",
         subtitle = sprintf("About %d%% fewer bat detections for every km away from the lot", round(-pct(bt$dist_km.trend)))) +
    theme_mgr

  # ── B: moths along the gradient, red vs white ─────────────────────────────
  mB <- as.data.frame(emmeans(moth, ~ color | dist_km, at = list(dist_km = grid), type = "response")) %>%
    mutate(light = factor(recode(as.character(color), R = "Red light", W = "White light"), c("Red light", "White light")))
  mt <- as.data.frame(test(emtrends(moth, ~ color, var = "dist_km")))
  ends <- mB %>% group_by(light) %>% slice_tail(n = 1) %>% ungroup()
  lc <- c("Red light" = red, "White light" = "#333333")
  pB <- ggplot(mB, aes(dist_km, response, color = light, fill = light)) +
    geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), alpha = 0.10, color = NA) +
    geom_line(aes(linetype = light), linewidth = 1.3) +
    geom_text(data = ends, aes(label = light), hjust = 0, nudge_x = 0.05, size = 4, fontface = "bold") +
    rug + km_axis +
    scale_color_manual(values = lc) + scale_fill_manual(values = lc) +
    scale_linetype_manual(values = c("Red light" = "solid", "White light" = "22")) +
    scale_y_continuous("Moths per trap night", limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(title = "White light depletes moths near the lot; red light does not",
         subtitle = sprintf("Under white light, moths rise %d%% per km away from the lot; under red, the gradient is flat",
                            round(pct(mt$dist_km.trend[mt$color == "W"])))) +
    theme_mgr

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
      geom_text(aes(y = hi, label = sprintf("%+.0f%%", est)), vjust = -0.7, size = 3.6, color = "#333333") +
      annotate("text", x = 0.45, y = Inf, label = up, hjust = 0, vjust = 1.3, size = 3.5, color = mute, fontface = "italic") +
      annotate("text", x = 0.45, y = -Inf, label = down, hjust = 0, vjust = -0.6, size = 3.5, color = mute, fontface = "italic") +
      scale_fill_manual(values = c(SB_INT_FILL, ns = "white")) +
      scale_x_discrete("Light intensity", labels = d$lab) +
      scale_y_continuous(ytitle, expand = expansion(mult = 0.16), labels = function(v) paste0(ifelse(v > 0, "+", ""), v, "%")) +
      labs(title = title, subtitle = subtitle) + theme_mgr
  }
  holm <- function(e) mutate(e, p = p.adjust(p.value, "holm"))

  # ── C: moth distance effect by intensity ──────────────────────────────────
  mc <- as.data.frame(test(emtrends(moth, ~ intensity, var = "dist_km"))) %>% holm() %>%
    mutate(est = pct(dist_km.trend), lo = pct(dist_km.trend - 1.96 * SE), hi = pct(dist_km.trend + 1.96 * SE))
  pC <- int_panel(mc, "Change in moths per km from the lot",
                  "↑ more moths farther from the lot", "↓ more moths near the lot",
                  "Dimming to 30% or less removed the moth gradient",
                  "Moths built up away from the lot only at 50% and 70% intensity")

  # ── D: bats, white vs red, by intensity ───────────────────────────────────
  bd <- as.data.frame(summary(contrast(emmeans(bat, ~ color | intensity), list(WvR = c(-1, 1))), infer = TRUE)) %>%
    holm() %>% mutate(est = pct(estimate), lo = pct(asymp.LCL), hi = pct(asymp.UCL))
  rng <- range(round(bd$est[bd$p < 0.05 & bd$est < 0]))
  pD <- int_panel(bd, "Bat activity, white vs red light",
                  "↑ more bats under white", "↓ more bats under red",
                  "Bats were less active under white than red light",
                  sprintf("%d to %d%% fewer bat detections under white at 30–70%% intensity", -rng[2], -rng[1]))

  note <- paste0(
    "Grand Teton National Park, Colter Bay, summer 2022; 11 acoustic detector and moth trap sites 0–2.1 km from the parking lot. ",
    "Lines and points are model predictions with 95% confidence intervals; ticks on the x-axis mark the sites.\n",
    "Filled points differ from zero; open points do not. Light color and intensity were set experimentally, night by night; ",
    "the distance pattern is observational, since sites farther from the lot also differ in habitat. ",
    "Results hold after accounting for overnight temperature and rain.\n",
    "† 100% (shaded): the white vs. red difference holds after accounting for moon, temperature and rain, but rests on four red-light nights ",
    "(two at the July full moon, the two wettest of the study) and reverses the pattern at all other intensities.")

  plate <- (pA | pB) / (pC | pD) +
    plot_annotation(
      tag_levels = "A",
      title = "How the Colter Bay parking lot affects bats and moths at a distance",
      caption = note,
      theme = theme(plot.title = element_text(face = "bold", size = 20, family = "Arial"),
                    plot.caption = element_text(size = 9, color = mute, hjust = 0, lineheight = 1.15, family = "Arial"))) &
    theme(plot.tag = element_text(face = "bold", size = 16))

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plate, width = 15, height = 11.2, dpi = 300, bg = "white")
  invisible(list(plate = plate, bat = bat, moth = moth))
}
