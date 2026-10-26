# 46_far_site_treatment_panels.R
# Bat response to light color and intensity at the far (~7 km) sites only,
# GRTE14 vs. GRTE13 + 15, for all bats and for each call-frequency guild.
#   Rows: all bats | low frequency | Myev | high frequency
#   Cols: GRTE14 | GRTE13 + 15 (separate models)
# Each panel: predicted detections per night (marginal means, response scale) at each
# intensity under red and white light, 95% CIs; asterisks = white vs red differs at that
# intensity (* p<0.05, ** p<0.01, *** p<0.001; unadjusted). Panel subtitles give the
# model's joint tests for color, intensity and color x intensity.
# Model per site set and response:
#   y ~ color x intensity [+ site, fixed, for 13 + 15]
# No season / moon terms: each color x intensity cell ran on only 1-4 blocks of nights, so
# at 1-2 sites those terms are confounded with the cells and adjusting for them extrapolates
# far beyond the data (e.g. GRTE14 red 50%: raw mean 184/night, season-adjusted 1,928).
# Estimates therefore equal the observed nightly means and include whatever differed between
# the treatment nights; each cell is labeled with its number of nights. No site random
# effect (1-2 sites); non-forest / brightness are constant within a site. All bats: nightly totals (NB); guilds: guild-nights (NB, zero-inflated where zero
# nights occur). Colors: storyboard crimson (red light), darkened gold (white light).
# Output: output/figures/fig_far_sites_bat_color_intensity.png

far_site_sets   <- list(g14 = "GRTE14", g1315 = c("GRTE13", "GRTE15"))
far_site_labels <- c(g14 = "GRTE14", g1315 = "GRTE13 + 15")

fit_far <- function(d, response, zi = FALSE) {
  library(glmmTMB)
  d   <- droplevels(d)
  rhs <- "color * intensity"
  if (n_distinct(d$site) > 1) rhs <- paste(rhs, "+ site")
  zi <- zi && any(d[[response]] == 0)   # zero-inflation only where zero nights occur
  glmmTMB(as.formula(paste(response, "~", rhs)), data = d,
          ziformula = if (zi) ~1 else ~0, family = nbinom2())
}

make_far_bat_color_intensity <- function(site_night, cluster_night, bat_theme,
                                         col_red = "#7A1C2E", col_white = "#A8842A",
                                         file = "output/figures/fig_far_sites_bat_color_intensity.png") {
  library(emmeans)
  library(tidyverse)
  library(patchwork)

  responses <- c("All bats" = "all", "Low frequency" = "LowFreq", "Myev" = "Myev",
                 "High frequency" = "HighFreq")
  int_levels <- c("10", "30", "50", "70", "100")
  stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
  fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p))
  ok    <- function(m) isTRUE(m$sdr$pdHess) && all(is.finite(sqrt(diag(vcov(m)$cond))))

  fits <- map(far_site_sets, function(s) imap(responses, function(r, nm) {
    if (r == "all") fit_far(filter(site_night, site %in% s), "total_det")
    else fit_far(filter(cluster_night, site %in% s, cluster == r), "detections", zi = TRUE)
  }))

  results <- imap(fits, function(ms, ds) imap(ms, function(m, nm) {
    if (!ok(m)) return(NULL)
    em  <- emmeans(m, ~ color | intensity, type = "response")
    est <- as.data.frame(em) %>%
      transmute(color = recode(color, R = "Red", W = "White"),
                intensity = factor(intensity, int_levels),
                y = response, lo = asymp.LCL, hi = asymp.UCL)
    ct  <- as.data.frame(contrast(em, "revpairwise", adjust = "none")) %>%
      transmute(intensity = factor(intensity, int_levels), p = p.value)
    jt  <- as.data.frame(joint_tests(m))
    jp  <- function(term) jt$p.value[jt$`model term` == term]
    n   <- model.frame(m) %>% count(color, intensity) %>%
      transmute(color = recode(color, R = "Red", W = "White"), intensity = factor(intensity, int_levels), n)
    list(est = est %>% left_join(n, by = c("color", "intensity")), ct = ct, joint = jt,
         sub = sprintf("Color p %s · Intensity p %s · Color × intensity p %s",
                       fmt_p(jp("color")), fmt_p(jp("intensity")), fmt_p(jp("color:intensity"))))
  }))

  dodge <- position_dodge(0.45)
  panel <- function(ds, nm) {
    r   <- results[[ds]][[nm]]
    ttl <- paste(nm, "·", far_site_labels[[ds]])
    if (is.null(r)) {
      return(ggplot() + annotate("text", x = 0, y = 0, label = "Model did not converge", color = "#666666") +
               labs(title = ttl) + theme_void() + theme(plot.title = element_text(size = 11, face = "bold")))
    }
    lab <- r$ct %>% left_join(r$est %>% group_by(intensity) %>% summarise(top = max(hi)), by = "intensity")
    ggplot(r$est, aes(intensity, y, color = color, group = color)) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.6, position = dodge) +
      geom_line(linewidth = 0.5, alpha = 0.5, position = dodge) +
      geom_point(size = 2.8, position = dodge) +
      geom_text(aes(y = 0, label = n), position = dodge, vjust = -0.4, size = 3,
                color = "#666666", show.legend = FALSE) +
      geom_text(data = lab, aes(intensity, top, label = stars(p)), inherit.aes = FALSE,
                vjust = -0.3, size = 4.2, color = "#222222") +
      scale_color_manual(values = c(Red = col_red, White = col_white), name = "Light color") +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
      labs(title = ttl, subtitle = r$sub, x = "Light intensity (%)", y = "Detections per night") +
      bat_theme +
      theme(plot.title = element_text(size = 11), plot.subtitle = element_text(size = 8.5),
            axis.text = element_text(size = 10), axis.title = element_text(size = 11))
  }

  panels <- unlist(map(names(responses), function(nm) map(names(far_site_sets), ~ panel(.x, nm))),
                   recursive = FALSE)
  plate <- wrap_plots(panels, ncol = 2, guides = "collect") +
    plot_annotation(
      tag_levels = "A",
      title = "Far sites: bat detections by light color and intensity",
      subtitle = paste0("GRTE14 and GRTE13 + 15 (~7 km from the lot) modeled separately; mean detections per night ",
                        "(color × intensity, + site for GRTE13 + 15), 95% CIs. Gray numbers: nights per cell.\n",
                        "No season/moon adjustment (confounded with the treatment nights at 1–2 sites), so differences include ",
                        "date effects as well as any light effect.\n",
                        "Low frequency = Epfu, Lano, Laci; High frequency = Mylu, Myvo, Myci, Myca, Myyu. ",
                        "Asterisks: white vs red differs at that intensity (* p < 0.05, ** p < 0.01, *** p < 0.001)."),
      theme = theme(plot.title = element_text(face = "bold", size = 15),
                    plot.subtitle = element_text(size = 10, color = "#555555"))) &
    theme(legend.position = "top", plot.tag = element_text(face = "bold", size = 14))

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plate, width = 13, height = 16, dpi = 300, bg = "white")

  list(plate = plate, results = results, models = fits,
       converged = map(fits, ~ map_lgl(.x, ok)))
}
