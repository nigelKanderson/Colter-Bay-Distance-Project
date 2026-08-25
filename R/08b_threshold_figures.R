library(ggplot2)
library(dplyr)
library(patchwork)
library(forcats)
library(emmeans)
library(ggrepel)

# ==============================================================================
# ── Setup — reuse shared theme and palette from 08_figures.R ──────────────────
# ==============================================================================

# data_30 is missing dist_km / dist_cat / pct_nonforest — join from data_env
message("data_env cols: ", paste(names(data_env), collapse = ", "))
message("data_30  cols: ", paste(names(data_30),  collapse = ", "))

site_meta <- data_env %>%
  select(site, any_of(c("dist_km", "dist_cat", "pct_nonforest"))) %>%
  distinct()

date_meta <- data_env %>%
  select(site, date, jd) %>%
  distinct()

d30 <- data_30 %>%
  select(-any_of(c("jd", "dist_km", "dist_cat", "pct_nonforest"))) %>%
  left_join(site_meta, by = "site") %>%
  left_join(date_meta, by = c("site", "date")) %>%
  mutate(
    intensity = factor(intensity, levels = c("10", "30", "50", "70", "100")),
    jd_c      = jd - mean(jd, na.rm = TRUE)
  )

message("d30 cols after join: ", paste(names(d30), collapse = ", "))

# Shared constants (defined in 08_figures.R; redefined here for standalone use)
bat_theme <- theme_classic(base_size = 12) +
  theme(
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "#EBEBEB", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    axis.line          = element_line(color = "#444444", linewidth = 0.5),
    axis.ticks         = element_line(color = "#444444"),
    axis.text          = element_text(color = "#333333", size = 10),
    axis.title         = element_text(color = "#222222", size = 11),
    plot.title         = element_text(color = "#111111", size = 13,
                                      face = "bold", margin = margin(b = 3)),
    plot.subtitle      = element_text(color = "#666666", size = 9,
                                      margin = margin(b = 8)),
    legend.background  = element_rect(fill = "white", color = NA),
    legend.key         = element_rect(fill = "white", color = NA),
    legend.text        = element_text(color = "#333333", size = 10),
    legend.title       = element_text(color = "#444444", size = 9, face = "bold"),
    strip.background   = element_rect(fill = "#F2F2F2", color = NA),
    strip.text         = element_text(color = "#222222", face = "bold", size = 10),
    plot.caption       = element_text(color = "#999999", size = 8,
                                      hjust = 0, margin = margin(t = 8)),
    plot.margin        = margin(12, 14, 10, 12)
  )

col_red    <- "#7A1C2E"
col_white  <- "#C9A84C"
col_navy   <- "#1B2A4A"
col_purple <- "#4A1259"
col_ivory  <- "#F5F0E3"

intensity_pal <- c("10"  = col_ivory,
                   "30"  = col_white,
                   "50"  = col_purple,
                   "70"  = col_red,
                   "100" = col_navy)

dist_order <- c("Close", "Medium", "Further", "Far")
# Inherit the qmd's 10-species focal set when sourced from the pipeline;
# fall back to the same list so the script still runs standalone.
if (!exists("focal_spp")) {
  focal_spp <- c("Epfu", "Laci", "Lano", "Myca", "Myci",
                 "Mylu", "Myvo", "Myyu", "Myth", "Myev")
}

dir.create("output/figures/threshold30", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# ── Fig T1: Light color effect ─────────────────────────────────────────────────
# ==============================================================================

color_summary <- d30 %>%
  group_by(color) %>%
  summarise(
    mean_occ = mean(occurrences, na.rm = TRUE),
    se       = sd(occurrences,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

fig_t1 <- ggplot(color_summary, aes(x = color, y = mean_occ, fill = color)) +
  geom_col(width = 0.5, color = NA) +
  geom_errorbar(aes(ymin = mean_occ - se, ymax = mean_occ + se),
                width = 0.12, linewidth = 0.9, color = "#333333") +
  scale_fill_manual(values = c(R = col_red, W = col_white), guide = "none") +
  scale_x_discrete(labels = c(R = "Red", W = "White")) +
  labs(
    title    = "Effect of Light Color on Bat Occurrences (30-sec threshold)",
    subtitle = "Mean ± SE",
    x        = "Light color",
    y        = "Mean occurrences per night"
  ) +
  bat_theme

# ==============================================================================
# ── Fig T2: Distance from light source ────────────────────────────────────────
# ==============================================================================

dist_summary <- d30 %>%
  group_by(site, dist_km) %>%
  summarise(
    mean_occ = mean(occurrences, na.rm = TRUE),
    se       = sd(occurrences,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  arrange(dist_km) %>%
  mutate(site = fct_reorder(site, dist_km))

fig_t2 <- ggplot(dist_summary, aes(x = site, y = mean_occ)) +
  geom_col(fill = col_ivory, color = col_navy, linewidth = 0.6, width = 0.7) +
  geom_errorbar(aes(ymin = mean_occ - se, ymax = mean_occ + se),
                width = 0.2, linewidth = 0.8, color = col_navy) +
  geom_text(aes(label = round(dist_km, 2), y = -0.4),
            size = 2.8, color = "#666666", vjust = 1) +
  labs(
    title    = "Bat Occurrences by Distance from Light Source (30-sec threshold)",
    subtitle = "Site-level means ± SE, ordered by increasing distance (km shown below bars)",
    x        = "Site",
    y        = "Mean occurrences per night"
  ) +
  bat_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ==============================================================================
# ── Fig T3: Color × Intensity interaction ─────────────────────────────────────
# ==============================================================================

intensity_summary <- d30 %>%
  group_by(intensity, color) %>%
  summarise(
    mean_occ = mean(occurrences, na.rm = TRUE),
    se       = sd(occurrences,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

fig_t3 <- ggplot(intensity_summary,
                 aes(x = intensity, y = mean_occ, color = color, group = color)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = mean_occ - se, ymax = mean_occ + se),
                width = 0.15, linewidth = 0.8,
                position = position_dodge(0.1)) +
  scale_color_manual(values = c(R = col_red, W = col_white),
                     labels = c(R = "Red", W = "White"),
                     name   = "Light color") +
  labs(
    title    = "Light Intensity Effect by Color (30-sec threshold)",
    subtitle = "Mean ± SE at each intensity level",
    x        = "Light intensity (%)",
    y        = "Mean occurrences per night"
  ) +
  bat_theme +
  theme(legend.position = c(0.88, 0.88))

# ==============================================================================
# ── Fig T4: Seasonal patterns by light color ──────────────────────────────────
# ==============================================================================

fig_t4 <- ggplot(d30, aes(x = jd, y = occurrences, color = color)) +
  geom_smooth(method = "loess", span = 0.4, se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = c(R = col_red, W = col_white),
                     labels = c(R = "Red", W = "White"),
                     name   = "Light color") +
  labs(
    title    = "Seasonal Bat Activity by Light Color (30-sec threshold)",
    subtitle = "LOESS smooth across Julian day",
    x        = "Julian day",
    y        = "Occurrences per night"
  ) +
  bat_theme +
  theme(legend.position = c(0.12, 0.88))

# ==============================================================================
# ── Fig T5: Species community heatmap by color ────────────────────────────────
# ==============================================================================

heat_data <- d30 %>%
  filter(species %in% focal_spp) %>%
  group_by(species, color) %>%
  summarise(mean_occ = mean(occurrences, na.rm = TRUE), .groups = "drop") %>%
  mutate(species = fct_reorder(species, mean_occ, .fun = sum))

fig_t5 <- ggplot(heat_data, aes(x = color, y = species, fill = mean_occ)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient(low = col_ivory, high = col_navy, name = "Mean\noccurrences") +
  scale_x_discrete(labels = c(R = "Red", W = "White")) +
  labs(
    title    = "Bat Community Response by Light Color (30-sec threshold)",
    subtitle = "Mean occurrences per species",
    x        = "Light color",
    y        = NULL
  ) +
  bat_theme +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

# ==============================================================================
# ── Fig T6: Moon phase effect by color ────────────────────────────────────────
# ==============================================================================

fig_t6 <- ggplot(d30 %>% filter(!is.na(mean_phase)),
                 aes(x = mean_phase, y = occurrences, color = color)) +
  geom_smooth(method = "loess", span = 0.8, se = FALSE, linewidth = 1.2) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("New", "1st Q", "Full", "3rd Q", "New")) +
  scale_color_manual(values = c(R = col_red, W = col_white),
                     labels = c(R = "Red", W = "White"),
                     name   = "Light color") +
  labs(
    title    = "Effect of Moon Phase on Bat Occurrences (30-sec threshold)",
    subtitle = "LOESS smooth by light color",
    x        = "Moon phase",
    y        = "Occurrences per night"
  ) +
  bat_theme +
  theme(legend.position = c(0.12, 0.88))

# ==============================================================================
# ── Fig T7: Occurrences vs continuous distance (site-level means) ──────────────
# ==============================================================================

site_means_t <- d30 %>%
  group_by(site, dist_km) %>%
  summarise(mean_occ = mean(occurrences, na.rm = TRUE), .groups = "drop")

fig_t7 <- ggplot(site_means_t, aes(x = dist_km, y = mean_occ)) +
  geom_smooth(method = "lm", se = TRUE,
              color = col_navy, fill = col_ivory, linewidth = 1) +
  geom_point(size = 3.5, color = col_red) +
  geom_text_repel(aes(label = site), size = 3, color = "#555555",
                  box.padding = 0.4) +
  labs(
    title    = "Bat Occurrences vs Distance from Light Source (30-sec threshold)",
    subtitle = "Site-level means with linear trend",
    x        = "Distance from light (km)",
    y        = "Mean occurrences per night"
  ) +
  bat_theme

message("Distance-occurrence Spearman correlation (30-sec threshold):")
print(cor.test(site_means_t$dist_km, site_means_t$mean_occ, method = "spearman"))

# ==============================================================================
# ── Fig T8: Species × Distance — community heatmap ────────────────────────────
# ==============================================================================

spp_dist_t <- d30 %>%
  filter(species %in% focal_spp) %>%
  mutate(dist_cat = factor(dist_cat, levels = dist_order)) %>%
  group_by(species, dist_cat) %>%
  summarise(mean_occ = mean(occurrences, na.rm = TRUE), .groups = "drop") %>%
  mutate(species = fct_reorder(species, mean_occ, .fun = sum))

fig_t8 <- ggplot(spp_dist_t, aes(x = dist_cat, y = species, fill = mean_occ)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient(low = col_ivory, high = col_navy, name = "Mean\noccurrences") +
  labs(
    title    = "Species Occurrences by Distance from Light Source (30-sec threshold)",
    subtitle = "Mean occurrences per species × distance category",
    x        = "Distance category",
    y        = NULL
  ) +
  bat_theme +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

# ==============================================================================
# ── Fig T9: Species × Site heatmap (z-scored within species) ──────────────────
# ==============================================================================

site_order_t <- d30 %>%
  distinct(site, dist_km) %>%
  arrange(dist_km) %>%
  pull(site)

spp_site_t <- d30 %>%
  filter(species %in% focal_spp) %>%
  group_by(species, site) %>%
  summarise(mean_occ = mean(occurrences, na.rm = TRUE), .groups = "drop") %>%
  group_by(species) %>%
  mutate(mean_occ_z = as.numeric(scale(mean_occ))) %>%
  ungroup() %>%
  mutate(
    site    = factor(site, levels = site_order_t),
    species = fct_reorder(species, mean_occ, .fun = sum)
  )

dist_labels_t <- d30 %>%
  distinct(site, dist_km) %>%
  arrange(dist_km) %>%
  mutate(site  = factor(site, levels = site_order_t),
         label = paste0(round(dist_km, 2), " km"))

fig_t9 <- ggplot(spp_site_t, aes(x = site, y = species, fill = mean_occ_z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = col_red, mid = col_ivory, high = col_navy,
                       midpoint = 0, name = "Relative\nactivity (z)") +
  scale_x_discrete(labels = setNames(dist_labels_t$label, dist_labels_t$site)) +
  labs(
    title    = "Species Activity by Site (30-sec threshold)",
    subtitle = "Z-scored within species · sites ordered by increasing distance from light source",
    x        = "Distance from light source",
    y        = NULL
  ) +
  bat_theme +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid      = element_blank(),
    axis.line       = element_blank(),
    axis.ticks      = element_blank(),
    legend.position = "right"
  )

# ==============================================================================
# ── Fig T10: Mean occurrences by intensity level ───────────────────────────────
# ==============================================================================

intensity_bar_t <- d30 %>%
  group_by(intensity) %>%
  summarise(
    mean_occ = mean(occurrences, na.rm = TRUE),
    se       = sd(occurrences,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

fig_t10 <- ggplot(intensity_bar_t, aes(x = intensity, y = mean_occ, fill = intensity)) +
  geom_col(width = 0.6, color = NA) +
  geom_errorbar(aes(ymin = mean_occ - se, ymax = mean_occ + se),
                width = 0.18, linewidth = 0.8, color = "#333333") +
  scale_fill_manual(values = intensity_pal, guide = "none") +
  labs(
    title    = "Bat Occurrences by Light Intensity (30-sec threshold)",
    subtitle = "Mean ± SE across all sites and colours",
    x        = "Light intensity (%)",
    y        = "Mean occurrences per night"
  ) +
  bat_theme

# ==============================================================================
# ── Fig T11: Occurrences by intensity × color per site ────────────────────────
# ==============================================================================

int_color_bar_t <- d30 %>%
  group_by(intensity, color, site, dist_km) %>%
  summarise(
    mean_occ = mean(occurrences, na.rm = TRUE),
    se       = sd(occurrences,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(site_label = paste0(site, "\n(", round(dist_km, 2), " km)"),
         site_label = fct_reorder(site_label, dist_km))

fig_t11 <- ggplot(int_color_bar_t,
                  aes(x = intensity, y = mean_occ, fill = color)) +
  geom_col(position = position_dodge(0.7), width = 0.6, color = NA) +
  geom_errorbar(aes(ymin = mean_occ - se, ymax = mean_occ + se),
                position = position_dodge(0.7),
                width = 0.18, linewidth = 0.7, color = "#333333") +
  scale_fill_manual(values = c(R = col_red, W = col_white),
                    labels  = c(R = "Red", W = "White"),
                    name    = "Light color") +
  facet_wrap(~ site_label, nrow = 2) +
  labs(
    title    = "Bat Occurrences by Intensity and Color at Each Site (30-sec threshold)",
    subtitle = "Mean ± SE · sites ordered by increasing distance from light source",
    x        = "Light intensity (%)",
    y        = "Mean occurrences per night"
  ) +
  bat_theme +
  theme(legend.position = "top", strip.text = element_text(size = 8))

# ==============================================================================
# ── Fig T12: Species × color × intensity (faceted by species) ─────────────────
# ==============================================================================

spp_int_color_t <- d30 %>%
  filter(species %in% focal_spp) %>%
  group_by(species, intensity, color) %>%
  summarise(
    mean_occ = mean(occurrences, na.rm = TRUE),
    se       = sd(occurrences,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

fig_t12 <- ggplot(spp_int_color_t,
                  aes(x = intensity, y = mean_occ, fill = color)) +
  geom_col(position = position_dodge(0.7), width = 0.6, color = NA) +
  geom_errorbar(aes(ymin = mean_occ - se, ymax = mean_occ + se),
                position = position_dodge(0.7),
                width = 0.18, linewidth = 0.6, color = "#333333") +
  scale_fill_manual(values = c(R = col_red, W = col_white),
                    labels  = c(R = "Red", W = "White"),
                    name    = "Light color") +
  facet_wrap(~ species, nrow = 2, scales = "free_y") +
  labs(
    title    = "Bat Occurrences by Light Color and Intensity per Species (30-sec threshold)",
    subtitle = "Mean ± SE · y-axes free within species",
    x        = "Light intensity (%)",
    y        = "Mean occurrences per night"
  ) +
  bat_theme +
  theme(legend.position = "top", strip.text = element_text(size = 9, face = "bold"))

# ==============================================================================
# ── Fig T13: Occurrences vs pct_nonforest ─────────────────────────────────────
# ==============================================================================

nonforest_means_t <- d30 %>%
  group_by(site, pct_nonforest) %>%
  summarise(mean_occ = mean(occurrences, na.rm = TRUE), .groups = "drop")

fig_t13 <- ggplot(nonforest_means_t, aes(x = pct_nonforest, y = mean_occ)) +
  geom_smooth(method = "lm", se = TRUE,
              color = col_navy, fill = col_ivory, linewidth = 1) +
  geom_point(size = 3.5, color = col_red) +
  geom_text_repel(aes(label = site), size = 3, color = "#555555",
                  box.padding = 0.4) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title    = "Bat Occurrences vs. Non-forest Cover (30-sec threshold)",
    subtitle = "Site-level means with linear trend",
    x        = "Non-forest cover within 50 m buffer (%)",
    y        = "Mean occurrences per night"
  ) +
  bat_theme

# ==============================================================================
# ── Fig T14: Predicted occurrences ~ distance by intensity (data_30_model) ────
# ==============================================================================

dist_grid_t <- seq(0, 2, by = 0.05)

pred_dist_t <- emmeans(
  data_30_model,
  ~ intensity | dist_km,
  at   = list(dist_km = dist_grid_t),
  type = "response"
) %>%
  as.data.frame() %>%
  rename(
    conf.low  = dplyr::any_of(c("asymp.LCL", "lower.CL", "lower.HPD"))[1],
    conf.high = dplyr::any_of(c("asymp.UCL", "upper.CL", "upper.HPD"))[1]
  ) %>%
  mutate(intensity = factor(intensity, levels = c("10", "30", "50", "70", "100")))

fig_t14 <- ggplot(pred_dist_t,
                  aes(x = dist_km, y = response,
                      color = intensity, fill = intensity)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = intensity_pal, name = "Intensity (%)") +
  scale_fill_manual(values  = intensity_pal, guide = "none") +
  scale_x_continuous(breaks = c(0, 0.3, 0.7, 1.0, 1.5, 2.0),
                     name   = "Distance from light (km)") +
  geom_vline(xintercept = c(0.3, 0.7, 2.0),
             linetype = "dashed", color = "#BBBBBB", linewidth = 0.4) +
  labs(
    title    = "Bat Occurrences vs. Distance from Light by Intensity (30-sec threshold)",
    subtitle = "Marginal predictions ± 95% CI from data_30_model",
    y        = "Predicted occurrences per night"
  ) +
  bat_theme +
  theme(legend.position = "right")

# ==============================================================================
# ── Fig T15: Species × color across distance gradient ─────────────────────────
# ==============================================================================

fig_t15_data <- d30 %>%
  filter(species %in% focal_spp) %>%
  group_by(species, color, site, dist_km) %>%
  summarise(
    mean_occ = mean(occurrences, na.rm = TRUE),
    se       = sd(occurrences,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

fig_t15 <- ggplot(fig_t15_data,
                  aes(x = dist_km, y = mean_occ,
                      color = color, fill = color, group = color)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.0, alpha = 0.15) +
  scale_color_manual(values = c(R = col_red, W = col_white),
                     labels  = c(R = "Red", W = "White"),
                     name    = "Light color") +
  scale_fill_manual(values  = c(R = col_red, W = col_white), guide = "none") +
  facet_wrap(~ species, scales = "free_y", nrow = 2) +
  labs(
    title    = "Bat Occurrences by Light Color Across the Distance Gradient (30-sec threshold)",
    subtitle = "Linear trend ± 95% CI · one observation per site",
    x        = "Distance from light source (km)",
    y        = "Mean occurrences per night"
  ) +
  bat_theme +
  theme(legend.position = "top", strip.text = element_text(size = 9, face = "bold"))

# ==============================================================================
# ── Fig T16: Species × intensity across distance gradient ─────────────────────
# ==============================================================================

fig_t16_data <- d30 %>%
  filter(species %in% focal_spp) %>%
  group_by(species, intensity, site, dist_km) %>%
  summarise(mean_occ = mean(occurrences, na.rm = TRUE), .groups = "drop")

fig_t16 <- ggplot(fig_t16_data,
                  aes(x = dist_km, y = mean_occ,
                      color = intensity, fill = intensity, group = intensity)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.0, alpha = 0.12) +
  scale_color_manual(values = intensity_pal, name = "Intensity (%)") +
  scale_fill_manual(values  = intensity_pal, guide = "none") +
  facet_wrap(~ species, scales = "free_y", nrow = 2) +
  labs(
    title    = "Bat Occurrences by Light Intensity Across the Distance Gradient (30-sec threshold)",
    subtitle = "Linear trend ± 95% CI · one observation per site",
    x        = "Distance from light source (km)",
    y        = "Mean occurrences per night"
  ) +
  bat_theme +
  theme(legend.position = "top", strip.text = element_text(size = 9, face = "bold"))

# ==============================================================================
# ── Fig T17: Mean occurrences by species ──────────────────────────────────────
# ==============================================================================

fig_t17_data <- d30 %>%
  filter(species %in% focal_spp) %>%
  group_by(species) %>%
  summarise(
    mean_occ = mean(occurrences, na.rm = TRUE),
    se       = sd(occurrences,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  ) %>%
  mutate(species = fct_reorder(species, mean_occ))

fig_t17 <- ggplot(fig_t17_data,
                  aes(x = mean_occ, y = species, fill = species)) +
  geom_col(width = 0.65, color = NA) +
  geom_errorbar(aes(xmin = mean_occ - se, xmax = mean_occ + se),
                height = 0.25, linewidth = 0.7, color = "#333333") +
  scale_fill_manual(values = colorRampPalette(c(col_ivory, col_navy))(length(focal_spp)),
                    guide  = "none") +
  labs(
    title    = "Mean Bat Occurrences by Species (30-sec threshold)",
    subtitle = "Mean ± SE across all sites, nights, and treatments",
    x        = "Mean occurrences per night",
    y        = NULL
  ) +
  bat_theme

# ==============================================================================
# ── Fig T18: Species × light color bar charts ─────────────────────────────────
# ==============================================================================

fig_t18_data <- d30 %>%
  filter(species %in% focal_spp) %>%
  group_by(species, color) %>%
  summarise(
    mean_occ = mean(occurrences, na.rm = TRUE),
    se       = sd(occurrences,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

fig_t18 <- ggplot(fig_t18_data, aes(x = color, y = mean_occ, fill = color)) +
  geom_col(width = 0.6, color = NA) +
  geom_errorbar(aes(ymin = mean_occ - se, ymax = mean_occ + se),
                width = 0.18, linewidth = 0.7, color = "#333333") +
  scale_fill_manual(values = c(R = col_red, W = col_white),
                    labels  = c(R = "Red", W = "White"),
                    name    = "Light color") +
  scale_x_discrete(labels = c(R = "Red", W = "White")) +
  facet_wrap(~ species, nrow = 2, scales = "free_y") +
  labs(
    title    = "Bat Occurrences by Light Color per Species (30-sec threshold)",
    subtitle = "Mean ± SE · y-axes free within species",
    x        = "Light color",
    y        = "Mean occurrences per night"
  ) +
  bat_theme +
  theme(legend.position = "none", strip.text = element_text(size = 9, face = "bold"))

# ==============================================================================
# ── Fig T19: Species × intensity bar charts ───────────────────────────────────
# ==============================================================================

fig_t19_data <- d30 %>%
  filter(species %in% focal_spp) %>%
  group_by(species, intensity) %>%
  summarise(
    mean_occ = mean(occurrences, na.rm = TRUE),
    se       = sd(occurrences,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

fig_t19 <- ggplot(fig_t19_data, aes(x = intensity, y = mean_occ, fill = intensity)) +
  geom_col(width = 0.65, color = NA) +
  geom_errorbar(aes(ymin = mean_occ - se, ymax = mean_occ + se),
                width = 0.18, linewidth = 0.7, color = "#333333") +
  scale_fill_manual(values = intensity_pal, guide = "none") +
  facet_wrap(~ species, nrow = 2, scales = "free_y") +
  labs(
    title    = "Bat Occurrences by Light Intensity per Species (30-sec threshold)",
    subtitle = "Mean ± SE · y-axes free within species",
    x        = "Light intensity (%)",
    y        = "Mean occurrences per night"
  ) +
  bat_theme +
  theme(legend.position = "none", strip.text = element_text(size = 9, face = "bold"))

# ==============================================================================
# ── Save all figures ──────────────────────────────────────────────────────────
# ==============================================================================

ggsave("output/figures/threshold30/figT1_color_effect.png",         fig_t1,  width = 5,  height = 5,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT2_distance_effect.png",      fig_t2,  width = 6,  height = 5,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT3_color_intensity.png",      fig_t3,  width = 7,  height = 5,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT4_seasonal.png",             fig_t4,  width = 9,  height = 5,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT5_community_heatmap.png",    fig_t5,  width = 5,  height = 5.5, dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT6_moon_phase.png",           fig_t6,  width = 7,  height = 5,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT7_distance_gradient.png",    fig_t7,  width = 7,  height = 5.5, dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT8_species_distance.png",     fig_t8,  width = 9,  height = 6,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT9_species_site.png",         fig_t9,  width = 11, height = 6,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT10_intensity_bar.png",       fig_t10, width = 6,  height = 5,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT11_intensity_color_bar.png", fig_t11, width = 16, height = 8,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT12_species_int_color.png",   fig_t12, width = 14, height = 7,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT13_nonforest.png",           fig_t13, width = 7,  height = 5,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT14_dist_intensity_pred.png", fig_t14, width = 10, height = 6,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT15_species_color_dist.png",  fig_t15, width = 14, height = 8,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT16_species_int_dist.png",    fig_t16, width = 14, height = 8,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT17_species_bar.png",         fig_t17, width = 7,  height = 5,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT18_species_color_bar.png",   fig_t18, width = 12, height = 7,   dpi = 300, bg = "white")
ggsave("output/figures/threshold30/figT19_species_intensity_bar.png", fig_t19, width = 12, height = 7, dpi = 300, bg = "white")

message("All threshold-30 figures saved to output/figures/threshold30/")
