# =============================================================================
# 35_moth_rarefaction.R
# Does moth SPECIES RICHNESS genuinely increase with distance, or is the
# observed increase just a by-product of catching more individuals?
# Uses coverage-based rarefaction/extrapolation (Hill q = 0; Chao & Jost) via
# iNEXT to standardize sampling effort.
#   (A) per-site richness: observed vs coverage-standardized, vs distance
#   (B) near / mid / far zone accumulation curves
# Output: output/figures/panel10_moth_rarefaction.{pdf,png} + console stats
# =============================================================================

suppressMessages({ library(tidyverse); library(iNEXT); library(patchwork) })
source("R/13_import_insects.R"); source("R/33_moth_taxonomy.R")

col_red <- "#7A1C2E"; col_navy <- "#1B2A4A"; col_purple <- "#4A1259"
panel_theme <- theme_classic(base_size = 15, base_family = "Arial") +
  theme(panel.grid.major.y = element_line(color = "#EBEBEB", linewidth = 0.4),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "#666666", size = 10),
        axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        legend.position = "top")

sp <- add_species_bin(import_insects("data/grte_distance_insectID.xlsx")$specimens) %>%
  filter(!is.na(Species_bin))
site_dist <- readRDS("data/insect_data_out.rds") %>% distinct(site, dist_km)

# ---- (A) per-site coverage-standardized richness ----------------------------
ab_list <- lapply(split(sp$Species_bin, sp$site), function(x) sort(as.numeric(table(x)), TRUE))
est <- estimateD(ab_list, q = 0, datatype = "abundance", base = "coverage") %>%
  as.data.frame()
names(est)[names(est) == "Assemblage"] <- "site"
site_std <- est %>% filter(Order.q == 0) %>%
  transmute(site, cov = round(SC, 3), std_rich = qD, lower = qD.LCL, upper = qD.UCL)

site_obs <- sp %>% group_by(site) %>%
  summarise(obs_rich = n_distinct(Species_bin), n_ind = n(), .groups = "drop")
site_tab <- site_dist %>% left_join(site_obs, by = "site") %>%
  left_join(site_std, by = "site")

cat("=== Per-site richness vs distance ===\n")
cat(sprintf("Standardized to common coverage SC = %.3f\n", unique(site_std$cov)[1]))
print(as.data.frame(site_tab %>% mutate(across(where(is.numeric), ~round(.,2)))), row.names = FALSE)
sp_obs  <- cor.test(site_tab$obs_rich,  site_tab$dist_km, method = "spearman")
sp_std  <- cor.test(site_tab$std_rich,  site_tab$dist_km, method = "spearman")
sp_ab   <- cor.test(site_tab$n_ind,     site_tab$dist_km, method = "spearman")
cat(sprintf("\nSpearman vs distance:  abundance rho=%.2f (p=%.3f) | observed richness rho=%.2f (p=%.3f) | STANDARDIZED richness rho=%.2f (p=%.3f)\n",
            sp_ab$estimate, sp_ab$p.value, sp_obs$estimate, sp_obs$p.value, sp_std$estimate, sp_std$p.value))

figA <- ggplot(site_tab, aes(dist_km)) +
  geom_smooth(aes(y = obs_rich), method = "lm", se = FALSE, color = col_purple, linetype = 2) +
  geom_point(aes(y = obs_rich), color = col_purple, size = 2.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, color = col_navy) +
  geom_smooth(aes(y = std_rich), method = "lm", se = FALSE, color = col_navy) +
  geom_point(aes(y = std_rich), color = col_navy, size = 2.5) +
  annotate("text", x = 4, y = max(site_tab$obs_rich), hjust = 0,
           label = sprintf("observed (dashed): rho=%.2f, p=%.3f\nstandardized (solid): rho=%.2f, p=%.2f",
                           sp_obs$estimate, sp_obs$p.value, sp_std$estimate, sp_std$p.value),
           size = 3.3, color = "#333333") +
  labs(title = "Per-site species richness vs. distance",
       subtitle = "Observed (purple, dashed) vs. coverage-standardized (navy, solid; ±95% CI)",
       x = "Distance from Colter Bay parking lot (km)", y = "Species richness") +
  panel_theme

# ---- (B) near / mid / far accumulation curves -------------------------------
zone <- sp %>% left_join(site_dist, by = "site") %>%
  mutate(zone = cut(dist_km, c(-Inf, 1, 3, Inf), labels = c("Near (<1 km)","Mid (1-3 km)","Far (>3 km)")))
zone_list <- lapply(split(zone$Species_bin, zone$zone), function(x) sort(as.numeric(table(x)), TRUE))
inx <- iNEXT(zone_list, q = 0, datatype = "abundance", endpoint = 2000)
curve <- fortify(inx, type = 1)
pt   <- curve %>% filter(Method == "Observed")
line <- curve %>% filter(Method != "Observed") %>%
  mutate(Method = factor(Method, c("Rarefaction","Extrapolation")))
zpal <- setNames(c(col_navy, col_purple, col_red),
                 c("Near (<1 km)","Mid (1-3 km)","Far (>3 km)"))
figB <- ggplot(line, aes(x, y, color = Assemblage)) +
  geom_ribbon(aes(ymin = y.lwr, ymax = y.upr, fill = Assemblage), alpha = 0.12, color = NA) +
  geom_line(aes(linetype = Method), linewidth = 1) +
  geom_point(data = pt, size = 2.5) +
  scale_color_manual(values = zpal, name = NULL) +
  scale_fill_manual(values = zpal, guide = "none") +
  scale_linetype_manual(values = c(Rarefaction = 1, Extrapolation = 2), name = NULL) +
  labs(title = "Species accumulation by distance zone",
       subtitle = "Coverage-based rarefaction (solid) / extrapolation (dashed) ± 95% CI",
       x = "Number of individuals", y = "Species richness") +
  panel_theme

p <- figA | figB
ggsave("output/figures/panels/panel10_moth_rarefaction.pdf", p, width = 14, height = 5.5, device = cairo_pdf)
ggsave("output/figures/panels/panel10_moth_rarefaction.png", p, width = 14, height = 5.5, dpi = 200,
       bg = "white", device = ragg::agg_png)

dir.create("output/moth_taxonomy", showWarnings = FALSE, recursive = TRUE)
write_csv(site_tab, "output/moth_taxonomy/rarefaction_per_site.csv")
cat("\nDone. Panel 10 written.\n")
