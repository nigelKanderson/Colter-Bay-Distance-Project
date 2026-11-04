# =============================================================================
# 57_white_light_heatmap.R
# Heat map of sky brightness with the parking-lot lights ON (white, 100%):
# brightness_white100 = mean DSU of the fisheye sky photos at each station
# (python/run_sky_brightness.py), interpolated across the stations around the lot.
# Interpolation: inverse-distance weighting (power 2) on log10 DSU, 20 m grid in
# UTM 12N, clipped to 300 m around the stations used. The far group (GRTE13-15,
# ~6.5 km south) is left out of the surface (nothing measured in between) and
# reported in the caption; GRTE09/GRTE10 have no white-light photos.
# Background: NLCD 2021 water / developed / other, greyed so the surface reads.
# Set HEAT_VAR <- "brightness_dark" before sourcing for the lights-off map (darkest lights-off photo per
# station, as used in the models), and HEAT_LIMS (log10 DSU) to share a color scale between maps.
# Output: output/figures/fig_white_light_heatmap.png (fig_dark_heatmap.png for lights off); leaves the plot
# as `p_white_heatmap`, plus `p_far_inset` (far stations GRTE13/14 as points on the same scale)
# =============================================================================
suppressMessages({ library(tidyverse); library(sf); library(terra); library(ggrepel) })
source("R/02_clean_data.R")          # `sites`
source("R/storyboard_sites.R")       # SB_FAR_DROP (GRTE13/14/15)

heat_var <- if (exists("HEAT_VAR")) HEAT_VAR else "brightness_white100"
heat_on  <- heat_var == "brightness_white100"
heat_lab <- if (heat_on) "white lights on" else "lights off"
# lights-on map: stations without white-light photos (GRTE09, GRTE10) take their lights-off value,
# since stations away from the lot barely change with the lights on (filled = TRUE, label marked *)
sky <- read_csv("data/sky_brightness_measurements.csv", show_col_types = FALSE) %>%
  transmute(site, val = .data[[heat_var]],
            filled = heat_on & is.na(val) & !is.na(brightness_dark),
            val = ifelse(filled, brightness_dark, val))
utm <- 32612
st <- sites %>% left_join(sky, by = "site") %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% st_transform(utm)
st <- bind_cols(st, as_tibble(st_coordinates(st)))
near <- filter(st, !site %in% SB_FAR_DROP)
obs  <- filter(near, !is.na(val))
far  <- filter(st, site %in% SB_FAR_DROP)
# the lit parking lot spans GRTE01, GRTE16 and GRTE17: outline their hull (60 m buffer)
lot_sites <- c("GRTE01", "GRTE16", "GRTE17")
lot  <- st_buffer(st_convex_hull(st_union(filter(st, site %in% lot_sites))), 60)
lot_lab <- as_tibble(st_coordinates(st_centroid(lot)))

# ---- IDW surface on log10 DSU ----------------------------------------------------
hull <- st_buffer(st_convex_hull(st_union(near)), 300)
bb   <- st_bbox(st_buffer(hull, 250))
grid <- expand_grid(X = seq(bb["xmin"], bb["xmax"], by = 20), Y = seq(bb["ymin"], bb["ymax"], by = 20))
d2   <- outer(grid$X, obs$X, "-")^2 + outer(grid$Y, obs$Y, "-")^2
w    <- 1 / pmax(d2, 1)                                    # power 2 on distance
grid$lg <- as.vector((w %*% log10(obs$val)) / rowSums(w))
inside  <- st_intersects(st_as_sf(grid, coords = c("X", "Y"), crs = utm), hull, sparse = FALSE)[, 1]
grid$lg[!inside] <- NA

# ---- background: NLCD (water / developed / other) ----------------------------------
nlcd <- rast("data/Annual_NLCD_LndCov_2021_CU_C1V1.tif")
e    <- st_bbox(st_transform(st_as_sfc(bb), crs(nlcd)))
rc   <- crop(nlcd, ext(e["xmin"] - 200, e["xmax"] + 200, e["ymin"] - 200, e["ymax"] + 200))
try(levels(rc) <- NULL, silent = TRUE)
rc   <- project(rc, paste0("EPSG:", utm), method = "near", res = 30)
bg   <- as.data.frame(rc, xy = TRUE); names(bg)[3] <- "code"
bg$fill <- case_when(bg$code == 11 ~ "#CFE0EC", bg$code %in% 21:24 ~ "#C9C9C9", TRUE ~ "#ECECEA")

# ---- plot ------------------------------------------------------------------------
brks <- c(3, 5, 10, 20, 40, 80)
lims <- if (exists("HEAT_LIMS")) HEAT_LIMS else log10(range(obs$val))
xr <- c(bb["xmin"], bb["xmax"]); yr <- c(bb["ymin"], bb["ymax"])
lab <- near %>% mutate(label = ifelse(is.na(val), paste0(site, "  no photo"),
                                      sprintf("%s  %.1f%s", site, val, ifelse(filled, "*", ""))))

p_white_heatmap <- ggplot() +
  geom_raster(data = bg, aes(x, y), fill = bg$fill) +
  geom_raster(data = filter(grid, !is.na(lg)), aes(X, Y, fill = lg), alpha = 0.85, interpolate = TRUE) +
  scale_fill_gradientn(colours = c("#FFF7D6", "#FEE391", "#FB9A29", "#CC4C02", "#7A2E05"),
                       limits = lims, breaks = log10(brks), labels = brks, oob = scales::squish,
                       name = paste0("Sky brightness,\n", heat_lab, "\n(DSU)"),
                       guide = guide_colourbar(barheight = unit(6, "cm"))) +
  geom_sf(data = lot, fill = NA, color = "black", linewidth = 0.8, linetype = "22") +
  annotate("text", x = lot_lab$X - 30, y = lot_lab$Y - 110, label = "Parking\nlot", size = 3.5, fontface = "bold",
           lineheight = 0.9) +
  geom_point(data = filter(near, is.na(val)), aes(X, Y), shape = 21, size = 3,
             fill = "white", color = "black", stroke = 0.8) +
  geom_point(data = obs, aes(X, Y), shape = 21, size = 3, fill = "black", color = "white", stroke = 0.6) +
  geom_label_repel(data = lab, aes(X, Y, label = label), size = 3, fill = alpha("white", 0.85), label.size = 0,
                   label.padding = unit(0.12, "lines"), box.padding = 0.5, min.segment.length = 0,
                   max.overlaps = Inf, seed = 1) +
  annotate("segment", x = xr[1] + 150, xend = xr[1] + 650, y = yr[1] + 120, yend = yr[1] + 120, linewidth = 1) +
  annotate("text", x = xr[1] + 400, y = yr[1] + 120, label = "500 m", vjust = -0.7, size = 3) +
  coord_sf(crs = utm, xlim = xr, ylim = yr, expand = FALSE, datum = NA) +
  labs(title = paste("Sky brightness around the Colter Bay lot,", if (heat_on) "white lights on (100%)" else "lights off"),
       subtitle = "Fisheye sky photos at each station, interpolated between stations (value after each name, DSU). Dashed outline: parking lot (GRTE01, GRTE16, GRTE17).",
       caption = paste0(
         "Inverse-distance-weighted interpolation of log DSU, clipped to 300 m around the stations; the surface is an estimate between measured points.\n",
         sprintf("Far stations (~6.5 km south, not interpolated): %s.  ",
                 paste(sprintf("%s %.1f", far$site, far$val), collapse = ", ")),
         if (any(is.na(near$val))) paste(paste(near$site[is.na(near$val)], collapse = " and "), "have no photos. ") else "",
         if (any(near$filled)) paste0("* ", paste(near$site[near$filled], collapse = " and "),
                                      ": no white-light photos; lights-off value used. ") else "",
         "Background: NLCD 2021 (blue = water, grey = developed)."),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(color = "#555555"),
        plot.caption = element_text(hjust = 0, color = "#666666", size = 8.5),
        panel.border = element_rect(fill = NA, color = "#888888"))

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
heat_file <- if (heat_on) "output/figures/fig_white_light_heatmap.png" else "output/figures/fig_dark_heatmap.png"
ggsave(heat_file, p_white_heatmap, width = 11, height = 9, dpi = 300, bg = "white")
message("saved ", heat_file, " (", nrow(obs), " stations in the surface)")

# ---- far stations inset (for figures that include them; GRTE15 dropped as in SB_FULL_DROP) ----
# points coloured on the same scale as the surface; NLCD background around the far group
far_in <- filter(far, !site %in% SB_FULL_DROP, !is.na(val))
fbb <- st_bbox(st_buffer(st_union(far_in), 450))
fe  <- st_bbox(st_transform(st_as_sfc(fbb), crs(nlcd)))
frc <- crop(nlcd, ext(fe["xmin"] - 200, fe["xmax"] + 200, fe["ymin"] - 200, fe["ymax"] + 200))
try(levels(frc) <- NULL, silent = TRUE)
frc <- project(frc, paste0("EPSG:", utm), method = "near", res = 30)
fbg <- as.data.frame(frc, xy = TRUE); names(fbg)[3] <- "code"
fbg$fill <- case_when(fbg$code == 11 ~ "#CFE0EC", fbg$code %in% 21:24 ~ "#C9C9C9", TRUE ~ "#ECECEA")

p_far_inset <- ggplot() +
  geom_raster(data = fbg, aes(x, y), fill = fbg$fill) +
  geom_point(data = far_in, aes(X, Y, fill = log10(val)), shape = 21, size = 4.5,
             color = "black", stroke = 0.8) +
  scale_fill_gradientn(colours = c("#FFF7D6", "#FEE391", "#FB9A29", "#CC4C02", "#7A2E05"),
                       limits = lims, oob = scales::squish, guide = "none") +
  geom_label_repel(data = far_in, aes(X, Y, label = sprintf("%s  %.1f", site, val)),
                   size = 3, fill = alpha("white", 0.85), label.size = 0, label.padding = unit(0.12, "lines"),
                   box.padding = 0.6, min.segment.length = 0, seed = 1) +
  coord_sf(crs = utm, xlim = fbb[c("xmin", "xmax")], ylim = fbb[c("ymin", "ymax")], expand = FALSE, datum = NA) +
  labs(title = "Far stations, ~6.5 km south", x = NULL, y = NULL) +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 9.5, hjust = 0.5, margin = margin(b = 3)),
        panel.border = element_rect(fill = NA, color = "#444444", linewidth = 0.6),
        plot.background = element_rect(fill = "white", color = NA))
