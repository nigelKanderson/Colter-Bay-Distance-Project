# =============================================================================
# 40_site_map.R
# Map of the 14 monitoring stations over NLCD 2021 land cover, coloured by
# distance from the Colter Bay parking lot. Uses only local data (no basemap
# download). Output: output/figures/site_map.png
# =============================================================================
suppressMessages({ library(tidyverse); library(sf); library(terra); library(ggrepel) })
source("R/02_clean_data.R")   # provides the `sites` table

# scenario: "full" (all 14 sites) or "nofar" (drop the three ~7 km sites)
.args    <- commandArgs(trailingOnly = TRUE)
SCENARIO <- if (length(.args) >= 1) .args[1] else "full"
DROP     <- if (SCENARIO == "nofar") c("GRTE13","GRTE14","GRTE15") else character(0)
SFX      <- if (SCENARIO == "nofar") "_nofar" else ""
if (length(DROP)) sites <- dplyr::filter(sites, !site %in% DROP)

ref_ll <- data.frame(lon = -110.645134, lat = 43.904879)  # Colter Bay parking lot (0 km)

# ---- distance from the lot (UTM 12N, metres -> km) --------------------------
ref_utm   <- st_transform(st_as_sf(ref_ll, coords = c("lon","lat"), crs = 4326), 32612)
sites_utm <- st_transform(st_as_sf(sites, coords = c("lon","lat"), crs = 4326, remove = FALSE), 32612)
sites$dist_km <- as.numeric(st_distance(sites_utm, ref_utm)) / 1000

# ---- NLCD land cover, cropped to the study area ----------------------------
r <- rast("data/Annual_NLCD_LndCov_2021_CU_C1V1.tif")
sites_sf <- st_transform(st_as_sf(sites, coords = c("lon","lat"), crs = 4326, remove = FALSE), crs(r))
ref_sf   <- st_transform(st_as_sf(ref_ll, coords = c("lon","lat"), crs = 4326), crs(r))

bb <- st_bbox(sites_sf); pad <- 1200
rc <- crop(r, ext(bb["xmin"]-pad, bb["xmax"]+pad, bb["ymin"]-pad, bb["ymax"]+pad))
try(levels(rc) <- NULL, silent = TRUE)   # strip category table -> raw integer codes

rdf <- as.data.frame(rc, xy = TRUE); names(rdf)[3] <- "code"
rdf$cover <- factor(dplyr::case_when(
  rdf$code == 11            ~ "Water",
  rdf$code %in% 21:24       ~ "Developed",
  rdf$code %in% 41:43       ~ "Forest",
  rdf$code %in% c(51,52)    ~ "Shrub",
  rdf$code %in% 71:74       ~ "Grassland",
  rdf$code %in% c(81,82)    ~ "Cropland",
  rdf$code %in% c(90,95)    ~ "Wetland",
  TRUE                      ~ "Other"),
  levels = c("Forest","Shrub","Grassland","Wetland","Water","Developed","Cropland","Other"))
pal <- c(Forest="#4E7A51", Shrub="#C7B77A", Grassland="#E9E2C6", Wetland="#79B0A6",
         Water="#A6CEE3", Developed="#BDBDBD", Cropland="#E8D8A0", Other="#F0F0F0")

sc <- st_coordinates(sites_sf); sites$X <- sc[,1]; sites$Y <- sc[,2]
rc_ref <- st_coordinates(ref_sf)

# scale bar (1 km) near bottom-left
xr <- range(rdf$x); yr <- range(rdf$y)
sbx <- xr[1] + 0.06*diff(xr); sby <- yr[1] + 0.05*diff(yr)

p <- ggplot() +
  geom_raster(data = rdf, aes(x, y, fill = cover)) +
  scale_fill_manual(values = pal, name = "Land cover", drop = TRUE) +
  geom_point(data = sites, aes(X, Y, color = dist_km), size = 3.4) +
  scale_color_viridis_c(option = "C", direction = -1, name = "Distance\nfrom lot (km)") +
  geom_point(data = data.frame(X = rc_ref[1], Y = rc_ref[2]), aes(X, Y),
             shape = 23, size = 4, fill = "white", color = "black", stroke = 1.1) +
  ggrepel::geom_text_repel(data = sites, aes(X, Y, label = site), size = 3,
                           color = "#111111", box.padding = 0.45, seed = 1, min.segment.length = 0) +
  annotate("text", x = rc_ref[1], y = rc_ref[2], label = "Colter Bay lot",
           hjust = -0.12, vjust = -1.1, size = 3, fontface = "italic") +
  annotate("segment", x = sbx, xend = sbx + 1000, y = sby, yend = sby, linewidth = 1.1) +
  annotate("text", x = sbx + 500, y = sby, label = "1 km", vjust = -0.6, size = 3) +
  annotate("text", x = xr[2]-0.03*diff(xr), y = yr[2]-0.04*diff(yr), label = "N", fontface = "bold", size = 5) +
  annotate("segment", x = xr[2]-0.03*diff(xr), xend = xr[2]-0.03*diff(xr),
           y = yr[2]-0.09*diff(yr), yend = yr[2]-0.055*diff(yr),
           arrow = arrow(length = unit(0.2,"cm"), type = "closed"), linewidth = 0.8) +
  coord_equal(expand = FALSE) +
  labs(title = paste0("Grand Teton ALAN study — ", nrow(sites), " monitoring stations",
                      if (length(DROP)) " (far sites removed)" else ""),
       subtitle = "NLCD 2021 land cover; points coloured by distance from the Colter Bay parking lot",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())

ggsave(sprintf("output/figures/site_map%s.png", SFX), p, width = 9, height = 9, dpi = 300, bg = "white")
message("saved output/figures/site_map", SFX, ".png  (", nrow(sites), " sites; distance 0-",
        round(max(sites$dist_km),2), " km)")
