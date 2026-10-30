# 55_moth_family_composition_distance.R
# Exploratory figures (not in the storyboards): community composition at each site, ordered by
# distance from the Colter Bay parking lot. One pie per site (all nights pooled); pie AREA is
# proportional to individuals per sampling night, so sites sampled on more nights don't look bigger.
# Sites are evenly spaced in distance order (real km under each pie) because the 11 near sites span
# only 0-2.1 km and the 3 far sites sit at 7.2-7.5 km.
#   make_moth_family_composition(): moth families per trap night (drop = families to leave out,
#     e.g. "Noctuidae"); families with < 1% of the moths shown are pooled as "Other"
#   make_bat_species_composition(): the 9 focal bat species per detector night (drop = species codes
#     to leave out, e.g. "Mylu"), from the
#     Distance_Project.qmd pipeline objects `data_env_focal` (single-species labels) and `site_night`
#   pie_distance_plot(): the shared drawing code
# Output: output/figures/moth_family_composition_distance[_no_<family>].png,
#         output/figures/bat_species_composition_distance[_no_<species>].png

# comp: site, group (factor, legend order), n;  sites: site, km, nights
pie_distance_plot <- function(comp, sites, pal, title, subtitle, size_title, legend_title,
                              legend_italic = TRUE) {
  library(tidyverse)
  comp <- comp %>% group_by(site) %>% mutate(tot = sum(n), prop = n / tot) %>% ungroup() %>%
    inner_join(sites, by = "site") %>% arrange(km, site, group)

  # pie geometry: x = rank in distance (+ a gap before the far sites), radius ~ sqrt(rate)
  sites <- comp %>% distinct(site, km, nights, tot) %>% arrange(km) %>% mutate(rate = tot / nights)
  ref   <- unique(signif(max(sites$rate) * c(0.25, 0.5, 1), 1))   # size key
  rmax  <- max(c(sites$rate, ref))
  sites <- sites %>% mutate(x = row_number() + ifelse(km > 5, 0.8, 0), r = 0.47 * sqrt(rate / rmax))
  wedge <- function(x0, r, a0, a1) {
    a <- seq(a0, a1, length.out = max(2, ceiling(80 * (a1 - a0))))
    tibble(px = c(x0, x0 + r * sin(2 * pi * a), x0), py = c(0, r * cos(2 * pi * a), 0))
  }
  poly <- comp %>% left_join(select(sites, site, x, r), by = "site") %>%
    group_by(site) %>% mutate(a1 = cumsum(prop), a0 = lag(a1, default = 0)) %>% ungroup() %>%
    filter(prop > 0) %>% mutate(id = row_number()) %>%
    rowwise() %>% mutate(g = list(wedge(x, r, a0, a1))) %>% ungroup() %>% unnest(g)

  # size key: reference circles, bottom-aligned, labels at their tops
  keyx <- max(sites$x) + 1.6
  key  <- tibble(n = ref, r = 0.47 * sqrt(n / rmax), y = -0.47 + r) %>%
    mutate(g = map(r, ~ tibble(a = seq(0, 2 * pi, length.out = 120)))) %>% unnest(g) %>%
    mutate(px = keyx + r * sin(a), py = y + r * cos(a))
  keylab <- tibble(n = ref, r = 0.47 * sqrt(ref / rmax)) %>% mutate(py = -0.47 + 2 * r)   # circle tops

  gap <- mean(c(max(sites$x[sites$km < 5]), min(sites$x[sites$km > 5])))
  ggplot() +
    annotate("segment", x = min(sites$x) - 0.5, xend = gap - 0.35, y = -0.86, yend = -0.86,
             arrow = arrow(length = unit(2.5, "mm"), type = "closed"), color = "#555555", linewidth = 0.6) +
    annotate("segment", x = gap + 0.35, xend = max(sites$x) + 0.5, y = -0.86, yend = -0.86,
             arrow = arrow(length = unit(2.5, "mm"), type = "closed"), color = "#555555", linewidth = 0.6) +
    annotate("text", x = gap, y = -0.86, label = "//", size = 5, color = "#555555") +
    annotate("text", x = min(sites$x) - 0.5, y = -0.99, label = "Farther from the parking lot  →",
             hjust = 0, size = 3.6, color = "#555555", fontface = "italic") +
    geom_polygon(data = poly, aes(px, py, group = id, fill = group), color = "white", linewidth = 0.25) +
    geom_text(data = sites, aes(x, -0.52, label = site), size = 3.2, fontface = "bold", vjust = 1) +
    geom_text(data = sites, aes(x, -0.68, label = sprintf("%.2f km", km)), size = 3, vjust = 1, color = "#333333") +
    geom_text(data = sites, aes(x, 0.52, label = sprintf("%.0f / night\nn = %d\n%d nights", rate, tot, nights)),
              size = 2.7, vjust = 0, color = "#555555", lineheight = 0.9) +
    geom_path(data = key, aes(px, py, group = n), color = "#777777", linewidth = 0.4) +
    geom_segment(data = keylab, aes(x = keyx, xend = keyx + 0.62, y = py, yend = py), color = "#999999", linewidth = 0.3, linetype = "22") +
    geom_text(data = keylab, aes(keyx + 0.66, py, label = n), size = 2.8, hjust = 0, color = "#555555") +
    annotate("text", x = keyx, y = 0.62, label = size_title, size = 3, fontface = "bold", lineheight = 0.9, vjust = 0) +
    scale_fill_manual(values = pal, name = legend_title) +
    coord_fixed(clip = "off", xlim = c(min(sites$x) - 0.55, keyx + 0.9), ylim = c(-1.02, 0.78)) +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_family = "Arial") +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 9.5, color = "#555555", margin = margin(b = 10)),
          legend.position = "bottom", legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(size = 9.5, face = if (legend_italic) "italic" else "plain"),
          plot.margin = margin(10, 15, 10, 15), plot.background = element_rect(fill = "white", color = NA))
}

make_moth_family_composition <- function(drop = NULL, file = NULL) {
  library(tidyverse)
  source("R/13_import_insects.R")
  fam  <- suppressWarnings(import_insects("data/grte_distance_insectID.xlsx"))$family %>%
    filter(!Family %in% drop)
  site <- readRDS("data/insect_data_out.rds") %>%
    group_by(site) %>% summarise(km = first(dist_km), nights = n_distinct(date), .groups = "drop")

  keep <- fam %>% count(Family, wt = detections, sort = TRUE) %>% filter(n / sum(n) >= 0.01) %>% pull(Family)
  comp <- fam %>%
    mutate(Family = ifelse(Family %in% keep, Family, "Other")) %>%
    count(site, Family, wt = detections, name = "n") %>%
    transmute(site, group = factor(Family, c(keep, "Other")), n)
  other <- setdiff(unique(fam$Family), keep)

  pal <- c(Noctuidae = "#332288", Lasiocampidae = "#88CCEE", Erebidae = "#44AA99", Geometridae = "#DDCC77",
           Sphingidae = "#CC6677", Notodontidae = "#AA4499", Euteliidae = "#117733", Drepanidae = "#999933",
           Crambidae = "#882255", Other = "#BBBBBB")   # Paul Tol muted (color-blind safe); fixed per family
  pal <- pal[levels(comp$group)]

  no  <- if (length(drop)) paste0(" — ", paste(drop, collapse = ", "), " excluded") else ""
  sub <- paste0("One pie per site, all nights pooled; pie area proportional to moths identified per trap night",
                if (length(drop)) paste0(" (excluding ", paste(drop, collapse = ", "), ")") else "",
                ". Sites evenly spaced in order of distance (gap = 2.1 to 7.2 km).",
                if (length(other)) paste0("\nOther = ", paste(other, collapse = ", "), " (",
                                          sum(comp$n[comp$group == "Other"]), " moths).") else "")
  p <- pie_distance_plot(comp, site, pal,
                         paste0("Moth family composition with distance from the Colter Bay parking lot", no),
                         sub, "Moths per night", "Moth family")
  file <- file %||% paste0("output/figures/moth_family_composition_distance",
                           if (length(drop)) paste0("_no_", paste(tolower(drop), collapse = "_")) else "", ".png")
  ggsave(file, p, width = 16, height = 4.3, dpi = 300, bg = "white")
  invisible(list(plot = p, composition = comp))
}

make_bat_species_composition <- function(data_env_focal, site_night, drop = NULL, file = NULL) {
  library(tidyverse)
  spp <- c(Epfu = "big brown", Lano = "silver-haired", Laci = "hoary", Myev = "long-eared myotis",
           Mylu = "little brown", Myvo = "long-legged myotis", Myci = "small-footed myotis",
           Myca = "California myotis", Myyu = "Yuma myotis")   # ordered by guild: low freq, Myev, high freq
  lab <- paste0(names(spp), " (", spp, ")")
  comp <- data_env_focal %>% filter(!species %in% drop) %>% count(site, species, wt = detections, name = "n") %>%
    transmute(site = as.character(site), group = factor(lab[match(species, names(spp))], lab), n)
  site <- site_night %>% group_by(site = as.character(site)) %>%
    summarise(km = first(dist_km), nights = n_distinct(date), .groups = "drop")

  pal <- setNames(c("#882255", "#CC6677", "#AA4499",   # low frequency: reds/purples
                    "#999933",                         # Myev
                    "#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77"), lab)   # high frequency
  pal <- pal[!names(spp) %in% drop]   # colors stay fixed per species
  no  <- if (length(drop)) paste0(" — ", paste(lab[names(spp) %in% drop], collapse = ", "), " excluded") else ""
  sub <- paste0("One pie per detector site, all treatment nights pooled; pie area proportional to detections per night of the ", sum(!names(spp) %in% drop), " focal species ",
                "(single-species calls only; compound multi-species labels excluded).\n",
                "Sites evenly spaced in order of distance (gap = 2.1 to 7.2 km). Legend ordered by call-frequency guild: ",
                "low frequency (Epfu, Lano, Laci), Myev, high frequency (Myotis).")
  p <- pie_distance_plot(comp, site, pal, paste0("Bat species composition with distance from the Colter Bay parking lot", no),
                         sub, "Detections per night", "Bat species", legend_italic = FALSE)
  file <- file %||% paste0("output/figures/bat_species_composition_distance",
                           if (length(drop)) paste0("_no_", paste(tolower(drop), collapse = "_")) else "", ".png")
  ggsave(file, p, width = 16, height = 4.5, dpi = 300, bg = "white")
  invisible(list(plot = p, composition = comp))
}

# Stack the four pie figures (moths, moths without owlet moths, bats, bats without little brown bats)
# into one plate, tagged A-D, for the end of distance_storyboard_nofar_distance.qmd
stack_composition_figures <- function(files = c("output/figures/moth_family_composition_distance.png",
                                                "output/figures/moth_family_composition_distance_no_noctuidae.png",
                                                "output/figures/bat_species_composition_distance.png",
                                                "output/figures/bat_species_composition_distance_no_mylu.png"),
                                      file = "output/figures/composition_distance_stacked.png") {
  library(magick)
  imgs <- lapply(seq_along(files), function(i) {
    im  <- image_trim(image_read(files[i])) %>% image_border("white", "40x40")
    tag <- image_blank(130, image_info(im)$height, "white") %>%   # left margin holding the panel letter
      image_annotate(LETTERS[i], size = 90, weight = 700, gravity = "northwest", location = "+20+20")
    image_append(c(tag, im))
  })
  w <- max(sapply(imgs, function(im) image_info(im)$width))
  imgs <- lapply(imgs, function(im) image_extent(im, sprintf("%dx%d", w, image_info(im)$height), gravity = "west", color = "white"))
  out <- image_append(do.call(c, imgs), stack = TRUE)
  image_write(out, file, density = 300)
  invisible(file)
}
