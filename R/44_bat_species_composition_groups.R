# 44_bat_species_composition_groups.R
# Bat species composition (% of single-species detections) for four site groups:
#   A  lot sites            GRTE01, 05, 16, 17
#   B  mid-gradient sites   GRTE02, 03, 06, 08, 09, 10, 11
#   C  GRTE14
#   D  far sites            GRTE13, 15
# Any other grouping can be passed as `groups` (named list of site vectors).
# Shares are pooled within a group (sum of detections per species / group total),
# so busier sites weigh more. Compound (multi-species) labels are excluded.
# Bars are colored by call-frequency guild (Okabe-Ito blue / pink + darkened gold,
# color-blind checked); species are grouped by guild on the axis as well.
# Output: output/figures/fig_bat_species_composition_groups.png

default_composition_groups <- list(
  "Lot sites (GRTE01, 05, 16, 17)"                      = c("GRTE01", "GRTE05", "GRTE16", "GRTE17"),
  "Mid-gradient sites (GRTE02, 03, 06, 08, 09, 10, 11)" = c("GRTE02", "GRTE03", "GRTE06", "GRTE08",
                                                           "GRTE09", "GRTE10", "GRTE11"),
  "GRTE14"                                              = "GRTE14",
  "Far sites (GRTE13, 15)"                              = c("GRTE13", "GRTE15")
)

make_bat_species_composition_groups <- function(data_env, bat_theme,
                                                groups = default_composition_groups,
                                                title  = "Bat species composition by site group",
                                                file   = "output/figures/fig_bat_species_composition_groups.png",
                                                height = 9) {
  library(tidyverse)
  library(patchwork)

  guild_of <- c(Epfu = "Low frequency", Lano = "Low frequency", Laci = "Low frequency",
                Myev = "Myev (medium)",
                Mylu = "High frequency", Myvo = "High frequency", Myci = "High frequency",
                Myca = "High frequency", Myyu = "High frequency")
  common <- c(Epfu = "big brown", Lano = "silver-haired", Laci = "hoary", Myev = "long-eared myotis",
              Mylu = "little brown", Myvo = "long-legged myotis", Myci = "small-footed myotis",
              Myca = "California myotis", Myyu = "Yuma myotis")
  guild_pal <- c("Low frequency" = "#0072B2", "Myev (medium)" = "#CC79A7", "High frequency" = "#9A7400")

  site_group <- enframe(groups, name = "group", value = "site") %>% unnest(site)

  comp <- data_env %>%
    mutate(site = as.character(site)) %>%
    filter(species %in% names(guild_of)) %>%
    inner_join(site_group, by = "site") %>%
    group_by(group, species) %>%
    summarise(det = sum(detections), .groups = "drop") %>%
    complete(group, species = names(guild_of), fill = list(det = 0)) %>%
    group_by(group) %>%
    mutate(n_det = sum(det), share = 100 * det / n_det) %>%
    ungroup() %>%
    mutate(guild = factor(guild_of[species], names(guild_pal)),
           label = paste0(species, " (", common[species], ")"))

  # species order: by guild, then by overall share (largest at top within guild)
  sp_order <- comp %>%
    group_by(species, label, guild) %>% summarise(s = sum(det), .groups = "drop") %>%
    arrange(desc(guild), s) %>% pull(label)
  comp <- comp %>% mutate(label = factor(label, sp_order))

  panel <- function(g) {
    d <- filter(comp, group == g)
    n_sites <- length(groups[[g]])
    ggplot(d, aes(share, label, fill = guild)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = case_when(det == 0 ~ "0%", share < 0.05 ~ "<0.1%",
                                      TRUE ~ sprintf("%.1f%%", share))),
                hjust = -0.15, size = 3.2, color = "#333333") +
      scale_fill_manual(values = guild_pal, name = "Call-frequency guild", drop = FALSE) +
      scale_x_continuous(limits = c(0, 105), breaks = seq(0, 100, 25),
                         labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0))) +
      labs(title = g,
           subtitle = sprintf("%s single-species detections · %d site%s",
                              format(d$n_det[1], big.mark = ","), n_sites, ifelse(n_sites > 1, "s", "")),
           x = "Share of detections", y = NULL) +
      bat_theme +
      theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 9.5),
            axis.text = element_text(size = 10), axis.title = element_text(size = 11),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(color = "#EBEBEB", linewidth = 0.4))
  }

  plate <- wrap_plots(lapply(names(groups), panel), ncol = 2, guides = "collect") +
    plot_annotation(
      tag_levels = "A",
      title = title,
      subtitle = "% of single-species detections, pooled within each group (compound multi-species calls excluded)",
      theme = theme(plot.title = element_text(face = "bold", size = 15),
                    plot.subtitle = element_text(size = 10, color = "#555555"))) &
    theme(legend.position = "top", plot.tag = element_text(face = "bold", size = 14))

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(file, plate, width = 13, height = height, dpi = 300, bg = "white")

  list(plate = plate, table = comp)
}
