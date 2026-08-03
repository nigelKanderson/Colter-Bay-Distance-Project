# =============================================================================
# 36_sonobat_compare_figs.R
# Better visualizations of a SonoBat version comparison than the single
# confusion heatmap. Reads the CSVs written by compare_sonobat_versions()
# (matched_recordings.csv, per_species_totals.csv) and builds:
#   (A) species-only confusion heatmap  -> where the two versions disagree
#        when BOTH commit to a species (NoID dropped so it's readable)
#   (B) commitment breakdown            -> the NoID asymmetry (who calls more)
#   (C) per-species detections dumbbell -> log-scale magnitude per species
# Usage: make_compare_figs("output/sonobat_compare_wide", "v4.4.5", "v30.1")
# =============================================================================

suppressMessages({ library(tidyverse); library(patchwork) })

col_a_hex <- "#1B2A4A"  # navy  = version A (v4.4.5)
col_b_hex <- "#7A1C2E"  # red   = version B (v30.1)

fig_theme <- theme_classic(base_size = 14, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "#666666", size = 10),
        axis.title = element_text(size = 14), axis.text = element_text(size = 12),
        legend.position = "top", legend.title = element_blank())

make_compare_figs <- function(dir, label_a = "v4.4.5", label_b = "v30.1") {
  san <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
  ca <- san(label_a); cb <- san(label_b)
  m       <- read_csv(file.path(dir, "matched_recordings.csv"), show_col_types = FALSE)
  per_spp <- read_csv(file.path(dir, "per_species_totals.csv"), show_col_types = FALSE)

  # ---- (A) species-only confusion (both committed) --------------------------
  both <- m %>% filter(id_a != "NoID", id_b != "NoID")
  spp_order <- both %>% count(id_a) %>% arrange(desc(n)) %>% pull(id_a)
  all_spp <- sort(union(both$id_a, both$id_b))
  confA <- both %>% count(id_a, id_b) %>%
    group_by(id_a) %>% mutate(row_pct = 100 * n / sum(n)) %>% ungroup() %>%
    mutate(id_a = factor(id_a, levels = rev(all_spp)),
           id_b = factor(id_b, levels = all_spp))
  agree_pct <- 100 * mean(both$id_a == both$id_b)
  figA <- ggplot(confA, aes(id_b, id_a, fill = row_pct)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n), size = 2.6) +
    scale_fill_gradient(low = "#F5F0E3", high = col_b_hex,
                        name = paste0("% of\n", label_a, " row"),
                        limits = c(0, 100)) +
    labs(title = "A. Species agreement (both versions committed)",
         subtitle = sprintf("n = %s recordings; %.0f%% on the diagonal. NoID excluded.",
                            format(nrow(both), big.mark = ","), agree_pct),
         x = paste0(label_b, " ID"), y = paste0(label_a, " ID")) +
    fig_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

  # ---- (B) commitment breakdown (the NoID asymmetry) ------------------------
  brk <- m %>% mutate(cat = case_when(
      id_a == "NoID" & id_b == "NoID" ~ "Both NoID",
      id_a != "NoID" & id_b == "NoID" ~ paste0(label_a, " only"),
      id_a == "NoID" & id_b != "NoID" ~ paste0(label_b, " only"),
      id_a == id_b                    ~ "Agree (species)",
      TRUE                            ~ "Disagree (species)")) %>%
    count(cat) %>%
    mutate(cat = factor(cat, levels = c("Both NoID", paste0(label_a, " only"),
                                        paste0(label_b, " only"),
                                        "Disagree (species)", "Agree (species)")))
  pal <- setNames(c("#BFBFBF", col_a_hex, col_b_hex, "#C98A3B", "#2E6B4F"),
                  levels(brk$cat))
  n_b_only <- brk$n[brk$cat == paste0(label_b, " only")]
  n_a_only <- brk$n[brk$cat == paste0(label_a, " only")]
  figB <- ggplot(brk, aes(n, "recordings", fill = cat)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = ifelse(n / sum(n) > 0.03, format(n, big.mark = ","), "")),
              position = position_stack(vjust = 0.5), color = "white", size = 3) +
    scale_fill_manual(values = pal) +
    scale_x_continuous(labels = scales::comma, expand = expansion(c(0, 0.02))) +
    labs(title = "B. What each version does on the same recordings",
         subtitle = sprintf("Of %s matched recordings, %s commits to a species on %s that %s leaves NoID (reverse: only %s)",
                            format(nrow(m), big.mark = ","), label_b,
                            format(n_b_only, big.mark = ","), label_a,
                            format(n_a_only, big.mark = ",")),
         x = "Matched recordings", y = NULL) +
    guides(fill = guide_legend(nrow = 1)) +
    fig_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  # ---- (C) per-species detections, log scale (dumbbell) ---------------------
  pl <- per_spp %>%
    transmute(id, a = .data[[ca]], b = .data[[cb]],
              ratio = .data[[cb]] / pmax(.data[[ca]], 1)) %>%
    mutate(id = fct_reorder(id, b))
  lab_at <- function(v) log10(v + 1)
  figC <- ggplot(pl, aes(y = id)) +
    geom_segment(aes(x = lab_at(a), xend = lab_at(b), yend = id),
                 color = "#BBBBBB", linewidth = 1) +
    geom_point(aes(x = lab_at(a), color = label_a), size = 3) +
    geom_point(aes(x = lab_at(b), color = label_b), size = 3) +
    geom_text(aes(x = pmax(lab_at(a), lab_at(b)),
                  label = ifelse(ratio >= 2, paste0(round(ratio), "×"),
                          ifelse(ratio <= 0.5, paste0(round(ratio, 2), "×"), ""))),
              hjust = -0.25, size = 3, color = "#333333") +
    scale_color_manual(values = setNames(c(col_a_hex, col_b_hex), c(label_a, label_b))) +
    scale_x_continuous(breaks = lab_at(c(0, 10, 100, 1000, 10000)),
                       labels = c("0", "10", "100", "1k", "10k"),
                       expand = expansion(c(0.02, 0.16))) +
    labs(title = "C. Confident detections per species (shared nights)",
         subtitle = "Log scale; label = v30.1 / v4.4.5 ratio. New species sit at the v4.4.5=0 floor.",
         x = "Detections (log)", y = NULL) +
    fig_theme

  fig <- (figA | figC) / figB + plot_layout(heights = c(2, 1))
  ggsave(file.path(dir, "comparison_figure.pdf"), fig, width = 14, height = 10,
         device = cairo_pdf)
  ggsave(file.path(dir, "comparison_figure.png"), fig, width = 14, height = 10,
         dpi = 200, bg = "white")
  cat("Wrote", file.path(dir, "comparison_figure.{pdf,png}"), "\n")
  invisible(fig)
}

# Run for the wide dataset (and the distance subset if present)
make_compare_figs("output/sonobat_compare_wide", "v4.4.5", "v30.1")
if (file.exists("output/sonobat_compare/matched_recordings.csv"))
  make_compare_figs("output/sonobat_compare", "v4.4.5", "v30.1")
