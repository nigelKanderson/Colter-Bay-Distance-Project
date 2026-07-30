# =============================================================================
# 32_sonobat_version_compare.R
# Compare bat species IDs between two SonoBat classifier versions (v4.4.5 vs
# v30.1) on the same acoustic recordings.
#
# Input : data/sonobat_compare/v445/*.xlsx  and  data/sonobat_compare/v30/*.xlsx
# Output: console summary + output/sonobat_compare/*.csv + confusion heatmap
#
# Key details discovered in the data:
#   * The OUTPUT filename suffix encodes the classification (…-Mylu.wav,
#     …-HiF.wav, …-noID.wav), so it differs by version. Recordings are matched
#     by the timestamp key SITE_YYYYMMDD_HHMMSS (identical raw input to both).
#   * v30.1 writes compound labels for ambiguous calls (SppAccp "Anpa/Epfu/Coto",
#     Prob "0.36/0.35/..."). We take each version's *confident single-species
#     accept*: a lone species in SppAccp with (first) Prob >= 0.9 — matching the
#     main pipeline's rule; everything else is "NoID".
# =============================================================================

suppressMessages({ library(tidyverse); library(readxl) })
dir.create("output/sonobat_compare", showWarnings = FALSE, recursive = TRUE)

PROB_CUT <- 0.9

read_version <- function(dir, label) {
  files <- list.files(dir, pattern = "\\.xlsx$", full.names = TRUE)
  message("Reading ", length(files), " ", label, " files ...")
  map_dfr(files, function(f) {
    d <- tryCatch(suppressWarnings(read_excel(f, col_types = "text")),
                  error = function(e) NULL)
    if (is.null(d) || !all(c("Filename", "SppAccp", "Prob") %in% names(d))) {
      message("  skipped (missing cols): ", basename(f)); return(NULL)
    }
    tibble(
      rec_key = str_extract(as.character(d$Filename), "^[A-Za-z0-9]+_[0-9]{8}_[0-9]{6}"),
      site    = toupper(str_extract(as.character(d$Filename), "^[A-Za-z0-9]+")),
      spp_raw = as.character(d$SppAccp),
      prob1   = suppressWarnings(as.numeric(str_extract(as.character(d$Prob), "^[0-9.]+"))))
  }) %>%
    filter(!is.na(rec_key)) %>%
    mutate(
      date     = str_extract(rec_key, "[0-9]{8}"),
      night    = paste(site, date),          # a site-night
      compound = str_detect(replace_na(spp_raw, ""), "/"),
      # confident single-species accept, else NoID
      id = if_else(!is.na(spp_raw) & !compound & spp_raw != "NoID" &
                     !is.na(prob1) & prob1 >= PROB_CUT, spp_raw, "NoID")) %>%
    # one row per recording: prefer a confident ID, then highest prob
    arrange(id == "NoID", desc(prob1)) %>%
    distinct(rec_key, .keep_all = TRUE)
}

v445 <- read_version("data/sonobat_compare/v445", "v4.4.5")
v30  <- read_version("data/sonobat_compare/v30",  "v30.1")

cat("\n=== Coverage: recordings ===\n")
cat(sprintf("v4.4.5 recordings: %d | v30.1 recordings: %d\n", nrow(v445), nrow(v30)))
common <- intersect(v445$rec_key, v30$rec_key)
cat(sprintf("Recordings in BOTH (matched by timestamp): %d\n", length(common)))
cat(sprintf("Only in v4.4.5: %d | Only in v30.1: %d\n",
            length(setdiff(v445$rec_key, v30$rec_key)),
            length(setdiff(v30$rec_key, v445$rec_key))))

# ---- Coverage: site-nights ---------------------------------------------------
nights <- full_join(v445 %>% count(night, name = "rec_v445"),
                    v30  %>% count(night, name = "rec_v30"), by = "night") %>%
  mutate(across(c(rec_v445, rec_v30), ~ replace_na(., 0)),
         in_both = rec_v445 > 0 & rec_v30 > 0) %>%
  arrange(desc(in_both), night)
shared_nights <- nights$night[nights$in_both]

cat("\n=== Coverage: site-nights ===\n")
cat(sprintf("site-nights in v4.4.5: %d (%d sites) | in v30.1: %d (%d sites) | in BOTH: %d\n",
            sum(nights$rec_v445 > 0), n_distinct(v445$site),
            sum(nights$rec_v30 > 0),  n_distinct(v30$site), length(shared_nights)))
cat(sprintf("Only v4.4.5: %d | Only v30.1: %d\n",
            sum(nights$rec_v445 > 0 & nights$rec_v30 == 0),
            sum(nights$rec_v30 > 0 & nights$rec_v445 == 0)))
cat("=> v4.4.5 is a partial re-run; every v4.4.5 night is also in v30.1.\n")
write_csv(nights, "output/sonobat_compare/site_night_coverage.csv")

# On the SHARED site-nights, do the versions output the same recordings?
v445s <- v445 %>% filter(night %in% shared_nights)
v30s  <- v30  %>% filter(night %in% shared_nights)
cat(sprintf("\nOn the %d shared nights: recordings v4.4.5 = %d, v30.1 = %d (ratio %.2fx)\n",
            length(shared_nights), nrow(v445s), nrow(v30s), nrow(v30s) / nrow(v445s)))
cat("=> nearly 1:1, so the larger v30.1 volume above is EXTRA NIGHTS, not more per night.\n")

# All the ID comparisons below use recordings matched by timestamp, which are by
# construction on shared nights.

# ---- Matched per-recording comparison ---------------------------------------
m <- inner_join(v445 %>% select(rec_key, site, id_445 = id),
                v30  %>% select(rec_key, id_30 = id), by = "rec_key")

cat("\n=== Species-ID agreement on matched recordings ===\n")
cat(sprintf("Overall agreement (incl. NoID):        %.1f%%  (n = %d)\n",
            100 * mean(m$id_445 == m$id_30), nrow(m)))
either_id <- m %>% filter(id_445 != "NoID" | id_30 != "NoID")
cat(sprintf("Agreement where either has an ID:      %.1f%%  (n = %d)\n",
            100 * mean(either_id$id_445 == either_id$id_30), nrow(either_id)))
both_id <- m %>% filter(id_445 != "NoID", id_30 != "NoID")
cat(sprintf("Agreement where BOTH have a species:   %.1f%%  (n = %d)\n",
            100 * mean(both_id$id_445 == both_id$id_30), nrow(both_id)))
cat(sprintf("v4.4.5 ID but v30.1 NoID: %d | v30.1 ID but v4.4.5 NoID: %d\n",
            sum(m$id_445 != "NoID" & m$id_30 == "NoID"),
            sum(m$id_30 != "NoID" & m$id_445 == "NoID")))

# Cohen's kappa on species (matched, incl. NoID)
lvls <- sort(union(m$id_445, m$id_30))
kap <- tryCatch(psych::cohen.kappa(cbind(match(m$id_445, lvls), match(m$id_30, lvls)))$kappa,
                error = function(e) NA)
cat(sprintf("Cohen's kappa (incl. NoID): %.3f\n", kap))

cat("\n=== Confusion matrix: v4.4.5 (rows) vs v30.1 (cols) ===\n")
conf <- table(v4.4.5 = m$id_445, v30.1 = m$id_30)
print(conf)

# ---- Per-species detection totals (SHARED nights = fair comparison) ----------
per_spp <- full_join(
  v445s %>% filter(id != "NoID") %>% count(id, name = "v4_4_5"),
  v30s  %>% filter(id != "NoID") %>% count(id, name = "v30_1"),
  by = "id") %>%
  mutate(across(c(v4_4_5, v30_1), ~ replace_na(., 0)),
         ratio = round(v30_1 / pmax(v4_4_5, 1), 1)) %>%
  arrange(desc(v4_4_5 + v30_1))
cat("\n=== Per-species confident detections on SHARED site-nights ===\n")
print(as.data.frame(per_spp), row.names = FALSE)

# ---- Save + figure ----------------------------------------------------------
write_csv(m, "output/sonobat_compare/matched_recordings.csv")
write_csv(as.data.frame(conf), "output/sonobat_compare/confusion_matrix.csv")
write_csv(per_spp, "output/sonobat_compare/per_species_totals.csv")

conf_df <- as.data.frame(conf) %>%
  group_by(v4.4.5) %>% mutate(row_pct = 100 * Freq / sum(Freq)) %>% ungroup()
p <- ggplot(conf_df, aes(v30.1, v4.4.5, fill = row_pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(Freq > 0, Freq, "")), size = 2.7) +
  scale_fill_gradient(low = "#F5F0E3", high = "#7A1C2E", name = "% of\nv4.4.5 row") +
  labs(title = "SonoBat species-ID agreement: v4.4.5 vs v30.1",
       subtitle = paste0("Matched recordings; confident single-species accept (Prob >= ",
                         PROB_CUT, "). Diagonal = agreement."),
       x = "v30.1 accepted ID", y = "v4.4.5 accepted ID") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/sonobat_compare/confusion_heatmap.png", p,
       width = 8, height = 7, dpi = 150, bg = "white")

cat("\nDone. Outputs in output/sonobat_compare/\n")
