# =============================================================================
# sonobat_compare_fns.R
# Shared functions for comparing two SonoBat classifier versions (e.g. v4.4.5
# vs v30.1) on the same acoustic recordings. Used by:
#   - 32_sonobat_version_compare.R   (distance-project subset)
#   - 30_vs_445.R                    (grandteton_wide dataset)
#
# Method (validated in 32_...):
#   * Recordings matched by timestamp key SITE_YYYYMMDD_HHMMSS (the raw .wav is
#     identical input to both versions; the output filename suffix encodes the
#     classification and so differs by version).
#   * Each version's "confident single-species accept" = a lone species in
#     SppAccp with (first) Prob >= prob_cut; compound labels ("Anpa/Epfu/...")
#     and low-prob calls become "NoID".
# =============================================================================

suppressMessages({ library(tidyverse); library(readxl) })

# Read a vector of SonoBatch .xlsx files into one-row-per-recording form.
read_sonobat_version <- function(files, label, prob_cut = 0.9) {
  message("Reading ", length(files), " ", label, " files ...")
  map_dfr(files, function(f) {
    d <- tryCatch(suppressWarnings(read_excel(f, col_types = "text")),
                  error = function(e) NULL)
    if (is.null(d) || !all(c("Filename", "SppAccp", "Prob") %in% names(d))) {
      message("  skipped (missing cols / unreadable): ", basename(f)); return(NULL)
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
      night    = paste(site, date),
      compound = str_detect(replace_na(spp_raw, ""), "/"),
      id = if_else(!is.na(spp_raw) & !compound & spp_raw != "NoID" &
                     !is.na(prob1) & prob1 >= prob_cut, spp_raw, "NoID")) %>%
    arrange(id == "NoID", desc(prob1)) %>%
    distinct(rec_key, .keep_all = TRUE)
}

# Full comparison of two versions given their file vectors. Writes CSVs + a
# confusion heatmap to out_dir and prints a summary. Returns results invisibly.
compare_sonobat_versions <- function(files_a, files_b, out_dir,
                                     label_a = "v4.4.5", label_b = "v30.1",
                                     prob_cut = 0.9) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  san <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
  ca <- san(label_a); cb <- san(label_b)

  va <- read_sonobat_version(files_a, label_a, prob_cut)
  vb <- read_sonobat_version(files_b, label_b, prob_cut)

  cat("\n=== Coverage: recordings ===\n")
  cat(sprintf("%s recordings: %d | %s recordings: %d\n", label_a, nrow(va), label_b, nrow(vb)))
  cat(sprintf("Recordings in BOTH (matched by timestamp): %d\n",
              length(intersect(va$rec_key, vb$rec_key))))
  cat(sprintf("Only in %s: %d | Only in %s: %d\n", label_a,
              length(setdiff(va$rec_key, vb$rec_key)), label_b,
              length(setdiff(vb$rec_key, va$rec_key))))

  # ---- site-night coverage ----
  nights <- full_join(va %>% count(night, name = "rec_a"),
                      vb %>% count(night, name = "rec_b"), by = "night") %>%
    mutate(across(c(rec_a, rec_b), ~ replace_na(., 0)),
           in_both = rec_a > 0 & rec_b > 0) %>%
    arrange(desc(in_both), night)
  shared_nights <- nights$night[nights$in_both]
  cat("\n=== Coverage: site-nights ===\n")
  cat(sprintf("site-nights in %s: %d (%d sites) | in %s: %d (%d sites) | in BOTH: %d\n",
              label_a, sum(nights$rec_a > 0), n_distinct(va$site),
              label_b, sum(nights$rec_b > 0), n_distinct(vb$site), length(shared_nights)))
  cat(sprintf("Only %s: %d | Only %s: %d\n", label_a,
              sum(nights$rec_a > 0 & nights$rec_b == 0), label_b,
              sum(nights$rec_b > 0 & nights$rec_a == 0)))
  write_csv(nights, file.path(out_dir, "site_night_coverage.csv"))

  vas <- va %>% filter(night %in% shared_nights)
  vbs <- vb %>% filter(night %in% shared_nights)
  cat(sprintf("\nOn the %d shared nights: recordings %s = %d, %s = %d (ratio %.2fx)\n",
              length(shared_nights), label_a, nrow(vas), label_b, nrow(vbs),
              if (nrow(vas) > 0) nrow(vbs) / nrow(vas) else NA))

  # ---- matched per-recording comparison ----
  m <- inner_join(va %>% select(rec_key, site, id_a = id),
                  vb %>% select(rec_key, id_b = id), by = "rec_key")
  cat("\n=== Species-ID agreement on matched recordings ===\n")
  cat(sprintf("Overall agreement (incl. NoID):        %.1f%%  (n = %d)\n",
              100 * mean(m$id_a == m$id_b), nrow(m)))
  either_id <- m %>% filter(id_a != "NoID" | id_b != "NoID")
  cat(sprintf("Agreement where either has an ID:      %.1f%%  (n = %d)\n",
              100 * mean(either_id$id_a == either_id$id_b), nrow(either_id)))
  both_id <- m %>% filter(id_a != "NoID", id_b != "NoID")
  cat(sprintf("Agreement where BOTH have a species:   %.1f%%  (n = %d)\n",
              100 * mean(both_id$id_a == both_id$id_b), nrow(both_id)))
  cat(sprintf("%s ID but %s NoID: %d | %s ID but %s NoID: %d\n",
              label_a, label_b, sum(m$id_a != "NoID" & m$id_b == "NoID"),
              label_b, label_a, sum(m$id_b != "NoID" & m$id_a == "NoID")))
  lvls <- sort(union(m$id_a, m$id_b))
  kap <- tryCatch(psych::cohen.kappa(cbind(match(m$id_a, lvls), match(m$id_b, lvls)))$kappa,
                  error = function(e) NA)
  cat(sprintf("Cohen's kappa (incl. NoID): %.3f\n", kap))

  conf <- table(a = m$id_a, b = m$id_b)
  cat(sprintf("\n=== Confusion matrix: %s (rows) vs %s (cols) ===\n", label_a, label_b))
  print(conf)

  per_spp <- full_join(
    vas %>% filter(id != "NoID") %>% count(id, name = ca),
    vbs %>% filter(id != "NoID") %>% count(id, name = cb), by = "id") %>%
    mutate(across(c(!!ca, !!cb), ~ replace_na(., 0)),
           ratio = round(.data[[cb]] / pmax(.data[[ca]], 1), 1)) %>%
    arrange(desc(.data[[ca]] + .data[[cb]]))
  cat("\n=== Per-species confident detections on SHARED site-nights ===\n")
  print(as.data.frame(per_spp), row.names = FALSE)

  # ---- save + figure ----
  write_csv(m, file.path(out_dir, "matched_recordings.csv"))
  write_csv(as.data.frame(conf), file.path(out_dir, "confusion_matrix.csv"))
  write_csv(per_spp, file.path(out_dir, "per_species_totals.csv"))

  conf_df <- as.data.frame(conf) %>%
    group_by(a) %>% mutate(row_pct = 100 * Freq / sum(Freq)) %>% ungroup()
  p <- ggplot(conf_df, aes(b, a, fill = row_pct)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(Freq > 0, Freq, "")), size = 2.7) +
    scale_fill_gradient(low = "#F5F0E3", high = "#7A1C2E",
                        name = paste0("% of\n", label_a, " row")) +
    labs(title = paste0("SonoBat species-ID agreement: ", label_a, " vs ", label_b),
         subtitle = paste0("Matched recordings; confident single-species accept (Prob >= ",
                          prob_cut, "). Diagonal = agreement."),
         x = paste0(label_b, " accepted ID"), y = paste0(label_a, " accepted ID")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out_dir, "confusion_heatmap.png"), p,
         width = 8, height = 7, dpi = 150, bg = "white")

  cat("\nDone. Outputs in ", out_dir, "\n")
  invisible(list(a = va, b = vb, matched = m, confusion = conf,
                 per_species = per_spp, coverage = nights))
}
