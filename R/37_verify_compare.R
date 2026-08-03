# =============================================================================
# 37_verify_compare.R
# Independent verification that a SonoBat version comparison (30_vs_445.R /
# 32_...) was done correctly. Two kinds of checks:
#   reconcile(dir)  : arithmetic self-consistency of the SAVED csv outputs
#                     (no file reads) -- fast, run on any comparison dir.
#   spotcheck(...)  : RE-READS the raw xlsx with a SEPARATE minimal parser and
#                     confirms, per sampled recording, that the assigned ID
#                     matches the raw SppAccp/Prob. Emits a human-readable CSV
#                     you can hand-check against the source files.
# =============================================================================

suppressMessages({ library(tidyverse); library(readxl) })

# ---- (1+2) reconcile saved outputs: totals + matching integrity -------------
reconcile <- function(dir) {
  cat("\n================ reconcile:", dir, "================\n")
  m    <- read_csv(file.path(dir, "matched_recordings.csv"), show_col_types = FALSE)
  conf <- read_csv(file.path(dir, "confusion_matrix.csv"),   show_col_types = FALSE)
  chk  <- function(name, ok) cat(sprintf("  [%s] %s\n", ifelse(ok, "PASS", "FAIL"), name))

  # matching integrity: one row per matched recording, no duplicate keys
  chk("matched rows are unique recordings (no duplicate rec_key)",
      n_distinct(m$rec_key) == nrow(m))
  chk("every rec_key matches the SITE_YYYYMMDD_HHMMSS format",
      all(str_detect(m$rec_key, "^[A-Za-z0-9]+_[0-9]{8}_[0-9]{6}$")))

  # confusion matrix must account for exactly the matched recordings
  chk("confusion-matrix cell counts sum to the matched total",
      sum(conf$Freq) == nrow(m))
  # rebuild the confusion table from matched and compare to the saved one
  rebuilt <- m %>% count(id_a, id_b, name = "Freq2")
  saved   <- conf %>% rename(id_a = a, id_b = b)
  cmp <- full_join(rebuilt, saved, by = c("id_a", "id_b")) %>%
    mutate(across(c(Freq2, Freq), ~replace_na(., 0)))
  chk("saved confusion matrix == recount of matched_recordings",
      all(cmp$Freq2 == cmp$Freq))

  # commitment breakdown categories partition the matched set
  brk <- m %>% mutate(cat = case_when(
      id_a == "NoID" & id_b == "NoID" ~ "both_noid",
      id_a != "NoID" & id_b == "NoID" ~ "a_only",
      id_a == "NoID" & id_b != "NoID" ~ "b_only",
      id_a == id_b                    ~ "agree_spp",
      TRUE                            ~ "disagree_spp")) %>% count(cat)
  chk("breakdown categories partition the matched set",
      sum(brk$n) == nrow(m))
  cat("   breakdown: ")
  cat(paste(brk$cat, brk$n, sep = "="), sep = "  "); cat("\n")

  # headline agreement numbers, recomputed straight from matched
  both <- m %>% filter(id_a != "NoID", id_b != "NoID")
  cat(sprintf("   agreement | overall %.1f%% (n=%d) | both-committed %.1f%% (n=%d)\n",
              100*mean(m$id_a == m$id_b), nrow(m),
              100*mean(both$id_a == both$id_b), nrow(both)))
  invisible(m)
}

# ---- minimal INDEPENDENT raw parser (deliberately not the pipeline's) -------
raw_read <- function(files, prob_cut = 0.9) {
  map_dfr(files, function(f) {
    d <- tryCatch(suppressMessages(suppressWarnings(read_excel(f, col_types = "text"))),
                  error = function(e) NULL)
    if (is.null(d) || !all(c("Filename","SppAccp","Prob") %in% names(d))) return(NULL)
    tibble(file = basename(f),
           filename = as.character(d$Filename),
           spp = as.character(d$SppAccp),
           prob = as.character(d$Prob))
  }) %>%
    mutate(rec_key = str_extract(filename, "^[A-Za-z0-9]+_[0-9]{8}_[0-9]{6}"),
           prob1 = suppressWarnings(as.numeric(str_extract(prob, "^[0-9.]+"))),
           compound = str_detect(replace_na(spp, ""), "/"),
           # the rule, re-implemented from scratch
           id_indep = if_else(!is.na(spp) & !compound & spp != "NoID" &
                                !is.na(prob1) & prob1 >= prob_cut, spp, "NoID")) %>%
    filter(!is.na(rec_key))
}

# ---- (3) spot-check assigned IDs against the raw source files ---------------
spotcheck <- function(dir, files_a, files_b, label_a = "v4.4.5", label_b = "v30.1",
                      n_per_cat = 3, seed = 1) {
  cat("\n================ spotcheck:", dir, "================\n")
  m  <- read_csv(file.path(dir, "matched_recordings.csv"), show_col_types = FALSE)
  ra <- raw_read(files_a); rb <- raw_read(files_b)

  # A) independent recompute of the winning id per rec_key, then compare to saved
  pick <- function(r) r %>% arrange(id_indep == "NoID", desc(prob1)) %>%
    distinct(rec_key, .keep_all = TRUE) %>% select(rec_key, id_indep)
  ia <- pick(ra); ib <- pick(rb)
  cmp <- m %>% select(rec_key, id_a, id_b) %>%
    inner_join(ia, by = "rec_key") %>% rename(id_a_indep = id_indep) %>%
    inner_join(ib, by = "rec_key") %>% rename(id_b_indep = id_indep)
  cat(sprintf("  matched keys re-derivable from raw: %d / %d\n", nrow(cmp), nrow(m)))
  cat(sprintf("  [%s] independent re-derivation of %s IDs matches saved (%d mismatches)\n",
              ifelse(all(cmp$id_a == cmp$id_a_indep), "PASS", "FAIL"), label_a,
              sum(cmp$id_a != cmp$id_a_indep)))
  cat(sprintf("  [%s] independent re-derivation of %s IDs matches saved (%d mismatches)\n",
              ifelse(all(cmp$id_b == cmp$id_b_indep), "PASS", "FAIL"), label_b,
              sum(cmp$id_b != cmp$id_b_indep)))

  # B) human-readable sample: raw rows from each version for a few recordings
  m2 <- m %>% mutate(cat = case_when(
      id_a == "NoID" & id_b == "NoID" ~ "both_noid",
      id_a != "NoID" & id_b == "NoID" ~ "a_only",
      id_a == "NoID" & id_b != "NoID" ~ "b_only",
      id_a == id_b                    ~ "agree_spp",
      TRUE                            ~ "disagree_spp"))
  set.seed(seed)
  samp <- m2 %>% group_by(cat) %>% slice_sample(n = n_per_cat) %>% ungroup()
  raw_of <- function(r, key) r %>% filter(rec_key == key) %>%
    summarise(spp = paste(spp, collapse = " | "), prob = paste(prob, collapse = " | ")) %>%
    unlist()
  detail <- samp %>% rowwise() %>% mutate(
      a_raw_spp = raw_of(ra, rec_key)["spp"], a_raw_prob = raw_of(ra, rec_key)["prob"],
      b_raw_spp = raw_of(rb, rec_key)["spp"], b_raw_prob = raw_of(rb, rec_key)["prob"]) %>%
    ungroup() %>%
    select(cat, rec_key, assigned_a = id_a, a_raw_spp, a_raw_prob,
           assigned_b = id_b, b_raw_spp, b_raw_prob)
  out <- file.path(dir, "verify_spotcheck.csv")
  write_csv(detail, out)
  cat("  wrote hand-checkable sample ->", out, "\n")
  print(as.data.frame(detail), row.names = FALSE)
  invisible(detail)
}

# --- run: reconcile both dirs (instant), spot-check the fast local subset ----
for (d in c("output/sonobat_compare_wide", "output/sonobat_compare"))
  if (file.exists(file.path(d, "matched_recordings.csv"))) reconcile(d)

spotcheck("output/sonobat_compare",
          list.files("data/sonobat_compare/v445", "\\.xlsx$", full.names = TRUE),
          list.files("data/sonobat_compare/v30",  "\\.xlsx$", full.names = TRUE))

# --- wide dataset raw re-derivation (re-reads the 176 grandteton_wide files) --
if (!identical(Sys.getenv("SKIP_WIDE"), "1")) {
  wb <- paste0("/Users/nanderson/Library/CloudStorage/",
               "GoogleDrive-nigel_anderson@brown.edu/.shortcut-targets-by-id/",
               "1QnPhsCxHagLxbcgwGcPhBMpzEksptTdA/barberlab_lightingprojects/grandteton_wide")
  spotcheck("output/sonobat_compare_wide",
            list.files(file.path(wb, "grte_wide_sonobatch_v4.4.5"), "\\.xlsx$",
                       full.names = TRUE, recursive = TRUE),
            list.files(file.path(wb, "grte_wide_sonobatch_v30.1"), "\\.xlsx$",
                       full.names = TRUE, recursive = TRUE))
}
