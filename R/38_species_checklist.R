# =============================================================================
# 38_species_checklist.R
# Flag each SonoBat-accepted taxon against the documented greater-Yellowstone /
# Grand Teton bat fauna (13 species; NPS Yellowstone inventory + GRTE list), and
# export a stratified sample of decision-relevant recordings for hand-vetting.
#
# In:  output/sonobat_compare_wide/{per_species_totals.csv, matched_recordings.csv}
# Out: output/sonobat_compare_wide/species_checklist.csv
#      output/sonobat_compare_wide/vetting_sample.csv
# =============================================================================

suppressMessages({ library(tidyverse) })
dir <- "output/sonobat_compare_wide"

# --- authoritative regional pool (NPS: 13 spp in the ecosystem) --------------
common <- c(
  Mylu="Little brown myotis", Epfu="Big brown bat", Myev="Long-eared myotis",
  Myvo="Long-legged myotis", Coto="Townsend's big-eared bat", Myth="Fringed myotis",
  Laci="Hoary bat", Lano="Silver-haired bat", Euma="Spotted bat",
  Anpa="Pallid bat", Myca="California myotis", Myci="Western small-footed myotis",
  Myyu="Yuma myotis",
  Labo="Eastern red bat", Myse="Northern long-eared bat",
  HiF="High-freq guild label (not a species)")
in_region    <- c("Mylu","Epfu","Myev","Myvo","Coto","Myth","Laci","Lano",
                  "Euma","Anpa","Myca","Myci","Myyu")
out_of_range <- c("Labo","Myse")
non_species  <- c("HiF")

# --- inventory: what each version calls, flagged --------------------------
per <- read_csv(file.path(dir, "per_species_totals.csv"), show_col_types = FALSE)
inv <- per %>%
  transmute(code = id, common = common[id],
            v4_4_5, v30_1, ratio,
            presence = case_when(code %in% out_of_range ~ "NOT in region",
                                 code %in% non_species  ~ "not a species",
                                 code %in% in_region    ~ "in region",
                                 TRUE                   ~ "unknown code"),
            called_by = case_when(v4_4_5 > 0 & v30_1 > 0 ~ "both",
                                  v4_4_5 == 0            ~ "v30.1 only",
                                  TRUE                   ~ "v4.4.5 only")) %>%
  arrange(match(presence, c("NOT in region","not a species","in region")), desc(v30_1))

cat("=== Species inventory (shared nights), flagged vs regional fauna ===\n")
print(as.data.frame(inv), row.names = FALSE)
write_csv(inv, file.path(dir, "species_checklist.csv"))

cat("\n--- v30.1 additions that are OUT OF REGION (suspected over-calls) ---\n")
print(as.data.frame(inv %>% filter(presence == "NOT in region")), row.names = FALSE)
cat("\n--- in-region species v30.1 newly resolves (legit gains) ---\n")
print(as.data.frame(inv %>% filter(called_by == "v30.1 only", presence == "in region")),
      row.names = FALSE)
cat("\n--- in-region species v30.1 nearly DROPS vs v4.4.5 (possible misses) ---\n")
print(as.data.frame(inv %>% filter(presence == "in region", ratio <= 0.2)), row.names = FALSE)

# --- stratified vetting sample -----------------------------------------------
m <- read_csv(file.path(dir, "matched_recordings.csv"), show_col_types = FALSE) %>%
  mutate(site = str_extract(rec_key, "^[A-Za-z0-9]+"),
         date = str_extract(rec_key, "[0-9]{8}"))
N <- 20; set.seed(42)
grab <- function(df, cat, k = N) df %>% slice_sample(n = min(k, nrow(df))) %>%
  transmute(stratum = cat, rec_key, site, date, v4.4.5 = id_a, v30.1 = id_b)

sample <- bind_rows(
  grab(m %>% filter(id_b == "Labo", id_a != "Labo"), "v30-only Labo (out of region)"),
  grab(m %>% filter(id_b == "Myse", id_a != "Myse"), "v30-only Myse (out of region)"),
  grab(m %>% filter(id_a == "Euma", id_b != "Euma"), "Euma lost by v30 (v4.4.5 called it)"),
  grab(m %>% filter(id_a == "Myev", id_b != "Myev"), "Myev lost by v30 (v4.4.5 called it)"),
  grab(m %>% filter(id_a != "NoID", id_b != "NoID", id_a != id_b), "species disagreement"),
  grab(m %>% filter(id_a == "NoID", id_b %in% in_region, id_b != "Labo"),
       "v30-only in-region ID (v4.4.5 NoID)"))

write_csv(sample, file.path(dir, "vetting_sample.csv"))
cat(sprintf("\nWrote %d recordings to vet -> %s\n", nrow(sample),
            file.path(dir, "vetting_sample.csv")))
cat("Strata:\n"); print(sample %>% count(stratum) %>% as.data.frame(), row.names = FALSE)
