# =============================================================================
# 26_people_pca.R
# Principal components analysis (PCA) of the streetlight perception items,
# following the standard motivation-segmentation recipe:
#
#   1. Assess suitability: Bartlett's test of sphericity (p < 0.05) and the
#      Kaiser-Meyer-Olkin statistic (KMO > 0.50).
#   2. Extract components with eigenvalue > 1 (Kaiser criterion).
#   3. An item loads on a component if |loading| >= 0.40.
#   4. Cross-loaded items (>= 0.40 on more than one component) are eliminated;
#      items that reach 0.40 on no component are also dropped.
#   5. Reliability of each component's items: Cronbach's alpha (>= 0.65 = ok).
#   6. Retained items are averaged into a single component score.
#
# NOTE: 
# Negatively-worded items are reverse-keyed (6 - x on the 1-5 scale) within a
# component before alpha and averaging so all items point the same way.
# =============================================================================

library(tidyverse)
library(readxl)
library(psych)

set.seed(42)
dir.create("output", showWarnings = FALSE)
dir.create("output/figures", showWarnings = FALSE)

LOAD_CUT  <- 0.40   # minimum loading to belong to a component
ALPHA_CUT <- 0.65   # minimum acceptable Cronbach's alpha

# ---- Load & clean ------------------------------------------------------------
raw <- read_excel("people_data.xlsx")
dat <- raw %>% mutate(across(everything(), ~ na_if(., 999)))

attitude_items <- c(
  "Streetlights_activitiesmorepleasurable",
  "Streetlights_eyetransition",
  "Streetlights_activitiesdifficult",
  "Streetlights_lesssafe",
  "Streetlights_easiernavigate",
  "Streetlights_wildlifebehavior",
  "Streetlights_wildlifebenefits",
  "Streetlights_reducehumanimpacts",
  "Streetlights_bright",
  "Streetlights_dark",
  "Streetlights_pointsofinterest",
  "Streetlights_affectGYE",
  "Streetlights_unlitareas"
)

# Median-impute so all 166 respondents keep a component score
items <- dat[attitude_items] %>%
  mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))

# ---- 1. Suitability ----------------------------------------------------------
R <- cor(items)
kmo  <- KMO(items)
bart <- cortest.bartlett(R, n = nrow(items))
cat("=== PCA suitability ===\n")
cat(sprintf("KMO overall MSA = %.3f  (want > 0.50)\n", kmo$MSA))
cat(sprintf("Bartlett chi-sq = %.1f, df = %d, p = %s  (want < 0.05)\n",
            bart$chisq, bart$df, format.pval(bart$p.value)))

# ---- 2. Number of components: eigenvalue > 1 ---------------------------------
eig <- eigen(R)$values
n_comp <- sum(eig > 1)
cat("\n=== Eigenvalues ===\n")
print(round(eig, 3))
cat(sprintf("Components with eigenvalue > 1: %d\n", n_comp))

# ---- 3. PCA with varimax rotation --------------------------------------------
pca <- principal(items, nfactors = n_comp, rotate = "varimax")
L <- unclass(pca$loadings)
cat("\n=== Rotated component loadings (varimax) ===\n")
print(round(L, 3))

# ---- 4. Assign items / flag cross-loadings -----------------------------------
assign_tbl <- tibble(
  item        = rownames(L),
  n_over_cut  = rowSums(abs(L) >= LOAD_CUT),
  best_comp   = apply(L, 1, function(r) which.max(abs(r))),
  best_load   = apply(L, 1, function(r) r[which.max(abs(r))])
) %>%
  mutate(status = case_when(
    n_over_cut == 0 ~ "dropped (loads < 0.40 everywhere)",
    n_over_cut  > 1 ~ "dropped (cross-loaded)",
    TRUE            ~ paste0("PC", best_comp)
  ))

cat("\n=== Item assignment (|loading| >= 0.40, no cross-loading) ===\n")
print(as.data.frame(assign_tbl %>%
                      mutate(best_load = round(best_load, 3))),
      row.names = FALSE)

retained <- assign_tbl %>% filter(str_starts(status, "PC"))

# ---- 5 & 6. Cronbach's alpha + averaged component scores ----------------------
comp_summary <- list()
scores <- tibble(.rows = nrow(items))

for (k in sort(unique(retained$best_comp))) {
  comp_items <- retained %>% filter(best_comp == k)
  itms <- comp_items$item
  if (length(itms) < 2) {
    cat(sprintf("\nPC%d has only %d item -> no alpha, skipped as a scale.\n",
                k, length(itms)))
    next
  }
  # reverse-key negatively-loaded items (6 - x on the 1-5 scale)
  sub <- items[itms]
  negs <- comp_items$item[comp_items$best_load < 0]
  for (nm in negs) sub[[nm]] <- 6 - sub[[nm]]

  a <- psych::alpha(sub, warnings = FALSE)
  alpha_raw <- a$total$raw_alpha

  cat(sprintf("\n=== Component PC%d ===\n", k))
  cat("Items:", paste(itms, collapse = ", "), "\n")
  if (length(negs)) cat("Reverse-keyed:", paste(negs, collapse = ", "), "\n")
  cat(sprintf("Cronbach's alpha = %.3f  (%s)\n",
              alpha_raw, ifelse(alpha_raw >= ALPHA_CUT, "OK", "BELOW 0.65")))

  score_name <- paste0("PC", k, "_score")
  scores[[score_name]] <- rowMeans(sub)
  comp_summary[[paste0("PC", k)]] <- tibble(
    component = paste0("PC", k),
    n_items   = length(itms),
    items     = paste(itms, collapse = "; "),
    reverse_keyed = paste(negs, collapse = "; "),
    cronbach_alpha = round(alpha_raw, 3),
    reliable  = alpha_raw >= ALPHA_CUT
  )
}

comp_summary <- bind_rows(comp_summary)
cat("\n=== Component reliability summary ===\n")
print(as.data.frame(comp_summary), row.names = FALSE)

# ---- Relate component scores to light condition + sky acceptability ----------
out <- dat %>% bind_cols(scores)
cat("\n=== Mean component score by light condition (1=red, 2=white) ===\n")
print(out %>%
        filter(!is.na(StreetlightCondition)) %>%
        mutate(light = factor(StreetlightCondition, labels = c("red", "white"))) %>%
        group_by(light) %>%
        summarise(across(ends_with("_score"), ~ round(mean(.), 2)),
                  n = n(), .groups = "drop"))

cat("\n=== Correlation of component scores with Acceptable_Sky_Conditions ===\n")
for (nm in names(scores)) {
  ct <- cor.test(out[[nm]], out$Acceptable_Sky_Conditions)
  cat(sprintf("%-12s r = %+.2f, p = %s\n", nm, ct$estimate,
              format.pval(ct$p.value, digits = 2)))
}

# ---- Save outputs ------------------------------------------------------------
loadings_out <- as.data.frame(round(L, 3)) %>%
  rownames_to_column("item") %>%
  left_join(assign_tbl %>% select(item, status), by = "item")
write_csv(loadings_out, "output/people_pca_loadings.csv")
write_csv(comp_summary, "output/people_pca_component_reliability.csv")
write_csv(out %>% mutate(row_id = row_number()) %>%
            select(row_id, ends_with("_score"),
                   StreetlightCondition, Acceptable_Sky_Conditions),
          "output/people_pca_component_scores.csv")

# Scree plot
png("output/figures/people_pca_scree.png", width = 1000, height = 700, res = 130)
plot(eig, type = "b", pch = 19, xlab = "Component", ylab = "Eigenvalue",
     main = "PCA scree - streetlight items")
abline(h = 1, lty = 2, col = "red")
dev.off()

cat("\nDone. Outputs written to output/ and output/figures/\n")
