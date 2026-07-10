# =============================================================================
# 25_people_factor_cluster.R
# Factor-then-cluster segmentation of the streetlight perception survey.
#
# Data: people_data.xlsx  (166 respondents x 15 variables)
#   StreetlightCondition       1 = red light, 2 = white light (experimental)
#   Streetlights_* (13 items)  1-5 Likert, 1 = "not at all true"
#   Acceptable_Sky_Conditions  1-7 item
#   999 is a missing-value code in every column.
#
# Pipeline:
#   1. Load, recode 999 -> NA, median-impute the 13 attitude items.
#   2. Exploratory factor analysis (EFA) on the 13 items.
#   3. Extract factor scores per respondent.
#   4. Cluster respondents on the factor scores (k-means, k chosen by silhouette;
#      Ward hierarchical as a cross-check).
#   5. Profile clusters and relate them to light condition + sky acceptability.
# =============================================================================

library(tidyverse)
library(readxl)
library(psych)
library(cluster)

set.seed(42)

dir.create("output", showWarnings = FALSE)
dir.create("output/figures", showWarnings = FALSE)

# ---- 1. Load & clean ---------------------------------------------------------
raw <- read_excel("people_data.xlsx")

# 999 is the missing code throughout
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

cat("\n--- Missingness (after recoding 999 -> NA) ---\n")
print(colSums(is.na(dat)))
cat("\nRespondents with >=1 missing attitude item:",
    sum(!complete.cases(dat[attitude_items])), "of", nrow(dat), "\n")

# Median-impute the Likert items so all 166 respondents keep a factor score.
items <- dat[attitude_items]
items_imp <- items %>%
  mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))

# ---- 2. Factor analysis suitability -----------------------------------------
cat("\n--- Sampling adequacy / sphericity ---\n")
kmo  <- KMO(items_imp)
bart <- cortest.bartlett(cor(items_imp), n = nrow(items_imp))
cat("KMO overall MSA:", round(kmo$MSA, 3), "\n")
cat("Bartlett chi-sq:", round(bart$chisq, 1),
    " df:", bart$df, " p:", format.pval(bart$p.value), "\n")

# ---- 2b. How many factors? (parallel analysis) -------------------------------
png("output/figures/people_parallel_scree.png", width = 1000, height = 750, res = 130)
pa <- fa.parallel(items_imp, fa = "fa", fm = "ml",
                  main = "Parallel analysis - streetlight attitude items")
dev.off()
n_factors <- max(pa$nfact, 2)   # at least 2 for an interpretable structure
cat("\nParallel analysis suggests", pa$nfact, "factor(s); using", n_factors, "\n")

# ---- 3. Run the EFA ----------------------------------------------------------
efa <- fa(items_imp, nfactors = n_factors, rotate = "oblimin", fm = "ml")

cat("\n--- Factor loadings (oblimin-rotated, |loading| >= 0.30) ---\n")
print(efa$loadings, cutoff = 0.30, sort = TRUE)
cat("\nProportion of variance explained:\n")
print(round(efa$Vaccounted, 3))

# Save a tidy loadings table
loadings_tbl <- unclass(efa$loadings) %>%
  as.data.frame() %>%
  rownames_to_column("item") %>%
  mutate(communality = efa$communality[item])
write_csv(loadings_tbl, "output/people_efa_loadings.csv")

# ---- 4. Factor scores + clustering ------------------------------------------
scores <- as.data.frame(efa$scores)          # standardized factor scores
colnames(scores) <- paste0("F", seq_len(ncol(scores)))

# Choose k by average silhouette width over k = 2..6 (k-means)
sil_width <- sapply(2:6, function(k) {
  km <- kmeans(scores, centers = k, nstart = 50, iter.max = 100)
  mean(silhouette(km$cluster, dist(scores))[, "sil_width"])
})
names(sil_width) <- 2:6
cat("\n--- Average silhouette width by k (k-means) ---\n")
print(round(sil_width, 3))
best_k <- as.integer(names(which.max(sil_width)))
cat("Selected k =", best_k, "\n")

km <- kmeans(scores, centers = best_k, nstart = 50, iter.max = 100)

# Ward hierarchical cross-check
hc  <- hclust(dist(scores), method = "ward.D2")
hc_clusters <- cutree(hc, k = best_k)
cat("\nk-means vs Ward agreement (cross-tab):\n")
print(table(kmeans = km$cluster, ward = hc_clusters))

# ---- 5. Profile the clusters -------------------------------------------------
result <- dat %>%
  mutate(cluster = factor(km$cluster)) %>%
  bind_cols(scores)

# Mean factor score by cluster
factor_profile <- result %>%
  group_by(cluster) %>%
  summarise(n = n(), across(all_of(colnames(scores)), ~ round(mean(.), 2)),
            .groups = "drop")
cat("\n--- Cluster sizes & mean factor scores ---\n")
print(factor_profile)

# Mean of each original item by cluster (on 1-5 scale, imputed)
item_profile <- items_imp %>%
  mutate(cluster = factor(km$cluster)) %>%
  group_by(cluster) %>%
  summarise(across(everything(), ~ round(mean(.), 2)), .groups = "drop")
cat("\n--- Mean attitude-item response by cluster (1-5) ---\n")
print(as.data.frame(item_profile))

# Relationship to light condition and sky acceptability
cat("\n--- Cluster x StreetlightCondition (1=red, 2=white) ---\n")
print(table(cluster = result$cluster,
            light   = factor(result$StreetlightCondition,
                             labels = c("red", "white"))))
cat("\n--- Mean Acceptable_Sky_Conditions (1-7) by cluster ---\n")
print(result %>% group_by(cluster) %>%
        summarise(mean_sky = round(mean(Acceptable_Sky_Conditions, na.rm = TRUE), 2),
                  n = n(), .groups = "drop"))

# ---- Save outputs ------------------------------------------------------------
write_csv(factor_profile, "output/people_cluster_factor_profile.csv")
write_csv(item_profile,   "output/people_cluster_item_profile.csv")
write_csv(result %>% mutate(row_id = row_number()) %>%
            select(row_id, cluster, all_of(colnames(scores)),
                   StreetlightCondition, Acceptable_Sky_Conditions),
          "output/people_cluster_assignments.csv")

# Cluster visualization in factor space (first two factors)
p <- ggplot(result, aes(F1, .data[[if (ncol(scores) >= 2) "F2" else "F1"]],
                        color = cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(type = "norm", linewidth = 0.4) +
  labs(title = "Respondent clusters in factor space",
       x = "Factor 1", y = if (ncol(scores) >= 2) "Factor 2" else "Factor 1") +
  theme_minimal()
ggsave("output/figures/people_clusters_factorspace.png", p,
       width = 7, height = 5, dpi = 130)

cat("\nDone. Outputs written to output/ and output/figures/\n")
