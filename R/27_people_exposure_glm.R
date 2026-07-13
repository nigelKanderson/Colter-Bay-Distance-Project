# =============================================================================
# 27_people_exposure_glm.R
# Logistic GLM: can survey attitudes discriminate red vs white light exposure?
#
# Outcome: StreetlightCondition (1 = red, 2 = white) -> white = 1, red = 0.
# Predictors: the reduced attitude scores from the two earlier analyses
#   - EFA factor scores   (F1 wildlife concern, F2 experience, F3 discomfort)
#   - PCA component scores (PC1, PC2) + Acceptable_Sky_Conditions
#
# This is a plain binomial GLM. NOTE ON DIRECTION: light condition was experimentally
# assigned and attitudes were measured after exposure, so this is a
# discriminant / association model ("which responses distinguish the groups"),
# not evidence that attitudes cause exposure.
#
# Requires the score files written by 25_people_factor_cluster.R and
# 26_people_pca.R.
# =============================================================================

library(tidyverse)

dir.create("output", showWarnings = FALSE)

efa <- read_csv("output/people_cluster_assignments.csv", show_col_types = FALSE)
pca <- read_csv("output/people_pca_component_scores.csv", show_col_types = FALSE)

d <- efa %>%
  left_join(pca %>% select(row_id, starts_with("PC")), by = "row_id") %>%
  filter(StreetlightCondition %in% c(1, 2)) %>%
  mutate(white = as.integer(StreetlightCondition == 2))   # 1 = white, 0 = red

cat(sprintf("n = %d  (white = %d, red = %d)\n\n",
            nrow(d), sum(d$white), sum(d$white == 0)))

# Area under the ROC curve (Mann-Whitney form, no extra packages)
auc <- function(m) {
  p <- predict(m, type = "response"); y <- m$model[[1]]
  r <- rank(p)
  (sum(r[y == 1]) - sum(y == 1) * (sum(y == 1) + 1) / 2) /
    (sum(y == 1) * sum(y == 0))
}

report <- function(m, label) {
  cat("===", label, "===\n")
  co <- summary(m)$coefficients
  tab <- data.frame(round(co, 4),
                    OR = round(exp(coef(m)), 3),
                    check.names = FALSE)
  print(tab)
  cat(sprintf("McFadden R2 = %.3f | AUC = %.3f | AIC = %.1f\n\n",
              1 - m$deviance / m$null.deviance, auc(m), AIC(m)))
  invisible(tab)
}

# ---- Models ------------------------------------------------------------------
mA <- glm(white ~ F1 + F2 + F3, data = d, family = binomial)
report(mA, "Model A: white ~ EFA factor scores (F1 wildlife, F2 experience, F3 discomfort)")

mB <- glm(white ~ PC1_score + PC2_score, data = d, family = binomial)
report(mB, "Model B: white ~ PCA component scores (PC1, PC2)")

mC <- glm(white ~ PC1_score + PC2_score + Acceptable_Sky_Conditions,
          data = d, family = binomial)
report(mC, "Model C: Model B + Acceptable_Sky_Conditions")

# ---- Save a tidy results table ----------------------------------------------
tidy_glm <- function(m, label) {
  co <- summary(m)$coefficients
  tibble(model = label,
         term = rownames(co),
         estimate = co[, 1], std_error = co[, 2],
         z = co[, 3], p = co[, 4],
         odds_ratio = exp(co[, 1]),
         mcfadden_r2 = 1 - m$deviance / m$null.deviance,
         auc = auc(m), aic = AIC(m))
}
bind_rows(tidy_glm(mA, "A: EFA factors"),
          tidy_glm(mB, "B: PCA components"),
          tidy_glm(mC, "C: PCA + sky")) %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  write_csv("output/people_exposure_glm_results.csv")

cat("Done. Results written to output/people_exposure_glm_results.csv\n")
