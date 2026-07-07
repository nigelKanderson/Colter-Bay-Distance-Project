# 18_bat_moth_predation.R
# Tests whether bat activity predicts (suppresses) moth activity at the same
# site-night, net of the shared light color/intensity/distance/season gradient
# both taxa already respond to independently. This is a same-night spatial
# co-occurrence test for the "predation sink" hypothesis -- it can show
# whether the association survives controlling for the shared drivers, but it
# is not proof of a direct predation mechanism (no diet, telemetry, or
# mortality data available to confirm bats are actually eating these moths).
#
# Input:  data_env (bat, one row per site/date/species), insect_total_env
#         (moth, one row per site/date/color/intensity with total detections),
#         bat_theme, col_navy, col_purple
# Output: list(model, data, n_total, n_matched, fig_partial)

run_bat_moth_predation_test <- function(data_env, insect_total_env, bat_theme, col_navy, col_purple) {
  library(dplyr)
  library(glmmTMB)
  library(emmeans)
  library(ggplot2)

  # Total bat activity per site-night, pooling across the 6 focal species --
  # matches the "general activity" concept already used for moths.
  bat_nightly <- data_env %>%
    mutate(site = as.character(site)) %>%
    group_by(site, date) %>%
    summarise(bat_detections = sum(detections, na.rm = TRUE), .groups = "drop")

  moth_bat_env <- insect_total_env %>%
    mutate(site = as.character(site)) %>%
    left_join(bat_nightly, by = c("site", "date")) %>%
    mutate(
      color     = factor(color),
      site      = factor(site),
      intensity = factor(intensity),
      # log2 so the coefficient reads directly as "% change in moth
      # detections per doubling of bat activity that same night"
      log_bat   = log2(bat_detections + 1)
    )

  n_total   <- nrow(insect_total_env)
  n_matched <- sum(!is.na(moth_bat_env$bat_detections))
  cat("Moth site-nights:", n_total, " | matched to bat activity that same night:", n_matched, "\n")

  moth_bat_env <- moth_bat_env %>% filter(!is.na(bat_detections))

  formula <- detections ~ log_bat + color * intensity + intensity * dist_km +
    mean_phase + pct_nonforest + jd_c + I(jd_c^2) + (1 | site)

  model <- tryCatch(
    glmmTMB(formula, data = moth_bat_env, ziformula = ~1, family = nbinom2()),
    error = function(e) NULL
  )
  if (is.null(model) || !isTRUE(model$sdr$pdHess)) {
    cat("  Retrying with BFGS optimizer...\n")
    model <- tryCatch(
      glmmTMB(formula, data = moth_bat_env, ziformula = ~1, family = nbinom2(),
              control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))),
      error = function(e) NULL
    )
  }
  if (is.null(model)) stop("Bat-moth predation model failed to fit.")
  if (!isTRUE(model$sdr$pdHess)) cat("  WARNING: convergence uncertain (non-PD Hessian)\n")

  # ── Partial-effect figure: predicted moth detections vs. bat activity ────
  bat_range <- range(moth_bat_env$bat_detections, na.rm = TRUE)
  log_grid  <- seq(log2(bat_range[1] + 1), log2(bat_range[2] + 1), length.out = 50)

  pred_df <- as.data.frame(
    emmeans(model, ~ log_bat, at = list(log_bat = log_grid), type = "response")
  ) %>%
    rename(lower_CL = asymp.LCL, upper_CL = asymp.UCL) %>%
    mutate(bat_detections = 2^log_bat - 1)

  fig_partial <- ggplot(pred_df, aes(x = bat_detections, y = response)) +
    geom_ribbon(aes(ymin = lower_CL, ymax = upper_CL), alpha = 0.15, fill = col_purple) +
    geom_line(linewidth = 1.1, color = col_purple) +
    geom_rug(data = moth_bat_env, aes(x = bat_detections), inherit.aes = FALSE,
             sides = "b", alpha = 0.3, color = col_navy) +
    labs(
      title    = "Moth Activity vs. Same-Night Bat Activity",
      subtitle = "Marginal prediction ± 95% CI, controlling for color, intensity, distance, moon phase, habitat, season",
      x        = "Total bat detections that night (same site)",
      y        = "Predicted moth detections"
    ) +
    bat_theme

  dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
  ggsave("output/figures/fig_bat_moth_predation.png", fig_partial,
         width = 8, height = 5, dpi = 300, bg = "white")

  list(
    model       = model,
    data        = moth_bat_env,
    n_total     = n_total,
    n_matched   = n_matched,
    fig_partial = fig_partial
  )
}
