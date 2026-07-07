# 20_bat_season_distance.R
# Bat-side counterpart to 19_moth_season_distance.R: tests whether the bat
# distance gradient strengthens over the course of the season. Same formula
# structure as the bat Primary Model, with a dist_km:jd_c term added.
#
# Input:  data_env, bat_theme, col_navy, col_red, col_purple
# Output: list(model, season_dist_df, fig_season_distance, fig_season_distance_relative)

run_bat_season_distance_test <- function(data_env, bat_theme, col_navy, col_red, col_purple) {
  library(dplyr)
  library(glmmTMB)
  library(emmeans)
  library(ggplot2)

  if (!"jd_c" %in% names(data_env)) {
    data_env$jd_c <- data_env$jd - mean(data_env$jd, na.rm = TRUE)
  }

  formula <- detections ~ color * intensity + intensity * dist_km +
    jd_c + I(jd_c^2) + dist_km:jd_c +
    mean_phase + pct_nonforest + brightness_dark + (1 | site)

  model <- tryCatch(
    glmmTMB(formula, data = data_env, ziformula = ~1, family = nbinom2()),
    error = function(e) NULL
  )
  if (is.null(model) || !isTRUE(model$sdr$pdHess)) {
    cat("  Retrying with BFGS optimizer...\n")
    model <- tryCatch(
      glmmTMB(formula, data = data_env, ziformula = ~1, family = nbinom2(),
              control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))),
      error = function(e) NULL
    )
  }
  if (is.null(model)) stop("Bat season x distance model failed to fit.")
  if (!isTRUE(model$sdr$pdHess)) cat("  WARNING: convergence uncertain (non-PD Hessian)\n")

  # ── Representative season stages: earliest / middle / latest sampled dates ─
  early <- min(data_env$jd_c, na.rm = TRUE)
  mid   <- median(data_env$jd_c, na.rm = TRUE)
  late  <- max(data_env$jd_c, na.rm = TRUE)
  season_levels <- c(Early = early, Mid = mid, Late = late)

  date_lookup <- data_env %>% distinct(jd_c, date)
  label_date  <- function(jdc) {
    closest <- date_lookup$date[which.min(abs(date_lookup$jd_c - jdc))]
    format(closest, "%b %d")
  }
  season_labels <- setNames(paste0(names(season_levels), " (", sapply(season_levels, label_date), ")"),
                            names(season_levels))

  dist_grid <- seq(min(data_env$dist_km, na.rm = TRUE),
                    max(data_env$dist_km, na.rm = TRUE), length.out = 60)

  season_dist_df <- as.data.frame(
    emmeans(model, ~ dist_km + jd_c,
            at = list(dist_km = dist_grid, jd_c = season_levels),
            type = "response")
  ) %>%
    rename(lower_CL = asymp.LCL, upper_CL = asymp.UCL) %>%
    mutate(
      season = factor(names(season_levels)[match(round(jd_c, 6), round(season_levels, 6))],
                      levels = names(season_levels))
    )

  fig_season_distance <- ggplot(season_dist_df,
                                aes(x = dist_km, y = response, color = season, fill = season)) +
    geom_ribbon(aes(ymin = lower_CL, ymax = upper_CL), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = c(Early = col_navy, Mid = col_red, Late = col_purple),
                       labels = season_labels, name = "Season stage") +
    scale_fill_manual(values  = c(Early = col_navy, Mid = col_red, Late = col_purple), guide = "none") +
    labs(
      title    = "Bat Distance Gradient Across the Season",
      subtitle = "Marginal predictions ± 95% CI · averaged over color and intensity",
      x        = "Distance from light source (km)",
      y        = "Predicted bat detections per night"
    ) +
    bat_theme +
    theme(legend.position = "top")

  # ── Rescaled version: % change from 0 km within each season stage ────────
  season_dist_rel_df <- season_dist_df %>%
    group_by(season) %>%
    mutate(pct_change = (response / response[which.min(dist_km)] - 1) * 100) %>%
    ungroup()

  fig_season_distance_relative <- ggplot(season_dist_rel_df,
                                         aes(x = dist_km, y = pct_change, color = season)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#999999") +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = c(Early = col_navy, Mid = col_red, Late = col_purple),
                       labels = season_labels, name = "Season stage") +
    labs(
      title    = "Bat Distance Gradient Shape Across the Season",
      subtitle = "% change from each stage's own 0 km baseline · isolates slope from the seasonal peak",
      x        = "Distance from light source (km)",
      y        = "% change from 0 km"
    ) +
    bat_theme +
    theme(legend.position = "top")

  dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
  ggsave("output/figures/fig_bat_season_distance.png", fig_season_distance,
         width = 8, height = 5, dpi = 300, bg = "white")
  ggsave("output/figures/fig_bat_season_distance_relative.png", fig_season_distance_relative,
         width = 8, height = 5, dpi = 300, bg = "white")

  list(
    model                        = model,
    season_dist_df               = season_dist_df,
    fig_season_distance          = fig_season_distance,
    fig_season_distance_relative = fig_season_distance_relative
  )
}
