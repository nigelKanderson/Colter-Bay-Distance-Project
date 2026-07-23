# 12_trait_models.R
# Per-trait GLMMs + emmeans summaries + figures — Distance Project
# Tests whether bat morphology (ear height, forearm length, wing loading,
# aspect ratio) moderates the effect of light intensity and distance,
# Input:  data_env (already joined to bat_traits), bat_theme, col_navy, col_red, col_purple
# Output: trait_models list, summary data frames, fig_trait_intensity, fig_trait_dist

run_trait_models <- function(data_env, bat_theme, col_navy, col_red, col_purple,
                             file_prefix = "") {

  library(glmmTMB)
  library(emmeans)
  library(tidyverse)

  traits <- c("ear_height", "forearm_length", "wing_loading", "aspect_ratio")
  traits <- traits[traits %in% names(data_env)]

  if (!"jd_c" %in% names(data_env)) {
    data_env$jd_c <- data_env$jd - mean(data_env$jd, na.rm = TRUE)
  }

  # Representative trait values (mean ± 1 SD) for simple-slopes style comparisons,
  trait_levels <- c(Low = -1, Mean = 0, High = 1)

  fit_trait <- function(trait) {
    z_col <- paste0(trait, "_z")
    dat <- data_env
    dat[[z_col]] <- as.numeric(scale(dat[[trait]]))

    form <- as.formula(paste0(
      "detections ~ ", z_col, " * intensity + ", z_col, " * dist_km + ",
      "mean_phase + pct_nonforest + jd_c + I(jd_c^2) + (1 | site)"
    ))

    mod <- tryCatch(
      glmmTMB(form, family = nbinom2, ziformula = ~1, data = dat),
      error = function(e) NULL
    )
    if (is.null(mod) || !isTRUE(mod$sdr$pdHess)) {
      mod <- tryCatch(
        glmmTMB(form, family = nbinom2, ziformula = ~1, data = dat,
                control = glmmTMBControl(optimizer = optim,
                                         optArgs   = list(method = "BFGS"))),
        error = function(e) NULL
      )
    }
    if (is.null(mod) || !isTRUE(mod$sdr$pdHess))
      cat("  WARNING:", trait, "— convergence uncertain after both optimizers\n")

    list(mod = mod, z_col = z_col)
  }

  trait_fits <- setNames(lapply(traits, fit_trait), traits)
  trait_fits <- trait_fits[!sapply(trait_fits, function(x) is.null(x$mod))]
  trait_models <- lapply(trait_fits, `[[`, "mod")

  if (length(trait_models) == 0) stop("No trait models converged.")
  cat("Successfully fit trait models:", paste(names(trait_models), collapse = ", "), "\n")

  norm_ci <- function(df) {
    lo <- intersect(c("asymp.LCL", "lower.CL", "lower.HPD", "lwr.CL", "lower"), names(df))[1]
    hi <- intersect(c("asymp.UCL", "upper.CL", "upper.HPD", "upr.CL", "upper"), names(df))[1]
    if (!is.na(lo) && !is.na(hi)) df <- rename(df, lower_CL = !!lo, upper_CL = !!hi)
    df
  }

  label_trait_level <- function(df, z_col) {
    df %>%
      mutate(
        trait_level = factor(names(trait_levels)[match(round(.data[[z_col]], 6),
                                                         round(trait_levels, 6))],
                             levels = names(trait_levels))
      ) %>%
      select(-all_of(z_col))
  }

  # ── Intensity effect at Low/Mean/High trait values ────────────────────────
  extract_trait_intensity <- function(trait) {
    mod   <- trait_fits[[trait]]$mod
    z_col <- trait_fits[[trait]]$z_col
    at_list <- setNames(list(trait_levels), z_col)
    tryCatch({
      emmeans(mod, as.formula(paste("~ intensity |", z_col)),
              at = at_list, type = "response") %>%
        as.data.frame() %>% norm_ci() %>%
        label_trait_level(z_col) %>%
        mutate(trait = trait)
    }, error = function(e) {
      cat("  emmeans (intensity|trait) failed for", trait, ":", conditionMessage(e), "\n")
      NULL
    })
  }

  df_trait_intensity <- map_dfr(names(trait_fits), extract_trait_intensity)

  # ── Distance slope at Low/Mean/High trait values ──────────────────────────
  extract_trait_slope <- function(trait) {
    mod   <- trait_fits[[trait]]$mod
    z_col <- trait_fits[[trait]]$z_col
    at_list <- setNames(list(trait_levels), z_col)
    tryCatch({
      emtrends(mod, as.formula(paste("~", z_col)), var = "dist_km",
               at = at_list, type = "response") %>%
        as.data.frame() %>% norm_ci() %>%
        label_trait_level(z_col) %>%
        mutate(trait = trait)
    }, error = function(e) {
      cat("  emtrends (dist_km|trait) failed for", trait, ":", conditionMessage(e), "\n")
      NULL
    })
  }

  df_trait_slope <- map_dfr(names(trait_fits), extract_trait_slope)

  # ── Full distance-gradient predictions at Low/Mean/High trait values ─────
  dist_grid <- seq(min(data_env$dist_km, na.rm = TRUE),
                    max(data_env$dist_km, na.rm = TRUE), length.out = 60)

  extract_trait_dist_pred <- function(trait) {
    mod   <- trait_fits[[trait]]$mod
    z_col <- trait_fits[[trait]]$z_col
    at_list <- setNames(list(trait_levels, dist_grid), c(z_col, "dist_km"))
    tryCatch({
      emmeans(mod, as.formula(paste("~", z_col, "+ dist_km")),
              at = at_list, type = "response") %>%
        as.data.frame() %>% norm_ci() %>%
        label_trait_level(z_col) %>%
        mutate(trait = trait)
    }, error = function(e) {
      cat("  dist emmeans failed for", trait, ":", conditionMessage(e), "\n")
      NULL
    })
  }

  df_trait_dist <- map_dfr(names(trait_fits), extract_trait_dist_pred)

  # ── Fig: intensity effect by trait level, faceted by trait ───────────────
  fig_trait_intensity <- ggplot(df_trait_intensity,
                                aes(x = trait_level, y = response,
                                    color = intensity, group = intensity)) +
    geom_line(linewidth = 0.8, position = position_dodge(0.4)) +
    geom_point(size = 3, position = position_dodge(0.4)) +
    geom_errorbar(aes(ymin = lower_CL, ymax = upper_CL),
                  width = 0.15, linewidth = 0.6,
                  position = position_dodge(0.4)) +
    scale_color_manual(values = c("10" = "#F5F0E3", "30" = "#C9A84C",
                                  "50" = "#4A1259", "70" = "#7A1C2E",
                                  "100" = "#1B2A4A"),
                       name = "Intensity (%)") +
    facet_wrap(~ trait, scales = "free_y") +
    labs(
      title    = "Light Intensity Effect Across Morphology Trait Values",
      subtitle = "Marginal means ± 95% CI at low/mean/high (±1 SD) trait values",
      x        = "Trait value",
      y        = "Predicted detections per night"
    ) +
    bat_theme +
    theme(legend.position = "top",
          strip.text = element_text(size = 9, face = "bold"))

  # ── Fig: distance gradient by trait level, faceted by trait ──────────────
  fig_trait_dist <- ggplot(df_trait_dist,
                           aes(x = dist_km, y = response,
                               color = trait_level, fill = trait_level)) +
    geom_ribbon(aes(ymin = lower_CL, ymax = upper_CL), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = c(Low = col_navy, Mean = col_red, High = col_purple),
                       name = "Trait value") +
    scale_fill_manual(values  = c(Low = col_navy, Mean = col_red, High = col_purple),
                      guide   = "none") +
    facet_wrap(~ trait, scales = "free_y") +
    labs(
      title    = "Bat Activity Across Distance Gradient by Morphology Trait Value",
      subtitle = "Marginal predictions ± 95% CI · low/mean/high (±1 SD) trait values",
      x        = "Distance from Colter Bay parking lot (km)",
      y        = "Predicted detections per night"
    ) +
    bat_theme +
    theme(legend.position = "top",
          strip.text = element_text(size = 9, face = "bold"))

  # ── Save figures + post-hoc tables ────────────────────────────────────────
  dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0("output/figures/", file_prefix, "fig_trait_intensity.png"), fig_trait_intensity,
         width = 10, height = 7, dpi = 300, bg = "white")
  ggsave(paste0("output/figures/", file_prefix, "fig_trait_dist.png"), fig_trait_dist,
         width = 10, height = 7, dpi = 300, bg = "white")

  write_csv(df_trait_intensity, "output/posthoc_trait_intensity.csv")
  write_csv(df_trait_slope,     "output/posthoc_trait_dist_slope.csv")
  write_csv(df_trait_dist,      "output/posthoc_trait_dist_pred.csv")
  cat("Saved trait figures and post-hoc CSVs to output/\n")

  list(
    trait_models        = trait_models,
    df_trait_intensity  = df_trait_intensity,
    df_trait_slope      = df_trait_slope,
    df_trait_dist       = df_trait_dist,
    fig_trait_intensity = fig_trait_intensity,
    fig_trait_dist      = fig_trait_dist
  )
}
