# 11_species_models.R
# Per-species GLMMs + emmeans summaries + figures — Distance Project
# Input:  data_env, bat_theme, col_red, col_white
# Output: sp_models list, summary data frames, figures S1–S4

run_species_models <- function(data_env, bat_theme, col_red, col_white) {

  library(glmmTMB)
  library(emmeans)
  library(tidyverse)
  library(patchwork)

  # Build formula — omit brightness_dark if Python script hasn't been run yet
  has_brightness <- "brightness_dark" %in% names(data_env)
  has_jd_c       <- "jd_c"           %in% names(data_env)

  if (!has_jd_c) {
    data_env$jd_c <- data_env$jd - mean(data_env$jd, na.rm = TRUE)
    cat("  jd_c computed on-the-fly\n")
  }

  fixed_terms <- c(
    "jd_c + I(jd_c^2)",
    "mean_phase",
    "pct_nonforest",
    if (has_brightness) "brightness_dark",
    "intensity * color",
    "dist_km"
  )
  sp_formula <- as.formula(
    paste("detections ~", paste(fixed_terms, collapse = " + "), "+ (1 | site)")
  )
  cat("Species formula:", deparse(sp_formula), "\n")

  fit_sp <- function(dat, sp) {
    mod <- tryCatch(
      glmmTMB(sp_formula, family = nbinom2, ziformula = ~1, data = dat),
      error = function(e) NULL
    )
    # Fallback: BFGS if Hessian not positive-definite
    if (is.null(mod) || !isTRUE(mod$sdr$pdHess)) {
      mod <- tryCatch(
        glmmTMB(sp_formula, family = nbinom2, ziformula = ~1, data = dat,
                control = glmmTMBControl(optimizer = optim,
                                         optArgs   = list(method = "BFGS"))),
        error = function(e) NULL
      )
    }
    if (is.null(mod) || !isTRUE(mod$sdr$pdHess))
      cat("  WARNING:", sp, "— convergence uncertain after both optimizers\n")
    mod
  }

  sp_models <- list()

  for (sp in unique(data_env$species)) {
    dat_sp <- data_env %>% filter(species == sp)

    if (sum(dat_sp$detections) < 30 || nrow(dat_sp) < 20) {
      cat("  Skipping", sp, "— insufficient data\n")
      next
    }

    cat("  Fitting:", sp, "\n")
    sp_models[[sp]] <- fit_sp(dat_sp, sp)
  }

  sp_models <- sp_models[!sapply(sp_models, is.null)]
  cat("Successfully fit:", paste(names(sp_models), collapse = ", "), "\n")

  if (length(sp_models) == 0) {
    stop("No species models converged. This usually means:\n",
         "  1. brightness_dark is missing (run the Python sky-brightness script first), OR\n",
         "  2. No species pass the 30-detection / 20-row threshold.\n",
         "  has_brightness = ", has_brightness, "\n",
         "  Species in data: ", paste(sort(unique(data_env$species)), collapse=", "))
  }

  # ── Extract emmeans summaries ──────────────────────────────────────────────
  # intensity is a FACTOR in data_env, so use factor-level emmeans (no `at`)
  # and pairwise contrasts instead of emtrends (which needs a continuous var)

  # Normalize CI column names — emmeans naming varies by version and model family.
  # Detect lower/upper CI columns dynamically and rename to lower_CL / upper_CL.
  norm_ci <- function(df) {
    lo <- intersect(c("asymp.LCL", "lower.CL", "lower.HPD", "lwr.CL", "lower"),
                    names(df))[1]
    hi <- intersect(c("asymp.UCL", "upper.CL", "upper.HPD", "upr.CL", "upper"),
                    names(df))[1]
    if (!is.na(lo) && !is.na(hi)) {
      df <- rename(df, lower_CL = !!lo, upper_CL = !!hi)
    } else {
      # Fallback: approximate ±1.96 SE on whatever the estimate column is
      est_col <- intersect(c("response", "ratio", "emmean", "estimate"), names(df))[1]
      if (!is.na(est_col) && "SE" %in% names(df)) {
        df$lower_CL <- df[[est_col]] - 1.96 * df$SE
        df$upper_CL <- df[[est_col]] + 1.96 * df$SE
        cat("  norm_ci: no CI cols found; approximated from SE. Cols:", paste(names(df), collapse=", "), "\n")
      }
    }
    df
  }

  extract_color_intensity <- function(mod, sp_name) {
    tryCatch({
      emmeans(mod, ~ color * intensity, type = "response") %>%
        as.data.frame() %>% norm_ci() %>% mutate(species = sp_name)
    }, error = function(e) {
      cat("  emmeans (color*intensity) failed for", sp_name, ":", conditionMessage(e), "\n")
      NULL
    })
  }

  # White vs Red contrast at each intensity level
  extract_color_contrast <- function(mod, sp_name) {
    tryCatch({
      emmeans(mod, ~ color | intensity, type = "response") %>%
        contrast(method = "revpairwise", adjust = "none") %>%
        as.data.frame() %>% norm_ci() %>% mutate(species = sp_name)
    }, error = function(e) {
      cat("  color contrast failed for", sp_name, ":", conditionMessage(e), "\n")
      NULL
    })
  }

  extract_color_main <- function(mod, sp_name) {
    tryCatch({
      emmeans(mod, ~ color, type = "response") %>%
        as.data.frame() %>% norm_ci() %>% mutate(species = sp_name)
    }, error = function(e) {
      cat("  emmeans (color) failed for", sp_name, ":", conditionMessage(e), "\n")
      NULL
    })
  }

  df_sp_ci      <- map2_dfr(sp_models, names(sp_models), extract_color_intensity)
  df_contrasts  <- map2_dfr(sp_models, names(sp_models), extract_color_contrast)
  df_sp_color   <- map2_dfr(sp_models, names(sp_models), extract_color_main)

  # ── Fig S1: Color effect by species ───────────────────────────────────────
  fig_s1 <- ggplot(df_sp_color,
                   aes(x = response,
                       y = fct_reorder(species, response),
                       color = color)) +
    geom_line(aes(group = species), color = "#CCCCCC", linewidth = 0.8) +
    geom_point(size = 3.5) +
    geom_errorbar(aes(xmin = lower_CL, xmax = upper_CL),
                  width = 0.3, linewidth = 0.7) +
    scale_color_manual(values = c(R = col_red, W = col_white),
                       labels = c(R = "Red", W = "White"),
                       name   = "Light color") +
    labs(
      title    = "Light Color Effect on Bat Detections by Species",
      subtitle = "Marginal means ± 95% CI from per-species negative binomial GLMMs",
      x        = "Predicted detections per night",
      y        = NULL,
      caption  = "Connected lines link red/white estimates for the same species"
    ) +
    bat_theme +
    theme(legend.position = "top")

  # ── Fig S2: White vs Red contrast by intensity level ──────────────────────
  # (replaces emtrends slope — intensity is a factor so no continuous slope)
  fig_s2 <- ggplot(df_contrasts,
                   aes(x    = ratio,
                       y    = fct_reorder(species, ratio, .fun = mean),
                       color = intensity)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "#AAAAAA", linewidth = 0.7) +
    geom_point(position = position_dodge(0.6), size = 2.8) +
    geom_errorbar(aes(xmin = lower_CL, xmax = upper_CL),
                  position = position_dodge(0.6),
                  width = 0.35, linewidth = 0.7) +
    scale_color_manual(values = c("10" = "#F5F0E3", "30" = "#C9A84C",
                                   "50" = "#4A1259", "70" = "#7A1C2E",
                                   "100" = "#1B2A4A"),
                       name = "Intensity") +
    labs(
      title    = "White vs. Red Light — Detection Ratio by Species & Intensity",
      subtitle = "Ratio > 1 = more detections under white; < 1 = more under red",
      x        = "White / Red detection ratio (response scale)",
      y        = NULL,
      caption  = "Pairwise EMMs contrast: white ÷ red at each intensity level"
    ) +
    bat_theme +
    theme(legend.position = "right")

  # ── Fig S4: Predicted detections ~ distance × color × intensity per species ─
  # Predict over a fine dist_km grid so the gradient is continuous.
  dist_grid <- seq(
    min(data_env$dist_km, na.rm = TRUE),
    max(data_env$dist_km, na.rm = TRUE),
    length.out = 40
  )

  extract_dist_pred <- function(mod, sp_name) {
    tryCatch({
      emmeans(mod, ~ color + dist_km,
              at   = list(dist_km = dist_grid),
              type = "response") %>%
        as.data.frame() %>% norm_ci() %>%
        mutate(species = sp_name)
    }, error = function(e) {
      cat("  dist emmeans failed for", sp_name, ":", conditionMessage(e), "\n")
      NULL
    })
  }

  df_dist_pred <- map2_dfr(sp_models, names(sp_models), extract_dist_pred)

  # Figure S4: one panel per species, x=dist_km, color=light color
  fig_s4a <- ggplot(df_dist_pred,
                    aes(x     = dist_km,
                        y     = response,
                        color = color,
                        fill  = color,
                        group = color)) +
    geom_ribbon(aes(ymin = lower_CL, ymax = upper_CL),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = c(R = col_red, W = col_white),
                       labels = c(R = "Red", W = "White"),
                       name   = "Light color") +
    scale_fill_manual(values  = c(R = col_red, W = col_white),
                      guide   = "none") +
    facet_wrap(~ species, scales = "free_y", nrow = 2) +
    labs(
      title    = "Predicted Bat Detections Across Distance Gradient by Light Color",
      subtitle = "Per-species GLMMs · marginal over intensity · ribbon = 95% CI",
      x        = "Distance from light source (km)",
      y        = "Predicted detections per night"
    ) +
    bat_theme +
    theme(
      legend.position = "top",
      strip.text      = element_text(size = 9, face = "bold")
    )

  fig_s4b <- NULL  # retired — single combined figure is cleaner

  # ── Fig S3: White:red log2 ratio heatmap ──────────────────────────────────
  df_heat <- df_sp_ci %>%
    select(any_of(c("species", "color", "intensity", "response"))) %>%
    pivot_wider(
      id_cols     = c(species, intensity),
      names_from  = color,
      values_from = response
    ) %>%
    mutate(
      log2_ratio = log2(W / R),
      intensity  = factor(intensity,
                          levels = c("0", "30", "50", "70", "100"),
                          labels = c("Dark (0%)", "Low (30%)", "Mid (50%)",
                                     "High (70%)", "Max (100%)"))
    )

  fig_s3 <- ggplot(df_heat,
                   aes(x    = intensity,
                       y    = fct_reorder(species, log2_ratio, .fun = mean),
                       fill = log2_ratio)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = round(log2_ratio, 2)),
              color = "white", size = 3.2, fontface = "bold") +
    scale_fill_gradient2(low      = col_red,
                         mid      = "#F5F0E3",
                         high     = col_white,
                         midpoint = 0,
                         name     = "log₂(White/Red)") +
    labs(
      title    = "White vs. Red Light Effect by Species and Intensity",
      subtitle = "Negative values (red) = white light suppresses relative to red",
      x        = "Light intensity",
      y        = NULL,
      caption  = "log₂ ratio of predicted detections: white ÷ red"
    ) +
    bat_theme +
    theme(panel.grid = element_blank(),
          axis.line  = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right")

  # ── Save figures ───────────────────────────────────────────────────────────
  dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

  ggsave("output/figures/species_color_effect.png",
         (fig_s1 | fig_s2) + plot_annotation(
           title = "Species-Level Light Effects — Distance Project",
           theme = theme(plot.title = element_text(size = 14, face = "bold"))
         ),
         width = 14, height = 7, dpi = 300, bg = "white")

  ggsave("output/figures/species_ratio_heatmap.png",
         fig_s3, width = 9, height = 6, dpi = 300, bg = "white")

  n_sp <- length(sp_models)
  ggsave("output/figures/fig_s4a_species_dist_color_intensity.png",
         fig_s4a,
         width  = max(14, n_sp * 1.8),
         height = 8,
         dpi    = 300, bg = "white")


  write_csv(df_contrasts, "output/posthoc_color_contrasts_by_species.csv")
  write_csv(df_sp_color,  "output/posthoc_color_main_by_species.csv")
  write_csv(df_dist_pred, "output/posthoc_dist_color_intensity_by_species.csv")
  cat("Saved species figures and post-hoc CSVs to output/\n")

  list(
    sp_models    = sp_models,
    df_sp_ci     = df_sp_ci,
    df_contrasts = df_contrasts,
    df_sp_color  = df_sp_color,
    df_dist_pred = df_dist_pred,
    fig_s1       = fig_s1,
    fig_s2       = fig_s2,
    fig_s3       = fig_s3,
    fig_s4a      = fig_s4a,
    fig_s4b      = fig_s4b
  )
}
