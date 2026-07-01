# 10_community_analysis.R
# Community-level bat activity analysis — Distance Project
# Input:     data_env  — cleaned detection data with all covariates
#            bat_theme — shared ggplot2 theme object
#            col_red / col_white — shared colour constants
# Returns:   named list of model objects, data frames, and figures

run_community_analysis <- function(data_env, bat_theme, col_red, col_white,
                                    taxon = "Bat", unit = "Species") {

  library(vegan)
  library(mvabund)
  library(tidyverse)
  library(ggrepel)

  # Drop sf geometry if present (left_join with sites_sf can attach geometry)
  if (inherits(data_env, "sf")) {
    data_env <- sf::st_drop_geometry(data_env)
  }
  data_env <- as.data.frame(data_env)

  cat("data_env columns:", paste(sort(names(data_env)), collapse = ", "), "\n")

  if (!"species" %in% names(data_env)) {
    # Look for a similarly named column and suggest a fix
    near <- grep("spec|spp|bat", names(data_env), ignore.case = TRUE, value = TRUE)
    stop("'species' column not found in data_env.\n",
         "  All columns: ", paste(names(data_env), collapse = ", "), "\n",
         if (length(near)) paste0("  Similar columns found: ", paste(near, collapse=", "), "\n"))
  }

  data_env <- data_env %>% filter(!is.na(dist_km))

  cat("Community analysis: ", nrow(data_env), "rows,",
      n_distinct(data_env$species), "species\n")

  # Optional covariates: handle case where sky_brightness hasn't been run yet
  has_brightness <- "brightness_dark" %in% names(data_env)
  has_year       <- "year"            %in% names(data_env)

  # ── Build effort-corrected community data ──────────────────────────────────
  grp_vars <- if (has_year) c("site","year","color","intensity") else c("site","color","intensity")

  effort <- data_env %>%
    group_by(across(all_of(grp_vars))) %>%
    summarise(nights = n_distinct(date), .groups = "drop")

  comm_dat <- data_env %>%
    group_by(across(all_of(c(grp_vars, "species")))) %>%
    summarise(detections = sum(detections), .groups = "drop") %>%
    left_join(effort, by = grp_vars) %>%
    mutate(det_per_night = detections / nights)

  # ── Build community matrix using id_cols (avoids select(... species ...)) ──
  community_matrix <- comm_dat %>%
    pivot_wider(
      id_cols     = all_of(grp_vars),
      names_from  = species,
      values_from = det_per_night,
      values_fill = 0
    ) %>%
    left_join(
      {
        env_sum <- data_env %>%
          group_by(across(all_of(grp_vars))) %>%
          summarise(
            mean_phase       = mean(mean_phase,    na.rm = TRUE),
            pct_nonforest    = first(pct_nonforest),
            dist_km          = first(dist_km[!is.na(dist_km)]),
            .groups = "drop"
          )
        if (has_brightness)
          env_sum <- env_sum %>%
            left_join(
              data_env %>% group_by(across(all_of(grp_vars))) %>%
                summarise(brightness_dark = first(brightness_dark), .groups = "drop"),
              by = grp_vars
            )
        env_sum
      },
      by = grp_vars
    )

  env_meta_cols <- c(grp_vars, "mean_phase", "pct_nonforest", "dist_km",
                     if (has_brightness) "brightness_dark")
  species_cols  <- setdiff(names(community_matrix), env_meta_cols)

  comm_full <- community_matrix %>% select(all_of(species_cols))
  env_full  <- community_matrix %>%
    select(all_of(env_meta_cols)) %>%
    mutate(color = factor(color))

  # Drop rows with any NA in env covariates (keeps comm and env aligned)
  complete_rows <- complete.cases(env_full)
  if (any(!complete_rows)) {
    cat("  Dropping", sum(!complete_rows), "rows with NA covariates from ordination\n")
  }
  comm <- comm_full[complete_rows, , drop = FALSE]
  env  <- env_full[complete_rows,  , drop = FALSE]

  # ── NMDS ──────────────────────────────────────────────────────────────────
  cat("Running NMDS...\n")
  nmds <- metaMDS(comm, distance = "bray", k = 2, trymax = 100, trace = FALSE)
  cat("  Stress:", round(nmds$stress, 3), "\n")

  multi_year <- has_year && n_distinct(env$year) > 1

  # Build covariate terms string (add brightness_dark only if available)
  cov_terms <- paste(
    c("color", "intensity", "dist_km", "mean_phase", "pct_nonforest",
      if (has_brightness) "brightness_dark"),
    collapse = " + "
  )

  # ── PERMANOVA ─────────────────────────────────────────────────────────────
  cat("Running PERMANOVA...\n")
  perm_formula <- as.formula(paste("comm ~", cov_terms))
  perm <- adonis2(
    perm_formula,
    data         = env,
    method       = "bray",
    permutations = 9999,
    strata       = if (multi_year) env$year else NULL
  )

  # ── Homogeneity of dispersion ──────────────────────────────────────────────
  bc        <- vegdist(comm, method = "bray")
  bd        <- betadisper(bc, group = env$color)
  disp_test <- anova(bd)

  # ── dbRDA (year partialled only if multiple years present) ────────────────
  cat("Running dbRDA...\n")
  db_rhs <- if (multi_year) {
    paste(cov_terms, "+ Condition(year)")
  } else {
    cov_terms
  }
  db_formula <- as.formula(paste("comm ~", db_rhs))
  db <- dbrda(db_formula, data = env, distance = "bray")
  db_global <- anova(db, permutations = 9999)
  db_margin <- anova(db, by = "margin", permutations = 9999)

  # ── mvabund ───────────────────────────────────────────────────────────────
  cat("Running mvabund (nBoot = 9999 — may take a while)...\n")
  Y      <- mvabund(comm)
  mv_rhs <- if (multi_year) {
    paste(cov_terms, "+ year")
  } else {
    cov_terms
  }
  mv_formula <- as.formula(paste("Y ~", mv_rhs))
  mv_fit <- manyglm(
    mv_formula,
    family = "negative.binomial",
    data   = env
  )
  mv_anova <- anova(mv_fit, p.uni = "adjusted", nBoot = 9999)

  # ── Fig C1: dbRDA ordination biplot ───────────────────────────────────────
  site_scores <- scores(db, display = "sites") %>%
    as.data.frame() %>%
    bind_cols(env)

  bp_mat     <- scores(db, display = "bp")
  env_scores <- as.data.frame(bp_mat) %>%
    mutate(
      variable = rownames(bp_mat),
      variable = dplyr::recode(variable,
        colorW            = "White light",
        intensity         = "Intensity",
        mean_phase        = "Moon phase",
        pct_nonforest     = "% Non-forest",
        brightness_dark   = "Sky brightness",
        dist_km           = "Distance (km)"
      )
    ) %>%
    filter(dbRDA1^2 + dbRDA2^2 > 1e-8)   # drop zero-length arrows

  p_global <- round(db_global$`Pr(>F)`[1], 3)

  fig_c1 <- ggplot(site_scores, aes(dbRDA1, dbRDA2, color = color)) +
    geom_point(size = 2.5, alpha = 0.75) +
    geom_segment(
      data = env_scores,
      aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
      inherit.aes = FALSE,
      arrow       = arrow(length = unit(0.18, "cm"), type = "closed"),
      color = "#444444", linewidth = 0.6
    ) +
    geom_text_repel(
      data = env_scores,
      aes(dbRDA1, dbRDA2, label = variable),
      inherit.aes = FALSE,
      size = 3, color = "#222222", seed = 42
    ) +
    scale_color_manual(
      values = c(R = col_red, W = col_white),
      labels = c(R = "Red", W = "White"),
      name   = "Light color"
    ) +
    labs(
      title    = paste0(taxon, " Community Composition — dbRDA"),
      subtitle = paste0("Year partialled (Condition) · Global p = ", p_global),
      x        = "dbRDA1", y = "dbRDA2",
      caption  = "Bray–Curtis distance · arrows = constrained predictors"
    ) +
    bat_theme +
    theme(legend.position = "top")

  # ── Fig C2: Activity heatmap (z-scored by species) ────────────────────────
  heat_dat <- data_env %>%
    group_by(species, intensity, color) %>%
    summarise(mean_det = mean(detections), .groups = "drop") %>%
    group_by(species) %>%
    mutate(mean_det_z = as.numeric(scale(mean_det))) %>%
    ungroup() %>%
    mutate(species = as.character(species))

  sp_order <- heat_dat %>%
    group_by(species) %>%
    summarise(mean_z = mean(mean_det_z, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_z) %>%
    pull(species)

  heat_dat <- heat_dat %>%
    mutate(species = factor(species, levels = sp_order))

  fig_c2 <- ggplot(heat_dat,
                   aes(x = factor(intensity), y = species, fill = mean_det_z)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient2(
      low      = col_red,
      mid      = "#F5F0E3",
      high     = col_white,
      midpoint = 0,
      name     = "Relative\nactivity (z)"
    ) +
    facet_wrap(~ color,
               labeller = labeller(color = c(R = "Red light", W = "White light"))) +
    labs(
      title   = paste0(unit, " Activity by Light Intensity and Color"),
      x       = "Light intensity (% of maximum)",
      y       = NULL,
      caption = paste0("Z-scored within ", tolower(unit), " across all treatment combinations")
    ) +
    bat_theme +
    theme(panel.grid = element_blank(),
          axis.line  = element_blank(),
          axis.ticks = element_blank())

  # ── Fig C3: Species × Distance heatmap (z-scored) ─────────────────────────
  dist_heat <- data_env %>%
    mutate(dist_cat = factor(dist_cat, levels = c("Close","Medium","Further","Far"))) %>%
    group_by(species, dist_cat) %>%
    summarise(mean_det = mean(detections), .groups = "drop") %>%
    group_by(species) %>%
    mutate(mean_det_z = as.numeric(scale(mean_det))) %>%
    ungroup() %>%
    mutate(species = factor(species, levels = sp_order))

  fig_c3 <- ggplot(dist_heat,
                   aes(x = dist_cat, y = species, fill = mean_det_z)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient2(
      low      = col_red,
      mid      = "#F5F0E3",
      high     = col_white,
      midpoint = 0,
      name     = "Relative\nactivity (z)"
    ) +
    labs(
      title   = paste0(unit, " Activity by Distance from Light Source"),
      x       = "Distance category",
      y       = NULL,
      caption = paste0("Z-scored within ", tolower(unit), " across distance categories")
    ) +
    bat_theme +
    theme(panel.grid = element_blank(),
          axis.line  = element_blank(),
          axis.ticks = element_blank())

  # ── Return ─────────────────────────────────────────────────────────────────
  list(
    comm_dat  = comm_dat,
    comm      = comm,
    env       = env,
    nmds      = nmds,
    perm      = perm,
    disp_test = disp_test,
    db        = db,
    db_global = db_global,
    db_margin = db_margin,
    mv_fit    = mv_fit,
    mv_anova  = mv_anova,
    fig_c1    = fig_c1,
    fig_c2    = fig_c2,
    fig_c3    = fig_c3
  )
}
