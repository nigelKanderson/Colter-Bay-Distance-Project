# 24_lamp_visual_model.R
# Models how a moth and a bat (genus Myotis, one of this project's own focal
# genera) perceive the red and white experimental lights, using measured
# lamp emission spectra and published photoreceptor sensitivity peaks.
#
# Visual system sources:
#   Moth (Manduca sexta rhodopsins P357/P450/P520, trichromatic):
#     UV 357 nm, Blue 450 nm, Green 520 nm
#     White RH, Xu H, Munch TA, Bennett RR, Grable EA (2003) "The retina of
#     Manduca sexta: rhodopsin expression, the mosaic of green-, blue- and
#     UV-sensitive photoreceptors, and regional specialization." Journal of
#     Experimental Biology 206:3337-3348.
#   Bat -- genus Myotis (this project's own focal genus), dichromatic:
#     S/UV ~360-365 nm, L/green 558 nm
#     L-opsin lambda-max of 558 nm for Myotis velifer is a sequence-based
#     spectral tuning prediction from:
#       Zhao H, Rossiter SJ, Teeling EC, Li C, Cotton JA, Zhang S (2009)
#       "The evolution of color vision in nocturnal mammals." PNAS
#       106(22):8980-8985.
#     UV-tuned S-cone sensitivity corroborated physiologically (though in
#     phyllostomid bats, not Myotis) by:
#       Muller B, Glosmann M, Peichl L, Knop GC, Hagemann C, Ammermuller J
#       (2009) "Bat Eyes Have Ultraviolet-Sensitive Cone Photoreceptors."
#       PLoS ONE 4(7):e6390 -- their own ERG action spectrum peaks at or
#       below 365 nm (UV) with a secondary maximum near 450 nm; L cone is
#       only bounded to 530-560 nm by their data, consistent with the
#       558 nm Myotis velifer estimate used here. Note: 558 nm is well
#       above the ~510 nm L-opsin value typical of UV-tuned *rodents* --
#       bats are longwave-shifted relative to that pattern, per Muller
#       et al.'s discussion.
#
# Photoreceptor sensitivity curves are generated with pavo::sensmodel(),
# which implements the Govardovskii et al. (2000) visual pigment template.
# Quantum catch is computed manually: measured radiance is converted to
# relative photon flux (photon energy is inversely proportional to
# wavelength, so photon flux is proportional to radiance x wavelength)
# before integrating against each receptor's sensitivity curve -- the
# standard approach in visual ecology, since photoreceptor response scales
# with photon count, not radiometric power.
#
# Input:  data/lamp_spectra/LightMeasurementsLong.csv
# Output: list(catch_long, fig_moth, fig_bat, fig_moth_triangle)

library(tidyverse)
suppressWarnings(library(pavo))

run_lamp_visual_model <- function(spectra_csv = "data/lamp_spectra/LightMeasurementsLong.csv",
                                   bat_theme, col_red, col_white) {

  # ── 1. Load + reshape lamp spectra into pavo's rspec format ─────────────
  # Noise floor: measured across three independent regions expected to carry
  # no real signal (Red light 400-550nm, White light 350-390nm, both colors'
  # 700-750nm tail), the residual sensor noise maxes out around
  # 0.0004-0.00044 W/sr/m2. At low intensity (esp. Red 10%) this noise is a
  # large enough fraction of the very weak true signal to visibly distort
  # relative quantum catch. 0.0005 sits just above that noise ceiling while
  # staying below every real spectral feature checked, including the
  # weakest ones (White 10%'s dimmest peak shoulders ~0.0004-0.0014,
  # Red 10%'s real peak ~0.0007-0.012).
  noise_floor <- 0.0005

  lamp_long <- read_csv(spectra_csv, show_col_types = FALSE) %>%
    mutate(`Watts/sr/m2` = if_else(`Watts/sr/m2` < noise_floor, 0, `Watts/sr/m2`))

  lamp_wide <- lamp_long %>%
    select(wl, LightSet, `Watts/sr/m2`) %>%
    pivot_wider(names_from = LightSet, values_from = `Watts/sr/m2`) %>%
    arrange(wl)

  lamp_rspec <- as.rspec(lamp_wide, whichwl = "wl")
  wl <- lamp_rspec$wl
  light_cols <- setdiff(names(lamp_rspec), "wl")

  # ── 2. Photoreceptor sensitivity curves (Govardovskii template) ─────────
  moth_sens <- sensmodel(peaksens = c(357, 450, 520), range = range(wl))
  names(moth_sens) <- c("wl", "moth_UV", "moth_Blue", "moth_Green")

  bat_sens <- sensmodel(peaksens = c(360, 558), range = range(wl))
  names(bat_sens) <- c("wl", "bat_UV", "bat_Green")

  # ── 3. Quantum catch per receptor per light condition ────────────────────
  # Convert radiance -> relative photon flux (photon energy ~ 1/wavelength)
  photon_flux <- as.matrix(lamp_rspec[, light_cols]) * wl

  catch_from <- function(sens_df, receptor_cols) {
    sens_mat <- as.matrix(sens_df[, receptor_cols])
    # trapezoidal-equivalent: 1 nm spacing, so a plain sum is proportional
    # to the integral; only relative catch matters here
    t(sens_mat) %*% photon_flux
  }

  moth_catch <- catch_from(moth_sens, c("moth_UV", "moth_Blue", "moth_Green"))
  bat_catch  <- catch_from(bat_sens,  c("bat_UV", "bat_Green"))

  to_long <- function(catch_mat, species) {
    as.data.frame(catch_mat) %>%
      rownames_to_column("receptor") %>%
      pivot_longer(-receptor, names_to = "light_set", values_to = "catch") %>%
      mutate(
        color     = if_else(str_starts(light_set, "Red"), "R", "W"),
        intensity = as.numeric(str_extract(light_set, "\\d+")),
        species   = species
      )
  }

  catch_long <- bind_rows(
    to_long(moth_catch, "Moth"),
    to_long(bat_catch,  "Bat")
  ) %>%
    group_by(species, light_set) %>%
    mutate(catch_rel = catch / sum(catch)) %>%  # relative catch within each light condition
    ungroup()

  # ── 4. Figures ────────────────────────────────────────────────────────────
  moth_df <- catch_long %>%
    filter(species == "Moth") %>%
    mutate(
      receptor  = factor(str_remove(receptor, "moth_"), levels = c("UV", "Blue", "Green")),
      intensity = factor(intensity, levels = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
    )

  fig_moth <- ggplot(moth_df, aes(x = intensity, y = catch_rel, color = color, group = color)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_wrap(~ receptor, scales = "free_y") +
    scale_color_manual(values = c(R = col_red, W = col_white),
                       labels = c(R = "Red", W = "White"), name = "Light color") +
    labs(
      title    = "Moth Photoreceptor Catch by Light Color and Intensity",
      subtitle = "Relative quantum catch per receptor (UV 357 / Blue 450 / Green 520 nm) — White et al. 2003",
      x        = "Light intensity (%)",
      y        = "Relative quantum catch"
    ) +
    bat_theme +
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

  bat_df <- catch_long %>%
    filter(species == "Bat") %>%
    mutate(
      receptor  = factor(str_remove(receptor, "bat_"), levels = c("UV", "Green")),
      intensity = factor(intensity, levels = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
    )

  fig_bat <- ggplot(bat_df, aes(x = intensity, y = catch_rel, color = color, group = color)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_wrap(~ receptor, scales = "free_y") +
    scale_color_manual(values = c(R = col_red, W = col_white),
                       labels = c(R = "Red", W = "White"), name = "Light color") +
    labs(
      title    = "Bat (Myotis) Photoreceptor Catch by Light Color and Intensity",
      subtitle = "Relative quantum catch per receptor (S/UV 360 / L/Green 558 nm) — Zhao et al. 2009",
      x        = "Light intensity (%)",
      y        = "Relative quantum catch"
    ) +
    bat_theme +
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

  # ── 5. Absolute (non-normalized) catch: the brightness dimension ─────────
  # catch_rel above shows chromatic balance (invariant to pure brightness
  # scaling); these companion figures show raw quantum catch, i.e. how much
  # total signal each receptor actually gets, which does scale with
  # intensity for every color including red. Units are arbitrary (relative
  # photon flux, uncalibrated) but consistent across all light conditions
  # and both species, so ratios between conditions are meaningful.
  fig_moth_abs <- ggplot(moth_df, aes(x = intensity, y = catch, color = color, group = color)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_wrap(~ receptor, scales = "free_y") +
    scale_y_log10() +
    scale_color_manual(values = c(R = col_red, W = col_white),
                       labels = c(R = "Red", W = "White"), name = "Light color") +
    labs(
      title    = "Moth Photoreceptor Catch (Absolute, Log Scale) by Intensity",
      subtitle = "Raw quantum catch per receptor -- the brightness dimension, not chromatic balance",
      x        = "Light intensity (%)",
      y        = "Quantum catch (log scale, arbitrary units)"
    ) +
    bat_theme +
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

  fig_bat_abs <- ggplot(bat_df, aes(x = intensity, y = catch, color = color, group = color)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_wrap(~ receptor, scales = "free_y") +
    scale_y_log10() +
    scale_color_manual(values = c(R = col_red, W = col_white),
                       labels = c(R = "Red", W = "White"), name = "Light color") +
    labs(
      title    = "Bat (Myotis) Photoreceptor Catch (Absolute, Log Scale) by Intensity",
      subtitle = "Raw quantum catch per receptor -- the brightness dimension, not chromatic balance",
      x        = "Light intensity (%)",
      y        = "Quantum catch (log scale, arbitrary units)"
    ) +
    bat_theme +
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

  dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
  ggsave("output/figures/fig_moth_photoreceptor_catch.png", fig_moth,
         width = 10, height = 5, dpi = 300, bg = "white")
  ggsave("output/figures/fig_bat_photoreceptor_catch.png", fig_bat,
         width = 8, height = 5, dpi = 300, bg = "white")
  ggsave("output/figures/fig_moth_photoreceptor_catch_absolute.png", fig_moth_abs,
         width = 10, height = 5, dpi = 300, bg = "white")
  ggsave("output/figures/fig_bat_photoreceptor_catch_absolute.png", fig_bat_abs,
         width = 8, height = 5, dpi = 300, bg = "white")

  list(
    lamp_rspec  = lamp_rspec,
    moth_sens   = moth_sens,
    bat_sens    = bat_sens,
    catch_long  = catch_long,
    fig_moth    = fig_moth,
    fig_bat     = fig_bat,
    fig_moth_abs = fig_moth_abs,
    fig_bat_abs  = fig_bat_abs
  )
}
