# Methods and Results — skeleton draft

> Statistical descriptions and
> reported values reflect the analyses in this repo; **[bracketed items]** are
> field/collection details or citations to fill in. Figure numbers refer to the
> eight composite panels (Panels 1–8).

---

# Methods

## Study design and data collection

We studied the effect of artificial light at night (ALAN) on bat and moth
activity along a distance gradient from the Colter Bay parking lot, Grand Teton
National Park, during the 2022 season (21 June – 16 August). Bats and
moths were sampled at **14 monitoring stations** spanning **< 0.52 km to
/> 2.1 km** (full range 0–7.46 km) from the lot. At each station the 
experimental light was displayed in a **red vs. white × randomized** 
**five-intensity (10, 30, 50, 70, 100 %)** design, rotated across nights for 
three day blocks. Three day white blocks were followed by three day red blocks 
at the same intensity. In this way, red and white lights were always 
alternating.

- **Bat activity** was recorded with **[sm4 microphones]** and classified with
  **[SonoBat v30]**; detections were retained at a **≥ 90 % classification
  probability** threshold **[confirm]**. Analyses focus on **eight focal
  species**, of which **six** (Epfu, Laci, Lano, Myev, Mylu, Myvo) met the
  minimum-data threshold for species-level models.
- **Moths (Lepidoptera)** were sampled with **UV bucket traps and passive 
  malaise traps** at the same 14 stations on the second and third day of the 
  light blocks; **9,864 specimens** across **241 identifications** 
  (225 site-nights) were identified to family/subfamily, spanning 
  **11 families and 31 subfamilies** (numerically dominated by Noctuidae, 
  Lasiocampidae, and Erebidae).
- **Visitor perceptions** were measured with a survey of **166 park visitors**,
  each exposed to a single light condition (**85 red, 79 white**), who rated 13
  streetlighting attitude statements on a 1–5 scale (1 = "not at all true") and
  a night-sky acceptability item (1–7). **[survey administration / sampling]**

## Environmental covariates

For every site-night we compiled: distance from the Colter Bay parking lot 
(`dist_km`, continuous), lunar illumination / phase (`mean_phase`, from 
**[moonlit / suncalc]**), local canopy openness (`pct_nonforest`, from 
**[NLCD 2021]** within a **50 m** buffer), site sky brightness 
(`brightness_dark`), and centered Julian day (`jd_c`, with a quadratic term for 
seasonality).

## Statistical analysis

### Bat activity models

Nightly per-species detections were modeled with a **zero-inflated negative
binomial GLMM** (`glmmTMB`, `nbinom2`, zero-inflation intercept), with a random
intercept for site:

```
detections ~ color * intensity + intensity * dist_km + mean_phase
           + pct_nonforest + brightness_dark + jd_c + I(jd_c^2) + (1 | site)
```

The two focal interactions were **color × intensity** and **intensity × **
**distance**. Estimated marginal means and contrasts obtained with `emmeans`
on the response scale, with **Holm** (color × intensity) or **Dunnett/Tukey**
(intensity) adjustment. Distance slopes per intensity were compared with
`emtrends`. (N ≈ **[3,457]** records, 14 sites.)

### Community composition

Site × species (or family) matrices were analyzed with **NMDS** (Bray–Curtis),
**PERMANOVA** (`vegan::adonis2`), **distance-based RDA** (`dbRDA`, marginal
tests), and a **multivariate GLM** (`mvabund`, negative binomial, PIT-trap,
`nBoot = 9999`) to attribute community turnover to the treatment and
environmental terms. Homogeneity of multivariate dispersion between colors was
checked with `betadisper`.

### Species-level and morphological analyses

Species meeting a **30-detection / 20-row** threshold were modeled individually
with the same GLMM structure. Univariate treatment effects per species were
also drawn from the `mvabund` fit (adjusted p-values). Four continuous
wing-morphology traits (**ear height, forearm length, wing loading, aspect
ratio**) were z-scored and interacted with intensity and distance, then compared
at low/mean/high (± 1 SD) trait values.

### Moth (Lepidoptera) analyses

Total nightly moth counts were modeled with a GLMM **mirroring the bat model**,
using the identical fixed-effect structure and random intercept for site:

```
detections ~ color * intensity + intensity * dist_km + mean_phase
           + pct_nonforest + brightness_dark + jd_c + I(jd_c^2) + (1 | site)
```

Six families passing the same threshold were modeled individually, and the
family community was analyzed with the same NMDS/PERMANOVA/dbRDA/`mvabund`
pipeline. (N = 226 site-nights, 14 sites.)

### Moth taxonomic resolution and richness

Specimens (9,864, identified to species where possible) were additionally
analyzed at **subfamily (30 taxa), genus (159), and species (253)** resolution.
At each level the community was tested for color, intensity, and distance
effects with the same NMDS/PERMANOVA/dbRDA/`mvabund` pipeline (`mvabund`
bootstraps reduced to 999 at the genus and species levels for tractability).
**Taxonomic richness** — the number of distinct species, genera, and subfamilies
per site-night — was modeled with negative binomial (Poisson where the negative
binomial was unstable) GLMMs of the same form as the abundance model, with
**log catch size** added as a covariate so treatment effects reflect richness
beyond that expected from catch volume.

### Bat–moth comparison

To compare the **shape** of each taxon's response independent of absolute
abundance, marginal predictions were rescaled to **% change** from each taxon's
own baseline (10 % intensity; nearest shared site for distance), and effect
sizes (rate ratios ± 95 % CI) were calculated on a shared scale.

### Visitor perception survey

The 13 attitude items were reduced with a **principal components analysis**
(varimax rotation). We assessed suitability with **Bartlett's test of
sphericity** and the **Kaiser–Meyer–Olkin (KMO)** statistic, retained components
with **eigenvalue > 1**, assigned items at a minimum loading of **≥ 0.40**,
eliminated **cross-loaded** items, and required **Cronbach's α ≥ 0.65** for a
component to be retained; reliable components were averaged into a single score
(negatively worded items reverse-keyed). Missing responses (coded 999) were
median-imputed. In parallel, an exploratory factor analysis with **k-means**
clustering segmented respondents, and a **binomial GLM** tested whether the
component/factor scores discriminated the (experimentally assigned) light
condition. 

### Software

All analyses used **R [4.6.0]** with `glmmTMB`, `emmeans`, `vegan`, `mvabund`,
`sf`/`exactextractr`, `psych`, `GPArotation`, and `cluster` **[+ versions]**.

---

# Results

## Bat activity: light color × intensity (Panels 1–2)

Bat activity depended on the **interaction of color and intensity**. Red light
produced **2.06× more detections than white at 30 %** (p < 0.0001) and **1.63×
at 50 %** (p < 0.0001); at 70 % and 100 % the colors were statistically
indistinguishable (Panel 2A). Within red light, activity peaked at 30 %
(1.83× the 10 % baseline, p = 0.0001) and collapsed at 100 % (30 % vs 100 %
ratio 1.77, p < 0.0001), whereas **white was flat across all intensities** (no
contrast survived correction). Pooled across color, only the 70 % level
differed significantly from the 10 % baseline (1.28×, p = 0.002; Panel 1A),
because averaging over the flat white response dilutes the red-specific effect.

## Bat activity: intensity × distance (Panel 2B)

The **intensity × distance interaction** was driven by the two brightest levels:
the distance slopes at 70 % (β = −0.128, p < 0.001) and 100 % (β = −0.090,
p = 0.001) were significantly shallower than at 10 %. At 0.3 km, 70 % produced
**1.66× more detections than 10 %**; this gap narrowed but persisted at 2.0 km
(ratio 0.748, p = 0.0003). High intensities thus **suppress bats uniformly
across distance**, not just near the Colter Bay parking lot.

## Community composition (Panel 3)

Light **color did not restructure the bat community** (dbRDA p = 0.193,
`mvabund` p = 0.415). **Distance was the dominant structuring variable**
(dbRDA F = 7.22, p < 0.001; `mvabund` Dev = 142.6, p < 0.001), with intensity
marginal and habitat/sky brightness contributing (PERMANOVA R² = 0.156,
F = 2.51, p < 0.001; NMDS stress 0.207). Dispersion was homogeneous between
colors (p = 0.938).

## Species-level and morphological results (Panel 3)

Species differed in sensitivity: **Epfu** responded to intensity, distance, and
moon phase; **Myev** and **Myvo** tracked habitat/sky brightness and distance;
color was non-significant for every species. Across the shared wing-morphology
gradient, **smaller, slower, more maneuverable species (Myotis)** increased in
activity with distance from the lot, whereas larger, faster, open-air fliers
(e.g., Laci) were comparatively flat (Panel 3A) — a single morphological
gradient viewed four ways rather than four independent effects.

## Moth responses (Panels 4–6)

Moths showed **no color preference at low-to-mid intensity**, but at **100 %,
white produced 3.1× more detections than red** (p < 0.001; Panel 5A) — the
mirror image of the bats, driven by red-avoidance at high intensity rather than
a white attraction peak. Moth activity **increased with distance** from the lot
(β = 0.096, p = 0.020), a positive slope at 10–70 % intensity that flattened to
essentially zero at 100 % (intensity100 × dist_km = −0.111, p = 0.035). Adding
**site sky brightness** to the model (matching the bat model) did **not** change
these conclusions: `brightness_dark` was a **positive but non-significant**
predictor of moth activity (β = 0.015, p = 0.18), and moon phase was retained as
a negative predictor (β = −0.47, p = 0.049). At the family level, Noctuidae
dominated numerically and, with Lasiocampidae and Geometridae, rose with
distance at mid-to-high intensity (Panel 6). The family community was strongly
structured by intensity, distance, and moon phase (all dbRDA p < 0.001;
`mvabund` p < 2e-16) but **not color** (p = 0.159), and treatment explained
proportionally **more** of the moth community's variation (R² = 0.305) than the
bats' (R² = 0.156).

## Moth taxonomic resolution and richness (Panel 9)

The family-level pattern held at finer resolution. At **subfamily** level the
community was significantly structured by **intensity** (dbRDA F = 6.05,
p = 0.01), **distance** (F = 6.23, p = 0.01), and **moon phase** (F = 3.97,
p = 0.01) but **not by light color** (F = 1.57, p = 0.17). **Distance remained a
strong structuring variable at the genus and species levels** (dbRDA marginal
F = 3.48, p = 0.001 and F = 3.49, p = 0.001, respectively). This reflects a
broad, community-wide increase in activity away from the parking lot: of the
genera occurring at ≥ 6 sites, **26 increased with distance vs. 10 decreased**
(species: 24 vs. 9), and every abundant genus rose with distance — *Malacosoma*
(Spearman ρ = 0.72), *Leucania* (0.63), Noctuinae (0.57), *Euxoa* (0.51),
*Egira* (0.45), *Apamea* (0.42) — as did the dominant species (*Malacosoma
californica* ρ = 0.74, *Egira curialis* 0.54, *Apamea unanimis* 0.43).

**Taxonomic richness tracked the abundance story.** Adjusting for catch size,
**red light supported higher richness than white at 30–50 % intensity**:
species richness was 1.36× higher at 30 % (p = 0.03) and 1.48× at 50 %
(p = 0.01), and genus richness 1.34× at 50 % (p = 0.02). At **100 % intensity
the pattern reversed** — white supported higher richness (species 0.75×,
p = 0.06; genus 0.77×, p = 0.03; subfamily 0.72×, p = 0.04) — mirroring the
100 % white-over-red abundance effect. So the color × intensity effect on moths
is one of **diversity as well as abundance**: dim-to-moderate red light both
attracts more moths and draws a richer assemblage, while the brightest white
light does so instead (Panel 9).

## Bat vs. moth (Panel 7)

The taxa **diverge most clearly on color**: averaged over intensity, bats show
**22 % fewer** detections under white (CI excludes 1) while moths show **24 %
more** (CI excludes 1). At 100 % specifically, bat activity is indistinguishable
from baseline (ratio 1.03, ns) while moth activity drops (ratio 0.68, CI
excludes 1). **Both taxa increase with distance** over the shared range, with
the relative increase larger for moths (+125 % vs +80 % at the far end).

## Visitor perceptions (Panel 8)

The attitude items were suitable for reduction (**KMO = 0.80**, Bartlett
p < 0.001) and yielded **two reliable components**: a broad **wildlife-plus-
experience benefit** dimension (PC1, α = 0.90; loading the three wildlife items
plus eye-transition comfort and activity enjoyment) and a **discomfort-in-the-
dark** dimension (PC2, α = 0.70). Attitudes tracked the light experienced:
**red-light visitors scored much higher on PC1** (3.82 vs 2.59) and were more
accepting of dark skies (5.43 vs 3.88 on the 1–7 scale), while white-light
visitors reported slightly more discomfort. A logistic model discriminated the
two exposure groups at **AUC = 0.80**, driven almost entirely by PC1 (OR = 0.36
per SD, p < 0.001). Because attitudes were measured after exposure, this is an
**association** — red-light exposure and pro-wildlife/dark-sky attitudes travel
together — rather than evidence that prior attitudes drove exposure.

---

## Reporting notes / to-do
- Confirm N for the bat model (documented **3,457**; current cached data ≈ 3,378).
- Fill field-collection details, instrument models, and citations in **[brackets]**.
- Panels 3A and 6A are **model-based** (per-unit GLMM marginal predictions ± 95 % CI);
  the supplemental `fig_s4a` figures are the color-based counterparts.
