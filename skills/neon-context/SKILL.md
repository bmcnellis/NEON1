---
name: neon-context
description: "Project-specific -- NEON PUUM only; do not apply to other projects unless explicitly asked. Project context for the PUUM (Pu'u Maka'ala, Hawai'i) NEON plant community analysis: data root layout, the NEON1 R package and how to run it, the NEON API token and the cached data_NEON downloads, the two Hmsc joint species distribution models refit on 2026-10-19 and how the refreshed data changed the trait-interaction result, manuscript state (MS2), and the package-drift traps that fail silently. Load before touching this project's data, models, figures, or manuscript text."
scope: project-specific
project: NEON PUUM
---

# NEON-context — PUUM plant community project

> **Scope — project-specific: the NEON PUUM project.** Everything here is specific to
> that project. Do not carry its paths, file names, conventions, model choices or
> findings into another project's work unless the user explicitly asks for it.

Working context for Brandon McNellis's analysis of NEON plant community data at
**PUUM** (Pu'u Maka'ala Natural Area Reserve, Hawai'i Island; NEON domain D20,
the only site in the Pacific Tropical Domain). Everything below was read off the
files on disk as of **2026-10-19**; the manuscript is `MS2` v0_5 (2026-10-19) and
the fitted models are from **2026-10-19**, fitted to NEON data current to
2025-06-24. The 2024-09-21/22 fits no longer exist on disk — see §7.

## 1. The scientific question

Spatial patterns of plant community abundance and composition at PUUM — a
first site-level analysis of NEON's PUUM data against the landscape factors
unique to the site. Four framing questions from the manuscript:

1. How are plant abundance, composition and structure affected by **disturbance
   legacy**? (pig fences, lava flow age, historic cattle, historic logging)
2. How are they affected by the **environmental gradient** within PUUM?
   (elevation, topography)
3. How are they affected by **biotic interactions**? (PAI, residual species
   correlations)
4. Can anything be explained by **traits**? (native vs. invasive status)

Target outlet: *Landscape Ecology*, Original Research, 8,500 words excl.
abstract + references. Authors: McNellis, Phillips, Long, Donoso, Uehana,
[PAI coauthor TBD], Cordell (USDA PSW/IPIF Hilo; USGS PIERC; NEON/Battelle).

Site facts as written in the manuscript: 7,580-ha reserve, set aside 1981,
managed by HI-DLNR; elevation 853–1899 m; rainfall 200–400 cm; overstory
*Metrosideros polymorpha*, understory *Cibotium* spp., *Carex alligata* in wet
areas, *Acacia koa* emergent. Kūlani Correctional Facility (113 ha) sits inside
the reserve boundary. Cattle grazed some tracts until 2005; pigs excluded
1991–2003 from 3,885 ha across 9 units.

## 2. Where things live

Data root (granted host path, rw): **`/media/bem/data/NEON_claude`**

| Path | What |
|---|---|
| `NEON1/` | the analysis R package + git repo (branch `main`) |
| `NEON1/R/` | package functions: `helpers.R`, `plot_functions.R`, `keys.R`, `external.R` |
| `NEON1/exec/` | 16 analysis scripts (the actual pipeline) |
| `NEON1/data/` | 15 `.rda` package data objects |
| `data_NEON/` | cached raw NEON downloads: `d0`-`d4.rda` plus `download_log.txt`, written by `import_NEON.R` |
| `rlib/` | R library for the CRAN-only packages (MCMCvis, neonUtilities, Rsagacmd) and the installed `NEON1`. Keep it; rebuilding costs ~20 min |
| `data/`, `meta/` | field data (DHP csv, log/graze xlsx, fence unit zips) and NEON site metadata |
| `results/model_results/` | fitted Hmsc models + diagnostics + trace PDFs |
| `results/figures/` | current manuscript figures (9 files) + `archive/` |
| `results/reports/` | DHARMa residual-check PDFs, fencing csv/zip |
| `spatial/` | AOP elevation mosaics, `AOP/DP3.30024.001/` LiDAR DSM+DTM tiles (416 each, 1-km UTM 5N grid), `NEON_Shapefiles/`, `fencelines_modified/` |
| `provided/` | NEON-supplied docs and shapefiles (raw copies of `spatial/NEON_Shapefiles`); `NEON Docs/PlotCharacterization/Docs/` holds the DLNR management plan and cultural study |
| `references/` | literature: Konig 2021 + Dryad zip, gllvm/niku, the DLNR management plan, Maly and Maly (2004) cultural study, Mecyclothorax site descriptions, substrate map |
| `NEON1/skills/` | copy of this skill, kept in sync by hand |
| `<skills>/neon-context/references/` | Konig et al. 2021 *J. Biogeogr.* JSDM scale-dependency paper + its Dryad appendix (doi:10.5061/dryad.qfttdz0g7) - the methodological template |
| `MS1/`, `MS2/` | manuscript drafts with `previous_versions/` |

`MS1` is a separate, earlier manuscript (NEON cross-site/citation framing, its
own `data/`); **`MS2/NEON1_MS2_0_4.docx` is the live project manuscript.**

## 3. The NEON1 package

`Package: NEON1`, v0.2.0, R >= 4.1.2. Imports dplyr, tidyr, ggplot2, viridis,
terra, Rsagacmd, neonUtilities, Hmsc, coda, corrplot. Scripts also use
`sf`, `MCMCvis`, `bayesplot`, `DHARMa`, `glmmTMB`, `vegan`, `gllvm`,
`neonPlantEcology`, `pairwiseAdonis`, `readxl`, `usethis`.

Exported helpers worth knowing:

- `logit_Warton(y, epsilon)` — Warton & Hui (2011) logit for proportional cover
- `div_to_long(div, type, add_zeros)` — NEON div product to long form, zero-fill
- `posterior_from_coda(mcmc, term, prob, average, drop_ns)` — posterior means +
  HPD intervals from a coda object; `drop_ns = T` drops terms whose interval
  spans zero. This is what generates the `post_beta_*`/`post_gamm_*`/`post_omeg_*`
  tables the figures are built from.
- `param_plot(post, ylab)`, `plot_vp(hM, VP)`, `ess_plot`, `plot_resid`,
  `pp_check(m_p, mp_pp, n_draws)` — figure functions
- `check_residuals(fit, dir, grp, nm)` — DHARMa wrapper that writes the
  `results/reports/check_residuals_*.pdf` files
- `key_age_range()`, `relevel_age()`, `key_field_domain_id()`, `key_spp()`,
  `sph()` (splits NEON "value +/- error" strings), `order_to_class()`,
  `c25(option)` (colour palette), `HighResTiff()`, `Multiplot()`

Package data objects (`NEON1::<name>`), 15 as of 2026-10-19: `div_one`,
`div_ten`, `spp`, `met`, `dhp`, `str`, `soil_init`, `soil_peri`, `lus`, `topo`,
`flow_meta`, `fence`, `mve`, `incl_plots`, `NEON_aop_crs`. (`str.rda` and
`soil_peri.rda` did not exist before the 2026-10-19 run of `import_data.R`.)

Not every declared Import is installable here: `gllvm`, `neonPlantEcology`,
`pairwiseAdonis` and `import` are unavailable from the configured channels, so
`analyze_NMDS.R` and `analyze_climate.R` cannot be run. The reproducible stack
is a conda R env (`r-base`, `r-hmsc`, `r-coda`, `r-dplyr`, `r-tidyr`,
`r-tibble`, `r-ggplot2`, `r-viridis`, `r-corrplot`, `r-terra`, `r-sf`,
`r-dharma`, `r-vegan`, `r-readxl`, `r-usethis`, `r-devtools`, `r-bayesplot`,
`r-rstan`, `r-arrow`) plus `install.packages(c('MCMCvis','neonUtilities','Rsagacmd'))`
into `rlib/`. MCMCvis needs rstan and neonUtilities needs arrow; both fail to
compile from source and must come from conda first.

## 4. Data sources

| Variable | Role | Source |
|---|---|---|
| Plant presence + percent cover | response | NEON **DP1.10058.001** |
| Vegetation structure | response (not yet modelled) | NEON **DP1.10098.001** |
| Elevation / LiDAR (DTM → elev, TWI) | fixed | NEON **DP3.30024.001** |
| Soil biogeochem / periodic soil chem | fixed (imported, unused in fits) | NEON **DP1.10047.001**, **DP1.10086.001** |
| Foliar traits | intended random effect (unresolved) | NEON **DP1.10026.001** |
| PAI from digital hemispherical photos | fixed | NEON DHP → `data/neon_dhp_puum.csv` |
| Lava flow age | fixed | USGS flow-age polygons (`spatial/NEON_Shapefiles/Flow Map/Geological_Units.shp`) |
| Pig exclosures / fence units | fixed | HI-DLNR `PuuMakaala_Fence_Units_20240517.zip` |
| Cover type, logging, grazing history | fixed | NEON-provided `data/puum plots log_graze.xlsx` |

**NEON requires an API token for every download as of June 2026.** The token is
at `NEON1/exec/NEON_TOKEN.R` (a single `NEON_TOKEN <- "..."` line), covered by
`.gitignore` via `exec/NEON_TOKEN.R` and `*_TOKEN.R`. It was never committed.
Only `import_NEON.R` reads it. A sandboxed agent cannot complete a download:
the metadata API at `data.neonscience.org` is reachable and the token
authenticates (rate limit 200 -> 2000), but the files come from signed
path-style URLs on `storage.googleapis.com`, which is permanently denylisted,
and the signature covers the host header so the bucket-qualified form does not
help. **BEM has to run `import_NEON.R`.**

## 5. Pipeline order

```
import_NEON.R     # the only script that touches the network. loadByProduct for
                  #   DP1.10058/10098/10047/10086/10026 -> data_NEON/d0-d4.rda
                  #   + download_log.txt. RUN FROM NEON1/exec/ (see §8)
import_data.R     # reads data_NEON/d0-d4.rda + DHP csv + log/graze xlsx
                  #   -> usethis::use_data(...). No network, runs in ~2 s
import_topo.R     # per-tile SAGA (Rsagacmd): flow accumulation -> SCA -> slope
                  #   -> TWI, extract at plot points -> topo.rda   [BLOCKED, §8]
import_fence.R    # fence unit shapefile -> enclosure name/ha/completion -> fence.rda
make_map_figures.R # Figures 1 and 2; ALSO builds flow_meta (plot -> flow age),
                  #   though its use_data call is commented out
make_models.R     # <- THE model script. Builds matrices, fits both Hmsc models,
                  #   writes m_p_mod{0,1}.rda, m_p_diag{0,1}.rda, rhat_*.csv, traces
eval_models.R     # diagnostics on model 0; writes incl_plots.rda and mve.rda
make_figures.R    # Figures 4-8 from m_p_mod0/m_p_diag0
do_model_xv.R     # k-fold CV (2, 5, 10) -> m_p_xv{0,1}_{1..3}.rda  [NEVER RUN]
```

**Reinstall `NEON1` after every step that rewrites package data**, or the next
script reads stale objects: after `import_data.R`, after `import_fence.R`, and
after `eval_models.R`. `R CMD INSTALL --library=/media/bem/data/NEON_claude/rlib NEON1`
takes ~20 s.

To run any of it outside RStudio:

```bash
export R_LIBS_USER=/media/bem/data/NEON_claude/rlib
export PROJ_DATA=<conda-env>/share/proj   # terra/sf fail with "empty srs" without
export PROJ_LIB=$PROJ_DATA                #   these; the conda PROJ data dir is not
export GDAL_DATA=<conda-env>/share/gdal   #   found automatically
cd /media/bem/data/NEON_claude/NEON1      # the package root, NOT exec/
Rscript exec/<script>.R
```

Side/earlier analyses, not part of the current manuscript path:
`analyze_Hmsc.R` (superseded by `make_models.R` - see §8),
`analyze_NMDS.R` (unconstrained ordination + PERMANOVA across PUUM/GUAN/WREF;
PUUM stress 0.1367 at k=3, no significant plot or year pairs),
`analyze_climate.R` (climatograph from DP1.00006.001 / DP1.00002.001 - becomes
Figure S1), `analyze_structure.R` (vegetation structure, incomplete),
`process_UNK_species.R` (errors: references `NEON1::PUUM_div`, which no longer
exists), `extra.R` (scratch; errors on an undefined `plant_df`),
`make_figures_old.R` (superseded).

## 6. The fitted models

Two Hmsc probit joint species distribution models, identical structure, differing
only in the response data. Both were refit on **2026-10-19** against NEON data
current to 2025-06-24; the numbers below are that fit.

- **Model 0** (`m_p_mod0.rda`, 4.23 GB) - 1-m2 percent-cover data (`div_one`),
  binarised at cover > 0.0001. **188 plot-dates x 35 species.** This is the
  model the manuscript figures and diagnostics use.
- **Model 1** (`m_p_mod1.rda`, 9.21 GB) - 10-m2/100-m2 presence data (`div_ten`,
  union with `div_one` records on shared dates). **116 plot-dates x 94 species.**
  Fitted and now converged, but still not evaluated or figured.

Specification (from `make_models.R`, unchanged since 2024):

```r
XFormula  = ~ age_median + pai + time_since_fence + elevation + log + cow + cover_type
TrFormula = ~ nativeStatusCode            # native ('N') vs invasive ('z_I')
ranLevels = list(plotDate = HmscRandomLevel(units = stu$plotDate))
distr = 'probit', XScale = F             # predictors pre-scaled by hand
sampleMcmc(thin = 10, samples = 8000, transient = 4000, nChains = 6, nParallel = 6)
set.seed(1)
```

Row unit is **plotDate** (`<plot number>_<endDate>`, `PUUM_` stripped); the
plot-date random level is what yields the residual species associations
(`Omega`). Hmsc parameter meanings as used here: `Beta` = species niches,
`Gamma` = trait effects on niches, `Omega` = residual species associations.

Objects inside each file:

- `m_p_mod{0,1}.rda` -> `m_p_{0,1}`, `mc_p_*` (coda), `ma_p_*` (associations),
  `mp_p_*` (predicted values), `mp_pp_*` (`predict()` draws)
- `m_p_diag{0,1}.rda` -> `mf_p_*` (model fit), `m_ca_*`, `m_vp_*` (variance
  partitioning), `mc_s_{beta,gamm,omeg}_*` (MCMCvis summaries),
  `post_{beta,gamm,omeg}_*` (90% HPD tables)

Wall time for `make_models.R` is ~56 min for both models on 16 cores.
`m_p_mod1.rda` needs several GB of RAM to load; `do_model_xv.R` works around
this by `load()`ing then immediately `rm()`ing everything except the model
object and calling `gc()`. Use the same pattern. For diagnostics, prefer the
`m_p_diag*.rda` files - they are under 600 KB and carry everything the figures
need.

## 7. Where the analysis stands

**Model selection was a 10-round convergence hunt** documented in the comment
block at the top of `make_models.R` - do not delete those comments, they are the
record of which species were dropped and why. Model 0 drops
`Hypochaeris_radicata`, `Hymenophyllum_recurvum` and `Axonopus_fissifolius` (all
weedy non-natives with bad Rhat) at occurrence threshold `oc >= 15`; model 1
drops nothing at `oc >= 5`. **Those decisions were made against 2023-era data
and have not been revisited** against the refreshed data.

### Manuscript

`MS2/NEON1_MS2_0_5.docx` (2026-10-19) is the live draft and the first with a
complete Introduction, Methods, Results and Discussion. The Abstract and
Figure S1 are still missing, some `REF` placeholders remain in the
Introduction, and the author list still reads `[PAI guy]`. **Before editing the
manuscript, renumbering a figure or adding a citation, read
`reference/manuscript_state.md`** - it holds the section-by-section state,
the figure-to-script table, the word counts and the list of citations added in
v0_5 that are not yet in `MS2/references.docx`.

Two things settled in v0_5 that should not be reopened: the Introduction's
numbered question list carries **no citations** (BEM's instruction), and the
statistical-analysis text now matches `make_models.R` rather than the generic
probit model v0_4 described.

### Site facts and the fence dates

The study-system numbers were checked against DLNR (2013), which is on disk.
**Before writing or checking a study-system sentence, citing DLNR (2013), or
touching the ungulate exclusion dates, read `reference/site_facts.md`.** The
headline corrections: the reserve is **7,570 ha**, not the 7,580 ha carried
through v0_4; the 2005 cattle removal applies to the adjoining former pasture,
not to tracts inside the reserve; and the 1991-2003 exclusion window is **not**
in conflict with the 1992-2019 completion years in the fence layer, because the
plan gives no per-unit dates and lists the two newest units as proposed.

### What the refreshed data changed

The models on disk were refit on 2026-10-19 against NEON data current to
2025-06-24, and the numbers moved. The one that matters: in the variance
explained by each covariate's interaction with the invasive trait, historic
grazing stays first but loses more than half its share, and years ungulate-free
roughly doubles and moves from seventh to second. **Any manuscript sentence
quoting 41.45 / 26.43 / 11.30 is stale, and the ordering it argues for has
changed.** **Before quoting any fit, variance-partitioning or trait-R2 number
that did not come from the current `m_p_diag*.rda`, read
`reference/model_history.md`** for the full 2024-vs-2026 comparison.

Also: the `cover_type` Beta coefficients are only partly identified. For a
species with no positive observations in a cover type the probit likelihood is
flat in the negative direction, so the posterior mean is set by the prior and
the chain's path. Across two fits of the same 2024 data, *Sticherus owhyensis*
in `koa_tall` moved from -1.10 to -15.20. Signs are interpretable, magnitudes
are not. Needs a caption note on `param_plot_cover.pdf`.

**Open gaps:**

- Manuscript: the Abstract is a heading only; Figure S1 has a caption and no
  figure; `REF`/`REFS` placeholders remain in Introduction paragraphs 1-3;
  `MS2/NEON1_MS2_0_4.docx` has not been moved to `previous_versions/`.
- Cross-validation (`do_model_xv.R`) has never been run. Paths are fixed, but it
  refits each model 17 times (2+5+10 folds) - budget ~11 h for both.
- Model 1 is converged and reported in Table 2, but is still unfigured and has
  no variance-partitioning or parameter plots of its own. `evaluateModelFit`
  returns `NA` for AUC and Tjur R2 for 1 of its 94 species, so those means need
  `na.rm = TRUE` (0.9622 and 0.5629).
- Trait data (DP1.10026.001) is imported but never connected - "unsure how to
  connect the samples to each other".
- Whether the three dropped species should still be dropped under the refreshed
  data.
- Soil chemistry is imported and still absent from `XFormula`; Table 1 now says
  so explicitly rather than implying it is a fixed effect.

## 8. Gotchas - read before running anything

### Working directory - two conventions

Every script except `import_NEON.R` is run from **`NEON1/`**, the package root.
`import_NEON.R` is run from **`NEON1/exec/`** (it reads `'NEON_TOKEN.R'` and
writes `'../../data_NEON'`). Both paths resolve to the same `data_NEON`
directory; `import_data.R` reads it as `'../data_NEON'` from the root. The
deviation is deliberate - it followed BEM's own path specification - and is
flagged in a `NOTE:` at the top of `import_NEON.R`. Normalising it to the root
(`'exec/NEON_TOKEN.R'` and `'../data_NEON'`) is a standing offer.

### Fixed on 2026-10-19 - do not reintroduce

The three stale absolute roots (`dir0 <- '/media/bem/data/NEON'` in
`make_models.R` and `do_model_xv.R`, `e1_dir` in `import_topo.R`) are now
`dir0 <- '..'` and `e1_dir <- '../spatial/AOP'`. The Windows OneDrive paths in
`import_fence.R` are now `'../results/reports'` and `'../data/...'`.
`make_map_figures.R`'s `'../spatial/NEON_shapefiles'` is now `NEON_Shapefiles`
(capital S, case-sensitive on Linux). On 2026-10-19 `make_map_figure.R` was
renamed `make_map_figures.R` (`git mv`, rename staged) and extended to build
Figure 2; the Figure 1 ggplot is now assigned to `map_flow` and passed to
`ggsave` explicitly instead of relying on `last_plot()`. `analyze_climate.R` and
`analyze_structure.R` still carry `C:/Users/BrandonMcNellis/...` paths - they
are off the manuscript path and were not touched.

### Package drift that bites silently

- **terra >= 1.8 returns shapefile date fields already parsed as `POSIXct`.**
  `import_fence.R`'s original `strptime(DateComp, '%Y/%m/%d')` therefore
  returned `NA` for every enclosure, with no error, collapsing
  `time_since_fence` to 0 for every plot. Now
  `as.integer(format(as.POSIXct(DateComp, tz = 'UTC'), '%Y'))`. With the fix,
  `time_since_fence` takes its documented values {0, 10, 22, 27, 29}.
- **MCMCvis `MCMCtrace()` calls `utils::browseURL()`** and dies in batch R
  (`'browser' must be a non-empty character string`). All six calls in
  `make_models.R` now pass `open_pdf = F`. Without it the script dies *after*
  both models are saved but *before* the Rhat tables are written.
- **ggplot2 4.0 errors on `geom_sf(show.legend = 'polygon')`** with ``` `x` must
  be a vector, not `NULL` ```. `TRUE` gives the same filled-polygon key.
  `NA` and `'line'` also build; only `'polygon'` fails.
- **terra/sf need `PROJ_DATA`/`PROJ_LIB` set** to the conda env's `share/proj`
  or every CRS operation fails with "empty srs" even though `proj.db` is on
  disk. See the run recipe in §5. In the `r` kernel, derive it with
  `normalizePath(file.path(R.home(), '..', '..'))`.
- **The fence unit shapefile cannot be transformed out of its own CRS.**
  `data/PuuMakaala_Fence_Units_20240517.zip` ships a `.prj` declaring NAD83 /
  UTM 5N (EPSG:26905), and `st_transform` out of it returns **all-NA
  coordinates** here because the NAD83 datum grid is not installed. The failure
  is silent: an empty map panel, no error. Declaring `EPSG:32605` (WGS 84 / UTM
  5N) before transforming is the fix, and is what `import_fence.R` and
  `make_map_figures.R` both do. **Do not "correct" that override.**
- **The NEON API returns rows in a different order now.** An unsorted
  `all.equal` on `met` reports elevation differing by 21 %; sorted by `plotID`
  it is `TRUE`. Same for `fence`. Sort before comparing.
- **NEON dropped `waterRetention15Bar`** from the soil biogeochem table, so
  `soil_init` has 15 columns, not 16.

### The R kernel

- **Never call `sink()` in the `neon1-r` kernel.** It kills the kernel with
  `Error in stdout() : invalid connection`. Use
  `out <- capture.output({...}); writeLines(out, 'file.txt')` instead.
- **Never call `rm(list = ls())`.** It deletes the harness helper
  `truncate_output()` from the global environment and kills the kernel too.
  Remove named objects only.
- The kernel does not see `NEON1` until
  `.libPaths(c('/media/bem/data/NEON_claude/rlib', .libPaths()))` is called.
  Run the pipeline scripts with `source('exec/<script>.R')` from the package
  root rather than shelling out to `Rscript`, which cannot find the conda R.

### Still blocked

- `import_topo.R`: `Rsagacmd::saga_gis()` fails with "SAGA-GIS installation not
  found". SAGA is not installable from the configured conda channels under any
  name tried (`saga`, `saga-gis`). The TWI it computes is not in either model
  formula. It also ends with `usethis::use_data(topo)` with **no**
  `overwrite = T`, so it would error on a rerun even with SAGA present.
- `import_NEON.R` cannot be run by a sandboxed agent - see §4.

### Superseded and broken code

- **`analyze_Hmsc.R` is superseded.** An earlier, larger draft (TWI, `pig`,
  lat/lon/date random level, logit-cover matrix `mat_f`, `m_p_mod.rda`) with a
  broken `stopifnot` on an undefined `re1`. `make_models.R` is the script that
  produced the models on disk.
- `process_UNK_species.R` errors: `'PUUM_div' is not an exported object from
  'namespace:NEON1'`. `extra.R` errors on an undefined `plant_df`.
- **Known-broken package functions.** `standardize_df()` writes to `df1` before
  it exists; `key_plots_by_type()` is an unfinished `data.frame` with a trailing
  comma. Neither is used by the current pipeline.

### Data caveats baked into the models

- **Imputed PAI.** `import_data.R` copies DHP values into plots `PUUM_032`,
  `PUUM_034`, `PUUM_042` from `PUUM_039`, `PUUM_031`, `PUUM_041` respectively.
  Three of the PAI values in the model are borrowed, not measured.
- **Fence-date fixes.** Two enclosures (Kata-Stein, Unencumbered Kulani Pasture)
  had no completion date and were set to 2014 (youngest neighbour); Army Road
  was listed as not ungulate-free so `time_since_fence = 0`. `time_since_fence`
  is `2024 - enclosure_completed` and is the **only** predictor not centred or
  scaled (`pai`, `elevation`, `age_median` are). That is why it dominates
  `param_plot_disturb.pdf` visually at a different magnitude from the others.
- **Flow age** is the median of each published USGS age range, hardcoded in
  `make_map_figure.R`; `key_age_range()` defaults to `current_year = 2023`
  while the fence variable uses 2024.
- **Taxonomy filters.** Records not identified to species are dropped
  (`_sp`/`_spp` patterns). `nativeStatusCode` of `NI?` and `UNK` are treated as
  native `N`; invasive is recoded `z_I` so it sorts last and becomes the Hmsc
  reference contrast. `cover_type` `ohia_tall` is recoded `AAA_ohia_tall` for
  the same reason.

### Figures and files

- **Figures on disk are 9 files.** `map_figure` and `disturbance_figure` as
  `.pdf` and `.png`, plus `vp_plot`, `param_plot_env`, `param_plot_disturb`,
  `param_plot_cover` and `ve_invasive` as `.pdf`. The manuscript numbers them
  1-7 in that order; `make_figures.R`'s own comments still say 4-8.
- **Figure labels are clipped at 35 species.** The first species name runs off
  the left edge of `vp_plot.pdf`, `param_plot_env.pdf` and
  `param_plot_cover.pdf`, and in the last the y-axis title collides with a
  species label. The figure code is unchanged - it just has more categories to
  fit. Unfixed; it is a styling decision on BEM's figures.
- **`make_figures.R` attaches no packages** yet calls `ggplot()`, `aes()`,
  `geom_bar()`, `theme_bw()` unqualified. It works only because
  `NEON1::param_plot()` runs `require(ggplot2)` earlier in the script. Fragile.
- **Trace PDFs are enormous** - `trace_1_omega.pdf` is 1.26 GB,
  `trace_0_omega.pdf` 174 MB. Don't open them casually; the Omega trace calls
  are deliberately at the end of `make_models.R`.
- **`delete_host_files` cannot remove directories under `/media/bem/data`** -
  the Trash is on the home disk and there is no per-filesystem Trash on the
  data disk, so it fails with "cannot trash a directory across filesystems".
  Individual files there trash fine. Ask BEM to remove directory trees.
- **Git state.** `NEON1` is on `main`. `git status` shows every tracked file
  modified, which is line-ending churn (the working tree is CRLF, the index
  LF) rather than real edits - run `git diff --numstat` and check whether
  insertions equal deletions before assuming uncommitted work. Preserve CRLF
  when editing. `exec/import_NEON.R` is untracked; `exec/NEON_TOKEN.R` is
  ignored.

## 9. Conventions specific to this project

General R house style — project shape, the function-file split, inline `pkg::`
qualification, `T`/`F`, `set.seed(1)`, `usethis::use_data()`, relative paths —
lives in the **`r-coding-conventions`** skill. Load that alongside this one; only
the NEON1-specific conventions are below.

- Predictors are centred/scaled in the data-prep step, and Hmsc is called with
  `XScale = F` — don't let Hmsc rescale.
- Posterior summaries are **90% HPD** via `posterior_from_coda(..., 0.9)`;
  `make_figures.R` has a commented-out 95% alternative.
- Species labels in figures use `binomialName` (`Genus_species`) with
  underscores replaced by spaces at plot time.
- Variance partitioning is reported multiplied by each species' explanatory R²
  ("amount of total variation explained by each covariate") — the note appears at
  the top of both model scripts.
- Scripts expect the package attached (`library(NEON1)`) and call its own
  functions fully qualified as `NEON1::`.

## 10. Helper functions

`kernel.R` ships with this skill and defines, in the R kernel:

- `neon_path(..., root = NULL)` — build a path under the data root
- `neon_load_model(which = 0, diag_only = F, root = NULL)` — load one of the
  four `.rda` model/diagnostic files into an isolated environment and return its
  contents as a named list (avoids clobbering the global namespace)

## 11. Reference files

Each is loaded only when its trigger fires; none of them repeats a rule that
must be followed unconditionally.

- `reference/manuscript_state.md` - **before editing `MS2/NEON1_MS2_0_5.docx`,
  renumbering or regenerating a figure, or adding a citation.** Section-by-section
  state, word counts, the figure-to-script table, equations, tables, and the
  citations not yet in `MS2/references.docx`.
- `reference/site_facts.md` - **before writing or checking a study-system
  sentence, citing DLNR (2013), or touching the ungulate exclusion dates.**
  Verified site facts with page numbers, the 13-unit fence table, and the
  reconciliation of the 1991-2003 window with the layer's completion years.
- `reference/model_history.md` - **before quoting any fit,
  variance-partitioning or trait-R2 number that did not come from the current
  `m_p_diag*.rda`.** The 2024-vs-2026 comparison tables.
