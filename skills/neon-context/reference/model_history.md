# Model fit history: 2024 data vs. refreshed data

Read this when a number in an older draft disagrees with the current fit, or
before quoting any variance-partitioning or trait-R2 figure that did not come
from `m_p_diag0.rda` / `m_p_diag1.rda` as they stand on disk. The current
numbers are in SKILL.md section 6; what follows is the record of what changed
when the models were refit on 2026-10-19 against NEON data current to
2025-06-24.

## What the refreshed data changed

NEON published two more field seasons after the 2024 snapshot the original
models were fitted to. `div_one` went 4,346 -> 5,837 rows, `div_ten` 5,748 ->
7,504, survey dates 85 -> 116, last `endDate` 2023-06-29 -> 2025-06-24, taxa in
`div_one` 125 -> 137, `spp` 187 -> 251. `met`, `dhp`, `lus` and `fence` are
unchanged in content (24 plots, same PAI, same land-use history).

Model 0 grew from 142 x 31 to 188 x 35. Four species cleared `oc >= 15` and none
dropped out: *Adenophorus tamariscinus*, *Dichanthelium hillebrandianum*,
*Lellingeria saffordii* and *Rubus ellipticus* - the last an invasive, so the
native/invasive contrast carries one more.

| | 2024 fit, 2024 data | 2026 fit, refreshed data |
|---|---|---|
| RMSE mean (range) | 0.2407 (0.082-0.372) | 0.2304 (0.070-0.370) |
| AUC mean | 0.9565 | 0.9599 |
| Tjur R2 mean | 0.5823 | 0.5846 |
| R2T mean (trait interaction) | 22.02 % | 17.78 % |
| Beta Rhat > 1.01 | 9 / 341 | 6 / 385 |
| Omega Rhat > 1.01 | 66 / 961 | 25 / 1225 |
| significant Beta terms | 129 | 165 |

Model 1's 2024 fit was badly non-converged (362/968 Beta terms over Rhat 1.01,
max 10.36). The current fit is clean: 10/1034 Beta, 0/22 Gamma, 51/8836 Omega,
mean RMSE 0.2092.

## The result that moved

Variance explained by each covariate's interaction with the invasive trait
(`m_vp_0$R2T`), which is what `ve_invasive.pdf` plots:

| covariate | 2024 data | refreshed | rank |
|---|---|---|---|
| historic grazing (`cowTRUE`) | 41.45 % | 17.35 % | 1 -> 1 |
| plant area index | 26.43 % | 13.51 % | 2 -> 3 |
| flow age (`age_median`) | 11.30 % | 6.03 % | 3 -> 6 |
| elevation | 9.37 % | 4.03 % | 4 -> 9 |
| **years ungulate-free** | **7.71 %** | **14.44 %** | **7 -> 2** |
| mature ohia-koa overstory | 7.12 % | 9.22 % | 8 -> 4 |

Grazing stays first but loses more than half its share; years ungulate-free
roughly doubles and moves from seventh to second. **Any manuscript sentence
quoting 41.45 / 26.43 / 11.30 is stale, and the ordering it argues for has
changed.** The ranking is the finding; the magnitudes all shrank because mean
R2T fell.

## Refit of 2026-11-08: reference year 2026

Same data, same species, same specification; the only change is
`time_since_fence = 2026 - enclosure_completed` (was 2024), giving values
{0, 12, 24, 29, 31} instead of {0, 10, 22, 27, 29}. Because unfenced plots stay
at 0, the shift is not a constant offset and the models had to be refit. The
2026-10-19 fit is archived in
`results/model_results/model_results_archive/fit_20261019/`. Current numbers
(model 0 unless stated), read from `m_p_diag0/1.rda` by `extract_06.R`-style
code:

| | 2026-10-19 (ref 2024) | 2026-11-08 (ref 2026) |
|---|---|---|
| VP mean share: elevation / cover / years-free / PAI | 23.1 / 26.8 / 18.6 / 12.9 % | **24.9 / 24.3** / 19.8 / 12.7 % |
| VP x TjurR2: cover / elev / years-free / PAI | 16.5 / 13.4 / 9.5 / 8.3 % | 15.0 / 14.4 / 10.1 / 8.2 % |
| R2T$Y model 0 / model 1 | 17.8 / 13.8 % | 17.2 / 14.5 % |
| R2T cow / PAI / years-free | 17.3 / 13.5 / 14.4 % | 17.5 / 13.6 / 13.3 % |
| credible Beta: elev / PAI / flow / years-free / cow / log | 17 / 23 / 13 / 16 / 8 / 6 | 17 / 25 / 13 / 17 / 8 / 7 |
| credible cover-type terms | 58 | 65 (39 for species never seen in that cover type) |
| Beta Rhat > 1.01 (model 0 / 1) | 6 / 10 | 6 / 8 |
| Omega unique pairs Rhat > 1.01 (model 0 / 1) | 10 / ? | 3 / 8 (one > 1.1 in model 1: 1.12) |
| credible association pairs (model 0 / 1) | 35 / 388 | 34 / 392 |
| credible Gamma interactions, model 1 | ytf x I -0.10, cow x I 0.86 | ytf x I -0.11, cow x I 0.85 |

Elevation now edges cover type for the largest mean VP share, but cover type
is still first once scaled by Tjur R2. Nothing changed direction.

## Counting conventions corrected in v0_6

- `mc_s_omeg_*` covers the full Omega matrix: diagonal (variances) and both
  triangles. v0_5's "25 of 1,225 association parameters (max 1.19)" counted
  both; the max came from a diagonal variance. Count unique off-diagonal pairs
  (ns*(ns-1)/2: 595 and 4,371).
- `post_omeg_*` values are **covariances**. v0_5 quoted them (0.97, -0.88) as
  associations. Residual correlations come from `m_ca_*[[1]]$mean`; for every
  credible pair they are |r| >= 0.94.
- Mean VP shares are `rowMeans(m_vp_0$vals)`; "scaled" shares multiply each
  species column by `mf_p_0$TjurR2` (v0_5 numbers reproduce exactly this way).


