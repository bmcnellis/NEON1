# MS2 v0_5: section-by-section state

Read this before editing `MS2/NEON1_MS2_0_5.docx`, before renumbering or
regenerating a figure, or before adding a citation. Written 2026-10-19 against
the file as saved.

## Build

v0_5 was generated with `python-docx` from `MS2/NEON1_MS2_0_4.docx`, which
supplied the styles: Times New Roman 12 pt, level-1 headings bold, level-2
headings italic, `Table Grid` tables. Every Results number was read from
`m_p_diag0.rda` / `m_p_diag1.rda`; none was carried over from v0_4 text.

## Word counts

| Section | Paragraphs | Words |
|---|---|---|
| Introduction | 7 prose + 4 numbered items | ~1,140 |
| Methods | 7 subsections | ~1,700 |
| Results | 8 | ~915 |
| Discussion | 9 | ~1,630 |

Body ~5,400 of the 8,500-word *Landscape Ecology* allowance. The Introduction
runs over the 800-word target BEM set in v0_4; the question-and-hypothesis list
accounts for ~200 words of the overage and the Hawaiian-forest paragraph for
~120. Cutting is BEM's editorial call.

## Structure as drafted

- **Introduction** ends on four numbered items, each a question then a
  hypothesis then a mechanism, per `manuscript-prose-voice`. Three are
  hypotheses; item 4 is an expectation. **No citations inside the numbered
  list** - this is deliberate, requested on 2026-10-19. Peyton et al. (2023)
  sits in the pigs paragraph and Medeiros et al. (2019) in the
  Hawaiian-forest-properties paragraph instead.
- **Methods** run in order of operations: study system, plant community data,
  explanatory variables, species traits, joint species distribution model,
  model evaluation, software and reproducibility.
- **Results** open with the model-evaluation paragraph, then variance
  partitioning, environment, disturbance, cover type, traits, residual
  associations, and the 10-m2 model.
- **Discussion** has one paragraph answering each Introduction paragraph, with
  the caveats last, pointed back at the paragraph that justified the NEON data.

## Equations

Three numbered equations, (1) to (3), with italic variables and real Word
subscripts, indented with the number on a right tab stop at 6.5 in. Written in
Hmsc notation: a summation over covariates rather than seven spelled-out beta
terms, and the trait model at vector level, `beta_j ~ Normal(Gamma t_j, V)`.
The v0_4 specification (generic probit, `Normal(0,10)` / `HalfCauchy(10)`
priors, `alpha_date` varying intercept) did not match the fitted code and is
gone.

## Figures and tables

Figures are embedded inline above their captions at 6.5 in wide.

| Fig | File | Script |
|---|---|---|
| 1 | `results/figures/map_figure.png` | `make_map_figures.R` |
| 2 | `results/figures/disturbance_figure.png` | `make_map_figures.R` |
| 3 | `results/figures/vp_plot.pdf` | `make_figures.R` |
| 4 | `results/figures/param_plot_env.pdf` | `make_figures.R` |
| 5 | `results/figures/param_plot_disturb.pdf` | `make_figures.R` |
| 6 | `results/figures/param_plot_cover.pdf` | `make_figures.R` |
| 7 | `results/figures/ve_invasive.pdf` | `make_figures.R` |
| S1 | not produced | `analyze_climate.R`, cannot be run |

**The comments in `make_figures.R` still number these 4-8 and are offset by
one.** Fix them when that script is next touched.

- **Table 1** data sources, rewritten with a `Used` column, since v0_4 listed
  soil chemistry as a fixed effect when it is imported but absent from
  `XFormula`.
- **Table 2** fit and convergence of both models side by side.

## Citations added in v0_5 and not yet in `MS2/references.docx`

Ovaskainen et al. (2017), Ovaskainen and Abrego (2020), Tikhonov et al. (2020,
and 2025 for the package), Konig et al. (2021), Gelman and Rubin (1992), Brooks
and Gelman (1998), Tjur (2009), Youngflesh (2018), Chadwick et al. (1999),
Vitousek (2004), Maly and Maly (2004), plus the package citations in the
software section (Plummer et al. 2006, Gabry and Mahr 2025, Gabry et al. 2019,
Lunch et al. 2026, Wickham et al. 2026 and 2025, Wickham 2016, Garnier et al.
2024, Pebesma 2018, Pebesma and Bivand 2023, Hijmans et al. 2026, R Core Team
2026). Sherrod et al. (2021) and Medeiros et al. (2019) were already there.

## Still open in v0_5

- **Abstract** is a heading with a flag. Write it last, against the final
  Results.
- **Figure S1** has a caption and no figure.
- `REF` / `REFS` placeholders remain in Introduction paragraphs 1, 2 and 3.
- The author list still carries `[PAI guy]`.
- `MS2/NEON1_MS2_0_4.docx` is still at the MS2 top level and has not been
  moved to `previous_versions/`; that is BEM's filing convention to apply.
