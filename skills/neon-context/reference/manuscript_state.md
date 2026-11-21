# Manuscript state: v0_6 (2026-11-08)

Read this before editing `ms/NEON1_0_6.docx`, before renumbering or
regenerating a figure, or before adding a citation. Paths are under
`/media/bem/data/NEON1`.

## Files and versioning

- `ms/NEON1_0_6.docx` - live draft, written 2026-11-08 from BEM's review copy.
- `ms/NEON1_0_5_BEM.docx` - BEM's tracked edits + 33 comments on 0_5 (the
  review copy; never edit). `ms/NEON1_0_5.docx` - the 0_5 as generated. BEM
  files older versions into `ms/previous_versions/` himself; leave them.
- Naming: `NEON1_<major>_<minor>.docx` (the `MS2_` infix was dropped with the
  directory change); BEM's review copy adds `_BEM`.
- `ms/references.docx` is retired: it was moved to
  `ms/previous_versions/references_08112026_0032.docx` (BEM's instruction). The
  Literature Cited now lives inside the manuscript (58 entries, Springer
  author-year style, DOIs from Crossref).

## How 0_6 was built (repeat this for 0_7)

Direct OOXML with lxml from the `_BEM` copy (not python-docx regeneration):
accept all tracked changes (drop `w:del`, unwrap `w:ins`, remove `w:ins`/`w:del`
markers in `w:pPr/w:rPr`; a `w:del` in a paragraph-mark rPr means BEM merged
that paragraph with the next one - honour it), strip addressed comments, then
replace paragraph content from a mini-markup (`*italic*`, `**bold**`,
`^{sup}`, `_{sub}`, leading `\t`). Body paragraphs begin with a tab and are
separated by blank paragraphs (BEM's 0_5 edits). BEM's settings keep
`trackRevisions` on. LibreOffice crashes in the sandbox, so verify by
re-parsing the XML (python-docx loads it) and dumping text with run markup.
Figures are swapped by overwriting `word/media/image1-7.png` (1 map, 2
disturbance, 3 vp_plot, 4 env, 5 disturb, 6 cover, 7 ve_invasive; PDFs
rendered with `pdftoppm -r 200`).

## Structure of 0_6

- **Notes for co-authors** sit at the top, above the Abstract, as a numbered
  list (BEM's instruction): changes in 0_6, answers to his comments, and
  open issues. Rewrite them each version; do not carry stale notes.
- Introduction paragraphs 1-2 are BEM's 0_4 **outline text, verbatim**, each
  followed by a bold bracket of five suggested references (all in Literature
  Cited). Do not re-draft them into prose unless asked.
- The four numbered questions carry no citations, and each ends with a bold
  sentence naming the model output and figures/tables that answer it.
- The rooting-and-grazing mechanism lives in the pig paragraph (moved out of
  Question 1).
- 'Species traits' is no longer a subsection; the trait text is appended to
  the disturbance paragraph of Explanatory variables (BEM: "combine with the
  above paragraph").
- Results headings: 'Model performance and evaluation', 'Factors affecting
  plant occurrence and abundance across PUUM', 'Scale effects in model
  results' (BEM's).
- Every DLNR (2013) citation carries a page number (BEM's rule).
- Data-handling/plotting packages are cited in `NEON1/DESCRIPTION`, not the
  text (BEM's rule); only Hmsc, coda, MCMCvis, bayesplot and neonUtilities are
  cited in Software.
- Table 1 has no 'Used' column; unused data rows say 'Obtained, not modelled'.

## Figures

| Fig | File | Script |
|---|---|---|
| 1 | `results/figures/map_figure.png` | `make_map_figures.R` |
| 2 | `results/figures/disturbance_figure.png` | `make_map_figures.R` (2026 reference year) |
| 3 | `vp_plot.pdf` | `make_figures.R` |
| 4-6 | `param_plot_env/disturb/cover.pdf` | `make_figures.R`; whiskers `err_width = 0.63` |
| 7 | `ve_invasive.pdf` | `make_figures.R`; 90% HPD from `NEON1::r2t_draws()` |
| S1 | not produced | `analyze_climate.R` - cannot run (see open issues) |

Figure 6 drops cover-type terms for species never recorded in that cover type
(39 of 65 credible terms, all negative) and uses `facet_grid(space =
'free_x')` with short labels (Koa, Ohia-ash, Ohia-koa, Ohia woodland). Figure
comments in `make_figures.R` now match the manuscript numbering.

## Open issues carried in the 0_6 notes

[PAI guy] author (HemiPy, Brown et al. 2023: Luke A. Brown or Courtney Meier?);
Meier et al. (2023) unidentified; Schimel et al. (2015) matched to Schimel and
Keller (2015) - confirm; a citable source for post-2013 fence dates (DLNR 2024
is the unpublished layer); Figure S1; no plot-level random effect; *Dodonaea
viscosa* coded introduced; grazing end dates per plot; acknowledgement for the
NEON land-use table; Abstract; 'WREF, JORN, REFS' placeholder; cross-validation;
10-m2 model figures; Introduction over 800 words. The convergence standard
(PSRF < 1.1, counts > 1.01 reported) was stated by Claude and awaits BEM's
confirmation.

## Reports written for BEM's comments

- `results/reports/foliar_traits/` - `foliar_traits.md`, `.R`, CSVs, PNG.
  DP1.10026.001 at PUUM: 49 samples (2019, 2025), 23 *M. polymorpha*; 10 of
  35 1-m2 species sampled, none introduced. Hmsc traits are per species, so
  plot-date linkage was never the obstacle.
- `results/reports/konig_2021/konig_2021.md` - what König et al. did and how
  NEON1 departs (no plot/year random effect; grains not nested; no
  out-of-sample test; 90% vs 95%).
