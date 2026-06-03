<!-- Parent: ../AGENTS.md -->
<!-- Generated: 2026-02-10 | Updated: 2026-06-03 -->

# vignettes/

## Purpose

R Markdown vignettes providing narrative documentation and tutorials for the fmrilatent package. Built during `R CMD check` and rendered to the pkgdown site.

## Key Files

| File | Description |
|------|-------------|
| `fmrilatent.Rmd` | Introduction to fmrilatent: core concepts, core workflows, and article map |
| `encode-factory.Rmd` | Guide to the `encode()` pipeline and `latent_factory()` convenience API |
| `working-with-latentneurovec.Rmd` | Object access, slicing, reconstruction, and matrix contracts |
| `choosing-basis-family.Rmd` | Basis-family trade-offs on a compact toy dataset |
| `shared-spatial-dictionaries.Rmd` | Shared parcel and spatial dictionary workflows |
| `transport-aware-encoding.Rmd` | Transport-aware shared-asset and subject decoder workflows |
| `boldzip.Rmd`, `standalone-codecs.Rmd`, `shared-structure-boldzip.Rmd`, `compression-diagnostics.Rmd` | BOLDZip-SR and codec documentation |
| `explicit-vs-implicit-latents.Rmd` | Explicit versus implicit latent object model |
| `albers-header.html` | Vignette header hook that loads local Albers theme assets |
| `albers.css` | Local CSS for the Albers theme |
| `albers.js` | Local JS for Albers theme interactivity |

## For AI Agents

### Working In This Directory

- Vignettes use `knitr` and `rmarkdown` for rendering
- Test one vignette locally: `Rscript -e "rmarkdown::render('vignettes/fmrilatent.Rmd')"`
- Keep code chunks small and fast (vignettes run during `R CMD check`)
- Use `eval = FALSE` for expensive or data-dependent examples
- Vignette YAML uses `includes: in_header: albers-header.html`; keep local `albers.css` and `albers.js` in sync with README and pkgdown theme notes
- `albersdown` is a vignette dependency and should remain declared in `DESCRIPTION` when vignettes call `albersdown::`

### Testing Requirements

- `Rscript -e "devtools::build_vignettes()"` builds all vignettes
- Vignettes must complete without error during `R CMD check`

<!-- MANUAL: Any manually added notes below this line are preserved on regeneration -->
