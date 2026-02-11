<!-- Parent: ../AGENTS.md -->
<!-- Generated: 2026-02-10 | Updated: 2026-02-10 -->

# vignettes/

## Purpose

R Markdown vignettes providing narrative documentation and tutorials for the fmrilatent package. Built during `R CMD check` and rendered to the pkgdown site.

## Key Files

| File | Description |
|------|-------------|
| `intro.Rmd` | Introduction to fmrilatent: core concepts, LatentNeuroVec construction, basic usage |
| `encode-factory.Rmd` | Guide to the `encode()` pipeline and `latent_factory()` convenience API |
| `albers.css` | Custom CSS for Albers color theme (used by pkgdown) |
| `albers.js` | Custom JS for Albers theme interactivity |

## For AI Agents

### Working In This Directory

- Vignettes use `knitr` and `rmarkdown` for rendering
- Test locally: `Rscript -e "rmarkdown::render('vignettes/intro.Rmd')"`
- Keep code chunks small and fast (vignettes run during `R CMD check`)
- Use `eval = FALSE` for expensive or data-dependent examples
- CSS/JS files are shared with the pkgdown site theme

### Testing Requirements

- `Rscript -e "devtools::build_vignettes()"` builds all vignettes
- Vignettes must complete without error during `R CMD check`

<!-- MANUAL: Any manually added notes below this line are preserved on regeneration -->
