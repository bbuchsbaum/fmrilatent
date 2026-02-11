<!-- Parent: ../AGENTS.md -->
<!-- Generated: 2026-02-10 | Updated: 2026-02-10 -->

# inst/

## Purpose

Files installed with the package. Contents are copied verbatim to the installed package directory and accessible at runtime via `system.file()`.

## Subdirectories

| Directory | Purpose |
|-----------|---------|
| `scripts/` | Standalone R scripts for building precomputed data (e.g., hierarchical atlas templates) |

## Key Files

| File | Description |
|------|-------------|
| `scripts/build_hier_templates.R` | Script to build hierarchical Schaefer atlas templates at multiple resolutions |

## For AI Agents

### Working In This Directory

- Files here are user-accessible after installation via `system.file("scripts", "build_hier_templates.R", package = "fmrilatent")`
- Do not place internal package code here — use `R/` for that
- Large data files should go in `inst/extdata/` (currently unused)

<!-- MANUAL: Any manually added notes below this line are preserved on regeneration -->
