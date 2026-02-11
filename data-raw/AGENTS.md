<!-- Parent: ../AGENTS.md -->
<!-- Generated: 2026-02-10 | Updated: 2026-02-10 -->

# data-raw/

## Purpose

Scripts that generate package datasets (saved to `data/` via `usethis::use_data()`). These scripts are not included in the installed package.

## Key Files

| File | Description |
|------|-------------|
| `hierarchical_template.R` | Generates precomputed hierarchical parcellation templates |

## For AI Agents

### Working In This Directory

- Scripts here are run manually during development, not during package installation
- Output should be saved with `usethis::use_data(obj, overwrite = TRUE)`
- Document generated datasets in `R/data.R` (if data files are added to `data/`)

<!-- MANUAL: Any manually added notes below this line are preserved on regeneration -->
