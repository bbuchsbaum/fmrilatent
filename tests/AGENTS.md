<!-- Parent: ../AGENTS.md -->
<!-- Generated: 2026-02-10 | Updated: 2026-02-10 -->

# tests/

## Purpose

Test suite for fmrilatent using testthat edition 3. Entry point is `testthat.R` which bootstraps the test runner.

## Key Files

| File | Description |
|------|-------------|
| `testthat.R` | Test runner entry point: loads `testthat` and `fmrilatent`, calls `test_check()` |

## Subdirectories

| Directory | Purpose |
|-----------|---------|
| `testthat/` | All test files (see `testthat/AGENTS.md`) |

## For AI Agents

### Working In This Directory

- Run tests: `Rscript -e "devtools::test()"` (from package root)
- Run a single test file: `Rscript -e "devtools::test(filter = 'haar_wavelet')"`
- Do not modify `testthat.R` unless changing the test framework version
- New test files go in `testthat/test-*.R`

### Testing Requirements

- Use `set.seed()` for reproducible random data
- Create small synthetic masks (e.g., 5x5x5) to keep tests fast
- Test round-trip fidelity: encode -> LatentNeuroVec -> reconstruct -> compare
- Use `testthat::expect_*` assertions (edition 3 style)

<!-- MANUAL: Any manually added notes below this line are preserved on regeneration -->
