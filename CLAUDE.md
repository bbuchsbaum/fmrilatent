# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**fmrilatent** is an R package for latent space representations of fMRI neuroimaging data. It stores 4D brain data as `basis × loadings + offset` factorizations instead of full arrays, enabling compact storage for PCA, ICA, wavelet, and other decompositions. Built on top of the `neuroim2` package.

## Quick Commands

```bash
Rscript -e "devtools::load_all()"           # Load package
Rscript -e "devtools::test()"               # Run all tests
Rscript -e "testthat::test_file('tests/testthat/test-encode.R')"  # Single test file
Rscript -e "devtools::document()"           # Regenerate roxygen2 docs + NAMESPACE
Rscript -e "devtools::check()"              # Full R CMD check
```

After modifying C++ files in `src/`, run `devtools::load_all()` to recompile. The package uses Rcpp/RcppEigen with optional OpenMP (see `src/Makevars`).

## Architecture

### Core Data Model

`LatentNeuroVec` (S4 class in `R/all_class.R`) is the central type. It inherits from `neuroim2::NeuroVec` and represents 4D neuroimaging data as:

```
data[voxel, time] = basis[time, k] × loadings[voxel, k]^T + offset[voxel]
```

Slots accept either dense matrices or **handles** (`BasisHandle`/`LoadingsHandle` in `R/latent_handles.R`) — lightweight references to dictionaries that materialize on demand. This enables implicit/lazy representations where the basis is never stored as a full matrix (e.g., DCT, Slepian, lifted wavelets).

### Handle System & Registry

- `R/latent_handles.R` — `BasisHandle`/`LoadingsHandle` S4 classes and the internal materialization cache (`.fmrilatent_cache_env`)
- `R/latent_neurovec_materialize.R` — `basis_mat()`/`loadings_mat()` internal generics that dispatch on Matrix vs Handle types for operator-aware access
- Handles use `kind` (e.g., "dct", "lifted", "explicit") and `spec` (parameter list) to reconstruct their matrix when needed

### Encoding Pipeline

The `encode()` function (`R/encode.R`) is the main entry point for creating latent representations from raw `NeuroVec` data. It works via a spec-dispatch pattern:

1. **Spec constructors** create typed config objects: `spec_time_dct()`, `spec_time_slepian()`, `spec_space_heat()`, `spec_space_pca()`, `spec_st()` (spatiotemporal), etc.
2. **`encode_spec()`** S3 generic dispatches on spec class to build the actual `LatentNeuroVec`
3. **Encoder registry** (`R/encoder_registry.R`) — `register_encoder()`/`list_encoders()` allow external packages to register new encoder families

### Spatial Encoder Families

Each spatial encoder produces a loadings dictionary for a different basis type:

| File | Basis Type |
|------|-----------|
| `R/slepian_spatial.R` / `R/slepian_temporal.R` | Slepian/DPSS concentration |
| `R/heat_wavelet.R` / `R/heat_wavelet_handle.R` | Heat-kernel wavelets on graphs |
| `R/diffusion_wavelet.R` / `R/diffusion_wavelet_handle.R` | Diffusion wavelets |
| `R/haar_wavelet.R` | Haar wavelets (+ aliases in `R/zzz_haar_aliases.R`) |
| `R/hrbf.R` | Hemodynamic response basis functions |
| `R/pca_spatial.R` | PCA spatial encoder |
| `R/wavelet_active.R` | Active-set wavelet selection |
| `R/bspline_basis.R` / `R/dct_basis.R` | B-spline and DCT temporal bases |

### Other Key Modules

- `R/reduction.R` — Abstract `GraphReduction` S4 classes (`ClusterReduction`, `CoarsenedReduction`) for voxel grouping/coarsening
- `R/hierarchical_template.R` + `R/hierarchical_helpers.R` — `HierarchicalBasisTemplate` for multi-level parcellation-based encoding
- `R/implicit_latent.R` — `ImplicitLatent` S3 class for decoder-only (basis-free) representations
- `R/compat_profile.R` — Compatibility profiles for external package integration (e.g., `neuroarchive`)
- `R/latent_indexing.R` + `R/latent_methods.R` — `[`, `[[`, `series()`, `linear_access()`, `concat()` methods on `LatentNeuroVec`
- `R/graph_bridge.R` — Bridges between graph structures and spatial encoders

### C++ (src/)

- `slepian_dpss_rcpp.cpp` — DPSS eigenvalue computation via tridiagonal solver (RcppEigen)
- `haar_wavelet_rcpp.cpp` — Fast Haar wavelet forward/inverse transforms
- `hrbf_atoms_rcpp.cpp` — HRBF atom generation
- `active_pencil_wavelet.cpp` — Active-set wavelet pencil operations

### File Load Order

The `Collate` field in DESCRIPTION controls load order. Key dependencies: `latent_handles.R` → `all_class.R` → `all_generic.R` → other files. When adding new files, update `Collate` if they define classes/generics needed by other files.

## Code Style

- 2-space indent, snake_case for functions/variables
- S4 classes and methods follow `setClass`/`setMethod` with `generic,Class` dispatch pattern
- S3 classes used for lightweight specs and implicit latent objects
- roxygen2 `@include` tags manage source ordering within the S4 system
- `Matrix` package types preferred over base `matrix` for sparse data

## Issue Tracking

This project uses **beads** (`bd`) for git-backed issue tracking. See AGENTS.md for workflow.

```bash
bd ready                               # Find next task
bd create "title" -p 1                 # Create issue (P0-P4)
bd close <id> --reason "text"          # Close completed task
bd sync                                # Sync to git
```
