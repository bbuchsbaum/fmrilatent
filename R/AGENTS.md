<!-- Parent: ../AGENTS.md -->
<!-- Generated: 2026-02-10 | Updated: 2026-02-10 -->

# R/

## Purpose

All R source code for the fmrilatent package. Contains S4 class definitions, generic functions, method implementations, basis constructors, encoding pipeline, and utility functions. Files are loaded in the collation order specified in `DESCRIPTION`.

## Key Files

### Core Class & Infrastructure

| File | Description |
|------|-------------|
| `latent_handles.R` | `BasisHandle` / `LoadingsHandle` S4 classes and registry cache (loaded first) |
| `all_class.R` | `LatentNeuroVec` S4 class definition with slots: basis, loadings, offset, map, meta |
| `all_generic.R` | S4 generics: `basis`, `loadings`, `offset`, `map`, `mask`, `lift`, `encode` |
| `latent_neurovector.R` | `LatentNeuroVec` constructor, `show()`, `series()`, `[`, `[[` methods |
| `latent_methods.R` | Additional methods: `concat`, `linear_access` |
| `latent_indexing.R` | Voxel and time-point indexing logic |
| `latent_neurovec_materialize.R` | Materialization of `BasisHandle`/`LoadingsHandle` to dense matrices |
| `fmrilatent-package.R` | Package-level roxygen2 documentation |
| `RcppExports.R` | Auto-generated Rcpp bindings. **Do not hand-edit.** |

### Encoding Pipeline

| File | Description |
|------|-------------|
| `encode.R` | `encode()` S3 generic, `spec_time_*`/`spec_space_*`/`spec_st` constructors, `encode_spec()` dispatch |
| `reduction.R` | `GraphReduction`, `ClusterReduction`, `CoarsenedReduction` classes + `lift()` methods |

### Temporal Bases

| File | Description |
|------|-------------|
| `dct_basis.R` | Discrete cosine transform basis builder + `dct_basis_handle()` |
| `bspline_basis.R` | B-spline basis builder + `bspline_basis_handle()` |
| `slepian_temporal.R` | DPSS/Slepian temporal basis via `slepian_temporal_latent()` |

### Spatial Bases

| File | Description |
|------|-------------|
| `hrbf.R` | Hierarchical radial basis functions: `hrbf_generate_basis()`, project/reconstruct, `hrbf_latent()` |
| `slepian_spatial.R` | Graph Slepian spatial basis: `slepian_spatial_latent()` |
| `slepian_handles.R` | Handle constructors for Slepian loadings |

### Wavelet Systems

| File | Description |
|------|-------------|
| `haar_wavelet.R` | Morton-ordered Haar lifting transform (forward/inverse), `haar_latent()` |
| `zzz_haar_aliases.R` | Backward-compatible aliases for Haar functions |
| `heat_wavelet.R` | Heat-kernel graph wavelets on cluster reductions |
| `heat_wavelet_handle.R` | `heat_wavelet_loadings_handle()` for lazy materialization |
| `diffusion_wavelet.R` | Diffusion wavelet bases on cluster graphs |
| `diffusion_wavelet_handle.R` | `diffusion_wavelet_loadings_handle()` for lazy materialization |
| `wavelet_active.R` | CDF 5/3 pencil wavelet on mask-active voxels |

### Spatiotemporal & Composite

| File | Description |
|------|-------------|
| `slepian_spatiotemporal.R` | Joint spatiotemporal Slepian decomposition |
| `latent_dct_heatwavelet.R` | Composite DCT temporal + heat-wavelet spatial pipeline |
| `implicit_latent.R` | `ImplicitLatent` — closure-based reconstruction without stored basis |

### Hierarchical Templates

| File | Description |
|------|-------------|
| `hierarchical_template.R` | `HierarchicalBasisTemplate` class, build/save/load/project functions |
| `hierarchical_helpers.R` | `spectral_ward_hclust`, `cut_hclust_nested`, `validate_nested_parcellations` |

### Utilities

| File | Description |
|------|-------------|
| `searchlight_utils.R` | `latent_searchlight()` for local ROI analysis |
| `spatial_plot.R` | `plot_spatial_atom()` for visualizing spatial basis components |
| `slepian_plot.R` | `plot_slepian_temporal()` for Slepian taper visualization |
| `graph_bridge.R` | Bridge utilities for converting between graph representations |
| `benchmark_roundtrip.R` | `benchmark_roundtrip()` encode-decode fidelity testing |

## For AI Agents

### Working In This Directory

- **Collation order matters**: `DESCRIPTION` specifies `Collate:` field. `latent_handles.R` must load before `all_class.R` (handles define class unions used in LatentNeuroVec slots). If adding a new file, add it to `Collate:` in `DESCRIPTION`.
- **roxygen2 `@include` tags**: Files use `@include` to declare source-level dependencies (e.g., `#' @include all_generic.R`). These must be consistent with the `Collate:` field.
- Run `Rscript -e "devtools::document()"` after modifying roxygen2 comments.
- Naming: snake_case for functions, `UpperCamelCase` for S4 classes.
- S4 methods: `setMethod("generic", signature("ClassName"), function(...) { ... })`.
- Keep `RcppExports.R` auto-generated — modify C++ source in `src/` and run `Rcpp::compileAttributes()`.

### Testing Requirements

- Each file should have a corresponding `tests/testthat/test-<name>.R`.
- Test round-trip fidelity: encode -> LatentNeuroVec -> reconstruct -> compare to original.
- Use small deterministic masks with `set.seed()` for reproducibility.

### Common Patterns

- **Null coalescing**: `%||%` defined locally in several files (not exported).
- **Handle pattern**: Create `BasisHandle`/`LoadingsHandle` with unique `id` (digest hash); store materialized matrix in registry on first access.
- **Spec + encode_spec**: Lightweight spec object -> `encode_spec.spec_*()` S3 method builds the actual basis/loadings.
- **Rcpp toggle**: Functions like `use_haar_rcpp()` check `getOption("fmrilatent.haar.use_rcpp")` and fall back to pure R.

## Dependencies

### Internal
- `neuroim2` — `NeuroVec`, `NeuroSpace`, `LogicalNeuroVol`, `IndexLookupVol`, `series`, `linear_access`, etc.
- `src/` — Rcpp functions via `RcppExports.R`

### External
- `Matrix` (sparse/dense matrices), `methods` (S4), `Rcpp`, `crayon`, `digest`, `utils`
- Suggested: `RSpectra` (fast SVD), `neuroatlas`, `neurosurf`, `rgsp`, `bench`, `ggplot2`

<!-- MANUAL: Any manually added notes below this line are preserved on regeneration -->
