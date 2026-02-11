<!-- Parent: ../AGENTS.md -->
<!-- Generated: 2026-02-10 | Updated: 2026-02-10 -->

# testthat/

## Purpose

Individual test files for each module in the fmrilatent package. Each file corresponds to one or more R source files and tests their public and key internal behavior.

## Key Files

### Core Tests

| File | Tests |
|------|-------|
| `test-latent_neurovector.R` | `LatentNeuroVec` constructor, dim, `[`, `[[`, `series()`, `show()` |
| `test-all_generic.R` | Generic function dispatch: `basis()`, `loadings()`, `offset()`, `map()`, `mask()` |
| `test-encode.R` | `encode()` pipeline with various spec combinations |
| `test-encode-factory.R` | `latent_factory()` convenience wrapper |
| `test-reduction.R` | `GraphReduction`, `ClusterReduction`, `lift()` methods |
| `test-handles.R` | `BasisHandle`/`LoadingsHandle` materialization and registry |
| `test-latent_handles.R` | Handle creation, caching, registry enable/disable |
| `test-registry.R` | Registry stats, clear, enable/disable operations |

### Wavelet Tests

| File | Tests |
|------|-------|
| `test-haar_wavelet.R` | Haar forward/inverse, Morton ordering, round-trip fidelity |
| `test-zzz_haar_aliases.R` | Backward-compatible Haar function aliases |
| `test-heat_wavelet.R` | Heat-kernel wavelet basis construction |
| `test-heat_wavelet_handle.R` | Heat wavelet handle materialization |
| `test-diffusion_wavelet.R` | Diffusion wavelet basis on cluster graphs |
| `test-wavelet_active.R` | CDF 5/3 pencil wavelet on active voxels |

### Basis Tests

| File | Tests |
|------|-------|
| `test-dpss.R` | DPSS/Slepian basis generation (dense and tridiagonal) |
| `test-hrbf.R` | HRBF basis generation, projection, reconstruction |
| `test-hrbf-atoms.R` | HRBF Rcpp atom assembly (Gaussian/Wendland kernels) |
| `test-bspline_hrbf_st.R` | B-spline temporal + HRBF spatial spatiotemporal encoding |
| `test-latent_dct_heatwavelet.R` | Composite DCT + heat-wavelet pipeline |
| `test-implicit_latent.R` | Implicit latent closure-based reconstruction |

### Slepian Tests

| File | Tests |
|------|-------|
| `test-slepian_recon.R` | Slepian reconstruction accuracy |
| `test-slepian_spatial.R` | Graph Slepian spatial basis |
| `test-slepian_spatiotemporal.R` | Joint spatiotemporal Slepian decomposition |
| `test-slepian_handles.R` | Slepian handle constructors |
| `test-slepian_plot.R` | Slepian taper plotting |

### Other Tests

| File | Tests |
|------|-------|
| `test-hierarchical_helpers.R` | Spectral Ward clustering, nested parcellation validation |
| `test-hierarchical_template.R` | `HierarchicalBasisTemplate` build/save/load/project |
| `test-graph_bridge.R` | Graph representation conversion utilities |
| `test-searchlight_utils.R` | `latent_searchlight()` local ROI analysis |
| `test-spatial_plot.R` | Spatial atom plotting |
| `test-benchmark_roundtrip.R` | Encode-decode benchmarking |
| `test-sim-fmri.R` | Simulated fMRI data for integration tests |

## For AI Agents

### Working In This Directory

- Name new test files `test-<module_name>.R` matching the R source file being tested
- Each test file should `library(testthat)` at the top (or rely on the testthat runner)
- Use `test_that("description", { ... })` blocks with descriptive names
- Keep tests independent — no reliance on execution order between files

### Common Patterns

- **Synthetic mask**: `mask_arr <- array(FALSE, dim = c(5, 5, 5)); mask_arr[2:4, 2:4, 2:4] <- TRUE`
- **NeuroSpace**: `spc <- neuroim2::NeuroSpace(c(5, 5, 5, n_time))`
- **LogicalNeuroVol**: `mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(5, 5, 5)))`
- **Round-trip check**: `expect_lt(max(abs(reconstructed - original)), tolerance)`
- **Rcpp toggle**: `withr::with_options(list(fmrilatent.haar.use_rcpp = FALSE), { ... })` to test R fallback

<!-- MANUAL: Any manually added notes below this line are preserved on regeneration -->
