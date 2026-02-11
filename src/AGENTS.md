<!-- Parent: ../AGENTS.md -->
<!-- Generated: 2026-02-10 | Updated: 2026-02-10 -->

# src/

## Purpose

C++ source code for performance-critical operations. Uses Rcpp and RcppEigen for R integration. Compiled to a shared library (`fmrilatent.so`) loaded at package startup. Optional OpenMP support for parallel computation.

## Key Files

| File | Description |
|------|-------------|
| `active_pencil_wavelet.cpp` | CDF 5/3 lifting wavelet: 1D temporal lift and 3D pencil wavelet on mask-active voxels with Morton ordering |
| `haar_wavelet_rcpp.cpp` | Haar lifting transform: Morton-ordered forward/inverse with 3D Z-curve encoding |
| `hrbf_atoms_rcpp.cpp` | HRBF atom assembly: computes sparse basis matrix from atom centers/sigmas using Gaussian or Wendland C6 kernels. OpenMP parallelized. |
| `slepian_dpss_rcpp.cpp` | DPSS (Slepian) basis generation: dense prolate matrix eigenvectors and tridiagonal solver via LAPACK `dstevr` |
| `RcppExports.cpp` | Auto-generated Rcpp glue code. **Do not hand-edit.** |
| `Makevars` | Compilation flags for Unix (OpenMP detection) |
| `Makevars.win` | Compilation flags for Windows |

## Exported C++ Functions

| Function | Source | Description |
|----------|--------|-------------|
| `cdf53_time_lift` | `active_pencil_wavelet.cpp` | 1D CDF 5/3 lifting on a matrix (rows=voxels, cols=time) |
| `active_pencil_wavelet` | `active_pencil_wavelet.cpp` | 3D pencil wavelet transform on mask-active voxels |
| `get_morton_ordered_indices_rcpp` | `haar_wavelet_rcpp.cpp` | Compute Morton Z-curve ordering for 3D voxel mask |
| `forward_lift_rcpp` | `haar_wavelet_rcpp.cpp` | Forward Haar lifting transform |
| `inverse_lift_rcpp` | `haar_wavelet_rcpp.cpp` | Inverse Haar lifting transform |
| `hrbf_atoms_rcpp` | `hrbf_atoms_rcpp.cpp` | Sparse HRBF basis matrix assembly |
| `generate_dpss_basis_rcpp` | `slepian_dpss_rcpp.cpp` | DPSS basis via dense prolate matrix |
| `generate_dpss_tridiag_rcpp` | `slepian_dpss_rcpp.cpp` | DPSS basis via tridiagonal LAPACK solver |

## For AI Agents

### Working In This Directory

- After modifying `.cpp` files, run `Rscript -e "Rcpp::compileAttributes('.')"` from the package root to regenerate `RcppExports.cpp` and `R/RcppExports.R`
- Then `Rscript -e "devtools::load_all()"` to recompile
- Use `// [[Rcpp::export]]` annotation above functions to expose to R
- Follow existing patterns: `#include <RcppEigen.h>` for Eigen types, `#include <Rcpp.h>` for basic Rcpp
- OpenMP guards: wrap parallel regions in `#ifdef _OPENMP` / `#endif`
- `.o`, `.gcda`, and `.so` files are build artifacts — do not commit

### Testing Requirements

- C++ functions are tested through their R wrappers in `tests/testthat/`
- Key tests: `test-haar_wavelet.R`, `test-hrbf-atoms.R`, `test-dpss.R`, `test-wavelet_active.R`
- Test both R fallback and Rcpp paths where applicable

### Common Patterns

- **RcppEigen**: Use `Eigen::Map<>` for zero-copy access to R matrices
- **Sparse output**: Build `std::vector<Eigen::Triplet<double>>` then construct `Eigen::SparseMatrix`
- **LAPACK direct**: `slepian_dpss_rcpp.cpp` calls LAPACK `dstevr` directly via `R_ext/Lapack.h`
- **Morton encoding**: `mortonEncode3D()` interleaves x/y/z bits for Z-curve spatial ordering

## Dependencies

### External
- `Rcpp` — R/C++ interface
- `RcppEigen` — Eigen linear algebra library (header-only)
- `R_ext/Lapack.h` — LAPACK routines bundled with R
- OpenMP — Optional thread parallelism (detected at compile time)

<!-- MANUAL: Any manually added notes below this line are preserved on regeneration -->
