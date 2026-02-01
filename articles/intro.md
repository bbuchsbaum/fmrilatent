# Introduction to fmrilatent

``` r
library(fmrilatent)
```

## The problem

An fMRI volume is mostly empty space—air, skull, tissue outside the
brain. Standard compression and denoising methods ignore this structure
and operate on the full dense array. This wastes memory and computation
on voxels that carry no signal.

fmrilatent represents brain data as a *latent factorization*: a low-rank
decomposition that respects the brain mask. You choose a basis family,
call [`encode()`](../reference/encode.md), and receive an object that
reconstructs the original data on demand—partially or fully—without ever
materializing the full array unless you ask for it.

## The interface

The API is declarative. You specify *what* you want, not *how* to
compute it:

1.  **Temporal bases** reduce the time dimension: Slepian/DPSS sequences
    (bandlimited), DCT (frequency), or B-splines (smooth functions).

2.  **Spatial bases** reduce the voxel dimension: graph Slepians
    (spectrally concentrated on regions), heat wavelets (multiscale
    diffusion), HRBF atoms (smooth radial functions), or CDF 5/3 lifting
    wavelets (fast, mask-aware).

3.  **Separable spatiotemporal bases** combine one temporal and one
    spatial family, storing only a small core tensor.

Under the hood, bases can be stored as explicit matrices or as lazy
“handles” that generate columns on demand. The user-facing call is
identical either way.

## A minimal example

``` r
# A toy 4x4x1 mask (all voxels in-mask)
mask <- array(TRUE, dim = c(4, 4, 1))
mask_vol <- neuroim2::LogicalNeuroVol(mask, neuroim2::NeuroSpace(dim(mask)))

# Random data: 6 time points, 16 voxels
X <- matrix(rnorm(6 * sum(mask)), nrow = 6)

# Encode with a Slepian temporal basis (4 components, TR = 2s, bandwidth = 0.08 Hz)
spec <- spec_time_slepian(tr = 2, bandwidth = 0.08, k = 4)
lat <- encode(X, spec, mask = mask_vol)

# Reconstruct and check error
recon <- as.matrix(lat)
sqrt(mean((recon - X)^2))
```

The latent object `lat` stores a temporal basis (6 × 4) and spatial
loadings (16 × 4). Reconstruction multiplies these factors; no 6 × 16
dense matrix is kept in memory.

## Choosing a basis family

The choice depends on your analysis goals:

| Family | Strengths | When to use |
|----|----|----|
| Slepian (time) | Optimal bandlimiting; controls spectral leakage | Resting-state, low-frequency fluctuations |
| DCT (time) | Fast; interpretable frequency components | General denoising, frequency-domain analyses |
| B-spline (time) | Smooth, localized; adjustable flexibility | HRF modeling, slow drifts |
| Slepian (space) | Graph-spectral concentration; sparse on regions | Parcellation-aware compression |
| Heat wavelet (space) | Multiscale; captures local and global structure | Exploratory spatial analysis |
| HRBF (space) | Analytic smooth atoms; tunable overlap | When you want explicit spatial basis functions |
| Wavelet active (space) | Fast lifting; touches only in-mask voxels | Large masks, speed-critical pipelines |
| Separable (time + space) | Compact core; partial decoding | Full 4D compression with selective reconstruction |

## Cheat sheet: storage, speed, accuracy

| Family | Stored | Implicit | Accuracy | Speed |
|----|----|----|----|----|
| slepian_space | basis (time × k), loadings (vox × k) | none (explicit) | exact on span; orthogonal per cluster | fast; sparse loadings |
| hrbf | basis (time × atoms), loadings (vox × atoms) | none (explicit) | improves as atoms ↑ | moderate; dense atoms |
| wavelet_active | coeff (time × vox) | decoder (active pencils) | exact (biorthogonal 5/3) | very fast; mask-aware |
| bspline_hrbf_st | core (kt × ks), B_t, L_s | decoder (separable) | approx → better with atoms | fast; separable |
| slepian_st | core (kt × ks), B_t, L_s | decoder (separable) | exact on span of bases | fast; sparse spatial |
| heat_space | basis, loadings | none (explicit) | approx (thresholded) | graph build cost; sparse |

## Reconstruction and partial access

Latent objects support several reconstruction paths:

- `as.matrix(lat)` returns the full time × voxels matrix
- `as.array(lat)` returns a 4D array matching the original volume
  dimensions
- For implicit (separable) representations, the decoder accepts
  `time_idx` and `roi_mask` arguments to reconstruct only what you need

Partial access is useful when you want, say, a single time point or a
region of interest, without paying the cost of full reconstruction.

## Checking performance

The helper
[`benchmark_roundtrip()`](../reference/benchmark_roundtrip.md) times the
encode-decode cycle and reports reconstruction error. This lets you
compare families on your own data and hardware:

``` r
res <- benchmark_roundtrip(
  mask_dims = c(16, 16, 8),
  n_time = 100,
  methods = c("slepian_time", "dct_time", "wavelet_active")
)
res
```

Visualization
helpers—[`plot_slepian_temporal()`](../reference/plot_slepian_temporal.md),
[`plot_basis_gram()`](../reference/plot_basis_gram.md),
[`plot_spatial_atom()`](../reference/plot_spatial_atom.md)—let you
inspect bases without manually extracting arrays.

## Summary

fmrilatent provides mask-aware latent representations for fMRI data. The
workflow is:

1.  Build a spec describing the basis family and parameters
2.  Call `encode(X, spec, mask = ...)`
3.  Use [`as.matrix()`](https://rdrr.io/r/base/matrix.html),
    [`as.array()`](https://rdrr.io/r/base/array.html), or partial
    decoding to recover data

The implementation handles the linear algebra; you describe what you
want.
