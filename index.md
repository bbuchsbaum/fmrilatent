# fmrilatent

Memory-efficient latent space representations for fMRI data in R.

## Overview

**fmrilatent** provides S4 classes and methods for representing
neuroimaging data using matrix factorization. The `LatentNeuroVec` class
stores 4D fMRI data as:

    X[v,t] = Σ basis[t,k] × loadings[v,k] + offset[v]

This enables compact storage and efficient computation for
dimensionality-reduced fMRI data (PCA, ICA, wavelet decompositions,
etc.).

## Installation

``` r
# Install from GitHub
remotes::install_github("bbuchsbaum/fmrilatent")
```

### Dependencies

The package requires [neuroim2](https://github.com/bbuchsbaum/neuroim2)
for neuroimaging data structures:

``` r
remotes::install_github("bbuchsbaum/neuroim2")
```

## Quick Start

### Basic Usage

``` r
library(fmrilatent)
library(neuroim2)
library(Matrix)

# Create example data
n_time <- 100
n_components <- 10
n_voxels <- 1000

# Temporal basis (time x components)
basis <- Matrix(rnorm(n_time * n_components), nrow = n_time)

# Spatial loadings (voxels x components)
loadings <- Matrix(rnorm(n_voxels * n_components), nrow = n_voxels, sparse = TRUE)

# Create 4D space and mask
space <- NeuroSpace(c(10, 10, 10, n_time))
mask <- LogicalNeuroVol(array(TRUE, c(10, 10, 10)), NeuroSpace(c(10, 10, 10)))

# Construct LatentNeuroVec
lvec <- LatentNeuroVec(
  basis = basis,
  loadings = loadings,
  space = space,
  mask = mask
)

# Access components
basis(lvec)      # Get temporal basis
loadings(lvec)   # Get spatial loadings
lvec[[1]]        # Extract first volume
```

### Using the Encoding API

The [`encode()`](reference/encode.md) function provides a high-level
interface for creating latent representations:

``` r
# Create data matrix (time x voxels)
X <- matrix(rnorm(n_time * n_voxels), nrow = n_time)

# Encode with temporal Slepian basis
lvec <- encode(X, spec_time_slepian(tr = 2, bandwidth = 0.1), mask = mask)

# Encode with DCT basis
lvec_dct <- encode(X, spec_time_dct(k = 15), mask = mask)

# Spatiotemporal encoding (separable)
lvec_st <- encode(X, 
  spec_st(
    time = spec_time_bspline(k = 10),
    space = spec_space_slepian(k = 3)
  ),
  mask = mask
)
```

### Available Basis Families

**Temporal bases:** -
[`spec_time_slepian()`](reference/spec_time_slepian.md) - DPSS/Slepian
sequences (bandlimited) -
[`spec_time_dct()`](reference/spec_time_dct.md) - Discrete cosine
transform - [`spec_time_bspline()`](reference/spec_time_bspline.md) -
B-spline basis

**Spatial bases:** -
[`spec_space_slepian()`](reference/spec_space_slepian.md) - Graph
Laplacian eigenvectors -
[`spec_space_heat()`](reference/spec_space_heat.md) - Heat/diffusion
wavelets - [`spec_space_hrbf()`](reference/spec_space_hrbf.md) -
Hierarchical radial basis functions

**Combined:** - [`spec_st()`](reference/spec_st.md) - Separable
spatiotemporal (time × space)

### Lazy Evaluation with Handles

For memory efficiency, bases can be computed lazily:

``` r
# Create with handle (lazy evaluation)
lvec <- encode(X, spec_time_slepian(tr = 2), mask = mask, materialize = "handle")

# Basis is computed on first access
b <- basis(lvec)  # Materializes and caches the basis
```

## Key Features

- **Memory efficient**: Store factorized representation instead of full
  4D array
- **Lazy evaluation**: Compute bases on-demand with caching
- **Multiple basis families**: Temporal, spatial, and spatiotemporal
  decompositions
- **neuroim2 integration**: Works with standard neuroimaging data
  structures
- **Sparse support**: Efficient handling of sparse loadings matrices

## Documentation

- [`vignette("intro", package = "fmrilatent")`](articles/intro.md) -
  Introduction and concepts
- [`vignette("encode-factory", package = "fmrilatent")`](articles/encode-factory.md) -
  Encoding API reference

## License

GPL (\>= 3)
