# LatentNeuroVec Class

A class that represents a 4-dimensional neuroimaging array using a
latent space decomposition. It stores the data as a set of basis
functions (dictionary) and a corresponding set of loadings
(coefficients), enabling efficient representation and manipulation of
high-dimensional data.

## Details

`LatentNeuroVec` inherits from
[`NeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html)
and
[`AbstractSparseNeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/AbstractSparseNeuroVec-class.html).
The original 4D data can be reconstructed as: \$\$data\[v,t\] = \sum_k
\bigl(basis\[t,k\] \times loadings\[v,k\]\bigr) + offset\[v\]\$\$ where
v indexes voxels within the mask and t indexes time points.

This approach is especially useful for large datasets where storing the
full 4D array is expensive. Common use cases include:

- PCA-reduced fMRI data

- ICA decompositions

- Dictionary learning representations

- Any matrix factorization of neuroimaging data

## Slots

- `basis`:

  A `Matrix` or `BasisHandle` where each column represents a basis
  vector in the latent space. Dimensions are (nTime x k) where k is the
  number of components.

- `loadings`:

  A `Matrix` or `LoadingsHandle` (often sparse) containing the
  coefficients for each basis vector across the spatial dimensions.
  Dimensions are (nVoxels x k).

- `offset`:

  A `numeric` vector representing a constant offset term for each voxel
  or spatial location. Length is nVoxels (within mask).

- `map`:

  A `IndexLookupVol` object representing the mapping from 3D coordinates
  to linear indices within the mask.

- `label`:

  A `character` string representing the label for the latent vector.

- `meta`:

  A `list` storing lightweight metadata (e.g., HRBF params, centers).

## Inheritance

`LatentNeuroVec` inherits from:

- [`NeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html) -
  base 4D neuroimaging class

- [`AbstractSparseNeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/AbstractSparseNeuroVec-class.html) -
  sparse representation framework

## See also

[`NeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html),
[`AbstractSparseNeuroVec-class`](https://bbuchsbaum.github.io/neuroim2/reference/AbstractSparseNeuroVec-class.html),
[`LatentNeuroVec`](LatentNeuroVec.md) constructor function.

## Examples

``` r
if (FALSE) { # \dontrun{
library(Matrix)
library(neuroim2)

# Example dimensions
n_timepoints <- 100
n_components <- 10
n_voxels <- 1000

# Create basis (temporal) & loadings (spatial)
basis <- Matrix(rnorm(n_timepoints * n_components),
  nrow = n_timepoints, ncol = n_components
)
loadings <- Matrix(rnorm(n_voxels * n_components),
  nrow = n_voxels, ncol = n_components, sparse = TRUE
)

# Create space (10x10x10 volume, 100 timepoints)
spc <- NeuroSpace(c(10, 10, 10, n_timepoints))

# Create mask
mask_array <- array(TRUE, dim = c(10, 10, 10))
mask_vol <- LogicalNeuroVol(mask_array, NeuroSpace(c(10, 10, 10)))

# Construct LatentNeuroVec
lvec <- LatentNeuroVec(
  basis = basis,
  loadings = loadings,
  space = spc,
  mask = mask_vol
)

print(lvec)
print(dim(lvec))
} # }
```
