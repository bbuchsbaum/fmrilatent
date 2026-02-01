# Create a Latent Space Representation of Neuroimaging Data

Constructs a [`LatentNeuroVec-class`](LatentNeuroVec-class.md) object,
which provides a memory-efficient representation of neuroimaging data
using matrix factorization. This is particularly useful for
dimensionality reduction techniques (e.g., PCA or ICA).

## Usage

``` r
LatentNeuroVec(
  basis,
  loadings,
  space,
  mask,
  offset = NULL,
  label = "",
  meta = list()
)
```

## Arguments

- basis:

  A numeric or `Matrix` object (\\n \times k\\) containing the temporal
  basis. Each row corresponds to a time point, each column to a
  component.

- loadings:

  A numeric or `Matrix` object (\\p \times k\\) containing spatial
  loadings. Each row corresponds to a voxel within the mask, each column
  to a component.

- space:

  A
  [`NeuroSpace-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroSpace-class.html)
  defining the spatial/temporal dimensions. Must be 4-dimensional.

- mask:

  A
  [`LogicalNeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.html)
  defining the brain mask. The number of TRUE values must equal
  nrow(loadings).

- offset:

  Optional numeric vector of length \\p\\ (voxel-wise offsets). If NULL,
  defaults to zero offset.

- label:

  Optional character label for the object.

- meta:

  Optional list of metadata (e.g., HRBF params or centres).

## Value

A new [`LatentNeuroVec-class`](LatentNeuroVec-class.md) instance.

## Details

Construct a LatentNeuroVec Object

The data is represented as the product: \$\$X = B \times L^T + c\$\$
where:

- B is the basis matrix (\\n \times k\\)

- L is the loadings matrix (\\p \times k\\)

- c is an optional offset vector (length p)

- n is the number of time points

- p is the number of voxels in the mask

- k is the number of components

## Examples

``` r
if (FALSE) { # \dontrun{
library(Matrix)
library(neuroim2)

# Example data
n_timepoints <- 100
n_components <- 10
n_voxels <- 1000

# Create basis & loadings
basis <- Matrix(rnorm(n_timepoints * n_components),
  nrow = n_timepoints,
  ncol = n_components
)
loadings <- Matrix(rnorm(n_voxels * n_components),
  nrow = n_voxels,
  ncol = n_components,
  sparse = TRUE
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
} # }
```
