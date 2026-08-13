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
  meta = list(),
  expect_dense = FALSE
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

- expect_dense:

  Logical; if \`TRUE\`, suppress the informational message emitted when
  a base-matrix \`basis\`/\`loadings\` is dense (\>50 families (e.g.
  diffusion wavelets) where a dense factor is by design. Defaults to
  \`FALSE\`, preserving the original message behavior.

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
# Example data
n_timepoints <- 4
n_components <- 2
mask_array <- array(TRUE, dim = c(2, 2, 1))
n_voxels <- sum(mask_array)

# Create basis & loadings
basis <- Matrix::Matrix(seq_len(n_timepoints * n_components),
  nrow = n_timepoints,
  ncol = n_components
)
loadings <- Matrix::Matrix(seq_len(n_voxels * n_components) / 10,
  nrow = n_voxels,
  ncol = n_components,
  sparse = TRUE
)

# Create space (2x2x1 volume, 4 timepoints)
spc <- neuroim2::NeuroSpace(c(2, 2, 1, n_timepoints))

# Create mask
mask_vol <- neuroim2::LogicalNeuroVol(mask_array, neuroim2::NeuroSpace(c(2, 2, 1)))

# Construct LatentNeuroVec
lvec <- LatentNeuroVec(
  basis = basis,
  loadings = loadings,
  space = spc,
  mask = mask_vol,
  expect_dense = TRUE
)
dim(lvec)
#> [1] 2 2 1 4
```
