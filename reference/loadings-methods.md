# Get the loadings matrix (spatial components)

Extract the loadings matrix from a latent space representation. For
`LatentNeuroVec` objects, this returns the spatial loadings matrix with
dimensions (nVoxels x k) where k is the number of components and nVoxels
is the number of voxels within the mask.

## Usage

``` r
loadings(x, ...)

# S4 method for class 'LatentNeuroSurfaceVector'
loadings(x)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
loadings(x)

# S4 method for class 'BlockLatentNeuroVector'
loadings(x)

# S4 method for class 'LatentNeuroVec'
loadings(x)
```

## Arguments

- x:

  An object containing loadings (e.g., `LatentNeuroVec`)

- ...:

  Additional arguments (currently unused)

## Value

The loadings matrix (typically voxels x components)

## Examples

``` r
mask <- neuroim2::LogicalNeuroVol(
  array(TRUE, dim = c(2, 2, 1)),
  neuroim2::NeuroSpace(c(2, 2, 1))
)
lvec <- LatentNeuroVec(
  basis = matrix(1:6, nrow = 3),
  loadings = matrix(seq_len(8) / 10, nrow = 4),
  space = neuroim2::NeuroSpace(c(2, 2, 1, 3)),
  mask = mask,
  expect_dense = TRUE
)
l_matrix <- loadings(lvec)
dim(l_matrix)
#> [1] 4 2
```
