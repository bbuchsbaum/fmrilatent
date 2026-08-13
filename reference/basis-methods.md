# Get the basis matrix (temporal components)

Extract the basis matrix from a latent space representation. For
`LatentNeuroVec` objects, this returns the temporal basis matrix with
dimensions (nTime x k) where k is the number of components.

## Usage

``` r
basis(x, ...)

# S4 method for class 'LatentNeuroSurfaceVector'
basis(x)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
basis(x)

# S4 method for class 'BlockLatentNeuroVector'
basis(x)

# S4 method for class 'LatentNeuroVec'
basis(x)
```

## Arguments

- x:

  An object containing a basis matrix (e.g., `LatentNeuroVec`)

- ...:

  Additional arguments (currently unused)

## Value

The basis matrix (typically time x components)

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
b_matrix <- basis(lvec)
dim(b_matrix)
#> [1] 3 2
```
