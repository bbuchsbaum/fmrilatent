# Get the offset vector

Extract the offset vector from a latent space representation. For
`LatentNeuroVec` objects, this returns the voxel-wise offset (mean or
intercept) that is added after the basis x loadings reconstruction.

## Usage

``` r
# S4 method for class 'ANY'
offset(object, ...)

# S4 method for class 'LatentNeuroSurfaceVector'
offset(object, ...)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
offset(object, ...)

# S4 method for class 'BlockLatentNeuroVector'
offset(object, ...)

# S4 method for class 'LatentNeuroVec'
offset(object, ...)
```

## Arguments

- object:

  An object containing an offset (e.g., `LatentNeuroVec`)

- ...:

  Additional arguments for methods. The generic keeps the `object` first
  argument used by
  [`stats::offset()`](https://rdrr.io/r/stats/offset.html) and adds
  `...` so offset accessors can follow the same extension pattern as
  sibling accessors such as [`basis()`](basis-methods.md) and
  [`loadings()`](loadings-methods.md).

## Value

The offset vector (length = number of voxels in mask)

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
off_vector <- offset(lvec)
length(off_vector)
#> [1] 0
```
