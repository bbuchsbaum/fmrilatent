# Concatenate LatentNeuroVec Objects

Concatenates two or more `LatentNeuroVec` objects along the temporal
dimension.

## Usage

``` r
# S4 method for class 'LatentNeuroVec,LatentNeuroVec'
concat(x, y, ...)
```

## Arguments

- x:

  First `LatentNeuroVec`.

- y:

  Second `LatentNeuroVec`.

- ...:

  Additional `LatentNeuroVec` objects to concatenate.

## Value

A new `LatentNeuroVec` if all objects are compatible, otherwise a
[`NeuroVecSeq-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVecSeq-class.html).

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
combined <- concat(lvec, lvec)
dim(combined)
#> [1] 2 2 1 6
```
