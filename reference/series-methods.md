# Extract time series from LatentNeuroVec

Extract time series data from a `LatentNeuroVec`. Reconstructs values
on-the-fly from the basis and loadings factorization.

Extract time series data from a `LatentNeuroVec` at specified spatial
coordinates or voxel indices.

## Usage

``` r
# S4 method for class 'LatentNeuroVec,integer'
series(x, i, j, k, ..., drop = TRUE)

# S4 method for class 'LatentNeuroVec,numeric'
series(x, i, j, k, ..., drop = TRUE)

# S4 method for class 'LatentNeuroVec,ANY'
series(x, i, j, k, ..., drop = TRUE)
```

## Arguments

- x:

  A `LatentNeuroVec` object

- i, j, k:

  Spatial indices (x, y, z coordinates)

- ...:

  Additional arguments passed to methods

- drop:

  Logical; drop dimensions of length 1 (default TRUE)

## Value

Matrix or vector of time series values

A matrix or vector of time series values
