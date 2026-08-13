# Reconstruct a matrix from a low-rank SVD baseline

Reconstruct a matrix from a low-rank SVD baseline

## Usage

``` r
boldzip_svd_reconstruct(X, rank, center = TRUE)
```

## Arguments

- X:

  Numeric matrix with dimensions \`voxels x time\`.

- rank:

  SVD rank.

- center:

  Whether to remove and restore row means.

## Value

Matrix with the same dimensions as \`X\`.
