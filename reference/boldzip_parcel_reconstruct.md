# Reconstruct a matrix from parcel-average time series

Reconstruct a matrix from parcel-average time series

## Usage

``` r
boldzip_parcel_reconstruct(X, parcels)
```

## Arguments

- X:

  Numeric matrix with dimensions \`voxels x time\`.

- parcels:

  Integer, factor, or character vector with one parcel label per row of
  \`X\`.

## Value

Matrix with the same dimensions as \`X\`, expanded from parcel means.
