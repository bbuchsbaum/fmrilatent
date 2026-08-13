# Compare BOLDZip-SR against simple reconstruction baselines

Compare BOLDZip-SR against simple reconstruction baselines

## Usage

``` r
compare_boldzip_sr(X, fit, parcels = NULL, svd_ranks = NULL)
```

## Arguments

- X:

  Original matrix with rows as voxels/grayordinates and columns as time
  points.

- fit:

  A \`BoldZipSR\` object.

- parcels:

  Optional parcel labels for a parcel-average baseline.

- svd_ranks:

  Optional integer vector of SVD ranks.

## Value

Data frame with method, MSE, RMSE, correlation, and scalar payload
estimates.
