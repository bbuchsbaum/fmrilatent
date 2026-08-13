# Create a PCA basis specification

\`basis_pca()\` creates a lightweight descriptor for parcel- or
cluster-local PCA bases. The descriptor is consumed by
\`lift(ClusterReduction, spec_pca, data = ...)\`, where \`data\`
supplies the time-by-voxel matrix used to estimate the components.

## Usage

``` r
basis_pca(k = 3, whiten = FALSE)
```

## Arguments

- k:

  Positive integer number of PCA components.

- whiten:

  Logical scalar. \`basis_pca()\` records this request, but
  \`lift(ClusterReduction, spec_pca)\` returns unwhitened loadings and
  emits a warning when \`whiten = TRUE\`; the higher-level
  \`encode(spec_space_pca())\` path applies whitening after projection.

## Value

A list with class \`spec_pca\` containing \`k\` and \`whiten\`.
