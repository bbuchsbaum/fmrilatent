# Build a spectral graph spatial basis for BOLDZip-SR

Constructs an experimental graph-adapted spatial basis from a symmetric
adjacency matrix. Low-frequency graph Laplacian eigenvectors become
coarse carrier support functions and the following higher-frequency
eigenvectors become detail atoms. This is a materialized small-to-medium
graph MVP, not a production graph-wavelet operator.

## Usage

``` r
boldzip_graph_spatial_basis(
  adjacency,
  n_coarse = 8L,
  n_detail = NULL,
  normalized = TRUE,
  label = NULL
)
```

## Arguments

- adjacency:

  Symmetric non-negative adjacency matrix.

- n_coarse:

  Number of low-frequency eigenvectors to use as \`phi_c\`.

- n_detail:

  Number of following eigenvectors to use as \`phi_d\`. Defaults to all
  remaining eigenvectors.

- normalized:

  Whether to use the normalized graph Laplacian.

- label:

  Optional label stored in metadata.

## Value

A \`BoldZipSRSpatialBasis\` object.
