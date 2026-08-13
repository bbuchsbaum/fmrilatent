# Compatibility wrapper for decoder-based covariance pushforward

Compatibility wrapper for decoder-based covariance pushforward

## Usage

``` r
project_vcov(
  x,
  Sigma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  diag_only = TRUE,
  ...
)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
project_vcov(
  x,
  Sigma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  diag_only = TRUE,
  ...
)

# S4 method for class 'ImplicitLatent'
project_vcov(
  x,
  Sigma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  diag_only = TRUE,
  ...
)

# S4 method for class 'LatentNeuroSurfaceVector'
project_vcov(
  x,
  Sigma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  diag_only = TRUE,
  ...
)

# S4 method for class 'LatentNeuroVec'
project_vcov(
  x,
  Sigma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  diag_only = TRUE,
  ...
)

# S4 method for class 'BlockLatentNeuroVector'
project_vcov(
  x,
  Sigma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  diag_only = TRUE,
  ...
)
```

## Arguments

- x:

  A latent object.

- Sigma:

  Coefficient covariance matrix.

- space:

  Output space to decode into.

- coordinates:

  Coordinate system used by `Sigma`.

- diag_only:

  Logical; if `TRUE`, return only the diagonal.

- ...:

  Additional arguments passed to methods.

## Value

Numeric vector or matrix in the requested output space.
