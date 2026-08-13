# Push coefficient covariance into an output space

Push coefficient covariance into an output space

## Usage

``` r
decode_covariance(
  x,
  Sigma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  diag_only = TRUE,
  ...
)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
decode_covariance(
  x,
  Sigma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  diag_only = TRUE,
  ...
)

# S4 method for class 'ImplicitLatent'
decode_covariance(
  x,
  Sigma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  diag_only = TRUE,
  ...
)

# S4 method for class 'LatentNeuroSurfaceVector'
decode_covariance(
  x,
  Sigma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  diag_only = TRUE,
  ...
)

# S4 method for class 'LatentNeuroVec'
decode_covariance(
  x,
  Sigma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  diag_only = TRUE,
  ...
)

# S4 method for class 'BlockLatentNeuroVector'
decode_covariance(
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
