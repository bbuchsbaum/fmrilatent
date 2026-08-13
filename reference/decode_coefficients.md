# Decode coefficient-space vectors into an output space

Decode coefficient-space vectors into an output space

## Usage

``` r
decode_coefficients(
  x,
  gamma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  wrap = c("none", "auto"),
  ...
)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
decode_coefficients(
  x,
  gamma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  wrap = c("none", "auto"),
  ...
)

# S4 method for class 'ImplicitLatent'
decode_coefficients(
  x,
  gamma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  wrap = c("none", "auto"),
  ...
)

# S4 method for class 'LatentNeuroSurfaceVector'
decode_coefficients(
  x,
  gamma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  wrap = c("none", "auto"),
  ...
)

# S4 method for class 'LatentNeuroVec'
decode_coefficients(
  x,
  gamma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  wrap = c("none", "auto"),
  ...
)

# S4 method for class 'BlockLatentNeuroVector'
decode_coefficients(
  x,
  gamma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  wrap = c("none", "auto"),
  ...
)
```

## Arguments

- x:

  A latent object.

- gamma:

  Coefficient-space vector or matrix.

- space:

  Output space to decode into.

- coordinates:

  Coordinate system used by `gamma`.

- wrap:

  If `"auto"`, wrap the decoded values with
  [`wrap_decoded()`](wrap_decoded.md) so the result is domain-aware (for
  example a `NeuroVol` for volumetric targets or a surface field for
  surface targets). Defaults to `"none"`, which returns the raw numeric
  vector or matrix.

- ...:

  Additional arguments passed to methods.

## Value

Numeric vector or matrix in the requested output space, or a
domain-aware wrapper when `wrap = "auto"`.
