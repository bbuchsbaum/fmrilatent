# Compatibility wrapper for decoder-based coefficient projection

Compatibility wrapper for decoder-based coefficient projection

## Usage

``` r
project_effect(
  x,
  gamma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  ...
)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
project_effect(
  x,
  gamma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  ...
)

# S4 method for class 'ImplicitLatent'
project_effect(
  x,
  gamma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  ...
)

# S4 method for class 'LatentNeuroSurfaceVector'
project_effect(
  x,
  gamma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  ...
)

# S4 method for class 'LatentNeuroVec'
project_effect(
  x,
  gamma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
  ...
)

# S4 method for class 'BlockLatentNeuroVector'
project_effect(
  x,
  gamma,
  space = c("native", "template"),
  coordinates = c("analysis", "raw"),
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

- ...:

  Additional arguments passed to methods.

## Value

Numeric vector or matrix in the requested output space.
