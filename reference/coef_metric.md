# Extract coefficient-space metric from a latent object

Extract coefficient-space metric from a latent object

## Usage

``` r
coef_metric(x, coordinates = c("analysis", "raw"), ...)

# S4 method for class 'LatentNeuroSurfaceVector'
coef_metric(x, coordinates = c("analysis", "raw"), ...)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
coef_metric(x, coordinates = c("analysis", "raw"), ...)

# S4 method for class 'ImplicitLatent'
coef_metric(x, coordinates = c("analysis", "raw"), ...)

# S4 method for class 'LatentNeuroVec'
coef_metric(x, coordinates = c("analysis", "raw"), ...)

# S4 method for class 'BlockLatentNeuroVector'
coef_metric(x, coordinates = c("analysis", "raw"), ...)
```

## Arguments

- x:

  A latent object.

- coordinates:

  Coordinate system for the requested metric.

- ...:

  Additional arguments passed to methods.

## Value

Matrix-like metric representation or `NULL`. For transport-backed latent
objects, analysis coordinates are Euclidean by contract in v1.
Raw-coordinate metrics are returned when the raw-to-analysis transform
exposes a linear matrix representation.
