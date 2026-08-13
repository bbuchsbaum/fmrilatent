# Extract coefficient time series from a latent object

Extract coefficient time series from a latent object

## Usage

``` r
coef_time(x, coordinates = c("analysis", "raw"), ...)

# S4 method for class 'LatentNeuroSurfaceVector'
coef_time(x, coordinates = c("analysis", "raw"), ...)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
coef_time(x, coordinates = c("analysis", "raw"), ...)

# S4 method for class 'ImplicitLatent'
coef_time(x, coordinates = c("analysis", "raw"), ...)

# S4 method for class 'LatentNeuroVec'
coef_time(x, coordinates = c("analysis", "raw"), ...)

# S4 method for class 'BlockLatentNeuroVector'
coef_time(x, coordinates = c("analysis", "raw"), ...)
```

## Arguments

- x:

  A latent object.

- coordinates:

  Coordinate system for the returned coefficients.

- ...:

  Additional arguments passed to methods.

## Value

Numeric matrix with rows = time and columns = latent coefficients.
