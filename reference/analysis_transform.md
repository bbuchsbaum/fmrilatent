# Describe the transform from raw to analysis coordinates

Describe the transform from raw to analysis coordinates

## Usage

``` r
analysis_transform(x, ...)

# S4 method for class 'LatentNeuroSurfaceVector'
analysis_transform(x, ...)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
analysis_transform(x, ...)

# S4 method for class 'ImplicitLatent'
analysis_transform(x, ...)

# S4 method for class 'LatentNeuroVec'
analysis_transform(x, ...)

# S4 method for class 'BlockLatentNeuroVector'
analysis_transform(x, ...)
```

## Arguments

- x:

  A latent object.

- ...:

  Additional arguments passed to methods.

## Value

Transform descriptor or `NULL`. Downstream model-fitting code should
ordinarily consume `coef_time(x, "analysis")` rather than reasoning
about raw coordinates directly.
