# Wrap flat decoded outputs into a domain-native representation

Wrap flat decoded outputs into a domain-native representation

## Usage

``` r
wrap_decoded(x, values, time_idx = NULL, space = c("native", "template"), ...)

# S4 method for class 'ImplicitLatent'
wrap_decoded(x, values, time_idx = NULL, space = c("native", "template"), ...)

# S4 method for class 'LatentNeuroSurfaceVector'
wrap_decoded(x, values, time_idx = NULL, space = c("native", "template"), ...)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
wrap_decoded(x, values, time_idx = NULL, space = c("native", "template"), ...)

# S4 method for class 'BlockLatentNeuroVector'
wrap_decoded(x, values, time_idx = NULL, space = c("native", "template"), ...)

# S4 method for class 'LatentNeuroVec'
wrap_decoded(x, values, time_idx = NULL, space = c("native", "template"), ...)
```

## Arguments

- x:

  A latent object.

- values:

  Flat decoded values, typically as a vector or matrix.

- time_idx:

  Optional integer time indices associated with `values`.

- space:

  Output space to wrap into.

- ...:

  Additional arguments passed to methods.

## Value

Domain-native wrapped representation.
