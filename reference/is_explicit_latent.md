# Test whether a latent object is explicit

Test whether a latent object is explicit

## Usage

``` r
is_explicit_latent(x, ...)

# S4 method for class 'ExplicitLatent'
is_explicit_latent(x, ...)

# S4 method for class 'ImplicitLatent'
is_explicit_latent(x, ...)
```

## Arguments

- x:

  A latent object.

- ...:

  Additional arguments (unused).

## Value

Logical scalar; `TRUE` when the object stores explicit basis and
loadings matrices.
