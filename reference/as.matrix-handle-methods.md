# Coerce latent handles to dense matrices

Coerce latent handles to dense matrices

## Usage

``` r
# S4 method for class 'BasisHandle'
as.matrix(x, ...)

# S4 method for class 'LoadingsHandle'
as.matrix(x, ...)
```

## Arguments

- x:

  A `BasisHandle` or `LoadingsHandle`.

- ...:

  Optional row/column subsetting arguments passed to `basis_mat()` or
  `loadings_mat()`.

## Value

A base dense matrix.
