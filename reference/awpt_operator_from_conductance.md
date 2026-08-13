# Build an AWPT field operator from a conductance matrix

Build an AWPT field operator from a conductance matrix

## Usage

``` r
awpt_operator_from_conductance(
  conductance,
  normalize = c("none", "sym", "rw"),
  tol = 1e-08
)
```

## Arguments

- conductance:

  Symmetric conductance matrix on the template graph.

- normalize:

  Laplacian normalization convention.

- tol:

  Numerical tolerance.

## Value

A field-space roughness operator.
