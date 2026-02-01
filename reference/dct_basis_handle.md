# Create a BasisHandle for a DCT temporal basis

Create a BasisHandle for a DCT temporal basis

## Usage

``` r
dct_basis_handle(n_time, k, norm = c("ortho", "none"), id = NULL, label = NULL)
```

## Arguments

- n_time:

  Number of time points.

- k:

  Number of components.

- norm:

  Normalization passed to \[build_dct_basis()\].

- id:

  Optional registry key; if NULL a deterministic id is derived from
  parameters.

- label:

  Optional human-readable label.

## Value

A `BasisHandle` object.
