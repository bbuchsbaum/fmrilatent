# Average subject conductance matrices on a shared template graph

Average subject conductance matrices on a shared template graph

## Usage

``` r
awpt_mean_conductance(
  conductances,
  method = c("log_euclidean", "arithmetic"),
  shrinkage = 0,
  enforce_psd = TRUE,
  tol = 1e-08
)
```

## Arguments

- conductances:

  List of symmetric conductance matrices on the same template graph.

- method:

  Averaging rule.

- shrinkage:

  Optional shrinkage toward the isotropic identity.

- enforce_psd:

  Logical; if `TRUE`, project the result to the PSD cone.

- tol:

  Numerical tolerance for SPD operations.

## Value

A symmetric averaged conductance matrix.
