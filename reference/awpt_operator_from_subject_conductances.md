# Build an AWPT field operator from subject conductance summaries

Build an AWPT field operator from subject conductance summaries

## Usage

``` r
awpt_operator_from_subject_conductances(
  conductances,
  mean_method = c("log_euclidean", "arithmetic"),
  normalize = c("none", "sym", "rw"),
  shrinkage = 0,
  enforce_psd = TRUE,
  tol = 1e-08
)
```

## Arguments

- conductances:

  List of subject conductance matrices on a shared template graph.

- mean_method:

  Averaging rule for the conductance mean.

- normalize:

  Laplacian normalization convention.

- shrinkage:

  Optional shrinkage toward an isotropic baseline before Laplacian
  construction.

- enforce_psd:

  Logical; if `TRUE`, project the mean conductance to the PSD cone.

- tol:

  Numerical tolerance.

## Value

A list with `conductance_mean` and `operator`.
