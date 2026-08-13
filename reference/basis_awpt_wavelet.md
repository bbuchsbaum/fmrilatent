# AWPT wave-packet basis specification

AWPT wave-packet basis specification

## Usage

``` r
basis_awpt_wavelet(
  scales = c(1, 2, 4, 8),
  order = 30L,
  threshold = 1e-06,
  k_neighbors = 6L,
  penalty_rule = c("inverse_scale", "inverse_scale_sq", "scale", "none", "custom"),
  custom_weights = NULL
)
```

## Arguments

- scales:

  Numeric vector of anatomical wave-packet scales.

- order:

  Polynomial approximation order for the underlying heat-wavelet
  construction.

- threshold:

  Threshold for small coefficients.

- k_neighbors:

  Graph neighborhood parameter.

- penalty_rule:

  Rule used to convert scales into roughness weights.

- custom_weights:

  Optional explicit weights matching `scales`.

## Value

A `spec_awpt_wavelet` object.
