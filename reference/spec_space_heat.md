# Spatial heat-wavelet spec (graph diffusion)

Spatial heat-wavelet spec (graph diffusion)

## Usage

``` r
spec_space_heat(
  scales = c(1, 2, 4, 8),
  order = 30L,
  threshold = 1e-06,
  k_neighbors = 6L
)
```

## Arguments

- scales:

  Heat scales.

- order:

  Polynomial order.

- threshold:

  Threshold for small coefficients.

- k_neighbors:

  k-NN graph parameter.

## Value

A \`spec_space_heat\` object.
