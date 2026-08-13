# Spatial heat-wavelet spec (graph diffusion)

Spatial heat-wavelet spec (graph diffusion)

## Usage

``` r
spec_space_heat(
  scales = c(1, 2, 4, 8),
  order = 30L,
  threshold = NULL,
  k_neighbors = 6L,
  sparsify_eps = NULL
)
```

## Arguments

- scales:

  Heat scales.

- order:

  Polynomial order.

- threshold:

  Deprecated alias for \`sparsify_eps\`.

- k_neighbors:

  k-NN graph parameter.

- sparsify_eps:

  Non-negative threshold for small heat-wavelet coefficients. This is
  stored as \`threshold\` for compatibility.

## Value

A \`spec_space_heat\` object.
