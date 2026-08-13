# Heat wavelet basis specification

Heat wavelet basis specification

## Usage

``` r
basis_heat_wavelet(scales = c(1, 2, 4, 8), order = 30, threshold = 1e-06)
```

## Arguments

- scales:

  Numeric vector of heat diffusion scales.

- order:

  Polynomial approximation order.

- threshold:

  Threshold for small coefficients.

## Value

A \`spec_heat_wavelet\` object for \`lift()\`.
