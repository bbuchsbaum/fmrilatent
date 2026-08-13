# Temporal B-spline spec

Temporal B-spline spec

## Usage

``` r
spec_time_bspline(
  k = NULL,
  degree = 3L,
  include_intercept = FALSE,
  orthonormalize = TRUE
)
```

## Arguments

- k:

  Optional number of components (df). If \`NULL\`, the encoder uses
  \`min(5, n_time)\` at encode time.

- degree:

  Spline degree (default 3).

- include_intercept:

  Logical include intercept.

- orthonormalize:

  Logical orthonormalize columns (default TRUE).

## Value

A \`spec_time_bspline\` object.
