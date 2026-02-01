# Build a B-spline temporal basis matrix

Build a B-spline temporal basis matrix

## Usage

``` r
build_bspline_basis(
  n_time,
  k,
  degree = 3L,
  knots = NULL,
  boundary_knots = NULL,
  include_intercept = FALSE,
  orthonormalize = TRUE
)
```

## Arguments

- n_time:

  Integer, number of time points.

- k:

  Integer, number of spline basis functions (df).

- degree:

  Degree of the spline (default 3 = cubic).

- knots:

  Optional numeric vector of interior knots in \[0, 1\].

- boundary_knots:

  Optional numeric vector length 2 giving boundary knots in \[0, 1\].

- include_intercept:

  Logical; passed to
  [`splines::bs()`](https://rdrr.io/r/splines/bs.html) (default FALSE).

## Value

Dense Matrix (n_time x k).
