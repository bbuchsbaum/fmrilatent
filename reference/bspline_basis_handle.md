# Create a BasisHandle for a B-spline temporal basis

Create a BasisHandle for a B-spline temporal basis

## Usage

``` r
bspline_basis_handle(
  n_time,
  k,
  degree = 3L,
  knots = NULL,
  boundary_knots = NULL,
  include_intercept = FALSE,
  id = NULL,
  label = NULL
)
```

## Arguments

- n_time:

  Integer, number of time points.

- k:

  Integer, number of spline basis functions (df).

- degree:

  Degree of the spline (default 3).

- knots:

  Optional interior knots (scaled 0-1).

- boundary_knots:

  Optional boundary knots (scaled 0-1).

- include_intercept:

  Logical; include intercept column (default FALSE).

- id:

  Optional registry key (generated if NULL).

- label:

  Optional human-readable label.

## Value

A `BasisHandle`.
