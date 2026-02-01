# DPSS temporal basis (Slepian sequences)

Generates Discrete Prolate Spheroidal Sequences (DPSS) for a given
series length and bandwidth. Uses an internal RcppEigen solver over the
prolate matrix; suitable for moderate \`n_time\`. For very long series,
a tridiagonal solver can replace the backend with the same interface.

## Usage

``` r
dpss_time_basis(
  n_time,
  tr,
  bandwidth,
  k = NULL,
  backend = c("tridiag", "dense")
)
```

## Arguments

- n_time:

  Integer length of the time dimension.

- tr:

  Repetition time in seconds.

- bandwidth:

  Half-bandwidth in Hz (typical BOLD range 0.008–0.1).

- k:

  Optional number of tapers; defaults to `floor(2 * NW) - 1` where
  `NW = n_time * bandwidth * tr`. Clamped to `[1, n_time]`.

- backend:

  Either `"tridiag"` (default, O(n^2)) or `"dense"` (O(n^3), debugging).

## Value

Numeric matrix (n_time x k) with orthonormal columns.
