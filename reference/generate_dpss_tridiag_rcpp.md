# Generate DPSS via tridiagonal prolate matrix (O(n^2))

Uses the classic tridiagonal formulation (Percival & Walden 1993;
Thomson 1982) and LAPACK dstev for efficiency on long series.

## Usage

``` r
generate_dpss_tridiag_rcpp(n, NW, k)
```

## Arguments

- n:

  Integer length of the time series.

- NW:

  Time-bandwidth product (N \* W), where W is normalized half-bandwidth.

- k:

  Number of tapers to return.

## Value

Matrix (n x k) with columns ordered by decreasing eigenvalue.
