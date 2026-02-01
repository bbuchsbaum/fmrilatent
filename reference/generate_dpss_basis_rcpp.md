# Generate DPSS (Slepian) basis via dense prolate matrix eigen-decomposition

Generate DPSS (Slepian) basis via dense prolate matrix
eigen-decomposition

## Usage

``` r
generate_dpss_basis_rcpp(n, NW, k)
```

## Arguments

- n:

  Integer length of the time series.

- NW:

  Time-bandwidth product (N \* W), where W is normalized half-bandwidth
  (cycles per sample).

- k:

  Number of tapers to return (columns).

## Value

Matrix (n x k) with columns ordered by decreasing eigenvalue (energy
concentration).

## Note

This implementation builds the dense prolate matrix A\[i,j\] = 2W (i==j)
else sin(2\*pi\*W\*(i-j)) / (pi\*(i-j)) and computes the top-k
eigenvectors. It is O(n^3) and intended as a simple, dependency-light
starting point; consider replacing with a tridiagonal solver for very
long n.
