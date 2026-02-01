# Build an orthonormal DCT-II basis matrix

Build an orthonormal DCT-II basis matrix

## Usage

``` r
build_dct_basis(n_time, k = n_time, norm = c("ortho", "none"))
```

## Arguments

- n_time:

  Number of time points.

- k:

  Number of components (columns); must satisfy k \<= n_time.

- norm:

  Normalization: "ortho" (default) or "none".

## Value

Dense Matrix (n_time x k).
