# Build a DCT-II temporal basis matrix

Build a DCT-II temporal basis matrix

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

  Normalization: \`"ortho"\` (default) gives orthonormal columns;
  \`"none"\` returns raw cosine columns with non-unit L2 norms.

## Value

Dense Matrix (\`n_time x k\`). Columns are orthonormal only when \`norm
= "ortho"\`; with \`norm = "none"\` they are unnormalized DCT-II cosine
columns.
