# Generate test data for encoder development

Creates a small synthetic dataset suitable for testing `encode_spec`
methods. Useful for extension packages that implement custom encoders.

## Usage

``` r
fmrilatent_test_data(dims = c(3L, 3L, 2L), n_time = 8L)
```

## Arguments

- dims:

  Integer vector of spatial dimensions (default `c(3, 3, 2)`).

- n_time:

  Number of time points (default 8).

## Value

A list with elements:

- X:

  Numeric matrix (n_time x n_voxels) of random data.

- mask:

  Logical array of dimensions `dims`, all `TRUE`.

- dims:

  The spatial dimensions used.

- n_time:

  The number of time points used.

## Examples

``` r
td <- fmrilatent_test_data()
dim(td$X)        # 8 x 18
#> [1]  8 18
dim(td$mask)     # 3 x 3 x 2
#> [1] 3 3 2
```
