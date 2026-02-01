# Benchmark encode/decode round-trips

Benchmark encode/decode round-trips

## Usage

``` r
benchmark_roundtrip(
  mask_dims = c(16, 16, 8),
  n_time = 10L,
  methods = c("slepian_space", "hrbf", "wavelet_active", "bspline_hrbf_st"),
  iterations = 1L
)
```

## Arguments

- mask_dims:

  Integer vector length 3 for spatial dims.

- n_time:

  Integer time points.

- methods:

  Character vector of families: "slepian_space", "hrbf",
  "wavelet_active", "bspline_hrbf_st".

- iterations:

  Integer iterations per method (default 1 to stay fast).

## Value

A data.frame/tibble with timings and reconstruction error.
