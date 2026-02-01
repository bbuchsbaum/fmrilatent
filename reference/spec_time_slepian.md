# Temporal Slepian/DPSS spec

Temporal Slepian/DPSS spec

## Usage

``` r
spec_time_slepian(
  tr,
  bandwidth = 0.1,
  k = NULL,
  backend = c("tridiag", "dense")
)
```

## Arguments

- tr:

  Repetition time (seconds).

- bandwidth:

  Half-bandwidth in Hz (default 0.1).

- k:

  Optional number of components (default floor(2\*NW)-1).

- backend:

  Backend to use ("tridiag" or "dense").

## Value

A \`spec_time_slepian\` object for \`encode()\` / \`spec_st()\`.
