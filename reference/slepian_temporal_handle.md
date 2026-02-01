# Create a BasisHandle for temporal Slepians (DPSS)

Create a BasisHandle for temporal Slepians (DPSS)

## Usage

``` r
slepian_temporal_handle(
  n_time,
  tr,
  bandwidth = 0.1,
  k = NULL,
  backend = c("tridiag", "dense"),
  id = NULL,
  label = NULL
)
```

## Arguments

- n_time:

  Integer number of time points.

- tr:

  Repetition time (seconds).

- bandwidth:

  Half-bandwidth in Hz (default 0.1).

- k:

  Optional number of tapers/components.

- backend:

  Backend passed to \`dpss_time_basis\` ("tridiag" or "dense").

- id:

  Optional registry key (generated if NULL).

- label:

  Optional human-readable label.

## Value

A `BasisHandle`.
