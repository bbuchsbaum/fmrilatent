# Construct a shared LoadingsHandle via heat-wavelet lifting

Wraps a heat-wavelet \`lift()\` call so multiple \`LatentNeuroVec\`
instances can share the same spatial dictionary without embedding the
full matrix in each object.

## Usage

``` r
heat_wavelet_loadings_handle(
  reduction,
  basis_spec = basis_heat_wavelet(),
  data = NULL,
  k_neighbors = 6L,
  id = NULL,
  label = "heat-wavelet"
)
```

## Arguments

- reduction:

  Graph/cluster reduction used by \`lift()\`.

- basis_spec:

  Basis specification; defaults to \`basis_heat_wavelet()\`.

- data:

  Ignored by heat-wavelet lifting; accepted only to keep the
  lifted-handle constructor signature aligned across families.

- k_neighbors:

  Number of neighbors used for local graph construction when
  materializing the lifted basis.

- id:

  Optional registry id; provide a stable string to reuse across
  sessions. If NULL, a deterministic id is derived from the spec and
  reduction.

- label:

  Optional human-readable label.

## Value

A `LoadingsHandle`.
