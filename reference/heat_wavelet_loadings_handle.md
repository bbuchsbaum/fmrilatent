# Construct a shared LoadingsHandle via heat-wavelet lifting

Wraps \`lift(reduction, basis_spec, data)\` so multiple
\`LatentNeuroVec\` instances can share the same spatial dictionary
without embedding the full matrix in each object.

## Usage

``` r
heat_wavelet_loadings_handle(
  reduction,
  basis_spec,
  data = NULL,
  id = NULL,
  label = "heat-wavelet"
)
```

## Arguments

- reduction:

  Graph/cluster reduction used by \`lift()\`.

- basis_spec:

  Basis specification, e.g., from \`basis_heat_wavelet()\`.

- data:

  Optional data passed through to \`lift()\` (often NULL).

- id:

  Optional registry id; provide a stable string to reuse across
  sessions. If NULL, a random id is generated.

- label:

  Optional human-readable label.

## Value

A `LoadingsHandle`.
