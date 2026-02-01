# Construct a shared LoadingsHandle via diffusion-wavelet lifting

Wraps \`lift(reduction, basis_spec, data)\` so multiple
\`LatentNeuroVec\` instances can share the same spatial dictionary
without embedding the full matrix in each object.

## Usage

``` r
diffusion_wavelet_loadings_handle(
  reduction,
  basis_spec,
  data = NULL,
  k_neighbors = 6L,
  id = NULL,
  label = "diffusion-wavelet"
)
```

## Arguments

- reduction:

  Graph/cluster reduction used by \`lift()\`.

- basis_spec:

  Basis specification, e.g., from \`basis_diffusion_wavelet()\`.

- data:

  Optional data passed through to \`lift()\` (often NULL).

- k_neighbors:

  k for graph building when lifting.

- id:

  Optional registry id; provide a stable string to reuse across
  sessions. If NULL, a random id is generated.

- label:

  Optional human-readable label.

## Value

A `LoadingsHandle`.
