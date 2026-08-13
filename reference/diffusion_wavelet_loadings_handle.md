# Construct a shared LoadingsHandle via diffusion-wavelet lifting

Wraps a diffusion-wavelet \`lift()\` call so multiple \`LatentNeuroVec\`
instances can share the same spatial dictionary without embedding the
full matrix in each object.

## Usage

``` r
diffusion_wavelet_loadings_handle(
  reduction,
  basis_spec = basis_diffusion_wavelet(),
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

  Basis specification; defaults to \`basis_diffusion_wavelet()\`.

- data:

  Ignored by diffusion-wavelet lifting; accepted only to keep the
  lifted-handle constructor signature aligned across families.

- k_neighbors:

  k for graph building when lifting.

- id:

  Optional registry id; provide a stable string to reuse across
  sessions. If NULL, a deterministic id is derived from the spec and
  reduction.

- label:

  Optional human-readable label.

## Value

A `LoadingsHandle`.
