# Diffusion wavelet latent constructor (explicit basis)

Diffusion wavelet latent constructor (explicit basis)

## Usage

``` r
diffusion_wavelet_latent(
  X,
  mask,
  reduction = NULL,
  spec = basis_diffusion_wavelet(),
  k_neighbors = 6L,
  label = ""
)
```

## Arguments

- X:

  Matrix time x voxels (mask order).

- mask:

  LogicalNeuroVol or 3D logical array.

- reduction:

  ClusterReduction; if NULL, defaults to one cluster per voxel.

- spec:

  diffusion wavelet spec (basis_diffusion_wavelet()).

- k_neighbors:

  k for graph building.

- label:

  Optional label.
