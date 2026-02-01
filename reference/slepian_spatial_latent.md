# Slepian spatial latent constructor (explicit basis)

Slepian spatial latent constructor (explicit basis)

## Usage

``` r
slepian_spatial_latent(
  X,
  mask,
  reduction = NULL,
  spec = basis_slepian(),
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

  Slepian basis spec (basis_slepian()).

- k_neighbors:

  k for local graph building.

- label:

  Optional label.
