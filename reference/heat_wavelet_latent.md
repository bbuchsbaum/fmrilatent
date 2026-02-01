# Heat wavelet latent constructor (explicit basis)

Heat wavelet latent constructor (explicit basis)

## Usage

``` r
heat_wavelet_latent(
  X,
  mask,
  reduction = NULL,
  spec = basis_heat_wavelet(),
  k_neighbors = 6L,
  label = ""
)
```

## Arguments

- X:

  Matrix time x voxels (mask order)

- mask:

  LogicalNeuroVol or 3D logical array

- reduction:

  ClusterReduction; if NULL, defaults to one cluster per voxel

- spec:

  heat wavelet spec (basis_heat_wavelet())

- k_neighbors:

  k for local graph building

- label:

  Optional label
