# Spatiotemporal Slepian latent (implicit, separable)

Spatiotemporal Slepian latent (implicit, separable)

## Usage

``` r
slepian_spatiotemporal_latent(
  X,
  mask,
  tr,
  bandwidth = 0.1,
  k_time = NULL,
  reduction = NULL,
  k_space = 3L,
  k_neighbors = 6L,
  label = ""
)
```

## Arguments

- X:

  Numeric matrix time x voxels (mask order).

- mask:

  LogicalNeuroVol or 3D logical array.

- tr:

  Repetition time (seconds).

- bandwidth:

  Half-bandwidth in Hz for temporal Slepians (default 0.1).

- k_time:

  Number of temporal Slepians; if NULL uses floor(2\*NW)-1.

- reduction:

  ClusterReduction for spatial graph; if NULL, one cluster per voxel.

- k_space:

  Number of spatial Slepians per cluster (default 3).

- k_neighbors:

  k-NN for spatial graph (default 6).

- label:

  Optional label.

## Value

An \`ImplicitLatent\` with decoder using separable Slepians.
