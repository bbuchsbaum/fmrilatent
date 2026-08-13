# Lift diffusion wavelets for clustered reduction

Lift diffusion wavelets for clustered reduction

## Usage

``` r
# S4 method for class 'ClusterReduction,spec_diffusion_wavelet'
lift(reduction, basis_spec, data = NULL, k_neighbors = 6L, ...)
```

## Arguments

- reduction:

  ClusterReduction describing voxel-to-cluster map.

- basis_spec:

  Diffusion wavelet spec (basis_diffusion_wavelet()).

- data:

  Ignored for this graph-only spatial dictionary; accepted only for the
  shared \`lift()\` method signature.

- k_neighbors:

  k for graph building (cluster centroids).

- ...:

  Additional arguments (unused).

## Value

Matrix of loadings (voxels x components) concatenating scaling bases
across scales.
