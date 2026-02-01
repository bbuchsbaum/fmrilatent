# Lift heat wavelets for clustered reduction

Lift heat wavelets for clustered reduction

## Usage

``` r
# S4 method for class 'ClusterReduction,spec_heat_wavelet'
lift(reduction, basis_spec, data = NULL, k_neighbors = 6L, ...)
```

## Arguments

- reduction:

  ClusterReduction describing voxel-to-cluster map.

- basis_spec:

  Heat wavelet spec (basis_heat_wavelet()).

- data:

  Optional (unused) for API symmetry.

- k_neighbors:

  k for local graph building.

- ...:

  Additional arguments (unused).

## Value

Sparse Matrix (voxels x components).
