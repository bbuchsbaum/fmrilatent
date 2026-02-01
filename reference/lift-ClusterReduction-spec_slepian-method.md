# Lift spatial Slepians for clustered reduction

Lift spatial Slepians for clustered reduction

## Usage

``` r
# S4 method for class 'ClusterReduction,spec_slepian'
lift(reduction, basis_spec, data = NULL, k_neighbors = 6L, ...)
```

## Arguments

- reduction:

  ClusterReduction describing voxel-to-cluster map.

- basis_spec:

  Slepian basis specification (from \`basis_slepian()\`).

- data:

  Optional (unused) for API symmetry.

- k_neighbors:

  k for local graph building.

- ...:

  Additional arguments (unused).

## Value

Sparse Matrix (voxels x components), block-concatenated over clusters.
