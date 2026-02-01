# Lift reduced bases back to voxel space (abstract generic)

Lift reduced bases back to voxel space (abstract generic)

## Usage

``` r
lift(reduction, basis_spec, data = NULL, ...)
```

## Arguments

- reduction:

  A \`GraphReduction\` subclass describing topology.

- basis_spec:

  A basis specification (e.g., \`basis_slepian()\`).

- data:

  Optional data for data-driven bases (e.g., PCA).

- ...:

  Additional arguments passed to methods (e.g., k_neighbors).

## Value

Typically a voxel x components \`Matrix\` (often sparse) or an
implementation-defined object for implicit decoders.
