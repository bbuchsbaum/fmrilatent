# Create a coarsened graph reduction

Create a coarsened graph reduction

## Usage

``` r
make_coarsened_reduction(mask, P_matrix, coarse_adj = NULL, info = list())
```

## Arguments

- mask:

  \`LogicalNeuroVol\` or array-like mask defining the fine domain.

- P_matrix:

  Fine-by-coarse sparse prolongation matrix.

- coarse_adj:

  Optional coarse-by-coarse sparse adjacency matrix.

- info:

  Optional metadata list.

## Value

A valid \`CoarsenedReduction\`.
