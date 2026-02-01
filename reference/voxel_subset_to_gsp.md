# Convert voxel subset to an rgsp graph

Convert voxel subset to an rgsp graph

## Usage

``` r
voxel_subset_to_gsp(mask, voxel_indices, k_neighbors = 6L)
```

## Arguments

- mask:

  LogicalNeuroVol or 3D logical array.

- voxel_indices:

  Integer vector of voxel indices (mask order).

- k_neighbors:

  Number of nearest neighbours for k-NN graph.

## Value

\`gsp_graph\` object (from rgsp).
