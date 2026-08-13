# Create a ClusterReduction from a mask and voxel-to-cluster map

Create a ClusterReduction from a mask and voxel-to-cluster map

## Usage

``` r
make_cluster_reduction(mask, map)
```

## Arguments

- mask:

  A `LogicalNeuroVol` or logical 3D array defining the brain mask.

- map:

  Integer vector (mask order) mapping each voxel to a cluster id.

## Value

A `ClusterReduction` object.
