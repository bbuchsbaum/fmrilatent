# Plot a spatial atom (loading vector) on a mask

Plot a spatial atom (loading vector) on a mask

## Usage

``` r
plot_spatial_atom(loadings, mask, idx = 1L, main = NULL)
```

## Arguments

- loadings:

  Matrix (voxels x k) or LoadingsHandle.

- mask:

  LogicalNeuroVol or logical array defining voxel order.

- idx:

  Component index (1-based).

- main:

  Optional title.

## Value

Invisibly, the 3D array plotted.
