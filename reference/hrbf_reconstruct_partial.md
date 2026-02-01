# Partially reconstruct selected voxels/timepoints

Computes only the requested voxels (and optional timepoints) without
materializing the full dense volume. Useful for fast indexed access.

## Usage

``` r
hrbf_reconstruct_partial(coeff, mask, params, voxel_idx, time_idx = NULL)
```

## Arguments

- coeff:

  Coefficient matrix (time x atoms).

- mask:

  \`LogicalNeuroVol\`.

- params:

  HRBF parameter list.

- voxel_idx:

  Integer vector of voxel linear indices (within the mask grid).

- time_idx:

  Optional integer vector of time indices; defaults to all rows.

## Value

Matrix of shape length(time_idx) x length(voxel_idx).
