# Low-level compatibility wrapper for matrix forward lifting

Mirrors neuroarchive matrix-batch forward-lift API.

## Usage

``` r
lna_forward_lift_matrix(
  data_morton,
  mask_flat_morton,
  mask_dims,
  levels,
  scalings,
  compat_profile = NULL
)
```

## Arguments

- data_morton:

  Numeric matrix (\`time x masked_voxels\`) in Morton order.

- mask_flat_morton:

  Logical full-grid mask in Morton order.

- mask_dims:

  Integer vector of length 3.

- levels:

  Integer decomposition levels.

- scalings:

  Precomputed scaling structure from \`precompute_haar_scalings\`.

- compat_profile:

  Optional compatibility profile override.

## Value

List with \`root_coeff\` and \`detail_coeffs_by_level\`.
