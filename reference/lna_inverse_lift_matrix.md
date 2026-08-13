# Low-level compatibility wrapper for matrix inverse lifting

Mirrors neuroarchive matrix-batch inverse-lift API.

## Usage

``` r
lna_inverse_lift_matrix(
  root_coeff,
  detail_coeffs_by_level,
  mask_flat_morton,
  mask_dims,
  levels,
  scalings,
  compat_profile = NULL
)
```

## Arguments

- root_coeff:

  Numeric matrix (\`time x n_root\`).

- detail_coeffs_by_level:

  List of detail coefficient matrices by level.

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

Numeric matrix (\`time x masked_voxels\`) in Morton order.
