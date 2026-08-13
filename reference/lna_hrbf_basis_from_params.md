# Build an HRBF Basis with Optional Neuroarchive Compatibility Semantics

Compatibility entry point for external callers that need explicit
control over centre/sigma inputs and output column-space semantics.

## Usage

``` r
lna_hrbf_basis_from_params(
  params,
  mask,
  centres = NULL,
  sigmas = NULL,
  full_grid = TRUE,
  compat_profile = NULL,
  mask_world_coords = NULL,
  mask_arr = NULL,
  mask_linear_indices = NULL
)
```

## Arguments

- params:

  HRBF parameter list.

- mask:

  \`LogicalNeuroVol\` mask defining voxel locations.

- centres:

  Optional numeric matrix (\`K x 3\`) of world-space centres.

- sigmas:

  Optional numeric vector (\`K\`) of sigma values for \`centres\`.

- full_grid:

  Logical. If \`TRUE\`, basis columns index the full mask grid
  (\`length(as.array(mask))\`); if \`FALSE\`, columns index active mask
  voxels.

- compat_profile:

  Optional explicit compatibility profile identifier.

- mask_world_coords:

  Optional precomputed active-voxel world coords.

- mask_arr:

  Optional precomputed logical mask array.

- mask_linear_indices:

  Optional precomputed active-voxel linear indices.

## Value

Sparse matrix with one row per HRBF atom.
