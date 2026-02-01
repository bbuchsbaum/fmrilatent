# Build hierarchical template from Schaefer surface atlas

Convenience wrapper that combines \`build_schaefer_levels()\` with
\`build_hierarchical_template()\` to produce a ready-to-use template
with geodesic-informed parcel clustering.

## Usage

``` r
build_schaefer_hierarchical_template(
  mask,
  parcels = 400,
  networks = 17,
  space = "fsaverage6",
  k_levels = c(16, 32, 64),
  k_per_level = c(8, 5, 3, 1),
  vol_atlas = NULL,
  alpha = 0.5,
  beta = 0.3,
  gamma = 0.2,
  d0 = 30,
  component_policy = "largest",
  cache = TRUE,
  k_neighbors = 6L,
  ridge = 1e-08,
  solver = c("chol", "qr")
)
```

## Arguments

- mask:

  LogicalNeuroVol defining the 3D brain mask (MNI space).

- parcels:

  Integer. Number of Schaefer parcels (100, 200, 300, 400, 500, 600,
  800, 1000).

- networks:

  Integer. Yeo network variant (7 or 17).

- space:

  Character. Surface space ("fsaverage6", "fsaverage", "fsaverage5").

- k_levels:

  Integer vector of cluster counts for coarser levels, ordered
  coarse→fine. The finest level is always the original Schaefer
  parcellation. Example: c(16, 32, 64) produces 4 levels: L0=16, L1=32,
  L2=64, L3=Schaefer-parcels.

- k_per_level:

  Integer vector of eigenmodes per parcel at each level. Length must
  match length(k_levels) + 1 (for the finest Schaefer level).

- vol_atlas:

  Optional. Pre-loaded volumetric Schaefer atlas (NeuroVol with parcel
  IDs). If NULL, attempts to load via neuroatlas::get_schaefer_atlas().

- alpha:

  Weight for boundary contact in similarity (default 0.5).

- beta:

  Weight for geodesic kernel in similarity (default 0.3).

- gamma:

  Weight for Yeo network agreement in similarity (default 0.2).

- d0:

  Scale for geodesic exponential kernel (default 30 mm).

- component_policy:

  How to handle fragmented parcels ("error", "largest", "each",
  "merge").

- cache:

  Logical. Use cached geodesic distances if available (default TRUE).

- k_neighbors:

  k for local graph construction inside parcels.

- ridge:

  Small diagonal ridge for Gram matrix stability.

- solver:

  Solver choice: "chol" or "qr".

## Value

HierarchicalBasisTemplate object.

## See also

[`build_hierarchical_template`](build_hierarchical_template.md),
[`build_schaefer_levels`](build_schaefer_levels.md)
