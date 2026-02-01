# Build hierarchical parcellation levels from Schaefer surface atlas

Constructs nested parcellation levels for the hierarchical template
by: 1. Loading Schaefer surface atlas from neuroatlas 2. Computing
geodesic distances and boundary contacts via neurosurf 3. Building a
parcel similarity matrix (geodesic + boundary + network) 4. Performing
spectral Ward clustering to create coarser levels 5. Mapping surface
parcels back to volumetric voxel labels

Real-data glue to finish: - Provide fsaverage Schaefer surfaces (via
\`neuroatlas::schaefer_surf\`) so we can compute geodesic/boundary
terms. - Provide the volumetric Schaefer atlas aligned to the MNI mask
(via \`neuroatlas::get_schaefer_atlas\`). - Ensure mask and volumetric
atlas share resolution/origin/spacing; this function assumes 2 mm MNI
unless you pass a custom \`vol_atlas\`. - Subcortex is \*\*not\*\*
handled here; add it separately when assembling the full template.

## Usage

``` r
build_schaefer_levels(
  mask,
  parcels = 400,
  networks = 17,
  space = "fsaverage6",
  k_levels = c(16, 32, 64),
  vol_atlas = NULL,
  alpha = 0.5,
  beta = 0.3,
  gamma = 0.2,
  d0 = 30,
  component_policy = c("largest", "error", "each", "merge"),
  cache = TRUE
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

## Value

A list with components:

- levels:

  List of integer vectors (length = \#voxels in mask), one per level,
  nested.

- network:

  Character vector of Yeo network labels per finest-level parcel.

- hemi:

  Character vector of hemisphere ("L"/"R") per finest-level parcel.

- geo_dist_lh, geo_dist_rh:

  Parcel geodesic distance matrices (for diagnostics).

- boundary_lh, boundary_rh:

  Boundary contact matrices (LH/RH).

## Details

Requires \`neuroatlas\` and \`neurosurf\` packages. The surface atlas
provides geodesic distances and boundary topology; the volumetric atlas
maps parcels to mask voxels.
