# HierarchicalBasisTemplate Class

A template-only, data-agnostic container for hierarchical Laplacian
frames. Stores the sparse spatial dictionary (primal basis), a cached
solver for the Gram matrix, and the parcellation hierarchy metadata.
Intended to be built offline for a fixed template (e.g., MNI) and reused
for encoding fMRI data.

## Slots

- `mask`:

  LogicalNeuroVol defining the domain (3D).

- `space`:

  NeuroSpace (typically 4D with a singleton time dim) matching the mask.

- `levels`:

  List of integer vectors (one per level) giving parcel ids per voxel
  (mask order).

- `parents`:

  List mapping child parcel ids to parent ids for each level \> 1.

- `loadings`:

  Sparse Matrix (voxels x atoms) containing concatenated atoms B.

- `gram_factor`:

  Cached factorization of \\G = t(B) \\\*\\ B\\ (e.g., dCHMsimpl from
  Matrix::Cholesky).

- `atoms`:

  data.frame describing each atom (col_id, level, parcel_id, parent_id,
  mode, label).

- `meta`:

  List for auxiliary metadata (atlas names, k_per_level, ridge,
  version).
