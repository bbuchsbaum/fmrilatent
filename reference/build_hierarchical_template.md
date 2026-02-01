# Build a hierarchical Laplacian template (offline)

Constructs a multi-level spatial basis from nested parcellations using
graph Laplacian eigenvectors. The resulting template can be reused to
efficiently encode multiple datasets that share the same mask geometry.

## Usage

``` r
build_hierarchical_template(
  mask,
  parcellations,
  k_per_level,
  k_neighbors = 6L,
  ridge = 1e-08,
  solver = c("chol", "qr"),
  label = "hierarchical_laplacian"
)
```

## Arguments

- mask:

  LogicalNeuroVol or logical array (3D) defining the domain.

- parcellations:

  List of integer vectors (one per level) of length = \#voxels in mask.
  Levels must be nested: each child parcel maps to exactly one parent
  parcel in the previous level.

- k_per_level:

  Integer vector giving \#modes per parcel at each level.

- k_neighbors:

  k for local graph construction inside parcels.

- ridge:

  Small diagonal ridge added to \\G = t(B) \\\*\\ B\\ for stability.

- solver:

  Solver choice: "chol" (default) or "qr" fallback.

- label:

  Optional label stored in meta.

## Value

HierarchicalBasisTemplate (primal basis B + cached solver and metadata).
