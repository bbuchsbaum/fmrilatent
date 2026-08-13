# Lift parcel/cluster-local PCA bases for ClusterReduction

Computes PCA eigenvectors within each cluster and assembles a global
block-sparse loadings matrix (voxels x components). This is typically
used with \`encode(..., spec_space_pca(...), reduction = ...)\`.

## Usage

``` r
# S4 method for class 'ClusterReduction,spec_pca'
lift(
  reduction,
  basis_spec,
  data = NULL,
  center = TRUE,
  scale = FALSE,
  offset = NULL,
  backend = c("auto", "svds", "svd"),
  ...
)
```

## Arguments

- reduction:

  ClusterReduction describing voxel-to-cluster map.

- basis_spec:

  PCA basis specification (from \`basis_pca()\`).

- data:

  Required numeric matrix (time x voxels in mask order). Unlike
  graph-only lift methods, PCA consumes \`data\` to estimate
  cluster-local components and aborts when it is \`NULL\`.

- center:

  Logical; center voxels before PCA (default TRUE).

- scale:

  Logical; scale voxels before PCA (default FALSE).

- offset:

  Optional numeric vector of voxel means (length n_vox). If provided and
  \`center = TRUE\`, this is used instead of recomputing
  \`colMeans(data)\`.

- backend:

  SVD backend: "auto" (default), "svds" (RSpectra), or "svd" (base).

- ...:

  Unused.

## Value

A sparse Matrix (voxels x components) with attribute
\`fmrilatent.singular_values\` giving per-component singular values.

## Details

\`lift(ClusterReduction, spec_pca)\` always returns unwhitened spatial
loadings. If \`basis_spec\$whiten\` is \`TRUE\`, the method emits a
classed warning because whitening requires the projected temporal scores
and is therefore the caller's responsibility. The standard \`encode(...,
spec_space_pca(whiten = TRUE))\` path handles this post-lift whitening
using the \`fmrilatent.singular_values\` attribute attached here.
