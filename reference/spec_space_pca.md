# Spatial PCA spec (cluster-local)

Computes PCA eigenvectors within each cluster/parcel specified by a
\`ClusterReduction\` and returns a block-sparse spatial dictionary.

## Usage

``` r
spec_space_pca(
  k = NULL,
  center = TRUE,
  whiten = FALSE,
  backend = c("auto", "svds", "svd"),
  scale = FALSE
)
```

## Arguments

- k:

  Optional components per cluster. If \`NULL\`, encoders use their
  family default (\`3L\`).

- center:

  Logical; center voxels before PCA (default TRUE). When TRUE, voxel
  means are stored in \`LatentNeuroVec@offset\`.

- whiten:

  Logical; if TRUE, return whitened scores (unit-variance) and rescaled
  loadings such that reconstruction is unchanged.

- backend:

  SVD backend: "auto" (default), "svds" (RSpectra), or "svd" (base).

- scale:

  Logical; scale voxels before PCA (default FALSE). Passed through to
  \`lift(ClusterReduction, spec_pca)\`.

## Value

A \`spec_space_pca\` object.
