# Build a shared parcel basis template

Computes a spatial dictionary within each parcel using
[`lift()`](lift.md) and caches the Gram factorization for efficient
projection. The resulting template can be reused across subjects via
[`spec_space_parcel()`](spec_space_parcel.md).

## Usage

``` r
parcel_basis_template(
  parcellation,
  basis_spec = basis_slepian(k = 5),
  data = NULL,
  center = TRUE,
  ridge = 1e-08,
  ...
)
```

## Arguments

- parcellation:

  A [`ClusterReduction`](ClusterReduction-class.md) or a
  `ClusteredNeuroVol` (coerced via
  [`as_cluster_reduction()`](as_cluster_reduction.md)).

- basis_spec:

  A basis specification for [`lift()`](lift.md). Default is
  `basis_slepian(k = 5)` which computes the k smallest Laplacian
  eigenvectors of the voxel adjacency graph within each parcel. Use
  `basis_pca(k)` with `data` for data-driven bases.

- data:

  Optional numeric matrix (time x voxels, mask order). Required for
  data-driven specs like [`basis_pca()`](basis_pca.md).

- center:

  Logical; if `TRUE` (default), center voxel time series both when
  building PCA templates and when projecting new subjects. The encoded
  object stores subject-specific voxel means in `LatentNeuroVec@offset`
  so reconstruction preserves voxel means. For geometric bases
  (Laplacian, Slepian) this controls whether subjects are centered prior
  to projection.

- ridge:

  Small positive scalar added to the Gram diagonal if Cholesky fails
  (default 1e-8).

- ...:

  Additional arguments passed to [`lift()`](lift.md) (e.g.,
  `k_neighbors` for graph-based specs). For PCA templates, projection
  must remain replayable at encode time, so preprocessing arguments such
  as `offset` and `scale` are not supported here.

## Value

A `"ParcelBasisTemplate"` object (S3) with components:

- loadings:

  Sparse `Matrix` (voxels x atoms), block-diagonal by parcel.

- gram_factor:

  Cached Cholesky factorization of \\L^T L\\.

- reduction:

  `ClusterReduction` used.

- basis_spec:

  The spec that produced the loadings.

- center:

  Whether centering is applied when building/projecting.

- meta:

  List with `family`, `k`, `ridge`, `label_map`, and `cluster_map`.

## Details

The default `basis_slepian(k)` computes the k smallest eigenvectors of
the graph Laplacian built from voxel coordinates within each parcel.
These smooth spatial functions form a data-independent dictionary
suitable for projecting any subject's data.

For data-driven bases, pass `basis_pca(k)` with `data =` a training
matrix (e.g., group-average data). The resulting PCA loadings are then
fixed and reused for each subject.

## See also

[`spec_space_parcel`](spec_space_parcel.md),
[`as_cluster_reduction`](as_cluster_reduction.md), [`lift`](lift.md),
[`basis_slepian`](basis_slepian.md), [`basis_pca`](basis_pca.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Geometric (Laplacian) basis -- no training data needed
atlas <- load_atlas("schaefer_400")
tmpl <- parcel_basis_template(atlas, basis_slepian(k = 8))
lvec <- encode(bold, spec_space_parcel(tmpl))

# Data-driven shared PCA basis
tmpl_pca <- parcel_basis_template(atlas, basis_pca(k = 5), data = group_bold)
lvec <- encode(subj_bold, spec_space_parcel(tmpl_pca))
} # }
```
