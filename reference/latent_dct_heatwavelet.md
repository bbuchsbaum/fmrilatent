# Create a template LatentNeuroVec with heat-wavelet spatial loadings

Builds a `LatentNeuroVec` whose loadings are a heat-wavelet
`LoadingsHandle` and whose basis is a placeholder zero matrix of
`n_time x k` (where `k` is determined by the heat-wavelet loadings, not
by the caller). The caller is expected to overwrite `lvec@basis` with
real coefficients (e.g. fitted from data) before the object represents a
valid factorization.

## Usage

``` r
latent_dct_heatwavelet(
  n_time,
  k_time = NULL,
  mask,
  cluster_map = NULL,
  reduction = NULL,
  hw_basis_spec = NULL,
  offset = numeric(0),
  label = "DCT + heat-wavelet"
)
```

## Arguments

- n_time:

  Number of time points.

- k_time:

  Optional ignored legacy argument. The number of components is
  determined by the heat-wavelet loadings. Supplying a non-NULL value
  warns with class \`fmrilatent_warning_deprecated\`.

- mask:

  LogicalNeuroVol or logical array mask (3D).

- cluster_map:

  Optional integer vector mapping voxels (mask order) to clusters.

- reduction:

  Graph reduction object; if NULL, built via
  \`make_cluster_reduction(mask, cluster_map)\` with default
  one-cluster-per-voxel map.

- hw_basis_spec:

  Heat-wavelet basis spec; defaults to \`basis_heat_wavelet()\`.

- offset:

  Optional voxel-wise offset (length n_vox).

- label:

  Optional label.

## Value

A `LatentNeuroVec` with placeholder basis matrix.

## Details

Despite the name, no DCT basis is constructed: the "dct" reference
predates the spec/encoder pipeline and is retained for API
compatibility. For an encoded DCT-temporal + heat-wavelet-spatial
pipeline see
`spec_st(time = spec_time_dct(...), space = spec_space_heat(...))`
passed to [`encode`](encode.md).
