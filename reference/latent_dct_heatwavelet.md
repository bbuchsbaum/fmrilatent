# Create a LatentNeuroVec with heat-wavelet spatial dictionary

Creates a template LatentNeuroVec with heat-wavelet spatial loadings.
The basis matrix is initialized to zeros and should be populated with
actual coefficients (e.g., via encoding fMRI data).

## Usage

``` r
latent_dct_heatwavelet(
  n_time,
  k_time,
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

  Ignored (kept for backwards compatibility). The number of components
  is determined by the heat-wavelet loadings.

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
