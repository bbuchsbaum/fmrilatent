# Inverse Haar wavelet transform

Inverse Haar wavelet transform

## Usage

``` r
haar_wavelet_inverse(
  coeff,
  mask,
  levels = NULL,
  z_seed = NULL,
  roi_mask = NULL,
  time_idx = NULL,
  levels_keep = NULL
)
```

## Arguments

- coeff:

  Output from \`haar_wavelet_forward()\` (list with root/detail
  matrices).

- mask:

  3D logical array or LogicalNeuroVol.

- levels:

  Integer levels (defaults to \`meta\$levels\` if present).

- z_seed:

  Integer Morton seed (defaults to \`meta\$z_seed\`).

- roi_mask:

  Optional ROI mask (same dims as mask) to subset voxels before
  returning.

- time_idx:

  Optional integer vector of time indices to subset rows.

- levels_keep:

  Optional integer vector of detail levels to include (1 = finest,
  \`levels\` = coarsest). Levels not listed are zeroed, enabling
  coarse/partial reconstruction without full detail.

## Value

Matrix (time x voxels) reconstructed (or subsetted if ROI/time
specified).
