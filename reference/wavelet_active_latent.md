# Active-pencil wavelet latent (CDF 5/3)

Active-pencil wavelet latent (CDF 5/3)

## Usage

``` r
wavelet_active_latent(
  X,
  mask,
  levels_space = 2L,
  levels_time = 0L,
  threshold = 0
)
```

## Arguments

- X:

  Numeric matrix (time x voxels in mask order) or 4D array.

- mask:

  LogicalNeuroVol or 3D logical array.

- levels_space:

  Integer spatial lifting levels (default 2).

- levels_time:

  Integer temporal lifting levels (default 0 = none).

- threshold:

  Optional hard threshold after space transform.

## Value

An \`ImplicitLatent\` with decoder supporting \`time_idx\` and
\`roi_mask\`.
