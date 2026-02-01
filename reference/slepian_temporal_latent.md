# LatentNeuroVec using a temporal DPSS basis

LatentNeuroVec using a temporal DPSS basis

## Usage

``` r
slepian_temporal_latent(
  X,
  mask,
  tr,
  bandwidth = 0.1,
  k = NULL,
  denoise = TRUE,
  backend = c("tridiag", "dense"),
  label = ""
)
```

## Arguments

- X:

  Numeric matrix (time x voxels within mask).

- mask:

  \`LogicalNeuroVol\`.

- tr:

  Repetition time in seconds.

- bandwidth:

  Half-bandwidth in Hz (default 0.1).

- k:

  Optional number of tapers; see \`dpss_time_basis\`.

- denoise:

  If TRUE, truncate to \`floor(2 \* NW) - 1\` (Shannon number).

- backend:

  DPSS computation backend passed to \`dpss_time_basis\`.

- label:

  Optional character label.

## Value

\`LatentNeuroVec\` with DPSS temporal basis and voxel loadings.
