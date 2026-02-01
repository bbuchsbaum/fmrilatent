# Build Haar latent representation

Build Haar latent representation

## Usage

``` r
haar_latent(
  X,
  mask,
  levels = NULL,
  z_seed = 42L,
  threshold = list(type = "none", value = 0)
)
```

## Arguments

- X:

  Numeric matrix (time x voxels within mask).

- mask:

  LogicalNeuroVol or 3D logical array.

- levels:

  Decomposition levels (optional).

- z_seed:

  Integer Morton seed (default 42).

- threshold:

  Threshold list passed to \`haar_wavelet_forward\`.

## Value

An object of class \`HaarLatent\` containing coefficients and metadata.
