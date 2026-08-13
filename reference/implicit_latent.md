# Construct an ImplicitLatent object

Construct an ImplicitLatent object

## Usage

``` r
implicit_latent(
  coeff,
  decoder,
  meta,
  mask = NULL,
  domain = NULL,
  support = NULL
)
```

## Arguments

- coeff:

  Arbitrary coefficient payload (list or matrix) needed by decoder.

- decoder:

  Function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL)
  returning matrix.

- meta:

  List metadata; must include \`family\` string.

- mask:

  Logical 3D array (or LogicalNeuroVol) describing volumetric support.

- domain:

  Optional decoded output domain. For non-volumetric latent objects,
  supply a domain such as a
  [`neurosurf::SurfaceGeometry`](https://bbuchsbaum.github.io/neurosurf/reference/SurfaceGeometry.html)
  together with `support`.

- support:

  Optional decoded output support. For volumetric latent objects this is
  usually derived from `mask`; for surface-like latent objects it should
  be a vector of vertex indices (or a logical mask over vertices).

## Value

An object of class \`ImplicitLatent\`.
