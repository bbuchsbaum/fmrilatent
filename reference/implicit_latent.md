# Construct an ImplicitLatent object

Construct an ImplicitLatent object

## Usage

``` r
implicit_latent(coeff, decoder, meta, mask)
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

  Logical 3D array (or LogicalNeuroVol) describing voxel support.

## Value

An object of class \`ImplicitLatent\`.
