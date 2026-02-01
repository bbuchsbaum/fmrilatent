# Encode data into a latent representation using a spec

Encode data into a latent representation using a spec

## Usage

``` r
encode(
  x,
  spec,
  mask,
  reduction = NULL,
  materialize = c("handle", "auto", "matrix"),
  label = "",
  ...
)
```

## Arguments

- x:

  Matrix (time x voxels in mask order).

- spec:

  Spec object created by \`spec_time\_\*\`, \`spec_space\_\*\`, or
  \`spec_st\`.

- mask:

  LogicalNeuroVol or logical array (required for spatial pieces).

- reduction:

  Optional GraphReduction (for spatial specs).

- materialize:

  "handle", "matrix", or "auto" (default "handle").

- label:

  Optional label.

- ...:

  Additional arguments passed to methods.

## Value

A \`LatentNeuroVec\` (explicit bases) or \`ImplicitLatent\` (separable
cases).
