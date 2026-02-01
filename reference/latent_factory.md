# Simple factory to build a spec and encode in one call

Simple factory to build a spec and encode in one call

## Usage

``` r
latent_factory(
  family,
  x,
  mask,
  reduction = NULL,
  ...,
  materialize = "handle",
  label = ""
)
```

## Arguments

- family:

  One of: "dct_time", "slepian_time", "slepian_space", "heat_space",
  "slepian_st".

- x:

  Data matrix (time x voxels).

- mask:

  Mask (required for spatial families).

- reduction:

  Optional GraphReduction for spatial specs.

- ...:

  Passed to spec constructors and encode().

- materialize:

  "handle", "matrix", or "auto" (default "handle").

- label:

  Optional label for the resulting object.

## Value

A LatentNeuroVec or ImplicitLatent object.
