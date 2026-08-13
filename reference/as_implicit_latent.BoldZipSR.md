# Coerce a BOLDZip-SR payload to an ImplicitLatent

Coerce a BOLDZip-SR payload to an ImplicitLatent

## Usage

``` r
# S3 method for class 'BoldZipSR'
as_implicit_latent(x, mask = NULL, domain = NULL, support = NULL, ...)
```

## Arguments

- x:

  A \`BoldZipSR\` object.

- mask:

  Optional logical 3D array or \`LogicalNeuroVol\` for volumetric
  wrapping.

- domain:

  Optional decoded output domain.

- support:

  Optional decoded support. Defaults to row indices.

- ...:

  Additional arguments ignored.

## Value

An \`ImplicitLatent\` whose decoder returns matrices as \`time x
voxels\`, matching the rest of the implicit latent API.
