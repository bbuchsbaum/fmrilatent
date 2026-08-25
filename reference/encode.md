# Encode data into a latent representation using a spec

Encode data into a latent representation using a spec

## Usage

``` r
encode(
  x,
  spec,
  mask,
  reduction = NULL,
  materialize = c("auto", "handle", "matrix"),
  label = "",
  units = NULL,
  ...
)
```

## Arguments

- x:

  Matrix (time x voxels in mask order).

- spec:

  Standard encode spec object created by \`spec_time\_\*\`,
  \`spec_space\_\*\`, or \`spec_st\`. AWPT specs created by
  \[basis_awpt_wavelet()\] describe shared templates and are
  intentionally not accepted by \`encode()\`; use \[encode_awpt()\] or
  \[encode_operator()\] for transport-backed AWPT fits.

- mask:

  LogicalNeuroVol or logical array (required for spatial pieces).

- reduction:

  Optional GraphReduction (for spatial specs).

- materialize:

  "handle", "matrix", or "auto" (default "handle").

- label:

  Optional label.

- units:

  Optional declared \[latent_units_record()\] captured on the encoded
  object. When omitted, \[latent_units()\] truthfully reports \`status =
  "undeclared"\`.

- ...:

  Additional arguments passed to methods.

## Value

The return class depends on the spec family:

- Explicit spatial families:

  \`spec_space_slepian\`, \`spec_space_heat\`, \`spec_space_hrbf\`,
  \`spec_space_pca\`, and \`spec_space_wavelet_active\` return a
  \[LatentNeuroVec\], which is a concrete \`ExplicitLatent\` (the
  virtual marker defined at \`R/all_class.R\`). It stores an explicit
  \`basis x loadings + offset\` factorization.

- Explicit temporal families:

  \`spec_time_slepian\`, \`spec_time_dct\`, and \`spec_time_bspline\`
  likewise return a \[LatentNeuroVec\] (\`ExplicitLatent\`).

- Spatiotemporal (\`spec_st\`):

  \*\*Always\*\* returns an \`ImplicitLatent\` (a decoder-only /
  separable representation), even when both the time and space factors
  are fully explicit dense bases. \`ImplicitLatent\` is \*not\* an
  \`ExplicitLatent\`.

- Transport (\`spec_transport\` / AWPT encoders):

  return a \`TransportLatent\`, which is also \*not\* an
  \`ExplicitLatent\`.

In short: \`ExplicitLatent\` is the virtual S4 marker inherited by
\[LatentNeuroVec\]; \`ImplicitLatent\` and \`TransportLatent\` are S3
classes that deliberately do not inherit it.

## Dispatch model

For standard in-mask matrix encoders, \`encode()\` routes to the S3
generic \[encode_spec()\], which dispatches on the spec class and builds
the latent object directly. Transport-backed AWPT is the explicit
exception: an AWPT basis spec is used to build a shared template, while
the subject fit also requires a \`basis_asset\` and \`field_operator\`.
Those assets are not part of the standard \`encode_spec()\` signature,
so AWPT enters through \[encode_awpt()\] or \[encode_operator()\]
instead of \[encode()\].

## Offset and centering contract

A \[LatentNeuroVec\] reconstructs its data as \`basis (length = number
of in-mask voxels) added back after the low-rank term. The offset is
owned by the encoder's \*lift\* step, which is the single place that
decides whether and how to center:

- Families that populate \`offset\`:

  Only \`spec_space_pca\` produces a non-empty offset, and only when
  \`center = TRUE\` (the default): the per-voxel column means are
  removed before the PCA and stored in the \`offset\` slot so
  reconstruction restores them. The mean removal is performed exactly
  once, inside \`lift.basis_pca\` (see the \`fmrilatent.mean_scores\`
  attribute it returns); the encode caller does not re-center. With
  \`center = FALSE\`, PCA stores \`offset = numeric(0)\`.

- Families that never center:

  All other explicit families (\`spec_space_slepian\`,
  \`spec_space_heat\`, \`spec_space_hrbf\`, \`spec_time_slepian\`,
  \`spec_time_dct\`, \`spec_time_bspline\`) store \`offset =
  numeric(0)\`, i.e. no offset.

By convention \`offset = numeric(0)\` means "no offset" and
reconstruction treats it as a zero vector; a populated \`offset\` always
has one entry per in-mask voxel. The shared \`.encode_center()\` helper
is the single implementation of column-mean centering used by the
offset-producing paths.
