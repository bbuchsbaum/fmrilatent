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
  materialize = "auto",
  label = ""
)
```

## Arguments

- family:

  Character scalar naming one of the standard \`encode()\` families. See
  \*\*Accepted family names\*\* for the canonical names and supported
  aliases.

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

The class follows the same per-family contract as \[encode()\]: explicit
spatial families and explicit temporal families return a
\[LatentNeuroVec\] (a concrete \`ExplicitLatent\`); the spatiotemporal
families (\`st_slepian\`, \`st_bspline_hrbf\`) build a \`spec_st\` and
therefore always return an \`ImplicitLatent\`. See the \`@return\`
section of \[encode()\] for the full taxonomy.

## Accepted family names

Canonical names are listed first; aliases in parentheses are accepted
for compatibility.

- Temporal:

  \`time_dct\` (\`dct_time\`), \`time_slepian\` (\`slepian_time\`).

- Spatial:

  \`space_slepian\` (\`slepian_space\`), \`space_pca\` (\`pca_space\`),
  \`space_parcel\` (\`parcel_space\`), \`space_heat\` (\`heat_space\`),
  \`space_hrbf\` (\`hrbf_space\`), \`space_wavelet_active\`
  (\`wavelet_active\`), and \`hierarchical\`.

- Spatiotemporal:

  \`st\` (requires explicit \`time\` and \`space\` specs),
  \`st_slepian\` (\`slepian_st\`), and \`st_bspline_hrbf\`
  (\`bspline_hrbf_st\`).

AWPT is intentionally not a \`latent_factory()\` family because it
requires a shared \`basis_asset\` and a subject \`field_operator\`; use
\[encode_awpt()\] or \[encode_operator()\] for AWPT subject fitting.
