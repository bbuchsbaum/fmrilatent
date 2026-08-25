# Identify latent coordinates and decoder domains

Returns a structured, deterministic identity record for the coefficient
coordinates and decoded output domain of a latent object. Coefficient
identity excludes observation values and observation count so objects
made from one fitted basis can be reused across independent data
collections.

## Usage

``` r
latent_space_id(
  x,
  coordinates = c("analysis", "raw"),
  space = c("native", "template"),
  ...
)

# S4 method for class 'LatentNeuroVec'
latent_space_id(
  x,
  coordinates = c("analysis", "raw"),
  space = c("native", "template"),
  ...
)

# S4 method for class 'LatentNeuroSurfaceVector'
latent_space_id(
  x,
  coordinates = c("analysis", "raw"),
  space = c("native", "template"),
  ...
)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
latent_space_id(
  x,
  coordinates = c("analysis", "raw"),
  space = c("native", "template"),
  ...
)

# S4 method for class 'BlockLatentNeuroVector'
latent_space_id(
  x,
  coordinates = c("analysis", "raw"),
  space = c("native", "template"),
  ...
)

# S4 method for class 'ImplicitLatent'
latent_space_id(
  x,
  coordinates = c("analysis", "raw"),
  space = c("native", "template"),
  ...
)
```

## Arguments

- x:

  A latent object.

- coordinates:

  Coefficient coordinate system to identify.

- space:

  Decoder output space to identify.

- ...:

  Additional arguments passed to methods.

## Value

A structured \`fmrilatent_latent_space_id\` record containing at least
\`coordinate_id\`, \`decoder_domain_id\`, and \`decoder_id\`.

## Details

The three identity layers have distinct meanings:

- \`coordinate_id\` identifies the ordered fitted latent axes and the
  requested raw-to-analysis transform.

- \`decoder_domain_id\` identifies the selected decoded domain, ordered
  support, and decoder contract, independently of coefficient values.

- \`decoder_id\` identifies the complete mapping from the selected
  coefficient coordinates to that decoded domain.

Identity is derived from stored receipts where possible, so
handle-backed basis or loading assets are not materialized. A plain
\`ImplicitLatent\` does not retain enough structured decoder provenance
and is rejected; transport-backed implicit objects satisfy the contract.
