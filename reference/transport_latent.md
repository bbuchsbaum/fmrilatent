# Construct a transport-backed implicit latent object

Construct a transport-backed implicit latent object

## Usage

``` r
transport_latent(
  coeff_raw,
  basis_asset,
  field_operator = NULL,
  observation_operator = NULL,
  mask = NULL,
  domain = NULL,
  support = NULL,
  coeff_analysis = NULL,
  analysis_transform = NULL,
  offset = numeric(0),
  run_info = NULL,
  meta = list()
)
```

## Arguments

- coeff_raw:

  Raw coefficient matrix (time x k).

- basis_asset:

  Shared basis asset.

- field_operator:

  Subject field operator satisfying the contract documented in
  [`encode_operator()`](encode_operator.md).

- observation_operator:

  Legacy alias for `field_operator`.

- mask:

  Target-domain mask for volumetric reconstruction.

- domain:

  Optional target domain for non-volumetric reconstruction.

- support:

  Optional target support for non-volumetric reconstruction.

- coeff_analysis:

  Optional analysis-coordinate coefficient matrix (time x k).

- analysis_transform:

  Optional transform descriptor from raw to analysis coordinates.

- offset:

  Optional voxel offset vector in target space.

- run_info:

  Optional run metadata.

- meta:

  Optional metadata list.

## Value

A `TransportLatent` object inheriting from `ImplicitLatent`.

## Details

\`transport_latent()\` stores the \`fmrilatent\` side of the
operator-backed workflow: shared basis asset, subject field operator,
raw and analysis-space coefficients, and decoder views for native or
template projection. It does not fit GLMs or carry statistical
summaries. Downstream model-fitting code should ordinarily consume
`coef_time(x, coordinates = "analysis")`, then return coefficient-space
effects and covariance matrices for projection with
[`decode_coefficients()`](decode_coefficients.md) and
[`decode_covariance()`](decode_covariance.md).
