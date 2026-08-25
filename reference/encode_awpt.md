# Encode data using the AWPT model

Encode data using the AWPT model

## Usage

``` r
encode_awpt(
  x,
  basis_asset,
  field_operator = NULL,
  observation_operator = NULL,
  mask = NULL,
  domain = NULL,
  support = NULL,
  spatial_lambda = 0,
  temporal_lambda = 0,
  temporal_order = 1L,
  sparse_lambda = 0,
  sparse_mode = c("none", "group_l2", "lasso"),
  max_iter = 200L,
  tol = 1e-06,
  center = TRUE,
  run_info = NULL,
  label = "",
  units = NULL,
  ...
)
```

## Arguments

- x:

  Numeric matrix (time x target samples) or a `NeuroVec`.

- basis_asset:

  An `AWPTBasisTemplate`.

- field_operator:

  Subject field operator. See [`encode_operator()`](encode_operator.md)
  for the required contract.

- observation_operator:

  Legacy alias for `field_operator`.

- mask:

  Optional volumetric target mask for the field-operator target domain.

- domain:

  Optional non-volumetric target domain.

- support:

  Optional target support aligned to `domain`.

- spatial_lambda:

  Strength of the anatomical roughness penalty.

- temporal_lambda:

  Strength of temporal smoothing.

- temporal_order:

  Temporal difference order for smoothing.

- sparse_lambda:

  Strength of optional sparse atom selection.

- sparse_mode:

  Sparse penalty mode. Use `"group_l2"` for atom-wise group shrinkage.

- max_iter:

  Maximum iterations for sparse AWPT optimization.

- tol:

  Relative convergence tolerance for sparse AWPT optimization.

- center:

  Logical; if `TRUE`, center target samples before solving.

- run_info:

  Optional run metadata; `run_lengths` control temporal blocks.

- label:

  Optional label stored in metadata.

- units:

  Optional declared \[latent_units_record()\] captured on the encoded
  object.

- ...:

  Reserved for future extensions.

## Value

A `TransportLatent` object with AWPT metadata.

## Details

AWPT does not expose a separate ridge penalty. The returned metadata
records `lambda = 0`; the anatomical roughness penalty is recorded
separately as `spatial_lambda`.
