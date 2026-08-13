# Encode data using a shared basis asset and subject field operator

Encode data using a shared basis asset and subject field operator

## Usage

``` r
encode_operator(
  x,
  template,
  field_operator = NULL,
  observation_operator = NULL,
  mask = NULL,
  domain = NULL,
  support = NULL,
  lambda = 0,
  center = TRUE,
  run_info = NULL,
  spatial_lambda = lambda,
  spatial_penalty = NULL,
  temporal_lambda = 0,
  temporal_order = 1L,
  sparse_lambda = 0,
  sparse_mode = c("none", "group_l2", "lasso"),
  max_iter = 200L,
  tol = 1e-06,
  label = "",
  ...
)
```

## Arguments

- x:

  Numeric matrix (time x target samples) or a `NeuroVec`.

- template:

  Shared basis asset providing [`basis_decoder()`](basis_decoder.md).

- field_operator:

  Subject field operator mapping template field samples to observed
  native samples.

- observation_operator:

  Legacy alias for `field_operator`.

- mask:

  Optional volumetric target mask for the field-operator target domain.

- domain:

  Optional non-volumetric target domain, for example a surface geometry.

- support:

  Optional target support aligned to `domain`, for example vertex
  indices on a surface. At least one of `mask` or `support` must be
  available either explicitly or via `field_operator$provenance`.

- lambda:

  Ridge penalty strength.

- center:

  Logical; if `TRUE`, center target samples before solving.

- run_info:

  Optional run metadata carried on the resulting latent object.

- spatial_lambda:

  Strength of the spatial coefficient penalty.

- spatial_penalty:

  Optional coefficient-space roughness matrix or diagonal weights.

- temporal_lambda:

  Strength of temporal smoothing.

- temporal_order:

  Difference order used for temporal smoothing.

- sparse_lambda:

  Strength of optional sparse coefficient shrinkage.

- sparse_mode:

  Sparse penalty mode. Use `"group_l2"` for atom-wise group shrinkage.

- max_iter:

  Maximum iterations for sparse AWPT optimization.

- tol:

  Relative convergence tolerance for sparse AWPT optimization.

- label:

  Optional label stored in metadata.

- ...:

  Reserved for future extensions.

## Value

A `TransportLatent` object.

## Details

The external field-operator contract is intentionally narrow.
\`fmrilatent\` consumes an operator-like object with:

- `n_source`, `n_target`:

  Source and target dimensions.

- `source_domain_id`, `target_domain_id`:

  Stable domain identifiers.

- `forward(z)`:

  Applies the operator to template field samples.

- `adjoint_apply(y)`:

  Applies the adjoint map.

- `provenance$target_mask`:

  Optional target-domain mask when the caller does not pass `mask`
  explicitly.

On the main quadratic and sparse transport/AWPT paths, coefficient
recovery is matrix-free: \`fmrilatent\` uses the operator's forward and
adjoint methods rather than materializing the full subject decoder.
