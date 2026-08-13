# Build an AWPT basis template

Build an AWPT basis template

## Usage

``` r
awpt_basis_template(
  parcellation,
  basis_spec = basis_awpt_wavelet(),
  loadings = NULL,
  anatomical_operator = NULL,
  conductance = NULL,
  coefficient_roughness = NULL,
  center = FALSE,
  ridge = 1e-08,
  label = "awpt_wavelet",
  ...
)
```

## Arguments

- parcellation:

  A `ClusterReduction` or `ClusteredNeuroVol`.

- basis_spec:

  An AWPT basis specification created by
  [`basis_awpt_wavelet()`](basis_awpt_wavelet.md).

- loadings:

  Optional explicit template loadings matrix. When supplied, fmrilatent
  skips wave-packet lifting and uses these loadings directly as the
  decoder basis \\B\\.

- anatomical_operator:

  Optional field-space roughness operator on the template domain. When
  supplied, the coefficient roughness is computed as \\Q = B^T L B\\. In
  the current v1 implementation this affects the roughness penalty, not
  the basis construction itself.

- conductance:

  Optional symmetric field-space conductance matrix. When supplied,
  fmrilatent builds the corresponding graph Laplacian and then forms \\Q
  = B^T L B\\. As with `anatomical_operator`, this shapes the v1
  roughness model rather than directly adapting the lifted basis.

- coefficient_roughness:

  Optional coefficient-space roughness matrix. This bypasses field-space
  construction and is stored directly as \\Q\\.

- center:

  Logical; if `TRUE`, center data before projection.

- ridge:

  Small positive ridge added to the Gram diagonal if needed.

- label:

  Optional label stored in metadata.

- ...:

  Additional arguments passed to the underlying lift path.

## Value

An `AWPTBasisTemplate`.
