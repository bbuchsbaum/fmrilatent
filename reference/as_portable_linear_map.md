# Coerce an operator to the portable linear-map contract

Normalizes an incoming object to the canonical
[`portable_linear_map`](portable_linear_map.md) representation used
internally by [`encode_operator()`](encode_operator.md),
[`decoder()`](decoder.md), and related routines. Matrices and
[`Matrix::Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html) objects
are wrapped with closure-based `forward`/`adjoint_apply` callbacks; list
inputs are validated and re-normalized with top-level `source_support` /
`target_support` preferred over `provenance$*` fallbacks.

## Usage

``` r
as_portable_linear_map(
  x,
  source_domain_id = "",
  target_domain_id = "",
  source_support = NULL,
  target_support = NULL,
  provenance = list(),
  context = "portable linear map"
)
```

## Arguments

- x:

  An input operator. Either a base `matrix`,
  [`Matrix::Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html), or a
  list satisfying the [`portable_linear_map`](portable_linear_map.md)
  contract.

- source_domain_id, target_domain_id:

  Optional character tags applied only when `x` is a matrix or
  [`Matrix::Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html).
  Ignored for list inputs (which carry their own).

- source_support, target_support:

  Optional support descriptors applied only when `x` is a matrix or
  [`Matrix::Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html).
  Typical values are a `LogicalNeuroVol`, a 3D logical mask, integer
  vertex indices, or `NULL`. Ignored for list inputs.

- provenance:

  Optional provenance list applied only when `x` is a matrix or
  [`Matrix::Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html).
  Ignored for list inputs.

- context:

  Optional label used in error messages.

## Value

A list in the canonical portable-linear-map form.

## See also

[`portable_linear_map`](portable_linear_map.md),
[`validate_portable_linear_map`](validate_portable_linear_map.md).

## Examples

``` r
m <- matrix(1:6, nrow = 2)
op <- as_portable_linear_map(m, target_domain_id = "native")
op$forward(c(1, 1, 1))
#> [1]  9 12
```
