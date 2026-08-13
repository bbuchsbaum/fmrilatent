# Validate an object against the portable linear-map contract

Returns the normalized canonical form invisibly when `x` satisfies the
[`portable_linear_map`](portable_linear_map.md) contract. The canonical
return value lets callers use the result directly rather than paying
normalization a second time via
[`as_portable_linear_map()`](as_portable_linear_map.md). Set
`error = FALSE` to receive `FALSE` on failure instead of a stop.

## Usage

``` r
validate_portable_linear_map(x, context = "portable linear map", error = TRUE)
```

## Arguments

- x:

  Object to validate.

- context:

  Optional label used in error messages.

- error:

  If `TRUE` (the default), signal an error on failure. If `FALSE`,
  return `FALSE` silently on failure.

## Value

Invisibly: the normalized portable-linear-map list on success, or
`FALSE` when `error = FALSE` and validation fails. On failure with
`error = TRUE` (the default), an error is raised.

## See also

[`portable_linear_map`](portable_linear_map.md),
[`as_portable_linear_map`](as_portable_linear_map.md).

## Examples

``` r
op <- validate_portable_linear_map(matrix(1:4, 2, 2))
op$n_source
#> [1] 2
op$n_target
#> [1] 2
validate_portable_linear_map(list(forward = function(x) x), error = FALSE)
```
