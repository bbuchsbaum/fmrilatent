# Default lift method (placeholder)

This method exists to provide a clear error when no concrete lift is
registered. External packages should implement methods for specific
(reduction, basis_spec) signatures.

## Usage

``` r
# S4 method for class 'GraphReduction,ANY'
lift(reduction, basis_spec, data = NULL, ...)
```

## Arguments

- reduction:

  A \`GraphReduction\` subclass.

- basis_spec:

  A basis specification object.

- data:

  Optional data for data-driven bases.

- ...:

  Additional arguments (unused in default method).
