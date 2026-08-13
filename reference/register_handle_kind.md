# Register a lazy handle materializer

Register a lazy handle materializer

## Usage

``` r
register_handle_kind(kind, materializer, type = c("basis", "loadings"))
```

## Arguments

- kind:

  Character scalar handle kind.

- materializer:

  Function that accepts a handle and returns a matrix-like object.

- type:

  Registry to update: `"basis"` or `"loadings"`.

## Value

Invisibly, the registered kind.
