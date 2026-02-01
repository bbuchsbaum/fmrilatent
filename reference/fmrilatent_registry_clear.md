# Clear the fmrilatent handle registry

Removes cached materialized matrices from the internal registry. Use
this to free memory or force re-materialization of handles.

## Usage

``` r
fmrilatent_registry_clear(type = c("all", "basis", "loadings"))
```

## Arguments

- type:

  Which registry to clear: "basis", "loadings", or "all" (default).

## Value

Invisibly, the number of entries removed.

## Examples

``` r
# Clear all cached matrices
fmrilatent_registry_clear()

# Clear only basis matrices
fmrilatent_registry_clear("basis")
```
