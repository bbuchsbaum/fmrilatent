# List entries in the fmrilatent handle registry

Returns the IDs of all cached matrices in the registry.

## Usage

``` r
fmrilatent_registry_list(type = c("all", "basis", "loadings"))
```

## Arguments

- type:

  Which registry to list: "basis", "loadings", or "all" (default).

## Value

Character vector of registered IDs.

## Examples

``` r
# See what's cached
fmrilatent_registry_list()
#> character(0)
```
