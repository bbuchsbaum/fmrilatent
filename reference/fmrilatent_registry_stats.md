# Get registry statistics

Returns count and approximate memory usage of cached matrices.

## Usage

``` r
fmrilatent_registry_stats(type = c("all", "basis", "loadings"))
```

## Arguments

- type:

  Which registry: "basis", "loadings", or "all" (default).

## Value

A list with count and bytes for each registry type.

## Examples

``` r
fmrilatent_registry_stats()
#> $basis
#> $basis$count
#> [1] 0
#> 
#> $basis$bytes
#> [1] 0
#> 
#> 
#> $loadings
#> $loadings$count
#> [1] 0
#> 
#> $loadings$bytes
#> [1] 0
#> 
#> 
#> $total
#> $total$count
#> [1] 0
#> 
#> $total$bytes
#> [1] 0
#> 
#> 
```
