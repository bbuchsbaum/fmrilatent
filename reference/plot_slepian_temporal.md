# Plot temporal Slepians (DPSS)

Plot temporal Slepians (DPSS)

## Usage

``` r
plot_slepian_temporal(basis, max_components = 6L)
```

## Arguments

- basis:

  Matrix (time x k) or BasisHandle (kind = "slepian_temporal").

- max_components:

  Maximum components to display (default 6).

## Value

A ggplot object (if ggplot2 available); otherwise invisibly plots with
base graphics.
