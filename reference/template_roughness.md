# Extract the spatial roughness operator for a template asset

Extract the spatial roughness operator for a template asset

## Usage

``` r
template_roughness(template, coordinates = c("analysis", "raw"), ...)

# S4 method for class 'SurfaceAWPTBasisTemplate'
template_roughness(template, coordinates = c("analysis", "raw"), ...)

# S4 method for class 'AWPTBasisTemplate'
template_roughness(template, coordinates = c("analysis", "raw"), ...)

# S4 method for class 'HierarchicalBasisTemplate'
template_roughness(template, coordinates = c("analysis", "raw"), ...)

# S4 method for class 'ParcelBasisTemplate'
template_roughness(template, coordinates = c("analysis", "raw"), ...)

# S4 method for class 'SurfaceBasisTemplate'
template_roughness(template, coordinates = c("analysis", "raw"), ...)
```

## Arguments

- template:

  A shared basis asset.

- coordinates:

  Coordinate system for the returned roughness operator.

- ...:

  Additional arguments passed to methods.

## Value

Matrix-like roughness operator or `NULL`.
