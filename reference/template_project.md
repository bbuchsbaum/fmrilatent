# Project data onto a template

Project data onto a template

## Usage

``` r
template_project(x, data, ...)

# S4 method for class 'SurfaceAWPTBasisTemplate'
template_project(x, data, ...)

# S4 method for class 'AWPTBasisTemplate'
template_project(x, data, ...)

# S4 method for class 'HierarchicalBasisTemplate'
template_project(x, data, ...)

# S4 method for class 'ParcelBasisTemplate'
template_project(x, data, ...)

# S4 method for class 'SurfaceBasisTemplate'
template_project(x, data, ...)
```

## Arguments

- x:

  A template object.

- data:

  Numeric matrix (time x voxels in mask order).

- ...:

  Additional arguments passed to methods.

## Value

A list with `coefficients` (time x atoms/components) and `offset` (voxel
means or `numeric(0)`).
