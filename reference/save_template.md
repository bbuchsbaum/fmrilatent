# Save a template object to disk

Save a template object to disk

## Usage

``` r
save_template(template, file, compress = "xz", ...)

# S4 method for class 'SurfaceAWPTBasisTemplate'
save_template(template, file, compress = "xz", ...)

# S4 method for class 'AWPTBasisTemplate'
save_template(template, file, compress = "xz", ...)

# S4 method for class 'HierarchicalBasisTemplate'
save_template(template, file, compress = "xz", ...)

# S4 method for class 'ParcelBasisTemplate'
save_template(template, file, compress = "xz", ...)

# S4 method for class 'SurfaceBasisTemplate'
save_template(template, file, compress = "xz", ...)
```

## Arguments

- template:

  A template object.

- file:

  Output path.

- compress:

  Compression passed to
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html).

- ...:

  Additional arguments passed to methods.

## Value

Normalized output path, invisibly.
