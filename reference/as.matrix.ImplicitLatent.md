# Reconstruct an ImplicitLatent as a matrix

Reconstruct an ImplicitLatent as a matrix

## Usage

``` r
# S3 method for class 'ImplicitLatent'
as.matrix(x, time_idx = NULL, roi_mask = NULL, ...)
```

## Arguments

- x:

  An `ImplicitLatent` object.

- time_idx:

  Optional integer time indices to keep.

- roi_mask:

  Optional logical ROI mask for spatial subsetting.

- ...:

  Additional arguments passed to the decoder.

## Value

A numeric matrix with rows = time and columns = voxels within the
requested mask support.
