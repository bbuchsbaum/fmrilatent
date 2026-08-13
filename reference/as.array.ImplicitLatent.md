# Reconstruct an ImplicitLatent as an array

Reconstruct an ImplicitLatent as an array

## Usage

``` r
# S3 method for class 'ImplicitLatent'
as.array(x, time_idx = NULL, roi_mask = NULL, ...)
```

## Arguments

- x:

  An `ImplicitLatent` object.

- time_idx:

  Optional integer time indices to keep.

- roi_mask:

  Optional logical ROI mask; voxels outside the ROI are zero.

- ...:

  Additional arguments passed to the decoder.

## Value

A numeric array with dimensions `c(x, y, z, time)`.
