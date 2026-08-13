# Reconstruct a latent object as a 4D array

Reconstructs a latent object into a 4D array over its spatial support.

## Usage

``` r
reconstruct_array(x, time_idx = NULL, roi_mask = NULL, ...)

# S4 method for class 'ImplicitLatent'
reconstruct_array(x, time_idx = NULL, roi_mask = NULL, ...)

# S4 method for class 'LatentNeuroSurfaceVector'
reconstruct_array(x, time_idx = NULL, roi_mask = NULL, ...)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
reconstruct_array(x, time_idx = NULL, roi_mask = NULL, ...)

# S4 method for class 'BlockLatentNeuroVector'
reconstruct_array(x, time_idx = NULL, roi_mask = NULL, ...)

# S4 method for class 'LatentNeuroVec'
reconstruct_array(x, time_idx = NULL, roi_mask = NULL, ...)
```

## Arguments

- x:

  A latent object.

- time_idx:

  Optional integer time indices to keep.

- roi_mask:

  Optional logical ROI mask; voxels outside the ROI are zero.

- ...:

  Additional arguments passed to methods.

## Value

Numeric 4D array.
