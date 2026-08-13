# Reconstruct a latent object as a matrix

Provides a common reconstruction interface for both explicit
`LatentNeuroVec` objects and implicit decoder-backed latent objects.

## Usage

``` r
reconstruct_matrix(x, time_idx = NULL, roi_mask = NULL, ...)

# S4 method for class 'ImplicitLatent'
reconstruct_matrix(x, time_idx = NULL, roi_mask = NULL, ...)

# S4 method for class 'LatentNeuroSurfaceVector'
reconstruct_matrix(x, time_idx = NULL, roi_mask = NULL, ...)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
reconstruct_matrix(x, time_idx = NULL, roi_mask = NULL, ...)

# S4 method for class 'BlockLatentNeuroVector'
reconstruct_matrix(x, time_idx = NULL, roi_mask = NULL, ...)

# S4 method for class 'LatentNeuroVec'
reconstruct_matrix(x, time_idx = NULL, roi_mask = NULL, ...)
```

## Arguments

- x:

  A latent object.

- time_idx:

  Optional integer time indices to keep.

- roi_mask:

  Optional logical ROI mask for spatial subsetting.

- ...:

  Additional arguments passed to methods.

## Value

Numeric matrix with rows = time and columns = voxels within the
requested mask support.
