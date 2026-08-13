# Reconstruct LatentNeuroVec as a matrix (time x voxels)

Reconstruct LatentNeuroVec as a matrix (time x voxels)

## Usage

``` r
# S4 method for class 'LatentNeuroSurfaceVector'
as.matrix(x, ...)

# S4 method for class 'BilatLatentNeuroSurfaceVector'
as.matrix(x, ...)

# S4 method for class 'BlockLatentNeuroVector'
as.matrix(x, ...)

# S4 method for class 'LatentNeuroVec'
as.matrix(x, ...)
```

## Arguments

- x:

  A `LatentNeuroVec` object.

- ...:

  Ignored.

## Value

A numeric matrix with rows = time points, columns = mask voxels.
