# Get the offset vector

Extract the offset vector from a latent space representation. For
`LatentNeuroVec` objects, this returns the voxel-wise offset (mean or
intercept) that is added after the basis x loadings reconstruction.

## Usage

``` r
# S4 method for class 'LatentNeuroVec'
offset(object)
```

## Arguments

- object:

  An object containing an offset (e.g., `LatentNeuroVec`)

## Value

The offset vector (length = number of voxels in mask)

## Examples

``` r
if (FALSE) { # \dontrun{
# For LatentNeuroVec:
off_vector <- offset(lvec)
print(length(off_vector)) # Should equal number of voxels in mask
} # }
```
