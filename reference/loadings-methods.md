# Get the loadings matrix (spatial components)

Extract the loadings matrix from a latent space representation. For
`LatentNeuroVec` objects, this returns the spatial loadings matrix with
dimensions (nVoxels x k) where k is the number of components and nVoxels
is the number of voxels within the mask.

## Usage

``` r
loadings(x, ...)

# S4 method for class 'LatentNeuroVec'
loadings(x)
```

## Arguments

- x:

  An object containing loadings (e.g., `LatentNeuroVec`)

- ...:

  Additional arguments (currently unused)

## Value

The loadings matrix (typically voxels x components)

## Examples

``` r
if (FALSE) { # \dontrun{
# For LatentNeuroVec:
l_matrix <- loadings(lvec)
print(dim(l_matrix)) # nVoxels x nComponents
} # }
```
