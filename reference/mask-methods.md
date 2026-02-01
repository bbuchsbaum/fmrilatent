# Get the mask

Extract the brain mask from a latent space representation. For
`LatentNeuroVec` objects, this returns the `LogicalNeuroVol` defining
which voxels are included in the representation.

## Usage

``` r
# S4 method for class 'LatentNeuroVec'
mask(x)
```

## Arguments

- x:

  An object containing a mask (e.g., `LatentNeuroVec`)

## Value

The logical mask volume
