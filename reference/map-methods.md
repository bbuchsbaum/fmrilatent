# Get the map object

Extract the index map from a latent space representation. For
`LatentNeuroVec` objects, this returns the `IndexLookupVol` that maps 3D
coordinates to linear mask indices.

## Usage

``` r
# S4 method for class 'LatentNeuroVec'
map(x)
```

## Arguments

- x:

  An object containing a map (e.g., `LatentNeuroVec`)

## Value

The index lookup object
