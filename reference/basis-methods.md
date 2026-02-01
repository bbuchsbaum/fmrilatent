# Get the basis matrix (temporal components)

Extract the basis matrix from a latent space representation. For
`LatentNeuroVec` objects, this returns the temporal basis matrix with
dimensions (nTime x k) where k is the number of components.

## Usage

``` r
basis(x, ...)

# S4 method for class 'LatentNeuroVec'
basis(x)
```

## Arguments

- x:

  An object containing a basis matrix (e.g., `LatentNeuroVec`)

- ...:

  Additional arguments (currently unused)

## Value

The basis matrix (typically time x components)

## Examples

``` r
if (FALSE) { # \dontrun{
# For LatentNeuroVec:
b_matrix <- basis(lvec)
print(dim(b_matrix)) # nTime x nComponents
} # }
```
