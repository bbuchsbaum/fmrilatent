# Linear access to LatentNeuroVec elements

Access elements of a `LatentNeuroVec` using linear (1D) indexing.
Reconstructs values on-the-fly from the basis and loadings
factorization.

Access elements of a `LatentNeuroVec` using linear (1D) indices into the
4D array representation.

## Usage

``` r
# S4 method for class 'LatentNeuroVec,numeric'
linear_access(x, i)

# S4 method for class 'LatentNeuroVec,integer'
linear_access(x, i)
```

## Arguments

- x:

  A `LatentNeuroVec` object

- i:

  Numeric index vector

## Value

Numeric vector of reconstructed values

The reconstructed values at the specified indices
