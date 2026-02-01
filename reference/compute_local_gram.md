# Compute local Gram matrices for neighborhoods

Given a loadings matrix L (voxels x k) and a list of voxel index sets
(neighborhoods), returns the per-neighborhood Gram matrices \`M_i =
t(L_Vi)

## Usage

``` r
compute_local_gram(loadings, neighborhoods, simplify = FALSE)
```

## Arguments

- loadings:

  Matrix or LoadingsHandle (voxels x k) from a LatentNeuroVec.

- neighborhoods:

  List of integer vectors of voxel indices (mask order).

- simplify:

  Logical; if TRUE and all neighborhoods same size, returns an array k x
  k x n_neighborhood; otherwise a list of matrices.

## Value

List (or array) of Gram matrices.
