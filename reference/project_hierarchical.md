# Project coefficients only (no LatentNeuroVec wrapper)

Project coefficients only (no LatentNeuroVec wrapper)

## Usage

``` r
project_hierarchical(template, X)
```

## Arguments

- template:

  HierarchicalBasisTemplate to use for projection

- X:

  matrix time x voxels (mask order) matching template mask

## Value

Matrix of coefficients (time x atoms)
