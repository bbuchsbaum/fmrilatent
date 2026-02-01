# Encode data using a hierarchical template

Encode data using a hierarchical template

## Usage

``` r
encode_hierarchical(X, template, label = NULL)
```

## Arguments

- X:

  matrix time x voxels (mask order) matching template mask

- template:

  HierarchicalBasisTemplate

- label:

  Optional label for the resulting LatentNeuroVec (defaults to template
  label)

## Value

LatentNeuroVec with basis = time x atoms coefficients, loadings =
template loadings
