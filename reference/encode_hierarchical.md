# Encode data using a hierarchical template

Note: unlike parcel and AWPT encoding, hierarchical encoding does not
center the data. The returned offset is always `numeric(0)`.

## Usage

``` r
encode_hierarchical(
  X,
  template,
  label = NULL,
  mask = NULL,
  materialize = c("handle", "auto", "matrix")
)
```

## Arguments

- X:

  matrix time x voxels (mask order) matching template mask

- template:

  HierarchicalBasisTemplate

- label:

  Optional label for the resulting LatentNeuroVec (defaults to template
  label)

- mask:

  Optional mask to validate against the template mask before encoding.
  When supplied, it must match the template exactly.

- materialize:

  "handle", "matrix", or "auto" for template loadings.

## Value

LatentNeuroVec with basis = time x atoms coefficients, loadings =
template loadings
