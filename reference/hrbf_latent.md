# Build a LatentNeuroVec using an HRBF basis

Build a LatentNeuroVec using an HRBF basis

## Usage

``` r
hrbf_latent(X, mask, params = list(), label = "")
```

## Arguments

- X:

  Numeric matrix (time x voxels within mask).

- mask:

  \`LogicalNeuroVol\`.

- params:

  HRBF parameter list (see \`hrbf_generate_basis\`).

- label:

  Optional character label.

## Value

A \`LatentNeuroVec\` with \`basis = coefficients\`, \`loadings =
t(HRBF_basis)\`.
