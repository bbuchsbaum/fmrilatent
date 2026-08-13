# Evaluate BOLDZip-SR reconstruction quality

Evaluate BOLDZip-SR reconstruction quality

## Usage

``` r
evaluate_boldzip_sr(X, object, reliability_weights = NULL)
```

## Arguments

- X:

  Original matrix with rows as voxels/grayordinates and columns as time
  points.

- object:

  A \`BoldZipSR\` object or a reconstructed matrix.

- reliability_weights:

  Optional matrix or vector of reliability weights.

## Value

Named numeric vector with reconstruction metrics.
