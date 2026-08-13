# Predict from a BOLDZip-SR codec payload

Predict from a BOLDZip-SR codec payload

## Usage

``` r
# S3 method for class 'BoldZipSR'
predict(object, time_idx = NULL, roi = NULL, ...)
```

## Arguments

- object:

  A \`BoldZipSR\` object.

- time_idx:

  Optional integer time indices to return.

- roi:

  Optional integer or logical row subset to return.

- ...:

  Additional arguments passed to \[boldzip_sr_decode()\].

## Value

Reconstructed matrix with rows as voxels/grayordinates and columns as
time points.
