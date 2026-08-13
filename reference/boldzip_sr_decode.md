# Decode an experimental BOLDZip-SR object

Decode an experimental BOLDZip-SR object

## Usage

``` r
boldzip_sr_decode(object, time_idx = NULL, roi = NULL)
```

## Arguments

- object:

  A \`BoldZipSR\` object.

- time_idx:

  Optional integer time indices to return.

- roi:

  Optional integer or logical row subset to return.

## Value

Reconstructed matrix with rows as voxels/grayordinates and columns as
time points.
