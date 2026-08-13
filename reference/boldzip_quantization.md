# Construct reliability-aware quantization settings for BOLDZip-SR

Construct reliability-aware quantization settings for BOLDZip-SR

## Usage

``` r
boldzip_quantization(base_step = 0, epsilon = 1e-06)
```

## Arguments

- base_step:

  Base quantization step. Set to 0 to disable quantization.

- epsilon:

  Small positive value used in reliability-shaped step sizes.

## Value

A \`BoldZipSRQuantization\` settings object.
