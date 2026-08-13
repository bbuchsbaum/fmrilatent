# Construct split-half reliability settings for BOLDZip-SR

Construct split-half reliability settings for BOLDZip-SR

## Usage

``` r
boldzip_reliability(
  split = c("odd_even", "halves"),
  min_texture_reliability = 0,
  min_temporal_reliability = 0
)
```

## Arguments

- split:

  Split strategy. \`"odd_even"\` pairs adjacent odd/even time points;
  \`"halves"\` pairs the first and second half of the run.

- min_texture_reliability:

  Minimum held-out reliability required to keep a fine-detail texture
  loading.

- min_temporal_reliability:

  Minimum carrier reliability required before temporal coefficients are
  kept. Carriers below this threshold are zeroed.

## Value

A \`BoldZipSRReliability\` settings object.
