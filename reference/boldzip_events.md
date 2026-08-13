# Construct sparse innovation event settings for BOLDZip-SR

Construct sparse innovation event settings for BOLDZip-SR

## Usage

``` r
boldzip_events(max_events = 256L, threshold_sd = 3, paired_fraction = 0.25)
```

## Arguments

- max_events:

  Maximum number of residual events to store.

- threshold_sd:

  Residual event threshold in robust standard deviations.

- paired_fraction:

  Minimum paired split amplitude agreement.

## Value

A \`BoldZipSREvents\` settings object.
