# Temporal DCT spec

Temporal DCT spec

## Usage

``` r
spec_time_dct(k = NULL, norm = c("ortho", "none"))
```

## Arguments

- k:

  Optional number of components. If \`NULL\`, the encoder uses all
  available DCT columns (\`k = n_time\`) at encode time.

- norm:

  Normalization ("ortho" or "none").

## Value

A \`spec_time_dct\` object.
