# Represent subject loadings as group shared loadings plus a delta

Represent subject loadings as group shared loadings plus a delta

## Usage

``` r
group_delta_loadings(group, delta = NULL, scale = 1, id = NULL, meta = list())
```

## Arguments

- group:

  Matrix-like shared/group loadings or a \`SharedReference\`.

- delta:

  Optional subject-specific loading delta. If \`NULL\`, a zero delta is
  used.

- scale:

  Numeric multiplier applied to \`delta\`.

- id:

  Optional id.

- meta:

  Optional advisory metadata.

## Value

A \`GroupDeltaLoadings\` object.
