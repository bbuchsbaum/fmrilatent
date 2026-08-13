# Construct a shared temporal basis descriptor

Construct a shared temporal basis descriptor

## Usage

``` r
shared_temporal_spec(
  kind = c("dct", "bspline", "slepian", "custom"),
  n_time = NULL,
  rank = NULL,
  basis = NULL,
  params = list(),
  id = NULL,
  meta = list()
)
```

## Arguments

- kind:

  Temporal basis kind: \`"dct"\`, \`"bspline"\`, \`"slepian"\`, or
  \`"custom"\`.

- n_time:

  Number of time points. Inferred from \`basis\` for custom specs.

- rank:

  Number of basis columns. Inferred from \`basis\` when omitted.

- basis:

  Matrix-like basis for \`kind = "custom"\`.

- params:

  Basis-specific parameters.

- id:

  Optional id.

- meta:

  Optional advisory metadata.

## Value

A \`SharedTemporalSpec\` descriptor.
