# Encode a matrix with experimental BOLDZip-SR compression

\`boldzip_sr_encode()\` is an experimental matrix-level implementation
of Split-Reliable Graph Carrier Compression. It expects a matrix with
rows as voxels or grayordinates and columns as time points. The
implementation stores temporally compressed carriers, sparse
high-resolution texture loadings, and sparse split-reliable residual
events. It is intended as a research prototype, not a production binary
codec.

## Usage

``` r
boldzip_sr_encode(
  X,
  spatial_basis = NULL,
  k_carriers = 96L,
  temporal_k = NULL,
  temporal_spec = NULL,
  q_texture = 2L,
  texture_lags = 0L,
  reliability = boldzip_reliability(),
  quantization = boldzip_quantization(),
  events = boldzip_events(),
  center = TRUE,
  label = NULL
)
```

## Arguments

- X:

  Numeric matrix with dimensions \`voxels x time\`.

- spatial_basis:

  Optional \`BoldZipSRSpatialBasis\` object or list with \`phi_c\` and
  \`phi_d\` matrices. If omitted, the detail basis is identity and
  carriers are learned from \`X\` directly.

- k_carriers:

  Number of carrier time courses to learn.

- temporal_k:

  Number of DCT coefficients per carrier. Defaults to \`ceiling(n_time /
  4)\` when \`temporal_spec\` is \`NULL\`.

- temporal_spec:

  Optional temporal basis descriptor. May be a \`SharedTemporalSpec\`,
  \`spec_time_dct\`, \`spec_time_bspline\`, or numeric basis matrix with
  rows equal to time points. When supplied, it determines the temporal
  basis and \`temporal_k\`.

- q_texture:

  Maximum number of carrier loadings per detail atom.

- texture_lags:

  Integer vector of allowed carrier lags for texture loading fits. A
  positive lag uses \`Z_k(t - lag)\`.

- reliability:

  Reliability settings from \[boldzip_reliability()\].

- quantization:

  Quantization settings from \[boldzip_quantization()\].

- events:

  Event settings from \[boldzip_events()\].

- center:

  Whether to store and remove a voxel-wise mean before fitting.

- label:

  Optional label stored in metadata.

## Value

A \`BoldZipSR\` object.
