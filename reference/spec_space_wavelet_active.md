# Spatial wavelet (active pencil) spec

Spatial wavelet (active pencil) spec

## Usage

``` r
spec_space_wavelet_active(
  levels_space = 2L,
  levels_time = 0L,
  select_threshold = NULL,
  threshold = NULL
)
```

## Arguments

- levels_space:

  Spatial lifting levels.

- levels_time:

  Optional time lifting levels.

- select_threshold:

  Non-negative active-set coefficient threshold after the spatial
  transform. This is stored as \`threshold\` for compatibility.

- threshold:

  Deprecated alias for \`select_threshold\`.

## Value

A \`spec_space_wavelet_active\` object.
