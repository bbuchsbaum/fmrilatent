# Spatiotemporal spec (separable)

Spatiotemporal spec (separable)

## Usage

``` r
spec_st(time, space, core_mode = c("auto", "explicit"))
```

## Arguments

- time:

  Temporal spec (\`spec_time\_\*\`).

- space:

  Spatial spec (\`spec_space\_\*\`).

- core_mode:

  Controls the returned representation. \`"auto"\` (the default) always
  returns an \`ImplicitLatent\` (decoder-only separable form).
  \`"explicit"\` returns an explicit \[LatentNeuroVec\] (the separable
  \`basis materializable, which holds for every currently supported
  \`spec_st\` combination; if the factors are not both materializable it
  falls back to \`ImplicitLatent\` with a
  \`fmrilatent_warning_explicit_core_unavailable\` warning. The
  reconstruction is identical either way.

## Value

A \`spec_st\` object.
