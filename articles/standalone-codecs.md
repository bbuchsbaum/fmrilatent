# Standalone Codecs

``` r

library(fmrilatent)
#> Warning in rgl.init(initValue, onlyNULL): RGL: unable to open X11 display
#> Warning: 'rgl.init' failed, will use the null device.
#> See '?rgl.useNULL' for ways to avoid this warning.
#> 
#> Attaching package: 'fmrilatent'
#> The following object is masked from 'package:stats':
#> 
#>     offset
```

The primary package workflow is `spec_*()` plus
[`encode()`](../reference/encode.md), which returns latent neuroimaging
representations such as `LatentNeuroVec` and transport-backed latents.

Standalone codecs are different. They work on ordinary matrices, own
their own encode/decode payloads, and are not dispatched through
[`encode()`](../reference/encode.md). BOLDZip-SR is currently the only
standalone codec in `fmrilatent`. For a full encode/decode/evaluate
workflow, see [`vignette("boldzip")`](../articles/boldzip.md).

## What sets a standalone codec apart?

BOLDZip-SR is currently the only standalone codec in `fmrilatent`, and
it makes the category concrete. Three properties distinguish it from a
`spec_*()` encoder:

- **Matrix in, matrix out.** It expects voxels (or grayordinates) in
  rows and time in columns — intentionally the transpose of the
  [`encode()`](../reference/encode.md) convention, which uses time by
  voxels.
- **Owns its payload.** The encoder returns a `BoldZipSR` object, not a
  `LatentNeuroVec`, and it carries its own decode logic.
- **Not in the encode-spec registry.** You call it directly; it is never
  dispatched through [`encode()`](../reference/encode.md).

A minimal fit shows all three at once:

``` r

set.seed(9)
n_voxels <- 10L
n_time <- 24L
carrier <- sin(seq(0, 2 * pi, length.out = n_time))
texture <- seq(-1, 1, length.out = n_voxels)
X <- outer(texture, carrier)               # voxels x time
X <- X + matrix(rnorm(n_voxels * n_time, sd = 0.01), nrow = n_voxels)

fit <- boldzip_sr_encode(X, k_carriers = 1L, q_texture = 1L,
                         temporal_k = 6L)

inherits(fit, "BoldZipSR")    # a codec object, not a LatentNeuroVec
#> [1] TRUE
inherits(fit, "LatentNeuroVec")
#> [1] FALSE
```

That is the whole distinction. For the full encode → evaluate → decode →
budget-tuning workflow (including the
[`as_implicit_latent()`](../reference/as_implicit_latent.md) bridge back
into the package’s decoder contract), see
[`vignette("boldzip")`](../articles/boldzip.md). For the lower-level
shared-structure contracts BOLDZip-SR can share with other
representations, see
[`vignette("shared-structure-boldzip")`](../articles/shared-structure-boldzip.md).
