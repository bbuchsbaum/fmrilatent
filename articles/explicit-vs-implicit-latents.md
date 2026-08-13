# Explicit vs decoder-backed latents

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

fmrilatent ships **two latent object types** that share a common
interface but store data very differently. Knowing which one you have —
and which one a given encoder returns — is the single most useful piece
of mental orientation when reading the docs.

## The two tiers

    ┌─────────────────────────────────────────┐
    │  Explicit:  basis × loadings + offset   │   LatentNeuroVec  (S4)
    │             matrices, on disk           │
    ├─────────────────────────────────────────┤
    │  Decoder-backed:  coeff + decoder()     │   ImplicitLatent  (S3)
    │             closure that materializes   │
    │             on demand                   │
    └─────────────────────────────────────────┘

| Property | Explicit (`LatentNeuroVec`) | Decoder-backed (`ImplicitLatent`) |
|----|----|----|
| Class system | S4, inherits [`neuroim2::NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html) | S3, plain list |
| Storage | `@basis`, `@loadings`, `@offset` matrices (or lazy handles) | `$coeff`, `$decoder`, `$meta`, `$mask` |
| Reconstruction | `as.matrix(x)`, `series(x, …)` | `predict(x, time_idx, roi_mask)` |
| Latent factors | `basis(x)`, `loadings(x)` | `x$coeff` (heterogeneous) |
| Saved to disk | Matrix bytes | Closure (captures its environment) |
| Typical use | Compact storage of pre-computed factorization | External solver / non-separable codec |

## Which encoders return which?

| Encoder family | Returns |
|----|----|
| `spec_time_dct` / `spec_time_slepian` / `spec_time_bspline` | Explicit `LatentNeuroVec` |
| `spec_space_slepian` / `spec_space_pca` / `spec_space_heat` / `spec_space_hrbf` / `spec_space_wavelet_active` | Explicit `LatentNeuroVec` |
| `spec_space_parcel` (with `parcel_basis_template`) | Explicit `LatentNeuroVec` |
| `spec_st(time = …, space = …)` (separable spatiotemporal) | Decoder-backed `ImplicitLatent` |
| `spec_hierarchical_template` | Explicit `LatentNeuroVec` |
| `encode_transport(...)` | Decoder-backed `ImplicitLatent` |
| `encode_awpt(...)` | Decoder-backed `ImplicitLatent` |
| `encode_operator(...)` | Decoder-backed `ImplicitLatent` |
| `haar_latent(...)` | Decoder-backed `ImplicitLatent` (subclass `HaarLatent`) |

The rule of thumb: if the basis can be written down as a matrix with
fewer rows than the time axis (or fewer columns than the voxel count),
the encoder produces an explicit object. If the underlying contract
requires a non-trivial decoder — separable Kronecker structure, operator
transport, lifted wavelets, learned codecs — the encoder produces a
decoder-backed object.

## Working with explicit latents

``` r

mask  <- array(TRUE, dim = c(4, 4, 4))
mask_vol <- neuroim2::LogicalNeuroVol(mask, neuroim2::NeuroSpace(dim(mask)))
set.seed(7)
X     <- matrix(rnorm(20 * sum(mask)), nrow = 20)

lvec  <- encode(X, spec_time_dct(k = 6), mask = mask_vol, materialize = "matrix")
class(lvec)
#> [1] "LatentNeuroVec"
#> attr(,"package")
#> [1] "fmrilatent"
isS4(lvec)
#> [1] TRUE

# Direct factor access:
dim(basis(lvec))     # 20 x 6
#> [1] 20  6
dim(loadings(lvec))  # 64 x 6
#> [1] 64  6

# Reconstruction:
recon <- as.matrix(lvec)
dim(recon)
#> [1] 20 64

# Slicing — same as a NeuroVec:
ts1 <- series(lvec, 1L)
length(ts1)
#> [1] 20
```

`LatentNeuroVec` is a subclass of
[`neuroim2::NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html),
so the standard neuroim2 operations work —
[`dim()`](https://rdrr.io/r/base/dim.html),
[`series()`](../reference/series-methods.md),
[`as.array()`](https://rdrr.io/r/base/array.html), `[`, `[[`.
[`basis()`](../reference/basis-methods.md) and
[`loadings()`](../reference/loadings-methods.md) give you the latent
matrices directly.

## Working with decoder-backed latents

``` r

spec_separable <- spec_st(
  time  = spec_time_dct(k = 4),
  space = spec_space_hrbf(params = list(sigma0 = 2, levels = 0,
                                        radius_factor = 2.5))
)
ilat <- encode(X, spec_separable, mask = mask_vol)
class(ilat)
#> [1] "ImplicitLatent"
isS4(ilat)
#> [1] FALSE

# Coefficients + decoder, not basis × loadings:
names(ilat)
#> [1] "coeff"   "decoder" "meta"    "mask"    "domain"  "support"
str(ilat$coeff, max.level = 1)
#> List of 3
#>  $ core: num [1:4, 1:64] 1.98 -1.23 -2.04 2.02 0.44 ...
#>  $ B_t : num [1:20, 1:4] 0.224 0.224 0.224 0.224 0.224 ...
#>  $ L_s : num [1:64, 1:64] 1 0 0 0 0 0 0 0 0 0 ...
ilat$meta$family
#> [1] "st_separable"

# Reconstruction goes through predict():
recon_full <- predict(ilat)
dim(recon_full)         # n_time x n_voxels
#> [1] 20 64

# Partial decode — only the first 5 time points:
recon_part <- predict(ilat, time_idx = 1:5)
dim(recon_part)
#> [1]  5 64
```

[`predict()`](https://rdrr.io/r/stats/predict.html) is the universal
decoder API for the implicit tier. It accepts `time_idx`, `roi_mask`,
and family-specific arguments (`levels_keep` for haar, etc.), and only
materializes the slice you ask for.

## Serialization implications

This is the most common gotcha. Both tiers can be
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html)’d, but the cost and
reproducibility characteristics differ.

``` r

# Explicit: matrices serialize natively. With handle-backed slots
# (e.g. dct_basis_handle), the @id + @spec are saved and the basis is
# rematerialized on first access in the new session.
saveRDS(lvec, "lvec.rds")
lvec2 <- readRDS("lvec.rds")
identical(as.matrix(basis(lvec)), as.matrix(basis(lvec2)))  # TRUE
```

``` r

# Decoder-backed: $decoder is a closure. saveRDS captures its
# environment — including any data the closure references. This means:
#   - Self-contained decoders (haar, st-separable) round-trip cleanly.
#   - Decoders that reference large external assets (subject field
#     operators) save a copy of the asset by default.
saveRDS(ilat, "ilat.rds")
ilat2 <- readRDS("ilat.rds")
identical(predict(ilat), predict(ilat2))                     # TRUE
```

When in doubt: round-trip through
[`tempfile()`](https://rdrr.io/r/base/tempfile.html) and check that
[`predict()`](https://rdrr.io/r/stats/predict.html) (or
[`as.matrix()`](https://rdrr.io/r/base/matrix.html)) returns the same
numbers. The package test suite has dedicated coverage for this on the
explicit side (`test-latent_serialization.R`) and the implicit decoders
are exercised indirectly through their family-specific tests.

## Shared structure is orthogonal

Both tiers can participate in the **shared structure protocol**
(`R/shared_structure.R`), which lets multiple objects reference the same
heavy data — a template basis, a parcel atlas, a precomputed graph —
instead of each carrying its own copy. The protocol works through
dictionary handles (`BasisHandle` / `LoadingsHandle` for the explicit
tier, decoder-side asset references for the implicit tier) plus an
in-session shared-reference registry. Use shared structures when:

- You’re encoding many subjects against a common template.
- Multiple `LatentNeuroVec`s in the same session would otherwise
  duplicate the same dictionary in memory.
- You’re building a benchmark and want strict equality of the basis
  across runs.

See
[`vignette("shared-spatial-dictionaries")`](../articles/shared-spatial-dictionaries.md)
for the parcel-template walkthrough, and
[`?fmrilatent_registry_enable`](../reference/fmrilatent_registry_enable.md)
for the in-session cache controls.

## Choosing between tiers

| Choose explicit if … | Choose decoder-backed if … |
|----|----|
| You want fast, predictable matrix access | You need partial decoding or operator transport |
| You’ll be slicing voxels or time often | The basis is non-separable or learned |
| You want to inspect / plot the basis directly | The decoder captures domain knowledge (Haar lifting, AWPT) |
| You’re storing many objects to disk and want straightforward bytes | You’re in a coefficient-space modeling pipeline |

In practice, most users start with explicit (DCT or B-spline temporal
encoding) and only reach for the implicit tier when they hit `spec_st`,
the transport pipeline, or a wavelet codec.

## Further reading

- [`?LatentNeuroVec`](../reference/LatentNeuroVec.md), `?ImplicitLatent`
  for the class contracts.
- [`vignette("transport-aware-encoding")`](../articles/transport-aware-encoding.md)
  — the implicit tier in depth, including the shared-asset + subject
  field-operator pipeline.
- [`vignette("shared-spatial-dictionaries")`](../articles/shared-spatial-dictionaries.md)
  — the shared-structure protocol applied to atlas-based encoders.
- [`vignette("compression-diagnostics")`](../articles/compression-diagnostics.md)
  — comparing tiers on the same data for compression vs. fidelity
  tradeoffs.
