# Enable or disable the fmrilatent handle registry

The handle registry caches materialized matrices for `BasisHandle` and
`LoadingsHandle` objects. This can improve performance and reduce memory
duplication when multiple `LatentNeuroVec` objects share the same handle
IDs.

## Usage

``` r
fmrilatent_registry_enable()

fmrilatent_registry_disable()

fmrilatent_registry_enabled()
```

## Value

Invisibly, `TRUE`.

## Details

Set `fmrilatent_registry_disable()` to turn off caching (useful for
deterministic benchmarking or to avoid retaining large matrices).

## Capacity

The cache is bounded. Each registry (basis and loadings) holds at most
`getOption("fmrilatent.registry.max_entries", 256)` materialized
matrices; once full, the least-recently-used entry is evicted on the
next registration. Both registration and retrieval count as a use. Set
the option to `Inf` (or a non-positive value) to restore an unbounded
cache, or call
[`fmrilatent_registry_clear()`](fmrilatent_registry_clear.md) to drop
everything immediately.

## Examples

``` r
fmrilatent_registry_disable()
fmrilatent_registry_enable()
```
