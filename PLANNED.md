# Planned extensions

Non-blocking work items identified during the transport-bridge audit
(see commits 0a4665b, 935b27a, b902a28). None of these block the
pipeline seam — fmrilatent already plays its role in the
ANTs/FNIRT/AFNI → neurotransform → neurofunctor → fmrilatent → fmrireg
pipeline. These are refinements to land when the downstream packages
(neurofunctor especially) are live enough to shape the requirements.

## Round 3 — Hybrid cortex+subcortex support type (deferred)

Full design captured in `dev/round3_cifti_design.md`. What we have
been calling a "hybrid support type" is essentially CIFTI-2
grayordinates: a flat length-N vector where each position is a
cortical surface vertex on L or R hemisphere or a subcortical voxel
in a per-structure volume grid. The design doc covers:

- Why component order is contract, not an implementation detail.
- The four-piece architectural split: `hybrid_support` primitive,
  `.wrap_decoded_hybrid()` helper, `HybridBasisTemplate` S3 class,
  `encode_transport()` matrix-input hybrid path.
- 13 gotchas ranked by how silently they break (component-order
  drift, medial-wall fills, per-structure subcortex NeuroSpaces,
  block-diagonal basis assumption, adjoint identity, wrap return
  shape, coefficient-space flatness, no-NeuroVec-input constraint,
  ROI subsetting dual modes, roughness penalties, sparse
  intermediates, serialization).
- A six-phase plan totalling ~700 LoC, ~40 tests, 3 new source files.
- Synthetic-fixture end-to-end test covering 20 assertions.
- Two open design questions (§6.2 producer-side hybrid_support
  strictness, §6.4 `template_loadings()` return type) that must be
  resolved before implementation starts.

Non-goals explicitly captured: CIFTI file I/O stays in `ciftiTools`,
no unified cross-domain basis functions, no `HybridNeuroObject`
wrapper class (wrap returns a plain named list), no nested subcortex
hierarchy for v1, no automatic neurofunctor composition.

Read `dev/round3_cifti_design.md` before starting implementation.

## Round 1 leftovers (deferred refinements)

### Producer-form vs canonical-form accessor

`transport_latent()` currently stores both the raw producer object at
`x$field_operator` / `x$observation_operator` AND the normalized
canonical map at `x$transport$field_operator`. The producer-form vs
canonical split is documented in `?portable_linear_map` but both forms
are still available on the top-level object. Consider a thin accessor
`transport_field_operator(x)` that always returns the canonical form,
and a separate `transport_raw_operator(x)` when debug access to the
producer is needed. Deprecate direct `$field_operator` reads over time.

### Flip `wrap = "auto"` default

`decode_coefficients()` currently defaults to `wrap = "none"` for
back-compat. Once downstream consumers (fmrireg engine adapters,
neurofunctor QA scripts) are all on the seam, flip the default to
`"auto"` in a minor-version NEWS cycle. This turns domain-aware
wrapping into the happy path and forces callers that want raw values
to opt in explicitly.

### Shared namespace for domain ids

The `strict` domain-id compatibility check in `.compose_linear_maps()`
is off by default because fmrilatent's own `basis_decoder()` methods
use content digests for `target_domain_id` while external field
operators use readable names. Once neurofunctor and fmrilatent agree
on a shared namespace convention (e.g., both sides use the same
content-addressed template id), enable `strict = TRUE` by default and
update the relevant tests.

## Critic's nits (low priority)

### Cargo-cult fields

`adjoint_convention` and `mode` are stored but not read at more than
one place each. `adjoint_convention` is now consumed by
`.require_euclidean_adjoint()` (round 1). `mode` is still advisory.
Either consume it in a specific code path or document it as advisory
in `?portable_linear_map`.

### `.compose_linear_maps()` provenance merge

Provenance lists are concatenated via base `c()` on named lists, which
lets later keys overwrite earlier ones on collision. For long
composition chains this can silently drop provenance entries. If this
becomes a real problem, introduce a `.merge_provenance()` helper that
preserves conflicts as a list of alternatives or warns on collision.

### `materialize()` fallback on large operators

`.materialize_linear_map()` falls back to `forward(diag(n_source))`
when no `materialize` callback is supplied. For full-brain operators
this can allocate gigabytes. Consider refusing the fallback above a
configurable `n_source` threshold and requiring the producer to supply
a real `materialize`.

### Contract version negotiation

The `contract_version` tag (round 1) currently just propagates the
max across composition. When the contract changes incompatibly (round
2+), add a real negotiation path: e.g., a version registry and a
`downgrade_portable_linear_map(op, target_version)` helper.
