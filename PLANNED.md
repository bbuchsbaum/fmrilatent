# Planned extensions

Non-blocking work items identified during the transport-bridge audit
(see commits 0a4665b, 935b27a, b902a28). None of these block the
pipeline seam — fmrilatent already plays its role in the
ANTs/FNIRT/AFNI → neurotransform → neurofunctor → fmrilatent → fmrireg
pipeline. These are refinements to land when the downstream packages
(neurofunctor especially) are live enough to shape the requirements.

## Round 3 — Hybrid cortex+subcortex support type (deferred)

CIFTI-style data flattens to a bare integer/logical support vector today,
which loses per-component semantics. To round-trip cortex-L + cortex-R +
subcortical volume cleanly:

- Introduce a `hybrid_support` type, probably as a named list of the
  form `list(cortex_l = lh_support, cortex_r = rh_support, subcortex = vol_mask)`
  with an explicit total cardinality and a components order attribute.
- Extend `.support_cardinality()` (R/transport_latent.R) to sum
  per-component cardinalities.
- Add a `.wrap_decoded_hybrid()` helper in R/implicit_latent.R that
  slices a flattened length-N value vector according to component
  boundaries and returns either a structured named list of
  NeuroVol/NeuroSurface objects or a new hybrid wrapper class.
- Route `wrap_decoded.ImplicitLatent` through `.wrap_decoded_hybrid()`
  when the resolved support is a hybrid_support.
- Decide whether hybrid_support should be a plain named list (minimal
  new API) or an S4 class with `cardinality()` / `slice()` methods
  (cleaner, more dispatch). Recommend waiting until neurofunctor has
  a concrete CIFTI consumer before committing to a shape.
- End-to-end test: encode a hybrid-target operator, decode back to
  the three components, verify per-component dimensions and content.

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
