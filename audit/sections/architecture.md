# Architecture

## Summary

The package is internally consistent and the layering it advertises (`latent_handles → all_class → all_generic → encoders`) is broadly honored, but the foundational layers carry more weight and more coupling than their stated scope warrants. No critical or high-severity structural defects were confirmed; every finding is a design-cohesion or load-order fragility issue rather than a correctness or security bug. The recurring theme is that the two "openness" claims of the architecture — the encoder registry (extensible families) and the handle/registry abstraction (self-contained lazy dictionaries) — are each leakier than documented: the lazy-handle path is effectively closed to external extension, and the bottom of the stack reaches upward into encoder/reduction contracts. A secondary theme is that the earliest-loaded file (`latent_handles.R`) has become a catch-all whose internal structure the rest of the package implicitly depends on.

## Findings (ordered by severity)

### 1. The lazy-handle mechanism is closed to extension and inverts the dependency direction

**Severity: 🟡 Medium** | `R/latent_neurovec_materialize.R:96-97,118-123` · `R/latent_handles.R:15,72`

*(Merged: "Registry/materialization reaches into encoder layer via lift()" + "Handle 'kind' is a closed enum hard-coded in three places" — shared root cause: handle kinds are a closed, hard-wired set rather than a registry.)*

**What's wrong.** Handle `kind` is a fixed enum enforced by `setValidity` in two places (`.known_basis_handle_kinds`, `.known_loadings_handle_kinds`) and re-encoded a third time as an `if/else` dispatch in `materialize_basis_from_spec()` / `materialize_loadings_from_spec()`. Two of those branches (`"lifted"`, `"slepian_spatial"`) call `lift()` — a generic from the reduction/contracts layer — so the bottom of the stack (handle materialization) depends upward on the encoder/reduction layer. The `@include` at `materialize.R:4` must name `bspline_basis.R`, `dct_basis.R`, `slepian_handles.R`, and `reduction.R` precisely because of this upward reach. A consequence is that a `BasisHandle` is not self-contained: its spec must carry a live, non-serializable `reduction` object to be materializable.

**Why it matters.** `encoder_registry.R` advertises `register_encoder()` + `encode_spec` S3 dispatch as the seam for external families. But any family needing a lazy/handle-backed basis cannot register a new `kind` without editing three locations in *this* package's source. The extensibility seam therefore stops at explicit dense/sparse output; the `"explicit"` kind (inline matrix) is the only escape hatch, and it forgoes lazy reconstruction. This is also the concrete manifestation of the load-order fragility CLAUDE.md already flags.

**Recommended fix.** Replace the closed enum + `if/else` with a kind registry: `register_handle_kind(kind, materializer)`. `setValidity` checks membership in the registry; `materialize_*` dispatches through it. Each encoder family then registers both its `encode_spec` method and its handle materializer, restoring the intended top-down dependency and making the lazy path open for extension.

### 2. Two parallel cache subsystems must be kept manually in sync

**Severity: 🟢 Low** | `R/latent_handles.R:107-159,180-195`

**What's wrong.** The LRU cache is split into `.fmrilatent_cache_env` (storage) and a separate `.fmrilatent_cache_order` tracker (parallel character vectors). The split exists (per the comment at 126-129) only so `.latent_get_registry_env()` can hand a plain environment to the stats/list/clear helpers. Every mutation path (register, get, clear, evict) must touch both structures, and `.latent_enforce_cap()` contains explicit drift-reconciliation logic ("any untracked id is treated as least-recently-used").

**Why it matters.** A permanent two-structure consistency obligation is bought to give three introspection helpers a bare environment. The reconciliation code is currently defensive — no live path actually produces drift — so the risk of a silent correctness bug is low, but the invariant is fragile under future edits.

**Recommended fix.** Store LRU metadata as a companion record inside the same environment (or use a single ordered structure) and adapt the three introspection helpers to read from it, eliminating the parallel tracker and its reconciliation branch.

### 3. The advertised `encode_spec` seam is not actually uniform

**Severity: 🟢 Low** | `R/encode.R:359-364,436-442`

**What's wrong.** `encode()` → `encode_spec` (S3 on spec class) is documented as the uniform extensibility seam, yet `encode_spec.spec_awpt_wavelet` exists only to abort and redirect to `encode_awpt()`/`encode_operator()`, and `latent_factory()` hard-codes an `awpt` guard that `stop()`s. The transport/AWPT family lives outside the advertised dispatch and needs its own entry points.

**Why it matters.** The single documented API surface (`encode`/`encode_spec`/`latent_factory`) is not uniform; callers must know which families are second-class. Impact is low because the abort is immediate and self-documenting (no silent misbehavior). *(The `spec_st` → `ImplicitLatent` return-type asymmetry, sometimes cited alongside this, is a working, intentional, documented method — not a broken contract — and is treated under Finding 4, not here.)*

**Recommended fix.** Either route `encode_awpt`/`encode_operator` through `encode_spec` like every other family (passing shared assets via the spec), or stop advertising `encode_spec` as universal and document transport as an explicit parallel API in the registry's dispatch-model section.

### 4. No unifying latent supertype; the type system is split across S4 and S3

**Severity: 🟢 Low** | `R/all_class.R:23-37,135`

**What's wrong.** Explicit latents are S4 (`ExplicitLatent` virtual parent; `LatentNeuroVec` et al. `contains` it), while `ImplicitLatent` and `TransportLatent` are S3, registered into S4 only via `setOldClass`. There is no `setClassUnion` or shared supertype spanning both branches. `is_explicit_latent()`/S4 dispatch covers one half; S3 `UseMethod` covers the other.

**Why it matters.** Any cross-cutting generic ("any latent": `concat`, `series`, `save`) must be written twice or guarded — and indeed `concat`/`series` are currently implemented only for the explicit S4 branch. The absence of a common dispatch target is what forces the long `@return` taxonomy at `encode.R:253-273`. The design is internally consistent and the split is documented, hence low severity.

**Recommended fix.** Introduce a single virtual/abstract `Latent` supertype (a `setClassUnion`, or a shared S4 virtual class both branches register under) so cross-cutting generics have one dispatch target — or explicitly document that no unifying type exists and standardize the dual-implementation pattern.

### 5. `latent_handles.R` is a mixed-concern, foundational catch-all (including `%||%`)

**Severity: 🟢 Low** | `R/latent_handles.R:457,459-583`

*(Merged: "latent_handles.R has grown into a mixed-concern utility module" + "%||% defined mid-file but used package-wide" — shared root cause: zero-dependency primitives and unrelated domain helpers accreted into the package's earliest-loaded file.)*

**What's wrong.** The file's stated purpose is "Handle classes and registry for shared/implicit dictionaries," but it now also defines mask/space/dimension utilities (`mask_to_array` :468, `.normalize_roi_mask` :494, `.space_with_time_from_mask` :517, `.assert_template_mask_match` :531) that have no relationship to the handle cache, plus the `%||%` null-coalesce operator (:457) consumed by 30+ files across the package. `latent_handles.R` is the second `Collate` entry, so everything transitively depends on it. (The `.latent_basis_dim`/`.latent_loadings_dim` helpers at 553-583 *do* dispatch on the handle classes and legitimately belong here.)

**Why it matters.** Foundational utilities whose availability is an implicit dependency on this file's internal structure amplify the load-order coupling CLAUDE.md already flags. There is no current broken load order, but any `Collate` reorder that moves a consumer ahead of this file — or any refactor that splits it — would silently break package load. The operator's placement is incidental to its package-wide importance.

**Recommended fix.** Move `%||%` and other zero-dependency primitives, plus the four mask/space helpers, into a dedicated `latent_utils.R` placed first in `Collate`. Keep `latent_handles.R` focused on the handle classes and cache (retaining the dim helpers, which belong with the classes).

## Systemic patterns

Three of the seven findings (1, 2, 5) trace to a single root cause: **the foundational base of the stack — handle classes, the cache, and `latent_handles.R` itself — carries responsibilities and couplings beyond its declared scope.** Specifically:

- **Closed sets where the architecture promises openness.** Handle `kind` (Finding 1) and the `encode_spec` seam (Finding 3) are both advertised as extension points but are in practice hard-wired enums / special-cased families. The fix pattern is identical in both cases: replace hard-coded dispatch with a registry, mirroring the existing `encoder_registry`.
- **Implicit dependencies on file/Collate structure.** Findings 1, 5 (and the spirit of CLAUDE.md's own load-order warning) all stem from foundational symbols whose availability depends on `@include` ordering and `Collate` position rather than on a clean dependency graph. Consolidating zero-dependency primitives into a first-loaded utilities file and inverting the handle-materialization dependency via a registry would remove most of this fragility at once.

None of these block correctness today; they are accumulating coupling that raises the cost and risk of future extension and refactoring.
