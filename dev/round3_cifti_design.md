# Round 3 design: hybrid support type (CIFTI-aware grayordinates)

**Status**: design, not yet implemented.
**Deferred from**: audit follow-ups after commits 0a4665b / 935b27a / b902a28.
**Entry point in PLANNED.md**: "Round 3 — Hybrid cortex+subcortex support type".

What we are calling a "hybrid cortex+subcortex support type" is exactly the
problem CIFTI-2 exists to solve. HCP and Workbench ship billions of
grayordinate-shaped analyses on this concept. That means three things:

1. We have a working model to steal from.
2. There are ordering and medial-wall conventions we cannot invent ourselves.
3. Keeping CIFTI file I/O out of scope is critical or this spirals.

This document captures the design decisions, gotchas, phased plan, and
open questions for the v1 hybrid support path in `fmrilatent`. It is
load-bearing for any implementation work — read it before editing
`R/hybrid_support.R`, `R/hybrid_template.R`, or the encode/decode/wrap
routing code in `R/transport_latent.R` and `R/implicit_latent.R`.

---

## 1. What the grayordinate problem actually looks like

A CIFTI-2 "grayordinate" array is a flat length-N vector where each
position is either a cortical surface vertex (on L or R hemisphere) or
a subcortical volume voxel in a specific reference grid (typically 2mm
MNI). The ordering is load-bearing: HCP's 91,282-grayordinate standard
at fs_LR 32k is

```
[0 .. 29695]     CORTEX_LEFT  (32,492 fs_LR vertices with the ~2,796 medial-wall verts removed)
[29696 .. 59411] CORTEX_RIGHT (same story)
[59412 .. 91281] 19 subcortical structures in a fixed order
                 (accumbens_l, accumbens_r, amygdala_l, amygdala_r,
                  brain_stem, caudate_l, caudate_r, cerebellum_l,
                  cerebellum_r, diencephalon_ventral_l,
                  diencephalon_ventral_r, hippocampus_l, hippocampus_r,
                  pallidum_l, pallidum_r, putamen_l, putamen_r,
                  thalamus_l, thalamus_r)
```

Every downstream tool, statistics library, visualization package, and
file format assumes this layout. If `fmrilatent` stores a hybrid
support and slices in the wrong order, the values land in the wrong
hemispheres and the bug is invisible to dimension checks.

> **Component order is the contract**, not an implementation detail.

The CIFTI metadata attached to a real file uses `BrainModel` blocks.
Each block has a `ModelType` (`CORTEX_LEFT`, `CORTEX_RIGHT`, `VOXELS`),
a `BrainStructure` label, a `SurfaceNumberOfVertices` or
`VoxelIndicesIJK`, a `VertexIndices` (which vertices are used after
medial-wall masking) or voxel IJK list, plus `IndexOffset` and
`IndexCount` locating the block in the flat grayordinate array.

fmrilatent's `hybrid_support` object is a faithful but minimal
abstraction over that metadata: a named list of per-component
sub-supports plus an explicit component order.

---

## 2. Architectural split

This is bigger than the support type alone. A full CIFTI-aware bridge
needs four pieces, in increasing blast radius:

### 2.1 `hybrid_support` (primitive, Phase 1)

A new value type that describes a flat grayordinate layout as a named
list of per-component sub-supports plus a component order. Knows its
own total cardinality. Small, isolated. New file `R/hybrid_support.R`.

Shape:

```r
hs <- hybrid_support(
  components = list(
    cortex_l = lh_vertex_indices,        # integer vector
    cortex_r = rh_vertex_indices,        # integer vector
    subcortex_accumbens_l = al_mask_vol, # LogicalNeuroVol in MNI grid
    subcortex_accumbens_r = ar_mask_vol,
    ...
    subcortex_thalamus_r = tr_mask_vol
  ),
  component_order = names(components),
  brain_structures = c(cortex_l = "CORTEX_LEFT", ...)  # optional CIFTI labels
)
```

Operations:
- `n_components(hs)` → integer scalar
- `component_order(hs)` → character vector
- `components(hs)` → named list in order
- `component_offsets(hs)` → named integer vector of start offsets
- `component_cardinalities(hs)` → named integer vector of per-component sizes
- `.support_cardinality(hs)` → sum of per-component cardinalities
- `print(hs)` → shows per-component breakdown and total
- `is_hybrid_support(x)` → predicate

### 2.2 `.wrap_decoded_hybrid()` (Phase 2)

A helper that slices a flat length-N value vector/matrix by component
offsets and dispatches each slice through the existing
`.wrap_decoded_volume()` / `.wrap_decoded_surface()` path. Returns a
named list of per-component neuroim2 / neurosurf objects in
`component_order`. Lives in `R/implicit_latent.R`.

Vector path (single-column gamma):

```r
.wrap_decoded_hybrid <- function(values, hybrid_support, hybrid_domain = NULL, ...) {
  stopifnot(length(values) == .support_cardinality(hybrid_support))
  offsets <- component_offsets(hybrid_support)
  cards   <- component_cardinalities(hybrid_support)
  order   <- component_order(hybrid_support)
  out <- setNames(vector("list", length(order)), order)
  for (nm in order) {
    slice <- values[seq.int(offsets[[nm]], length.out = cards[[nm]])]
    comp_support <- components(hybrid_support)[[nm]]
    comp_domain  <- if (!is.null(hybrid_domain)) hybrid_domain[[nm]] else NULL
    out[[nm]] <- .wrap_decoded_by_support(
      slice, comp_support, comp_domain,
      context = paste0("wrap_decoded.hybrid[", nm, "]")
    )
  }
  out
}
```

Matrix path (multi-column gamma): slice rows, dispatch each component
slice through the matrix branch of `.wrap_decoded_volume` /
`.wrap_decoded_surface`, returning a named list of `NeuroVec` /
`NeuroSurfaceVector` objects.

### 2.3 `HybridBasisTemplate` (Phase 3)

A new basis asset class that holds:
- `components`: list of per-component sub-templates (e.g., a
  `SurfaceBasisTemplate` for each hemisphere + a `ParcelBasisTemplate`
  or a small-rank PCA template for each subcortical structure)
- `component_order`: character vector
- `coeff_layout`: named integer list mapping `name → integer offsets`
  into the total coefficient vector (so cortex_l uses coefficients
  1..k1, cortex_r uses k1+1..k1+k2, etc.)
- `meta`: list with `family = "hybrid"`, `label`, `k = sum(k_components)`,
  and per-component `components` label

Implemented as an S3 class (matches `SurfaceBasisTemplate` and
`ParcelBasisTemplate` which are both S3).

Template protocol methods:

| Method | Returns |
|---|---|
| `template_rank(x)` | `sum(k_components)` |
| `template_support(x)` | a `hybrid_support` |
| `template_domain(x)` | named list of per-component domains |
| `template_measure(x)` | named list of per-component measures (or NULL) |
| `template_roughness(x)` | block-diagonal sparse penalty (`Matrix::bdiag`) over per-component roughness matrices |
| `template_loadings(x)` | block-diagonal sparse matrix (`Matrix::bdiag`) over per-component loadings, so downstream code expecting a single loadings matrix keeps working |
| `template_mask(x)` | `stop()` — not meaningful for hybrid; use `template_support()` |
| `template_meta(x)` | meta list |
| `basis_decoder(x)` | a portable linear map (see below) |

`basis_decoder(hybrid_template)` returns a portable linear map whose
`forward` is block-diagonal:

```r
forward(z) {
  # z is k_total x ncols
  blocks <- vector("list", n_components)
  for (i in seq_along(components)) {
    zi <- z[coeff_layout[[i]], , drop = FALSE]
    blocks[[i]] <- basis_decoder(components[[i]])$forward(zi)
  }
  rbind.matrix(blocks)  # in component_order
}
adjoint_apply(y) {
  # y is N_total x ncols
  blocks <- vector("list", n_components)
  for (i in seq_along(components)) {
    yi <- y[support_offsets[i]:(support_offsets[i] + cards[i] - 1L), , drop = FALSE]
    blocks[[i]] <- basis_decoder(components[[i]])$adjoint_apply(yi)
  }
  rbind.matrix(blocks)  # in coeff_layout order
}
```

with `target_support = hybrid_support(...)`,
`target_domain_id = digest::digest(list(component_order, components))`,
`adjoint_convention = "euclidean_discrete"`, and an optional
`materialize()` that builds `Matrix::bdiag` of the sub-template
materializations (refuses above a configurable N threshold).

### 2.4 `encode_transport()` hybrid path (Phase 4)

The bare-matrix input branch of `encode_operator()` accepts
`support = hybrid_support(...)` and dispatches naturally;
`.compose_linear_maps(basis_map, obs_map)` composes a block-diagonal
`basis_map` with any hybrid-compatible `obs_map` without modification.

The NeuroVec input branch of `encode_operator()` is unchanged — it
stays volumetric-only. Hybrid targets require bare-matrix input.

Plumbing updates:

- `.support_cardinality()` gains a `hybrid_support` dispatch branch.
- `.subset_reconstruction_by_support()` gains a hybrid branch that
  accepts either a per-component named list ROI or a flat logical
  ROI vector of length N_total.
- `.resolve_transport_target_support()` passes through hybrid supports
  without trying to coerce them to a LogicalNeuroVol.
- `wrap_decoded.ImplicitLatent` detects hybrid support on the latent
  and routes to `.wrap_decoded_hybrid()`.

Everything else (coefficient recovery, decoder composition, covariance
pushforward) works unchanged because the portable linear-map contract
is already support-agnostic — it only cares about `forward`,
`adjoint_apply`, and dimensions.

---

## 3. Gotchas, ranked by how silently they break

### 3.1 Component order drift between operator and template

The field operator's `target_support` and the basis template's
component layout MUST agree on the grayordinate order. Dimensions
match trivially (total N is the same) so mismatches will not error —
they will just paint cortex_r values onto cortex_l and produce
plausible-looking nonsense.

**Mitigation**: at `encode_transport` time, assert that

```r
identical(
  component_order(op$target_support),
  component_order(template_support(basis_asset))
)
```

and name the mismatched components in the error. If the operator
carries a flat support (not hybrid), warn explicitly that
`wrap = "auto"` will fall back to flat wrapping.

### 3.2 Medial-wall fills

CIFTI cortex uses ~29,696 of the 32,492 fs_LR vertices per hemisphere;
the rest are medial wall. When wrapping back to a full-vertex
`neurosurf::NeuroSurface` the unused positions must be filled with 0
(or NA, per `ciftiTools` convention).

**Mitigation**: `.wrap_decoded_surface()` already does this correctly
via `neurosurf::NeuroSurface(domain, support, values)`. Reuse, do not
re-implement. Regression-test with a support that skips ~10% of
vertices to prove the medial-wall path works.

### 3.3 Subcortex has its own NeuroSpace

Subcortical voxels live in the MNI 2mm reference grid, not in any
other volume in the analysis. The subcortex component's
`LogicalNeuroVol` carries its own `NeuroSpace` and that's the right
place for it — do not try to promote it to a shared reference for
the whole latent object.

**Mitigation**: each subcortical component stores its own
`LogicalNeuroVol`. Do not put a single global `target_mask` on the
hybrid latent; rely on per-component supports.

### 3.4 Subcortex is 19 structures, not one

Real CIFTI stores ~19 subcortical BrainModels as peer components.
Two shapes are available:

- **Flat**: `list(cortex_l = ..., cortex_r = ..., accumbens_l = ..., amygdala_l = ..., ..., thalamus_r = ...)` — every structure is a peer.
- **Nested**: `list(cortex_l = ..., cortex_r = ..., subcortex = <one LogicalNeuroVol over all 19>)` — simpler, loses labels.

**Decision**: **flat**. Matches CIFTI fidelity, lets each subcortical
structure carry its own `NeuroSpace`, and costs nothing compared to a
monolithic volume. The constructor accepts any number of peer
components keyed by arbitrary strings with an optional
`brain_structures` attribute carrying CIFTI `BrainStructure` labels
for producers who care.

### 3.5 Block-diagonal basis is the v1 assumption

For `HybridBasisTemplate`, the composite `B` is block-diagonal: cortex-L
coefficients only decode to cortex-L samples, etc. This is correct for
concatenated surface bases (Slepian-L, Slepian-R) plus a subcortex
volume basis. It is NOT correct for a unified basis that crosses the
cortex/subcortex boundary — that is an open research problem (joint
graph-Laplacian eigenmodes, etc.) that `fmrilatent` should not pretend
to solve in v1.

**Mitigation**: document the assumption in `?hybrid_basis_template`.
Sub-templates produce their own non-block-diagonal bases within their
own component; hybrid composition is block-diagonal at the hybrid level.

### 3.6 Adjoint identity

If each sub-template's `basis_decoder` satisfies
`<Bx, y> = <x, B^T y>` in its Euclidean inner product, and the hybrid
composition is block-diagonal, the composite inherits the identity.

**Mitigation**: test this explicitly on a synthetic fixture — do not
just trust it. Include a `<Bx, y> = <x, B^T y>` probe at 1e-10
tolerance in `test-hybrid_support.R`.

### 3.7 Wrap return shape

`decode_coefficients(lat, gamma, space = "native", wrap = "auto")`
should return a named list in `component_order`:

```r
list(
  cortex_l = neurosurf::NeuroSurface(...),         # or NeuroSurfaceVector for multi-col
  cortex_r = neurosurf::NeuroSurface(...),
  subcortex_accumbens_l = neuroim2::DenseNeuroVol(...),  # or DenseNeuroVec
  ...
)
```

**Decision**: plain named list. **Do not** invent a new
`HybridNeuroObject` wrapper class. R users consume per-domain objects
directly; downstream code can wrap in a custom class if needed.
Keeping the return type as a plain list avoids forcing a new
dispatchable class on the ecosystem for v1.

### 3.8 Coefficient space stays flat

`k_total = sum(k_components)`. `coef_time(lat, "analysis")` returns a
`time × k_total` matrix. `fmrireg` fits in this flat coefficient space
and never needs to know the hybrid structure.
`decode_coefficients()` handles the demultiplexing on the way out.
This is actually a feature — it keeps the GLM-facing interface totally
unchanged, so the hybrid path is invisible to downstream modeling code.

### 3.9 No NeuroVec input path

`encode_operator()` currently accepts a `NeuroVec` input when the
target mask is volumetric. For hybrid targets, the only input path is
a bare `time × N_total` matrix. The user (or a helper in `ciftiTools`
/ `neurofunctor`) pre-extracts grayordinate-ordered samples.
`fmrilatent` does not read CIFTI files.

**Mitigation**: in `encode_operator()`, detect hybrid support and
error with a clear message if `x` is a NeuroVec; tell the user to
pre-extract into a matrix.

### 3.10 ROI subsetting

`.subset_reconstruction_by_support()` needs a third mode for hybrid.
Accept both forms:

- **Per-component ROI**: a named list matching component names
  (`list(cortex_l = <lh_roi>, cortex_r = NULL, ...)`; NULL means
  "keep all"); slice per-component.
- **Flat ROI**: a logical or integer vector of length N_total; apply
  directly to the flat reconstruction.

Error otherwise. Both are genuinely useful (cortex-only vs top-K
grayordinates).

### 3.11 Roughness and penalties

`template_roughness(hybrid_template, coordinates = "analysis")`
returns a block-diagonal penalty via `Matrix::bdiag` over per-component
roughness matrices. `encode_operator(spatial_penalty = ...)` accepts
either:

- a single `k_total × k_total` matrix (pre-stacked), or
- a named list of per-component penalties (stacked internally).

For v1, require the user to pass a single block-diagonal penalty or
accept the default identity — do not auto-compose per-component
penalties with different weights.

### 3.12 Sparse intermediate allocations

Block-diagonal composition temptingly allocates an
`N_total × k_total` dense decoder. For whole-brain grayordinates
(N ~ 91k, k ~ a few thousand) this is hundreds of MB.

**But**: the existing `.solve_transport_coefficients_matrix_free()`
path never materializes the full decoder — it calls
`forward` / `adjoint_apply` in blocks. The block-diagonal hybrid
decoder's `forward` / `adjoint_apply` are cheap because they just
dispatch per-component.

**Mitigation**: v1 keeps the encode path matrix-free throughout.
The `materialize` fallback in `.materialize_linear_map()` refuses at a
configurable N threshold (see Round 1 leftovers in `PLANNED.md`).

### 3.13 Serialization

`HybridBasisTemplate` stores a list of sub-templates in an S3 slot.
Each sub-template serializes fine via `saveRDS`. Nesting them in a
list is safe. **No closures in slots** — the `basis_decoder()`
closures live on the result of `basis_decoder(template)`, not on the
template object itself.

---

## 4. Key things to get right (the top five)

1. **Component order is contract.** Enforce
   `identical(component_order(op$target_support), component_order(template))`
   at encode time. Name mismatched components in the error.
2. **Block-diagonal is the v1 composite basis shape.** Do not promise
   cross-domain unified bases. Document the assumption.
3. **Wrap returns a plain named list** of per-component neuroim2 /
   neurosurf objects in `component_order`. No new wrapper class.
4. **Bare matrix is the only input path.** `fmrilatent` stays
   CIFTI-file-naive; file I/O lives in `ciftiTools` / `neurofunctor`.
5. **Subcortex peers, not a monolith.** Match CIFTI's per-structure
   `BrainModel` fidelity so downstream tools can consume `fmrilatent`
   outputs without losing labels.

---

## 5. Phased plan

| Phase | Scope | ~LoC | Tests | New files |
|---|---|---|---|---|
| **1** | `hybrid_support` primitive: constructor, accessors, `.support_cardinality()` dispatch | ~100 | 10 | `R/hybrid_support.R` |
| **2** | `.wrap_decoded_hybrid()` helper + `wrap_decoded.ImplicitLatent` routing | ~60 | 8 | — |
| **3** | `HybridBasisTemplate` S3 class + `basis_decoder` method + all `template_*` methods | ~250 | 15 | `R/hybrid_template.R` |
| **4** | `encode_operator()` matrix path accepts hybrid support + `.subset_reconstruction_by_support()` routing | ~50 | 6 | — |
| **5** | End-to-end integration test on a synthetic 2-hemisphere + tiny volume fixture (~30 grayordinates total) | ~150 | 1 big test | `tests/testthat/test-hybrid_support.R` |
| **6** | Docs: `?hybrid_support`, `?hybrid_basis_template`, vignette section, `PLANNED.md` cleanup | ~100 | — | — |

**Total: ~700 LoC, ~40 tests, 3 new source files + 1 new test file + 5 modified files.**

Phases 1 and 2 are independently valuable (they make flat surface
latents wrappable block-wise) and can land before Phase 3. Phase 3 is
the core and the biggest. Phase 4 is small plumbing. Phase 5 is the
integration test that proves the whole thing works end-to-end. Phase
6 is docs cleanup.

---

## 6. Open design questions to resolve before coding

### 6.1 Sub-template class system

`SurfaceBasisTemplate` and `ParcelBasisTemplate` are S3 (list with
class attribute). Match that — `HybridBasisTemplate` as S3.
**Resolved.**

### 6.2 Producer contract: operator-side hybrid_support

Must the field operator's `target_support` be a `hybrid_support`, or
just carry enough metadata to reconstruct hybrid structure from a flat
descriptor?

**Recommended**: require `hybrid_support` on the operator when the
producer knows its data is hybrid. If the operator only has a flat
support, encoding still works but `wrap = "auto"` falls back to flat
wrapping (dispatch by total cardinality). Assert strict agreement when
both sides are hybrid — named components, same order, same
cardinalities.

**Still open**: what does "strict agreement" mean concretely — do we
compare component names, per-component cardinalities, or full
per-component supports (e.g. the LogicalNeuroVol / vertex index
vectors)? Full support comparison is safest but may be expensive.
Cardinality + order is the minimum acceptable.

### 6.3 Coefficient `analysis_transform` under hybrid

Block-diagonal identity for v1. Per-component custom transforms can
be added later if needed; this preserves
`coef_metric(lat, "analysis") == diag(k_total)`. **Resolved.**

### 6.4 What exactly does `template_loadings(HybridBasisTemplate)` return?

Options:

- a block-diagonal `Matrix::bdiag(per_component_loadings)` — returns a
  single sparse matrix; downstream code treating `template_loadings()`
  as a matrix keeps working
- a list of per-component loadings — cleaner semantically but breaks
  any caller that expects a matrix
- `stop()` — most conservative, forces callers to go through
  `basis_decoder()` instead

**Recommended**: `Matrix::bdiag()`. Preserves the caller contract,
makes sparsity visible, and any code doing matrix math on hybrid
loadings has to deal with the block structure one way or another.

### 6.5 Where does this design doc live?

Committed as `dev/round3_cifti_design.md` with `dev/` added to
`.Rbuildignore`. Versioned alongside the code it describes, reviewable
on GitHub, survives session resets. **Resolved.**

---

## 7. Non-goals (explicit)

- CIFTI file I/O (read or write). Stay in `ciftiTools`.
- Unified cross-domain basis functions (joint Laplacian eigenmodes,
  graph heat kernels that span cortex and subcortex, etc.).
- Automatic neurofunctor operator composition for hybrid targets.
  That is neurofunctor's job; `fmrilatent` consumes the compiled
  operator through the portable linear-map contract.
- Per-component weighting in `encode_transport()` penalties beyond a
  single block-diagonal penalty for v1.
- Hybrid `analysis_transform` matrices that mix components.
  v1 is block-diagonal identity.
- Automatic promotion of subcortex components to a shared reference
  NeuroSpace. Each subcortical component carries its own.
- A new `HybridNeuroObject` wrapper class. Decoded output is a plain
  named list of per-component neuroim2 / neurosurf objects.
- A new `hybrid_support` S4 class. It is an S3 class (named list with
  class attribute) for ergonomic consistency with `LogicalNeuroVol`
  support descriptors in the existing code.

---

## 8. Interaction with the rest of the bridge

### 8.1 `?portable_linear_map` contract

`hybrid_support` is added to the list of accepted values for
`source_support` and `target_support` in the portable linear-map
contract. The concept page (`R/transport_latent.R` top) is updated
accordingly.

### 8.2 `.project_covariance_diag()` gate

The Round 1 `.require_euclidean_adjoint()` check is unchanged — the
hybrid composite inherits the euclidean_discrete convention from its
sub-templates when all of them are euclidean. If any sub-template is
non-euclidean, `HybridBasisTemplate`'s `basis_decoder()` propagates
that convention to the composite and `.project_covariance_diag()`
will correctly refuse at use time.

### 8.3 `.compose_linear_maps()` strict mode

The `strict_compose` domain-id check from Round 1 is still opt-in.
For a hybrid target, both sides can set a content-digest
`target_domain_id` derived from the ordered components list; when
both sides use the same digest convention, users can set
`options(fmrilatent.strict_compose = TRUE)` to enforce composition
safety. Not required for v1 correctness.

### 8.4 Round 2 template-space wrap

The `wrap_decoded.ImplicitLatent` method landed in Round 2 already
routes through `basis_asset(x)` → `template_support()` /
`template_domain()`. For a hybrid template, `template_support()`
returns a `hybrid_support` and `template_domain()` returns a named
list of per-component domains. The new `.wrap_decoded_hybrid()` path
handles both native and template wrapping without further changes to
the routing logic.

---

## 9. Test plan summary

`tests/testthat/test-hybrid_support.R` covers the full stack on a
synthetic fixture:

**Fixture**: `~30 grayordinates` total
- 10 "left hemisphere" vertices (a tiny fake `SurfaceGeometry` with
  12 total vertices, 2 of them medial wall; support = 1:10)
- 10 "right hemisphere" vertices (similar)
- 2x2x2 "subcortex_mni" volume (8 voxels, all in the mask)

**Per-component basis**: 2-atom surface basis per hemisphere + a
2-parcel volume basis for subcortex. `k_total = 6`.

**Tests**:
1. `hybrid_support()` constructor validates components, rejects empty,
   preserves component_order, computes offsets and cardinalities
   correctly
2. `.support_cardinality(hs)` returns 28
3. `component_offsets(hs)` returns `c(cortex_l = 1, cortex_r = 11, subcortex_mni = 21)`
4. `hybrid_basis_template()` assembles sub-templates, validates k_total,
   builds coeff_layout
5. `basis_decoder(tpl)` returns a portable linear map with
   `n_source = 6, n_target = 28`
6. `forward(diag(6))` returns a 28x6 block-diagonal matrix
7. Adjoint identity: `<Bx, y> == <x, B^T y>` to 1e-10
8. `template_loadings(tpl)` returns a sparse block-diagonal matrix
9. `template_roughness(tpl)` returns a block-diagonal penalty
10. `encode_transport(Y, basis_asset = tpl, field_operator = identity_op,
    support = hs)` returns a `TransportLatent`
11. `coef_time(lat, "analysis")` returns time x 6
12. `decode_coefficients(lat, gamma, space = "native", wrap = "none")`
    returns a length-28 numeric vector
13. `decode_coefficients(lat, gamma, space = "native", wrap = "auto")`
    returns a named list with `cortex_l` (NeuroSurface), `cortex_r`
    (NeuroSurface), `subcortex_mni` (DenseNeuroVol)
14. `decode_coefficients(lat, gamma, space = "template", wrap = "auto")`
    returns the same shape of named list for the template support
15. Multi-column gamma returns a named list of NeuroSurfaceVector /
    DenseNeuroVec objects with the correct number of time points
16. ROI mask: per-component named list (drop `subcortex_mni`) returns
    only `cortex_l` and `cortex_r` in the reconstruction
17. ROI mask: flat logical vector of length 28 (keep first 15)
    returns a sliced reconstruction with correct length
18. Component order drift: an operator whose `target_support` has
    reversed component order errors at encode time with a clear
    message naming the mismatch
19. `decode_covariance(lat, Sigma, space = "native", diag_only = TRUE)`
    returns a length-28 vector matching the block-diagonal reference
    to 1e-10
20. Round-trip: reference block-diagonal dense matrix D equals
    `.materialize_linear_map(decoder(lat, "native"))` to 1e-10, and
    `decode_coefficients(lat, gamma) == D %*% gamma` to 1e-10

---

## 10. Estimated scope (repeat from Phase table for easy reference)

- **Phase 1** (hybrid_support primitive): ~100 LoC, 10 tests, 1 new R file
- **Phase 2** (.wrap_decoded_hybrid + routing): ~60 LoC, 8 tests, 0 new files
- **Phase 3** (HybridBasisTemplate): ~250 LoC, 15 tests, 1 new R file
- **Phase 4** (encode_transport + subset routing): ~50 LoC, 6 tests, 0 new files
- **Phase 5** (end-to-end integration test): ~150 LoC, 1 big test, 1 new test file
- **Phase 6** (docs + vignette + PLANNED.md cleanup): ~100 LoC, 0 tests

**Total: ~700 LoC, ~40 tests, 3 new source files + 1 new test file + 5 modified files.**

---

## 11. Things to decide before coding starts

When this moves from design to implementation, the open questions in
§6 need resolution:

1. **§6.2** — producer-side `hybrid_support` requirement: strict
   comparison shape (cardinality only? full support?) and fallback
   semantics when the operator carries a flat support.
2. **§6.4** — `template_loadings(HybridBasisTemplate)` return type:
   confirm `Matrix::bdiag` vs per-component list.

Everything else (§6.1 S3 class, §6.3 block-diagonal analysis
transform, §6.5 doc location) is resolved.

---

## 12. Follow-up items this design explicitly defers

- A cross-domain unified basis (joint graph-Laplacian eigenmodes, etc.)
- Nested subcortex hierarchy (e.g., `subcortex$basal_ganglia$putamen_l`).
  Flat peer structures for v1; nesting can be added via a helper
  constructor later if anyone asks.
- CIFTI metadata round-trip (BrainStructure labels, IndexOffset,
  IndexCount). The `hybrid_support` object can carry an attribute for
  these but `fmrilatent` does not parse or emit CIFTI files.
- Performance tuning for the block-diagonal materialize path. v1
  relies on the matrix-free solver and refuses to materialize above a
  configurable N threshold.
- Hybrid `analysis_transform` matrices that mix components (e.g., a
  cross-domain whitening). v1 is block-diagonal identity.

These land on `PLANNED.md` as Round 4+ items once Round 3 is done.
