# Documentation

## Summary

The package's documentation is structurally sound and largely complete, but carries a cluster of CRAN-readiness gaps and accuracy drift that would surface on submission rather than in normal use. No documentation defect causes a runtime failure today. The recurring issues are: (1) exported topics missing a `\value` section, (2) self-contained examples needlessly hidden behind `\dontrun{}`, and (3) roxygen contracts that have drifted out of sync with the implementation (stale comments, under-enumerated argument choices, inconsistent `@param` wording, and undocumented return-shape caveats). All findings below are documentation-grade; none requires a behavioral code change to remediate, though several reflect genuine API-contract ambiguity worth resolving in code or docs.

## Findings (by severity)

### 1. Ten exported man pages lack a `\value` section (CRAN NOTE)
**Severity: medium** — `man/` (generated from `R/reduction.R`, `R/slepian_spatial.R:44-50`, `R/heat_wavelet.R:68-74`, `R/diffusion_wavelet.R:161-167`, `R/parcel_basis.R:185`, `R/surface_template.R:117`, `R/hierarchical_template.R:4-6,240-244`, `R/reduction.R:87-91`)

**What's wrong:** Ten non-internal Rd topics have `\usage` but no `\value`, which triggers the R CMD check / CRAN incoming "Functions with \usage but no \value" NOTE. Affected pages: `basis_slepian.Rd`, `basis_diffusion_wavelet.Rd`, `diffusion_wavelet_latent.Rd`, `heat_wavelet_latent.Rd`, `is_hierarchical_template.Rd`, `lift-GraphReduction-ANY-method.Rd`, `print.ParcelBasisTemplate.Rd`, `print.SurfaceBasisTemplate.Rd`, `save_hierarchical_template.Rd`, `slepian_spatial_latent.Rd`. Root causes by group: spec constructors documented with `@param` but no `@return`; `*_latent` constructors that return a `LatentNeuroVec` but omit `@return`; two `print` methods with no `@return`; and the hierarchical/lift helpers.

**Why it matters:** This is a guaranteed CRAN NOTE on submission and degrades the help pages for ten user-facing topics (users cannot see what each function returns).

**Fix:** Add `@return` to each function's roxygen. For `*_latent` constructors use `@return A \code{\link{LatentNeuroVec-class}}`; for the two `print` methods use `@return The input \code{x}, invisibly.`; for `basis_slepian`/`basis_pca`/`basis_flat` document the returned `spec_*` descriptor object. Re-run `devtools::document()` and `R CMD check` until the NOTE clears.

### 2. Runnable examples needlessly hidden behind `\dontrun{}` (CRAN policy)
**Severity: medium** — `R/all_generic.R:16,39,60`; `R/encoder_registry.R:105`; `R/latent_methods.R:195`; `R/latent_indexing.R:250`; `R/all_class.R:88`; `R/latent_neurovector.R:91`

**What's wrong:** The package has 10 `\dontrun{}` blocks and zero `\donttest{}` blocks. CRAN reserves `\dontrun{}` for examples that genuinely cannot run (credentials, network, interactive input, unavailable data). Most of these wrap trivially runnable code: `basis()`/`loadings()`/`offset()` only need an `lvec` the constructor builds in a few lines; the full `LatentNeuroVec(...)` constructor examples are self-contained (built from `Matrix` + `neuroim2`); `concat(...)` and `lvec[[1]]` are self-contained. The inconsistency is visible in `encoder_registry.R` itself: `get_encoder("time_dct")` is wrapped at line 105 while the equivalent `register_encoder`/`list_encoders`/`get_encoder` example at lines 57-60 runs unwrapped. Only `parcel_basis.R:64` and `:315` are legitimately non-runnable (they call `load_atlas()`).

**Why it matters:** Blanket `\dontrun{}` hides examples from R CMD check (so they are never exercised) and is a routinely-flagged CRAN reviewer objection; it also denies users runnable, copy-pasteable examples.

**Fix:** Convert self-contained examples to live code (construct `lvec` inline via the constructor). Use `\donttest{}` for slow-but-valid cases. Keep `\dontrun{}` only for the atlas-dependent `parcel_basis.R` examples.

### 3. Spec-constructor docs are thin and merged: `basis_slepian`/`basis_pca`/`basis_flat`
**Severity: low** — `R/reduction.R:47-66` (generates `man/basis_slepian.Rd`)

**What's wrong:** These three exported, user-facing spec constructors are aliased into a single Rd ("Basis specifications (lightweight descriptors)") whose `@description` merely repeats the title, with `@param k/type/whiten` but no `@return` and nothing distinguishing the three constructors. This is also the root cause of `basis_slepian.Rd` appearing in finding #1's missing-`\value` set. (Note: `R/slepian_spatial.R` itself is clean; `basis_slepian` is defined in `R/reduction.R`.)

**Why it matters:** Users cannot tell the three constructors apart or learn what object they return.

**Fix:** Expand the roxygen with a real per-spec `@description`, an `@return` describing the returned `spec_*` descriptor, and ideally split into per-function blocks. This subsumes the `basis_slepian` entry of finding #1.

### 4. `latent_factory(family=)` documents only 5 of ~19 accepted names
**Severity: low** — `R/encode.R:344` (vs. `choices`/`family_alias` at `R/encode.R:365-381`)

**What's wrong:** The `@param family` doc lists only `"dct_time"`, `"slepian_time"`, `"slepian_space"`, `"heat_space"`, `"slepian_st"`, but `match.arg` accepts ~19 values via the canonical `choices` and `family_alias` map (e.g. `space_pca`/`pca_space`, `space_parcel`/`parcel_space`, `space_wavelet_active`/`wavelet_active`, `st_bspline_hrbf`/`bspline_hrbf_st`, `st_slepian`, `hierarchical`).

**Why it matters:** Reachable families are undiscoverable from the docs; users learn the full set only by reading source or triggering a `match.arg` error.

**Fix:** Enumerate the full accepted set (canonical + aliases) in `@param family` and keep it in sync with `choices`/`family_alias`.

### 5. Stale comment + incomplete `@include`: validator dependencies misattributed
**Severity: low** — `R/encode_spec.R:5-10`

**What's wrong:** The file-level comment claims the `.validate_*` validators "live in `R/encode.R`" and are available because `encode.R` loads first. They actually live exclusively in `R/encoder_validators.R` (Collate line 83), and `.encoder_cli_abort` lives in `R/encoder_utilities.R` (Collate line 82) — both loaded *after* `encode_spec.R` (Collate line 78). The `@include encode.R` tag does not declare these true dependencies. (`encode.R:234-237` already has the correct pointer, making the `encode_spec.R` comment doubly wrong.)

**Why it matters:** No present failure — R resolves functions at call time, not parse time, and the hand-maintained Collate order happens to load dependencies first. But if Collate were ever regenerated solely from `@include` tags, `encode_spec.R` would order before its dependencies and the first spec-constructor call would hit "could not find function".

**Fix:** Correct the comment to point to `R/encoder_validators.R`, and add `#' @include encoder_utilities.R` and `#' @include encoder_validators.R`.

### 6. Inconsistent `@param data` contract across lift methods / handle constructors
**Severity: low** — `R/slepian_spatial.R:7`; `R/heat_wavelet.R:21`; `R/diffusion_wavelet.R:140`; `R/heat_wavelet_handle.R:11`; `R/diffusion_wavelet_handle.R:11`; `R/pca_spatial.R:57`

**What's wrong:** The same-position `data` argument is documented three contradictory ways: slepian/heat/diffusion lift methods say `Optional (unused) for API symmetry`; the heat/diffusion handle constructors say `Optional data passed through to lift() (often NULL)` — but for those families `lift()` ignores `data` entirely, so the "passed through" claim is misleading; and `pca_spatial.R:57` says `Required` (and indeed aborts on `NULL` at lines 72-76).

**Why it matters:** Readers of a shared positional API cannot tell which families consume `data` versus ignore it.

**Fix:** Standardize `@param data` wording per actual behavior — state explicitly which families ignore it, fix the handle-constructor "passed through" claim, and mark `spec_pca`'s `data` as required consistently in the lift method and any wrapper.

### 7. `basis_spec$whiten` silently ignored inside `lift()`
**Severity: low** — `R/pca_spatial.R:67-179` (esp. line 90); whitening done in `R/encode_methods_space.R:126-143`

**What's wrong:** `basis_pca(k, whiten)` correctly forwards both fields (`R/reduction.R:59`), but `lift.ClusterReduction,spec_pca` reads only `basis_spec$k` and never references `basis_spec$whiten`. Whitening is instead applied post-lift by the caller `encode_spec.spec_space_pca` (keyed on the *encode* spec's `whiten`, using the `fmrilatent.singular_values` attribute `lift()` attaches). A user calling `lift(reduction, basis_pca(k=5, whiten=TRUE), data=X)` directly gets un-whitened loadings with no warning, and the `lift()` roxygen (`R/pca_spatial.R:49-66`) never mentions the caller's responsibility. (The original evidence's claim that `basis_pca` hardcodes `whiten=FALSE` is incorrect; it forwards correctly — but the silent-ignore in `lift()` stands.)

**Why it matters:** Direct callers of `lift()` (including external packages) silently get the wrong result. Only affects the direct-`lift()` path, not the normal `encode()` pipeline.

**Fix:** Either document in the `lift,spec_pca` roxygen that `whiten` is not honoured by `lift()` (whitening is the caller's job), or emit a classed warning when `isTRUE(basis_spec$whiten)` is passed directly.

### 8. `build_dct_basis()` doc claims unconditional orthonormality
**Severity: low** — `R/dct_basis.R:3-8`

**What's wrong:** The title ("Build an orthonormal DCT-II basis matrix") and bare `@return Dense Matrix (n_time x k)` are unconditional, but with `norm="none"` the columns are raw cosines with non-unit L2 norms (column 0 ≈ `sqrt(n_time)`, others ≈ `sqrt(n_time/2)`). Note: the originally suspected "silent reconstruction breakage" does **not** occur — `encode_spec.spec_time_dct` uses the same `norm` to build the basis, compute loadings, and store the handle spec, and `materialize_basis_from_spec` rematerializes with the same `norm`, so the encode/decode pair is always internally consistent (including across cache eviction). The defect is purely documentation; the suggested runtime warning is not warranted.

**Why it matters:** Misleading return-shape contract for an internal (`@keywords internal`) helper; a caller relying on orthonormality with `norm="none"` would be surprised.

**Fix:** Make the title and `@return` conditional on `norm` (orthonormal only for `norm="ortho"`; raw, non-unit-norm columns for `norm="none"`). No code change needed.

## Systemic patterns

Three root causes recur across these findings:

- **CRAN-readiness gaps (findings #1, #2, #3):** missing `\value` sections and blanket `\dontrun{}` are the two most common CRAN-reviewer objections, and both are pervasive here rather than isolated. They should be fixed together as a single "CRAN doc pass" with a follow-up `R CMD check`.
- **Contract drift between roxygen and implementation (findings #4, #5, #6, #7):** documented argument choices, file-load comments, `@param` wording, and honoured spec fields have all drifted from what the code actually does. The shared signature seam (`data` across the lift/handle family; `family` choices in `latent_factory`; `whiten` across `basis_pca`/`lift`/`encode_spec`) is where drift concentrates, suggesting these multi-method families need a single authoritative contract doc.
- **Under-specified return/shape contracts (findings #1, #3, #8):** several topics omit or over-claim what they return. Adopting a convention (always document `@return`, and qualify it when behavior is parameter-dependent) would close this class.

None of these are functional defects today; all are remediable by documentation edits plus, for #5/#7, optional defensive hardening (`@include` declarations, a classed warning).
