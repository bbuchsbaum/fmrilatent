# Encode API Contracts — Normalization Proposal

Mote: bd-01KSGC9SAWBRH0VHDKJACK32WS "Normalize encode API contracts".
Type: investigation + proposal. No source files were modified.
Date: 2026-05-30.

All file:line citations were read directly from source. Line numbers for
the dense-message site (Finding 5) were confirmed via `grep -n` after the
read tool returned duplicated lines for that region; the load-bearing
facts (gating option, threshold, constructor path) are verified.

---

## 1. Summary (executive)

- **Return type is tri-modal and under-documented.** `encode()`'s
  `@return` (R/encode.R:207) says "`LatentNeuroVec` (explicit bases) or
  `ImplicitLatent` (separable cases)", but a third class `ExplicitLatent`
  exists (virtual marker, R/all_class.R:23) and `spec_st` *always*
  returns `ImplicitLatent` (R/encode_methods_st.R:103) even for fully
  explicit separable bases (e.g. DCT × HRBF). Callers cannot predict the
  class from the spec family without reading source. Highest-value doc
  fix; matters for `neuroarchive` compat.

- **Argument-name drift is narrow.** The component count is uniformly
  `k` across every spec constructor (R/encode_spec.R:21,43,57,73,92) —
  there is NO `rank`/`n_components`/`n_basis` drift; that consistency
  should be locked in by docs. Real drift is concentrated in
  `latent_factory()`, which invents `k_time`/`k_space` and applies the
  `k` fallback unevenly (R/encode.R:300 vs :312-316).

- **Offset/centering contract is inconsistent.** Only `spec_space_pca`
  computes/stores an offset (R/encode_methods_space.R:104-123); every
  other explicit family hard-codes `offset = numeric(0)`. The PCA path
  applies centering in two places (inside `lift.basis_pca` via
  `center=`/`offset=`, R/pca_spatial.R:69/128-130, AND again at the
  caller via `mu_scores`, R/encode_methods_space.R:120-123). A third
  independent centering implementation lives in
  `.template_projection_payload` (R/encode_template.R:159-165). No
  documented "who owns the offset" contract.

- **Duplicate encode kernels exist.** The shared Gram/least-squares
  helper `.basis_coefficients` (R/encode.R:62-68) is used by the temporal
  methods but bypassed by the spatial atom paths, which re-implement the
  same `gram/rhs/.robust_gram_solve` kernel by hand: HRBF
  (R/encode_methods_space.R:74-76) and spec_st atoms
  (R/encode_methods_st.R:71-73). The two sparse transport solvers
  duplicate an entire FISTA scaffold and objective tail
  (R/encode_transport_solve.R:279 vs :347; :43 vs :128).

- **The "dense, 100% non-zero" message is a constructor-level message,
  not a diffusion-specific warning, and the diffusion dense path triggers
  it legitimately.** The `message()` lives in the shared
  `LatentNeuroVec`/`.make_latent_neurovector` constructor at
  R/latent_neurovector.R:204 (basis) and :214 (loadings), gated by a
  `density > 0.5` threshold (R/latent_neurovector.R:203, :213). Diffusion
  reaches it because `diffusion_wavelet_loadings()` returns a *dense*
  Matrix by design (R/diffusion_wavelet.R:129) and
  `diffusion_wavelet_latent()` feeds it to `.make_latent_neurovector`
  (R/diffusion_wavelet.R:179). For the global diffusion basis the dense
  output is expected, so the message is a false positive there.

- **Highest-value fixes, in order:** (1) document the return-type
  contract per family + mention `ExplicitLatent`; (2) consolidate the
  spatial Gram projection onto `.basis_coefficients`; (3) document +
  unify the offset contract; (4) suppress the constructor dense message
  for known-dense families (diffusion) via a per-call flag; (5) make
  `latent_factory()` arg aliasing uniform.

---

## 2. Findings

### Finding 1 — Return-type contract

**Evidence.**

- `encode()` `@return`: R/encode.R:207
  `#' @return A LatentNeuroVec (explicit bases) or ImplicitLatent (separable cases).`
- `latent_factory()` `@return`: R/encode.R:261
  `#' @return A LatentNeuroVec or ImplicitLatent object.`
- Third class `ExplicitLatent` (virtual S4 marker): R/all_class.R:23
  `setClass("ExplicitLatent", representation("VIRTUAL"))`; `LatentNeuroVec`
  and surface/block latents inherit it (R/all_class.R:135,166,187,207).
  `ImplicitLatent`/`TransportLatent` are S3 and explicitly NOT explicit
  (R/all_class.R:32-33; predicate FALSE at R/implicit_latent.R:301). Not
  named in any `@return`.
- Per-family actual returns:
  - slepian/heat/hrbf/pca space -> `LatentNeuroVec`
    (R/encode_methods_space.R:31, :63, :89, :159)
  - `spec_space_wavelet_active` -> whatever `wavelet_active_latent()`
    returns; class not asserted at the seam
    (R/encode_methods_space.R:179-185)
  - time slepian/dct/bspline -> `LatentNeuroVec`
    (R/encode_methods_time.R:30, :55, :96)
  - `spec_st` -> **always** `ImplicitLatent` via `implicit_latent(...)`
    (R/encode_methods_st.R:103-108), even when both time and space are
    fully explicit dense bases.

**Problem.** Three return classes, no per-family documentation, and a
family (`spec_st`) whose return class is independent of whether an
explicit factorization exists. Downstream `neuroarchive` compat must
branch on class with no documented contract. The wavelet_active seam
does not even guarantee a known return class.

**Proposed fix.**
1. Add a `@return` table to `?encode`/`?encode_spec` enumerating the
   class per spec family (doc-only).
2. Document the class taxonomy: explicit families -> `LatentNeuroVec`
   (an `ExplicitLatent`); decoder-only/separable families ->
   `ImplicitLatent`; explain `ExplicitLatent` as the virtual marker.
3. Assert the return class of `wavelet_active_latent()` at
   R/encode_methods_space.R:179 (cheap classed `stopifnot`).
4. **DECISION NEEDED:** Should fully-explicit `spec_st` (e.g. DCT×HRBF)
   return `LatentNeuroVec` instead of `ImplicitLatent`? The `core_mode`
   slot already exists (`c("auto","explicit")`, R/encode_spec.R:159) but
   is currently inert. Activating `core_mode = "explicit"` to emit an
   explicit pair would be more consistent, but changing the *default*
   return class is a behavior change that breaks callers relying on
   `predict()`/`as.matrix.ImplicitLatent`. Recommended interim: keep
   `ImplicitLatent` as the documented contract for all `spec_st`, add the
   opt-in `core_mode = "explicit"`.

**Risk.** Doc: low. wavelet_active assert: low. `spec_st` default change:
high. `core_mode` opt-in: medium.
**Blast radius.** Doc: none functional. `spec_st` default change: every
spec_st caller + tests/testthat/test-encode*.R + neuroarchive compat.
**Deprecation needed?** Yes if `spec_st` default return changes — gate
behind `core_mode`, keep `"auto"` == current behavior.

---

### Finding 2 — Argument-name drift

**Evidence (consistent — lock in, do not "fix").**
- Component count is uniformly `k`:
  `spec_time_slepian(tr, bandwidth, k=NULL, ...)` R/encode_spec.R:21;
  `spec_time_dct(k, norm)` R/encode_spec.R:43;
  `spec_time_bspline(k, ...)` R/encode_spec.R:57;
  `spec_space_slepian(k=3L, k_neighbors=6L)` R/encode_spec.R:73;
  `spec_space_pca(k=3L, ...)` R/encode_spec.R:92;
  `basis_pca` reads `basis_spec$k` (R/pca_spatial.R:90).
  No `rank`/`n_components`/`n_basis` spelling appears anywhere in the
  spec layer — genuine non-finding worth recording.
- `tol` is used consistently in solvers
  (R/encode_transport_solve.R:203,288,448;
  `.symmetric_matrix_factor(tol=1e-10)` R/encode_template.R:74). No
  `tolerance` spelling in the encode layer — non-finding.
- `mask` is the uniform public arg (R/encode.R:209 and all methods);
  internal arrays are `mask_arr` (R/encode_methods_space.R:9). No public
  `mask_3d_array` in the encode layer (that name appears only in
  lower-level Haar code, e.g. R/haar_wavelet.R:60, out of scope here).

**Evidence (real drift).**
1. **`latent_factory()` invents `k_time`/`k_space`** read from `...`:
   - `st_bspline_hrbf` branch: `args$k_time %||% args$k %||% 5L`
     (R/encode.R:300) — has the `k` fallback.
   - `st_slepian` branch: `k = args$k_time %||% NULL` (R/encode.R:312-314)
     and `k = args$k_space %||% 3L` (R/encode.R:315) — does NOT fall back
     to `k`. So `latent_factory("slepian_st", k = 5)` silently ignores
     `k`. Inconsistent within the same function.
2. **`seed` placement is inconsistent across families.** Top-level
   `basis_diffusion_wavelet(..., seed = 1L)` (R/diffusion_wavelet.R:17),
   but HRBF buries `seed` inside the `params` bag of
   `spec_space_hrbf(params = list())` (R/encode_spec.R:130,133;
   validated at R/encoder_validators.R:81-82). Same concept, two
   surfaces (named arg vs list key).

**Problem.** Public-surface drift is limited to (a) `latent_factory`'s
uneven `k`/`k_time`/`k_space` handling and (b) HRBF `seed`/params hidden
in a bag while diffusion exposes `seed` as a named arg.

**Proposed fix.**
- Make the `latent_factory` fallback uniform:
  `k_time = args$k_time %||% args$k`,
  `k_space = args$k_space %||% args$k` in the `st_slepian` branch
  (R/encode.R:312-316), matching the bspline_hrbf branch.
- Document HRBF `params` keys (`sigma0`, `levels`, `radius_factor`,
  `kernel_type`, `seed`) in `?spec_space_hrbf` (R/encode_spec.R:128); the
  doc line already lists them (R/encode_spec.R:130) — promote to a
  `\describe` block. Optionally accept top-level `seed`/`sigma0` via `...`
  merged into `params` for surface parity with diffusion (additive).

**Risk.** Low (additive fallbacks + docs).
**Blast radius.** `latent_factory` callers; HRBF additive only.
**Deprecation needed?** No — additive/back-compatible.

---

### Finding 3 — Offset / centering contract

**Evidence.**
- **Only PCA stores an offset.** `encode_spec.spec_space_pca`:
  `offset <- numeric(0); if (isTRUE(spec$center)) offset <- colMeans(x)`
  (R/encode_methods_space.R:104-107); stored on the object via
  `offset = offset` (R/encode_methods_space.R:165).
- **PCA centering is applied twice.**
  (a) Inside `lift.basis_pca`: it accepts `center=`/`offset=`
  (R/pca_spatial.R:69) and subtracts cluster-local means
  `X_loc <- sweep(X_loc, 2, mu[vox_idx], "-")` (R/pca_spatial.R:128-130).
  (b) Again at the caller: after `basis <- x %*% loadings`
  (R/encode_methods_space.R:119), it subtracts the offset's projection
  ```
  mu_scores <- as.matrix(crossprod(offset, loadings))
  basis <- basis - matrix(1, n_time, 1) %*% mu_scores
  ```
  (R/encode_methods_space.R:120-123). Correctness depends on the
  undocumented invariant that `lift.basis_pca` does NOT itself fold the
  mean into the returned basis scores.
- **Every other explicit family hard-codes no offset:** slepian
  R/encode_methods_space.R:32, heat :64, hrbf :90; time-slepian
  R/encode_methods_time.R:31, dct :56, bspline :97.
- **`center` exists only on `spec_space_pca`** (R/encode_spec.R:92;
  doc at R/encode_spec.R:85-87 states means are stored in `@offset`). No
  `center`/`centering` arg on any other spec.
- **A third centering implementation** lives in
  `.template_projection_payload(..., center=FALSE, offset=NULL, ...)`:
  `offset_out <- offset %||% colMeans(X); X_proj <- sweep(X,2,offset_out,"-")`
  (R/encode_template.R:159-165).
- **Reconstruction contract** is fixed by the class doc:
  `data[v,t] = sum_k basis[t,k]*loadings[v,k] + offset[v]`
  (R/all_class.R:63). Families with `offset = numeric(0)` rely on the
  reconstruction treating empty offset as zero.

**Problem.**
1. No single helper or doc for "who computes the offset, when, how decode
   applies it". Each family re-decides.
2. PCA centering applied in two coupled places; fragile invariant.
3. `center = FALSE` is only meaningful for PCA/template; other families
   cannot center at all, and `numeric(0)` vs a zero vector is an implicit
   convention.

**Proposed fix.**
- Document the offset contract in `?encode`/`?LatentNeuroVec`: encoder
  computes any per-voxel offset, stores it in `@offset`, reconstruction
  adds it back; non-centering families store `offset = numeric(0)`.
- Extract one helper `.encode_center(x, center, offset = NULL)` returning
  `list(offset, x_centered)`, used by the PCA caller and
  `.template_projection_payload`, deduping the `colMeans`/`sweep` logic.
- **DECISION NEEDED:** Should the PCA caller-side `mu_scores` correction
  (R/encode_methods_space.R:120-123) move into `lift.basis_pca` so
  centering is owned in exactly one place? This changes the documented
  `lift` generic contract (it would return offset-corrected basis
  scores) and `lift` may have external implementers — needs human
  sign-off.

**Risk.** Helper extraction: low-med (numerics must stay bit-identical;
guarded by existing PCA round-trip tests). Moving correction into `lift`:
medium-high (public generic contract change).
**Blast radius.** PCA + template encoders; the `lift` change crosses the
`lift` generic API.
**Deprecation needed?** Helper: no. `lift` contract change: yes (NEWS +
behavior-change note).

---

### Finding 4 — Duplicate encode kernels

**Evidence (Gram / least-squares projection).**
Shared helper already exists and is the intended path:
```
.basis_coefficients(x, basis, ridge, context)   # R/encode.R:62-68
.temporal_loadings_from_basis(x, basis, ...)     # R/encode.R:70-73
```
Temporal methods use it: R/encode_methods_time.R:17, :44, :83. The
spec_st temporal projection also uses it (R/encode_methods_st.R:75-79).
But two spatial/atom paths re-implement the identical kernel by hand:
- **HRBF** R/encode_methods_space.R:74-76:
  ```
  gram <- as.matrix(B_atoms %*% Matrix::t(B_atoms))
  rhs  <- t(as.matrix(x) %*% Matrix::t(B_atoms))
  coeff <- t(.robust_gram_solve(gram, rhs))
  ```
- **spec_st space atoms** R/encode_methods_st.R:71-73 (same shape):
  ```
  gram <- B_atoms %*% t(B_atoms)
  rhs_st <- t(as.matrix(x) %*% t(B_atoms))
  C_atoms <- t(.robust_gram_solve(as.matrix(gram), rhs_st))
  ```

**Evidence (transport solver near-duplication).**
`.solve_transport_coefficients_sparse` (R/encode_transport_solve.R:279)
and `.solve_transport_coefficients_sparse_matrix_free`
(R/encode_transport_solve.R:347) share the whole FISTA scaffold:
t_k update :327 vs :410; prox `.prox_sparse_awpt` :321 vs :404;
divergence guard :332 vs :424; the only real difference is dense `A`/`Lt`
matmuls vs matrix-free `.transport_apply_quadratic_system`. The two
objectives `.awpt_objective` (R/encode_transport_solve.R:43) and
`.awpt_objective_matrix_free` (:128) duplicate the penalty tail
character-for-character (lasso/group_l2 at :54-60 vs :143-149).

**Evidence (default-k re-derivation).**
The slepian default `floor(2*n_time*bw*tr) - 1` is computed independently
in `encode_spec.spec_time_slepian` (R/encode_methods_time.R:10) and
`encode_spec.spec_st` (R/encode_methods_st.R:18).

**Problem.** Three copies of the Gram projection, two copies of the FISTA
loop + objective tail, two copies of the slepian-k default. Bug fixes
must be replicated N times.

**Proposed fix.**
1. Route HRBF (R/encode_methods_space.R:74-76) and spec_st atoms
   (R/encode_methods_st.R:71-73) through one `.atom_coefficients()`
   wrapper around `.robust_gram_solve` (or generalize
   `.basis_coefficients` to accept a pre-transposed atom matrix).
2. Factor the FISTA scaffold into `.fista(grad_fn, obj_fn, prox, step,
   ...)`; both sparse solvers supply closures.
3. Factor the penalty tail into `.awpt_penalty_term(Z, sparse_lambda,
   sparse_mode)`.
4. Extract `.slepian_default_k(n_time, tr, bandwidth)` used by both
   slepian methods.

**Risk.** (1)/(4) low (mechanical, covered by encode round-trip tests).
(2)/(3) medium (solver refactor; needs transport tests green).
**Blast radius.** (1) hrbf + st encoders; (4) slepian time + st; (2)/(3)
AWPT/transport only.
**Deprecation needed?** No — all `.`-prefixed internal helpers.

---

### Finding 5 — Diffusion dense-loadings "false warning"

**Correction to the framing.** The "dense, 100% non-zero" string is NOT
emitted from `R/diffusion_wavelet.R`. It is a generic `message()` in the
shared `LatentNeuroVec` constructor (`.make_latent_neurovector`):

- R/latent_neurovector.R:204 (basis):
  `message("Input 'basis' is dense (", round(density_basis * 100), "% non-zero); storing as dense dgeMatrix.")`
- R/latent_neurovector.R:214 (loadings): analogous `message()`.
- Both are gated only by a `density > 0.5` threshold
  (R/latent_neurovector.R:203, :213) — there is currently NO option-based
  or per-call suppression, only the implicit "store sparse if <=50%
  dense" branch.

**Why diffusion trips it.** `diffusion_wavelet_loadings()` returns a
*dense* matrix by design:
`Matrix::Matrix(do.call(cbind, loadings_list), sparse = FALSE)`
(R/diffusion_wavelet.R:129). `diffusion_wavelet_latent()` passes that
through `lift()` (R/diffusion_wavelet.R:177) and into
`.make_latent_neurovector(X, mask, loadings, ...)`
(R/diffusion_wavelet.R:179), which runs the density check and emits the
message. For the global `phi_global` diffusion basis the dense result is
the expected output, so the message is a **false positive** in that path.
(The constructor message only fires when the input arrives as a base
`matrix`; if loadings are already a `Matrix`, the `is.matrix() && !is(
,"Matrix")` guard at R/latent_neurovector.R:201/211 skips it — but the
diffusion path produces a `Matrix`, so confirm at fix time whether the
message actually fires in the dense-Matrix case or only when a base
matrix is supplied. Either way the no-suppression design is the gap.)

**Problem.** A by-design dense result emits a user-facing `message()`
that the diffusion family cannot suppress, and there is no per-call or
per-family flag to opt out. That pollutes `neuroarchive` logs and trains
callers to ignore a diagnostic that IS useful for slepian/heat/hrbf where
unexpected post-threshold density signals a bug.

**Proposed fix.**
- Add a `quiet_dense`/`expect_dense` argument (default `FALSE`) to
  `.make_latent_neurovector()` / `LatentNeuroVec()`
  (R/latent_neurovector.R:126) that, when `TRUE`, skips the density
  `message()` at R/latent_neurovector.R:203-205 and :213-215.
- Have `diffusion_wavelet_latent()` (R/diffusion_wavelet.R:179) pass
  `expect_dense = TRUE`. Keep the message active for slepian/heat/hrbf.
- Lower-effort alternative: gate the message behind
  `getOption("fmrilatent.quiet_dense", FALSE)` so callers can silence it
  globally. Inferior — global, not family-scoped — but trivial.

**Risk.** Low (message-suppression only; no numeric change).
**Blast radius.** `.make_latent_neurovector`/`LatentNeuroVec` signature
(additive arg) + diffusion encoder; check tests/testthat for any
`expect_message` asserting on this string before changing it.
**Deprecation needed?** No — additive arg with back-compatible default.

---

## 3. Prioritized plan (incrementally committable)

### Bucket A — Low-risk mechanical (each ~1 commit + targeted test)

**Status (2026-05-30):** A1, A2, A3, A4, A5, A6, A7 are **DONE**
(implemented + tested + committed). B3 is also **DONE**. B1, B2, C1, C2
remain open.

- A1. **DONE.** Doc-only: added per-family `@return` table to `?encode`
  and `?latent_factory`, documenting the `ExplicitLatent` /
  `ImplicitLatent` / `TransportLatent` taxonomy. (F1.1/1.2)
- A2. **DONE.** Assert return class of `wavelet_active_latent()` at the
  `spec_space_wavelet_active` seam via a classed
  `fmrilatent_error_invalid_type` abort; class-asserting test added. (F1.3)
- A3. **DONE.** Routed HRBF (R/encode_methods_space.R) and spec_st
  atoms (R/encode_methods_st.R) through one shared `.atom_coefficients()`
  helper (R/encode.R), which forms the Gram over the voxel axis and
  returns time x atoms coefficients. Numerically identical to the prior
  inline `.robust_gram_solve` expressions. (F4.1)
- A4. **DONE.** Extracted `.slepian_default_k(n_time, tr, bandwidth)`
  (R/encode.R) and routed the temporal (R/encode_methods_time.R) and
  spatiotemporal (R/encode_methods_st.R) slepian encoders through it.
  Formula unchanged. (F4.4)
- A5. **DONE.** Added `expect_dense` arg (default `FALSE`) to
  `LatentNeuroVec()` / `.make_latent_neurovector()` suppressing the
  by-design dense message; `diffusion_wavelet_latent()` passes
  `expect_dense = TRUE`. Tests cover both the suppressed and the
  genuinely-unexpected dense paths. (F5)
- A6. **DONE.** `latent_factory()` bare `k` now falls back uniformly to
  `k_time`/`k_space` in the `st_slepian` branch, matching
  `st_bspline_hrbf`; test added. (F2.1)
- A7. **DONE.** Documented HRBF `params` keys in `?spec_space_hrbf` as a
  `\describe` block (incl. `seed`). (F2.2)

### Bucket B — Medium (aliases/shims/helper extraction + tests)
- B1. Extract `.encode_center()` and use it in the PCA caller and
  `.template_projection_payload`; preserve numerics (PCA round-trip test
  stays green). (F3 fix-1)
- B2. Document the single offset contract in `?encode`/`?LatentNeuroVec`
  and add an `@section Offset` to family docs. (F3)
- B3. **DONE.** Factored `.fista()` (shared FISTA scaffold via
  grad/obj/prox closures) and `.awpt_penalty_term()` (shared sparsity
  penalty tail) out of the two sparse transport solvers in
  R/encode_transport_solve.R. Step sizes, prox operator, momentum
  update, divergence guard, and convergence checks are bit-identical to
  the prior inline loops; behavior-preserving. (F4.2/4.3)

### Bucket C — DECISION NEEDED (human sign-off before implementing)
- C1. **`spec_st` return class.** Keep always-`ImplicitLatent` (document
  as the contract) vs activate `core_mode = "explicit"`
  (R/encode_spec.R:159) to emit `LatentNeuroVec` for fully-explicit
  separable specs. Behavior change; affects neuroarchive compat.
- C2. **PCA centering ownership.** Move the caller-side `mu_scores`
  correction (R/encode_methods_space.R:120-123) into `lift.basis_pca`
  so centering lives in one place. Changes the public `lift` generic
  contract; possible external implementers.

---

## 4. Out of scope (belongs to other motes)

- **codec_boldzip** -> bd-01KSQQPCEQSHQ6WTVMCXZR1VME.
  `R/codec_boldzip.R` is not part of the encode-spec contract surface
  (its `stop(...)` at R/codec_boldzip.R:152 etc. belongs there).
- **Pure `stop()` -> classed condition conversions** ->
  bd-01KSWB3Y55C6XAM44AHHXQB51M. Plain `stop()` calls remain in the
  encode layer (R/encode.R:146,165,170; R/encode_template.R:38,45,97;
  R/encode_transport_solve.R:73,84,234,333,425,491;
  R/latent_neurovector.R:186) and are tracked there, not here.
- Lower-level `lift()`/handle materialization internals beyond the offset
  and dense-message seams (e.g. the `mask_3d_array` arg name in Haar code,
  R/haar_wavelet.R:60) are not part of the encode-spec public contract.
