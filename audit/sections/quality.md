# Code Quality

## Summary

The package is in good quality health. The deep audit surfaced **no correctness bugs that affect normal usage** — every finding is either a defensive-validation gap, a diagnostic/error-message quality issue, or dead/duplicated code. Production paths are guarded (constructors validate, public defaults are safe), and the codebase is otherwise stylistically clean (consistent `seq_len`/`seq_along`, zero `1:length()`/`1:nrow()` occurrences, `setValidity` on most S4 classes). The defects that remain cluster into a handful of recurring patterns: missing input validation that defers failures to cryptic downstream errors, dead variables left by refactoring, and local variables that shadow generics. None are blockers; all are low-cost hardening opportunities that improve debuggability and maintainability.

Total findings: **2 medium, 17 low** (deduplicated to 13 themed entries below).

---

## Findings (ordered by severity)

### 1. `HierarchicalBasisTemplate` has no `setValidity` despite complex cross-slot invariants
**Severity: 🟡 Medium** · `R/all_class.R:227-238`

**What's wrong:** The class declares 8 slots with non-trivial cross-slot invariants (`nrow(loadings)` must equal `sum(mask)`; `nrow(atoms)` must equal `ncol(loadings)`; `gram_factor` must be a valid Cholesky factor of `t(loadings) %*% loadings`; `levels`/`parents` must be mutually consistent), yet no `setValidity("HierarchicalBasisTemplate", ...)` exists anywhere in the codebase. Six sibling S4 classes (`LatentNeuroVec`, `BasisHandle`, `LoadingsHandle`, `BlockLatentNeuroVector`, `LatentNeuroSurfaceVector`, `BilatLatentNeuroSurfaceVector`) all define `setValidity`, making this a clear deviation from the established pattern.

**Why it matters:** The `build_hierarchical_template()` constructor guards the primary path, but any path that bypasses it — `readRDS()`, `load()`, or a direct `new()` in tests — produces a silently invalid object. The first symptom is a cryptic error deep inside `Matrix::solve()` when `gram_factor` is used.

**Fix:** Add a `setValidity` immediately after the `setClass`, checking at minimum: `nrow(loadings) == sum(mask)`; `nrow(atoms) == ncol(loadings)`; `length(levels) >= 1L` with each element integer of length `nrow(loadings)`; `!is.null(gram_factor)`.

### 2. Silent rank reduction during orthonormalization emits no diagnostic
**Severity: 🟡 Medium** · `R/codec_boldzip_spatial.R:101-103` (root cause in `R/codec_boldzip_validate.R:126-138`)

**What's wrong:** `.boldzip_orthonormalize_columns()` returns `qr.Q(q)[, seq_len(rank), drop = FALSE]`; when the QR rank is less than `ncol(phi)`, it silently returns fewer columns than were passed in. The only hard stop is `rank < 1`. The call site assigns the result back to `phi` with no column-count check, and downstream `boldzip_sr_encode` validates row count and orthonormality but never the expected column count.

**Why it matters:** A caller passing a template with `k` columns expecting `k` basis vectors gets fewer, and the mismatch surfaces much later as a cryptic dimension error inside the encoder — hard to trace back to a dropped linearly-dependent column.

**Fix:** Emit a `warning()` in the rank-drop branch: `if (rank < ncol(x)) warning("orthonormalize: dropped ", ncol(x) - rank, " linearly dependent column(s); returning ", rank, " columns.", call. = FALSE)`.

### 3. Missing input validation that defers failure to opaque downstream errors
**Severity: 🟢 Low** · multiple files (themed)

Three findings share the root cause of *accepting an argument that is documented/typed as valid but has no working code path or no positivity check, then failing late with a misleading message*:

- **`R/shared_structure.R:510-541`** — `shared_temporal_spec()` accepts `kind = "slepian"` via `match.arg` (first-class, documented), but `materialize_shared_temporal_spec()` has branches only for `"dct"` and `"bspline"`; `"slepian"` falls through to an abort saying it is "a descriptor only and cannot be materialized without a custom basis." A real materialization path exists (`dpss_time_basis()` in `R/slepian_temporal.R`). **Fix:** add a `slepian` branch calling `dpss_time_basis()` with `params$n_time`/`half_bandwidth`/`n_tapers`, or validate at construction time that `kind="slepian"` requires an explicit `basis`.
- **`R/heat_wavelet.R:13-16`** — `basis_heat_wavelet()` stores `scales` with no positivity check; the value flows unvalidated into `rgsp::wavelet_heat_transform`. Heat kernels are `exp(-t·λ)`, so `t=0` yields a flat all-pass "filter" and `t<0` an amplifier — silent incorrectness, not a crash. The default `c(1,2,4,8)` is safe, so only erroneous user input triggers it. **Fix:** `if (!is.numeric(scales) || any(scales <= 0)) .encoder_cli_abort("scales must be a positive numeric vector", class = "fmrilatent_error_invalid_argument")`; also guard `order >= 1`, `threshold >= 0`.
- **`R/codec_boldzip_validate.R:42-44`** — `.boldzip_validate_matrix` rejects single-column input with `ncol(x) < 2L` and the message "must have at least one row and two columns," using column terminology for a voxels×time matrix and embedding a pairing-specific constraint (already enforced by `.boldzip_pair_indices` with a clearer `min=2L`) in a general validator. **Fix:** relax to `ncol(x) < 1L` (matching the sibling basis validators) and let `encode` enforce `n_time >= 2` with a domain-appropriate message.

### 4. Misleading error message when a list has none of the required contract fields
**Severity: 🟢 Low** · `R/transport_latent.R:337-358`

**What's wrong:** When a list has none of the four required fields (`forward`, `adjoint_apply`, `n_source`, `n_target`), `any(present)` is `FALSE` so the partial-field error branch is skipped, and execution falls through to the catch-all `.encoder_cli_abort("Unsupported ... contract.")`. An error is always raised (no silent wrong result), but the message does not tell the user which fields are missing — almost certainly a wrong-type mistake.

**Fix:** Change the guard to `if (length(missing) > 0L)` so partial and fully-absent lists both get the actionable "missing required fields: ..." message; removes the dead fall-through path.

### 5. `is_explicit_latent` dispatch failure swallows the informative error
**Severity: 🟢 Low** · `R/latent_block_vector.R:29-36`

**What's wrong:** `is_explicit_latent` has methods only for `ExplicitLatent` and `ImplicitLatent` (no `ANY` fallback). For a block of an unregistered class, the generic signals "unable to find an inherited method"; `try()` catches it and sets `ok_explicit <- FALSE`, producing "Block 'x' must be an explicit latent object" instead of surfacing the dispatch failure. All first-party block classes dispatch correctly, so only a foreign object hits this path.

**Fix:** Inspect the condition before discarding: when `explicit_try` is a `try-error`, include `conditionMessage(attr(explicit_try, "condition"))` in the reported error.

### 6. Validity check short-circuits dimension cross-checks on unrelated early errors
**Severity: 🟢 Low** · `R/latent_neurovector.R:353-407`

**What's wrong:** The cross-dimension checks are gated by `if (length(errors) == 0)`. The guard correctly avoids calling `dim(object@space)` etc. on a bad `@space`, but it also suppresses all dimension diagnostics when an *unrelated* earlier check fails (e.g. a bad `@label`). A `new()` call with both a bad label and wrong dimensions reports only the label error.

**Why it matters:** Incomplete (not misleading) error messages, only on edge-case direct `new()` with multiple simultaneous errors. Note: the blunt guard also protects `sum(object@mask)` / `object@map@indices` from type failures, so a naive fix that gates only on `space_ok` is incomplete.

**Fix:** Use per-prerequisite flags (`space_ok`, `mask_ok`, `map_ok`) so only the checks that genuinely require a given valid slot are skipped; label/offset errors should not suppress dimension diagnostics.

### 7. Redundant first validation pass discarded then repeated
**Severity: 🟢 Low** · `R/shared_structure.R:664-665`

**What's wrong:** Lines 664-665 call `validate_shared_component_contract` and `validate_template_protocol` via `invisible(lapply(...))`, discarding return values; lines 690-691 then store components as-is and re-run `validate_template_protocol`. The first pass is pure wasted work — any error it would catch is also caught by the second pass.

**Fix:** Remove the redundant first pass; rely on the validation done at store time (and store originals vs. manifests deliberately — see the separately-tracked manifests-vs-templates concern).

### 8. Local variables shadowing generics / base functions
**Severity: 🟢 Low** · multiple files (themed)

Three findings share the root cause of *binding a local variable to the same name as a function it shadows*:

- **`R/parcel_basis.R:340`** — `template_mask <- template_mask(tmpl)` shadows the exported S4 generic `template_mask` (`R/all_generic.R:349`) for the rest of `encode_spec.spec_space_parcel`. Works today only because no further call to the generic is made. The sibling `encode_hierarchical` already uses `tmpl_mask`. **Fix:** rename to `tmpl_mask`.
- **`R/codec_boldzip_diagnostics.R:152-153`** — `summary <- boldzip_sr_payload_summary(object)` shadows `base::summary` inside the `if (is_boldzip)` block. **Fix:** rename to `psum`.

Both are latent footguns: any future edit calling the shadowed name in scope silently dispatches on the local object.

### 9. Dead variables left by refactoring
**Severity: 🟢 Low** · multiple files (themed)

Five findings share the root cause of *an assigned variable that is never read*, misleading maintainers into thinking side-channel state or an optimization exists:

- **`R/latent_indexing.R:431`** — `idx <- 1` before the triple loop; never read (loop uses `ii`/`jj`/`kk`). **Fix:** delete the line.
- **`R/wavelet_active.R:47`** — `need_full_time <- meta$levels_time > 0`; never read (decoder always reconstructs all time points; the branch re-evaluates `meta$levels_time > 0` directly). The comment implies a partial-time-decode optimization that was never implemented. **Fix:** delete the variable and stale comment.
- **`R/encode_methods_space.R:70`** — `mask_arr <- .extract_mask_array(...)` in `encode_spec.spec_space_hrbf`; never read (the function uses raw `mask`). The only effect is the validation side-effect. Contrast `encode_spec.spec_space_pca:96` which uses `sum(mask_arr)`.
- **`R/encode_methods_space.R:179`** — same dead `mask_arr` in `encode_spec.spec_space_wavelet_active`; the side-effect is even redundant because `wavelet_active_latent()` re-validates the mask itself.
- **`R/parcel_basis.R:335`** — `reduction` accepted in the `encode_spec.spec_space_parcel` signature (per the generic contract) but never referenced; the template's own reduction is always used, silently ignoring a caller-supplied `reduction`.

**Fix:** For the dead `mask_arr` cases, replace the assignment with a bare validation call (or document the intent). For the ignored `reduction`, emit a warning when non-`NULL` (e.g. `class = "fmrilatent_warning_unused_reduction"`) so the silent discard becomes visible.

### 10. `k_time` is a required argument that is silently ignored
**Severity: 🟢 Low** · `R/latent_dct_heatwavelet.R:32`

**What's wrong:** `k_time` has no default (it precedes defaulted args), so callers must supply it, but the body never reads it — the component count comes from `k_components` derived from the heat-wavelet handle (line 70). The roxygen `@param` honestly documents it as "Ignored (kept for backwards compatibility)," which prevents silent runtime surprises but still forces callers to invent a meaningless value.

**Fix:** Give `k_time = NULL` and warn (`class = "fmrilatent_warning_deprecated"`) when a non-`NULL` value is supplied, signalling eventual removal.

### 11. `sapply()` used where type-stable `vapply()` is the package convention
**Severity: 🟢 Low** · `R/latent_methods.R:209, 283, 296`

**What's wrong:** These are the only three `sapply()` call sites in the entire `R/` tree. The return type/length is known at all three, and lines 209 & 296 feed the result straight into `all()`. On degenerate empty-list input `sapply()` silently returns `list()`, so the wrapping `all()` returns a vacuous `TRUE` rather than signalling — a silent logic inversion.

**Fix:** `all(vapply(all_objects, is, logical(1), "LatentNeuroVec"))` (209); `vapply(all_objects, function(obj) dim(obj@space)[4], integer(1))` (283); `all(vapply(all_basis, is, logical(1), "sparseMatrix"))` (296).

### 12. Duplicated identical `basis_mat`/`loadings_mat` method bodies
**Severity: 🟢 Low** · `R/latent_neurovec_materialize.R:25-61`

**What's wrong:** Four `setMethod` bodies (`basis_mat`/`loadings_mat` × `Matrix`/`matrix`) are byte-identical: `i <- i %||% seq_len(nrow(x)); j <- j %||% seq_len(ncol(x)); x[i, j, drop = FALSE]`. Real duplication in a core hot-path file.

**Fix:** Extract one private `.subset_mat(x, i, j)` and register a single method on a `setClassUnion(c("Matrix", "matrix"))`, or share one implementation across the two generics.

---

## Systemic patterns

Four root causes account for the bulk of findings and are worth addressing as classes rather than one-off:

1. **Deferred-failure validation gaps** (Findings 1, 2, 3) — arguments are accepted without checking invariants, so failures surface late and far from the cause (cryptic `Matrix::solve()` / dimension / `rgsp` errors). The package already has a validation idiom (`setValidity`, `.encoder_cli_abort`, `.extract_mask_array`); the fix is to apply it consistently at construction/spec boundaries.
2. **Misleading diagnostics on legitimate aborts** (Findings 2, 3, 4, 5) — the code correctly refuses bad input but emits a message that names the wrong problem (wrong terminology, "unsupported" instead of "missing fields," swallowed dispatch conditions). Low risk, high debuggability payoff.
3. **Refactoring residue** (Findings 8, 9, 12) — dead variables, generic-shadowing locals, and copy-pasted method bodies. Pure cleanup; a linter rule for unused locals and shadowed names would prevent recurrence.
4. **Generic-contract drift** (Findings 9 `reduction`, 10 `k_time`) — arguments retained for signature/back-compat conformance but silently ignored. Each should either do something or warn that it does not.
