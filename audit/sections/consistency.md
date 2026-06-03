# Consistency Audit

## Summary

The package is in the middle of a deliberate, incomplete migration toward classed
conditions (`.encoder_cli_abort()` / `.encoder_cli_warn()` with stable
`fmrilatent_error_*` classes), and almost every consistency defect found is a
straggler from that migration or a related API-surface drift. The encode
pipeline, spec constructors, and accessor generics are internally well-documented
and mathematically correct, but they exhibit recurring surface-level
inconsistencies: an entire codec family never adopted the classed-condition
convention; ~20 individual bare `stop()`/`message()` call sites remain across
validators, constructors, optional-dependency guards, plot helpers, and the C++
boundary; and sibling spec constructors disagree on argument names, defaults,
positions, return types, and signatures. None of these produce wrong numerical
output today — the one near-correctness item (the ROI-mask decoder path) fails
silently only under a caller mistake. The dominant risk is a fractured
error-handling contract: callers cannot reliably catch package errors by class,
and users moving between sibling encoder families cannot assume a stable API.

Severity legend: `[HIGH]` `[MEDIUM]` `[LOW]`

---

## Findings (by severity)

### 1. `[HIGH]` The entire `codec_boldzip` family bypasses the classed-condition convention

**Where:** `R/codec_boldzip.R:84-333`, `R/codec_boldzip_validate.R:20-134`,
`R/codec_boldzip_diagnostics.R:73-296`, `R/codec_boldzip_spatial.R:33-151`,
`R/codec_boldzip_decompose.R:39-114` (plus `codec_boldzip_st`).

**What's wrong:** The package standardized on `.encoder_cli_abort()`
(`cli::cli_abort` with `class = "fmrilatent_error_*"`), but the whole BOLDZip
codec family raises ~40 bare, unclassed `stop(..., call. = FALSE)` calls
(`simpleError` only). Per-file counts: `codec_boldzip.R` 17,
`codec_boldzip_validate.R` 11, `codec_boldzip_diagnostics.R` 7,
`codec_boldzip_spatial.R` 4, `codec_boldzip_decompose.R` 4. The codec tests
(`test-boldzip-sr.R`, `test-boldzip-sr-contracts.R`) call `expect_error()` ~18
times with **zero** `class =` assertions, and there is no
`test-codec_boldzip_error_classes.R` — every other migrated family has a
corresponding `test-*_error_classes.R` (awpt, encoder, haar_wavelet,
hierarchical_helpers, hierarchical_template, implicit_latent, latent_indexing,
shared_structure, transport).

**Why it matters:** This is the single largest unmigrated cluster. Callers cannot
catch BOLDZip failures by a stable `fmrilatent_error_*` class, and the test suite
locks in message-string matching, which is brittle and blocks future message
edits.

**Fix:** Add a codec-scoped classed abort helper mirroring `.encoder_cli_abort`
(e.g. `.boldzip_cli_abort` with `class = "fmrilatent_error_boldzip*"`), convert
the ~40 bare `stop()` calls, and add `tests/testthat/test-codec_boldzip_error_classes.R`
asserting `expect_error(..., class = "fmrilatent_error_*")`.

---

### 2. `[MEDIUM]` Decoder's `roi_mask` branch bypasses `.normalize_roi_mask`

**Where:** `R/wavelet_active.R:59-63`.

**What's wrong:** The decoder applies ROI subsetting with raw
`which(as.logical(roi_mask))` instead of routing through `.normalize_roi_mask()`,
which the package utility `roi_subset_columns()` (`implicit_latent.R:407-413`) and
`reconstruct_array.ImplicitLatent` both use. The raw path skips three guards:
(1) the dimension-identity check that aborts with `fmrilatent_error_dim` on a
shape mismatch, (2) `.mask_to_array` normalization of `LogicalNeuroVol` objects,
and (3) the out-of-mask voxel validation.

**Why it matters:** This is the only near-correctness item in the section. When
`roi_mask` has dimensions differing from `mask_arr` (e.g. a padded or sub-brain
mask), `which(as.logical(roi_mask))` enumerates flat indices in the wrong
coordinate frame, so the `global_idx %in% roi_global` comparison silently returns
wrong column indices and wrong reconstructed output — rather than an informative
error. Requires a caller mistake, but fails silently.

**Fix:** Replace the inline block with the existing validated drop-in:
```r
if (!is.null(roi_mask)) {
  rec_mat <- roi_subset_columns(rec_mat, mask_arr, roi_mask)
}
```

---

### 3. `[MEDIUM]` User-facing `encode()` and core-constructor validation use bare `stop()`

**Where:** `R/encode.R:191, 211, 217` (`.run_lengths_from_info`,
`.normalize_penalty_matrix`); `R/latent_dct_heatwavelet.R:76-78` (offset length
check). The same file mixes both mechanisms — `.encoder_cli_abort` with classed
conditions appears at `encode.R:108, 130, 309, 430, 438`, bare `stop()` at
`191, 211, 217, 361`.

**What's wrong:** Top-level `encode()` and a `LatentNeuroVec` constructor emit
unclassed `stop()` for dimension/length validation, while the surrounding pipeline
uses classed conditions. Two divergent error-emission mechanisms coexist inside
the very layer that defines the convention.

**Why it matters:** These are the highest-value catch points (user-facing entry
points and a constructor), yet downstream `tryCatch`/`withCallingHandlers` keyed
on `fmrilatent_error_*` cannot intercept them.

**Fix:** Route through `.encoder_cli_abort` with appropriate classes (e.g.
`fmrilatent_error_dim`, `fmrilatent_error_value`).

---

### 4. `[MEDIUM]` Optional-dependency guards use bare `stop()`

**Where:** `R/graph_bridge.R:12` (rgsp/heat wavelets);
`R/benchmark_roundtrip.R:14` (`bench` Suggests).

**What's wrong:** Missing-Suggests guards raise unclassed `stop()`, even though
the package already has a `fmrilatent_error_missing_dependency` class pattern
(`R/hierarchical_helpers.R:71-76, 401-402`). The cccde9c sweep converted 12 files
but missed these two.

**Why it matters:** Missing-dependency errors are exactly what callers want to
catch programmatically (to skip or fall back), yet they are not classed.

**Fix:** Raise
`.encoder_cli_abort(..., class = c("fmrilatent_error_missing_dependency", "fmrilatent_error"))`
at both sites.

---

### 5. `[LOW]` Bare `stop()` stragglers across validators, encoders, and plot helpers

A cluster of individual unmigrated call sites, all the same root cause as #1/#3/#4:

- **`R/encode_spec.R:119`** — `spec_space_heat` `scales` guard is the lone bare
  `stop()` in a file that otherwise calls `.validate_*`/`.encoder_cli_abort` 25
  times. Fix: `.encoder_cli_abort(..., class = "fmrilatent_error_invalid_scalar")`.
- **`R/slepian_spatiotemporal.R:25`** — `stop("X must have time in rows")` is
  unclassed **and** omits `call. = FALSE`, leaking call context; the sibling
  `slepian_temporal.R:63-65` uses `.encoder_cli_abort(..., class = "fmrilatent_error_invalid_argument")`
  for the identical check.
- **`R/slepian_plot.R:14, 48`** and **`R/spatial_plot.R:15-16`** — interactive
  plot helpers validate `basis`/`loadings`/`idx` with bare `stop()`. Lower
  priority (interactive), but inconsistent.

**Why it matters:** Each path is uncatchable by class; the spatiotemporal one
also produces a messier error message than its sibling.

**Fix:** Convert all to `.encoder_cli_abort` with appropriate
`fmrilatent_error_*` classes.

---

### 6. `[LOW]` Cholesky fallback uses bare `message()` instead of classed `.encoder_cli_warn()`

**Where:** `R/hierarchical_template.R:120` (also bare `message()` at
`R/latent_neurovector.R:231, 243`).

**What's wrong:** On Cholesky failure the QR-fallback diagnostic uses base
`message()` rather than the package-standard `.encoder_cli_warn()`. Bare
`message()` is not suppressible with `suppressWarnings()` and carries no class
for selective handling.

**Why it matters:** Inconsistent suppressibility and no classed handling; flagged
as a priority in commit cccde9c.

**Fix:**
`.encoder_cli_warn("Cholesky factorisation of Gram matrix failed; falling back to sparse QR.", class = "fmrilatent_warning_cholesky_fallback")`.

---

### 7. `[LOW]` C++ boundary aborts are unclassed while the R layer moved to classed conditions

**Where:** `src/haar_wavelet_rcpp.cpp:18-48, 61-66`;
`src/active_pencil_wavelet.cpp:106-126`; `src/hrbf_atoms_rcpp.cpp:19-31`;
`src/slepian_dpss_rcpp.cpp:26-30`.

**What's wrong:** All C++ validation failures raise plain `Rcpp::stop()`
(`simpleError`, no class). The matching R Morton guard uses
`class = "fmrilatent_error_morton_overflow"` (`R/haar_wavelet.R:81-86`, asserted
in `test-morton-bitwidth.R:46`), but the C++ aborts on identical code paths are
unclassed. `test-cpp_safety.R` therefore matches on message text only.

**Why it matters:** R-side `tryCatch` handlers keyed on a `fmrilatent_*` class
catch the R Morton overflow but silently miss the identical C++ abort
(dimension, level, Morton-cap, sigma). Consistency gap with the stated direction
— not a memory bug.

**Fix:** Translate known C++ message prefixes into classed `fmrilatent_*`
conditions at the R boundary (or re-raise in the thin RcppExports wrappers).

---

### 8. `[LOW]` `spec_*` sibling constructors disagree on argument contracts (API drift)

A family of related drifts where sibling spec constructors expose the same
concepts inconsistently. Root cause: organic growth without a shared signature
convention.

- **Component count `k` has three contracts** (`R/encode_spec.R:21, 43, 57, 73, 92`):
  required positional no-default (`spec_time_dct`, `spec_time_bspline`), optional
  `k = 3L` (`spec_space_slepian`, `spec_space_pca`), and optional `k = NULL`
  meaning auto-derive (`spec_time_slepian`). Heat/wavelet_active express it under
  entirely different names (`scales`, `levels_space`).
- **`threshold` has three meanings and defaults** (`R/encode_spec.R:117, 158`,
  `R/diffusion_wavelet.R:16`): `1e-6` sparsification floor (heat), `0` active-set
  cutoff (wavelet_active), `1e-5` compressed-operator sparsity (diffusion).
- **`k_neighbors` position drifts** (`R/encode_spec.R:73, 117`): position 2 in
  `spec_space_slepian` vs final position (4) in `spec_space_heat`.
- **`spec_space_hrbf` uses an opaque `params = list()` bag** (`R/encode_spec.R:147`)
  while every other spec uses explicit named, validated formals — hiding the
  parameter surface from help/autocomplete.
- **`spec_space_pca` omits the `scale` knob** its own `lift.spec_pca` implements
  (`R/pca_spatial.R:67-70`); the capability is unreachable through the public
  spec.

**Why it matters:** A user moving between encoder families cannot rely on `k`
being present/positional/same-defaulted, on `threshold` meaning the same thing,
on positional `k_neighbors`, or on a discoverable HRBF parameter surface.
Positional calls do not transfer. No wrong output — pure ergonomics/discoverability.

**Fix:** Adopt one convention per shared concept (uniform `k` with a documented
`k = NULL` auto rule; disambiguate `threshold` by name — `coef_eps` /
`select_threshold` / `sparsify_eps`; fix shared-arg positions; promote HRBF keys
to explicit formals with `params`/`...` as escape hatch; surface or document
`scale` in `spec_space_pca`).

---

### 9. `[LOW]` `encode_spec.*` methods and accessor generics have divergent signatures

- **`encode_spec` method formals diverge** (`R/encode_methods_time.R:7, 35, 60`;
  `R/encode_methods_hierarchical.R:22`): the three temporal methods and the
  hierarchical method omit `reduction` and rely on `...` to swallow it, while all
  spatial methods, `spec_st`, parcel, and awpt declare it explicitly. Named
  dispatch (`encode.R:323`) prevents runtime failure today, but sibling methods of
  one generic should share a signature; a positional call would misbind.
- **`offset` generic lacks `...`** (`R/all_generic.R:71-76`): declared as
  `function(object)` (to shadow `stats::offset`), forcing every method to use
  `object` and forbidding extra args, unlike sibling accessors `basis(x, ...)`
  and `loadings(x, ...)`. A method author following the `basis()`/`loadings()`
  pattern would mismatch the signature (S4 warns at load, not silent).

**Why it matters:** Latent API hazards — non-uniform formals invite
positional-call drift and method-registration mistakes.

**Fix:** Give every `encode_spec` method the uniform
`(x, spec, mask, reduction, materialize, label, ...)` signature; declare
`offset` as `function(object, ...)` (still `stats::offset`-compatible) and add a
roxygen note on the promotion behavior.

---

### 10. `[LOW]` Public name surfaces and return types are not aligned

- **Registry vs `latent_factory` family names** (`R/zzz.R:2-32` vs
  `R/encode.R:367-380`): `register_encoder("st", ...)` and
  `register_encoder("space_hrbf", ...)` produce `list_encoders()` names that are
  **not** accepted by `latent_factory()` — the factory only accepts `st_slepian`
  / `st_bspline_hrbf` (so `latent_factory("st", ...)` errors via `match.arg`),
  and `space_hrbf` is reachable only via the composite `st_bspline_hrbf`. Two
  registry entries break the `list_encoders()` -> `latent_factory()`
  discoverability contract. (Note: the broader "noun-last vs noun-first" framing
  is inaccurate — both surfaces use the same `domain_type` canonical convention;
  the real gap is these two unmapped entries.)
- **`spec_st` return-type asymmetry** (`R/encode_methods_st.R:101-145`): under the
  default `core_mode = "auto"`, `encode_spec.spec_st` returns an `ImplicitLatent`
  (S3, decoder-only) even when both factors are dense, whereas every explicit
  sibling returns a `LatentNeuroVec` (S4). Documented in `encode()`'s `@return`,
  and `core_mode = "explicit"` exists as opt-in, so the asymmetry is known — but
  callers iterating over specs cannot assume a uniform return class, and
  `LatentNeuroVec` methods silently do not apply to default `spec_st` output.

**Why it matters:** `list_encoders()` output cannot be fed directly into
`latent_factory()`; spec iteration cannot assume return class.

**Fix:** Map/align the two registry entries (or document the mapping). For
`spec_st`, keep the `@return` taxonomy authoritative and consider defaulting
`core_mode = "explicit"` when both factors are dense.

---

### 11. `[LOW]` `*_loadings_handle` constructors: only slepian gives `basis_spec` a default

**Where:** `R/slepian_handles.R:87-92` vs `R/heat_wavelet_handle.R:21-26` and
`R/diffusion_wavelet_handle.R:20-25`.

**What's wrong:** The three lifting-handle constructors share an identical body
shape (spec_payload + id + lift + register), but `slepian_spatial_loadings_handle`
defaults `basis_spec = basis_slepian()` while the heat and diffusion siblings
require `basis_spec` with no default.

**Fix:** Make the three symmetric — give all a sensible `basis_spec` default
(`basis_heat_wavelet()` / `basis_diffusion_wavelet()`) or none.

---

### 12. `[LOW]` `save_template` SurfaceBasisTemplate method ignores the generic's `xz` default

**Where:** `R/surface_template.R:246` vs generic `R/all_generic.R:428` and
`R/hierarchical_template.R:363`.

**What's wrong:** The generic declares `compress = "xz"` and the
`HierarchicalBasisTemplate` method mirrors it, but the `SurfaceBasisTemplate`
method declares `compress = TRUE`, so `saveRDS` uses gzip. A caller relying on the
documented default silently gets gzip output (~30-50% larger, format-inconsistent).

**Fix:** Change the method signature to `function(template, file, compress = "xz", ...)`.

---

### 13. `[LOW]` `encode_awpt` records `lambda = spatial_lambda` in meta, conflating two parameters

**Where:** `R/encode_operator.R:309`.

**What's wrong:** `encode_awpt()` passes `lambda = spatial_lambda`, so
`meta$lambda` (recorded as "ridge penalty strength", line 187) and
`meta$spatial_lambda` (line 188) are always equal. `encode_awpt()` exposes no
independent `lambda`, so callers cannot set them apart. Note: `lambda` is never
passed to any solver (the math at lines 144-170 uses only `spatial_lambda`,
`temporal_lambda`, `sparse_lambda`), so there is **no numerical corruption** —
the `meta$lambda` field merely carries a misleading label.

**Fix:** Either expose a separate `lambda` (default 0) and pass it independently,
or document that AWPT does not support a separate ridge and pass `lambda = 0`.

---

### 14. `[LOW]` DESCRIPTION says "offsets" (plural), diverging from the single per-voxel offset model

**Where:** `DESCRIPTION:12-18`.

**What's wrong:** DESCRIPTION reads "matrix factorization (basis x loadings plus
offsets)". Every other doc surface uses the singular per-voxel offset:
`README.md:14`, `R/all_class.R:11, 63`, `R/latent_neurovector.R:79-84`,
`vignettes/fmrilatent.Rmd:204`. The model has exactly one numeric vector of
length nVoxels.

**Why it matters:** The plural slightly misrepresents the model on the most
user-visible surface (the CRAN package page). No functional impact.

**Fix:** Reword to "matrix factorization (basis x loadings plus a per-voxel offset)".

---

## Systemic patterns

1. **Incomplete classed-condition migration (root cause of #1, #3, #4, #5, #6,
   #7).** Commit cccde9c began converting bare `stop()`/`warning()`/`message()`
   to `fmrilatent_error_*`/`fmrilatent_warning_*` classed conditions but covered
   only 12 files. The stragglers cluster predictably: an entire untouched family
   (BOLDZip, ~40 sites), user-facing entry points, optional-dependency guards,
   the C++ boundary, and scattered single-site validators/plot helpers. The fix
   is mechanical and uniform — route everything through `.encoder_cli_abort` /
   `.encoder_cli_warn` (adding a codec-scoped helper for BOLDZip and an R-boundary
   translation layer for C++), and add the missing
   `test-codec_boldzip_error_classes.R`.

2. **Sibling-constructor API drift (root cause of #8, #9, #10, #11).** Spec
   constructors and `encode_spec` methods grew family-by-family without a shared
   signature convention, producing divergent argument names, defaults, positions,
   formals, return types, and registry-vs-factory naming. None breaks current
   numeric output (named dispatch and per-constructor validation hold the line),
   but the surface is unpredictable for users moving between families. A single
   "spec constructor & method signature convention" doc plus a uniform
   `encode_spec` signature would resolve most of this cluster.

3. **One near-correctness outlier (#2).** The only finding with potential for
   wrong output is the ROI-mask decoder bypassing `.normalize_roi_mask`, and even
   that requires a caller-supplied dimension mismatch. Everything else is a
   catchability, suppressibility, ergonomics, or documentation concern — the
   package is numerically sound but contractually fractured.
