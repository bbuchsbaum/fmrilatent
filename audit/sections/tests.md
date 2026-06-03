# Tests

## Summary

The package's test suite is healthy at the end-to-end level: the headline pipelines (BOLDZip super-resolution encode/decode, the transport/AWPT solver path, and the encoder families) all have happy-path round-trip tests with tight reconstruction tolerances (e.g. `boldzip_sr_decode == manual` at 1e-12, `reconstruct_matrix == Y` at 1e-8). The gaps are not absences of coverage but deficits in *assertion depth*: several large, numerically load-bearing internal helpers — transport solvers, shared encoder validators, shared encoder utilities, and BOLDZip texture-fit internals — are reached only transitively through these end-to-end tests and have zero direct, path-isolated assertions. A regression inside any of them that still leaves the end-to-end output within tolerance would pass undetected. None of the verified findings rise above medium severity, and one prior basename-based "missing coverage" flag was retracted after function-name verification.

## Findings

### 1. Transport solvers (519-line numerical core) have no path-isolated equivalence or conditioning tests

`![medium]` &nbsp; **`R/encode_transport_solve.R:295-519`** (also `:191`, `:235`)

**What's wrong.** The largest encode-pipeline file is the numerical core of the transport/AWPT path — `.solve_transport_coefficients_matrix_free` (295), `.solve_transport_coefficients_sparse` (319), `.solve_transport_coefficients_sparse_matrix_free` (382), the dispatcher `.solve_transport_coefficients` (459), the CG solver `.cg_transport_quadratic` (235), the FISTA driver `.fista` (77), and `.estimate_transport_lipschitz` (191). A grep of the entire `tests/testthat/` tree for any of these symbols returns zero matches. They are exercised only indirectly through `encode_operator()` / `encode_transport()` in the AWPT/transport suites. The wrapper round trip is asserted (`test-transport_latent.R:145`, `reconstruct_matrix == Y` at tol 1e-8), and there are two path-labelled "does not materialize" tests, but nothing pins the relationships that distinguish the internal paths from each other.

**Why it matters.** A single happy-path test over a 519-line solver cannot catch the regressions that separate these code paths: (a) matrix-free vs. explicit dense normal-equation solve diverging, (b) the sparse path drifting from the dense path at `sparse_lambda=0`, (c) CG failing to converge to the direct solution on an ill-conditioned Gram, or (d) a Lipschitz under-estimate that destabilizes FISTA steps.

**Fix.** Add `test-encode_transport_solve.R` with metamorphic equivalence tests on a small fixed decoder map: (1) matrix-free result == explicit dense ridge solve within 1e-8; (2) sparse path == dense path at `sparse_lambda=0`; (3) CG result == direct solve on a deliberately ill-conditioned Gram; (4) `.estimate_transport_lipschitz >=` true operator spectral norm on a known matrix.

### 2. Shared encoder validators are nearly untested; NA/Inf/non-integer and scalar/flag branches have zero assertions

`![medium]` &nbsp; **`R/encoder_validators.R:11-99`**

**What's wrong.** Six shared validators gate every encoder family. Only two assertions touch them — `test-encoder_error_classes.R` calls `.validate_positive_count(0, "k")` and `.validate_hrbf_params(list(kernel_type = "bad"))`. Untested entirely: `.validate_nonnegative_count`, `.validate_positive_scalar`, `.validate_nonnegative_scalar`, `.validate_flag_scalar`. The error classes `fmrilatent_error_invalid_scalar` and `fmrilatent_error_invalid_flag` appear nowhere in the suite. The NA / non-finite / `length != 1` / non-integer rejection branches that every validator contains are also unexercised — e.g. line 12 (`length(x) != 1L || is.na(x) || !is.finite(x) || x < 1 || !isTRUE(all.equal(x, round(x)))`) is hit only via the value `0`, the `x < 1` branch. Most of `.validate_hrbf_params` (sigma0 / levels / radius_factor / seed / kernel_type_fine_levels, lines 66-97, and non-list input) is unasserted.

**Why it matters.** A regression in a shared validator silently weakens input checking across *all* encoders at once — exactly the cross-cutting failure that is hardest to localize from downstream symptoms.

**Fix.** Add `test-encoder_validators.R` calling each validator directly with NA, Inf, `length != 1`, negative, and non-integer (for counts) inputs, asserting the exact condition class (`fmrilatent_error_invalid_count`/`_scalar`/`_flag`), plus one positive case per validator confirming the coerced return value. Add `.validate_hrbf_params` cases for bad sigma0/levels/seed/kernel_type_fine_levels.

### 3. Shared encoder utilities: empty-cluster sparse-dictionary branches and the basis-synthesis contract are untested

`![medium]` &nbsp; **`R/encoder_utilities.R:39-117`**

**What's wrong.** These helpers are reused across encoder families (`.build_sparse_dictionary_from_clusters` is called by both `slepian_spatial.R` and `heat_wavelet.R`; `.make_latent_neurovector` by `encode_methods_space.R` and others). A grep for their names across `tests/testthat/` matched only `test-hrbf.R`, and that match is for the unrelated public wrapper `extract_mask_array_logical` — neither helper is invoked directly anywhere. Two correctness-load-bearing branches have no assertion: (1) `.build_sparse_dictionary_from_clusters` skips empty clusters (line 49) and returns a 0-column sparse matrix when all clusters are empty (line 62); (2) `.make_latent_neurovector` synthesizes the basis as `as.matrix(X) %*% as.matrix(loadings)` when `basis` is NULL (line 104) — the exact `basis × loadings` factorization contract the package is built on.

**Why it matters.** An off-by-one in the column accounting of the sparse-dictionary path would silently corrupt the output of *every* cluster-based encoder, surfacing only as wrong-shape outputs or downstream dimension mismatches far from the cause.

**Fix.** Add `test-encoder_utilities.R`: (1) `.build_sparse_dictionary_from_clusters` with a reduction containing one empty cluster (assert column count) and an all-empty reduction (assert 0-column result); (2) `.make_latent_neurovector` with `basis=NULL` asserting `basis == X %*% loadings` plus a reconstruction round trip; (3) `.extract_mask_array` list-mask fallback (lines 12-21) and the missing-dim error branch (line 22).

### 4. BOLDZip texture-fit internals are exercised only through the end-to-end RMSE round trip

`![low]` &nbsp; **`R/codec_boldzip_decompose.R:138-172`**

**What's wrong.** The rest of the BOLDZip codec is well covered — `boldzip_sr_encode`/`boldzip_sr_decode` round trips assert reconstruction equality (`test-boldzip-sr-contracts.R:182`, tol 1e-12) and `rmse < 0.15` (`:243`), `.boldzip_quantize_values` has 7 direct invocations, and `.boldzip_encode_events`/`.boldzip_pair_indices` are tested directly (`:493-495`). But the carrier/texture prediction core — `.boldzip_lagged_carrier_bank` (138), `.boldzip_predict_texture` (159), `.boldzip_fit_sparse_texture` (172) — has zero direct grep matches and is reached only via the encode/decode RMSE round trip.

**Why it matters.** A regression that still keeps end-to-end RMSE under 0.15 — e.g. wrong lag selection or a suboptimal sparse texture fit — would pass undetected. This is an assertion-strength gap, not an absence of coverage, hence low severity.

**Fix.** Add focused tests: feed `.boldzip_lagged_carrier_bank` a known signal and assert the lagged columns equal the hand-shifted carriers; feed `.boldzip_fit_sparse_texture` a synthetic detail signal that is an exact sparse combination of the carrier bank and assert it recovers the known sparse support and coefficients.

### 5. Retraction: codec_boldzip diagnostics and decompose are NOT missing-coverage (methodological note)

`![low]` &nbsp; **`tests/testthat/test-boldzip-sr-contracts.R:182-495`** — *informational, no action*

The `codec_boldzip_*.R` files have no same-named test file, which a pure filename cross-reference flags as untested. Verifying by *function name* shows full coverage: every exported diagnostics function is referenced (`evaluate_boldzip_sr` 2 files, `compare_boldzip_sr` 1, `boldzip_sr_payload_summary` 2, `boldzip_parcel_reconstruct` 1, `boldzip_svd_reconstruct` 1, `boldzip_reliability` 3, `boldzip_events` 2, `boldzip_sr_simulate` 3), and the decompose pipeline is exercised end-to-end with tight tolerances. **Lesson:** this package names tests by feature/symbol, not by source filename, so coverage gaps must be confirmed by function-name grep, not basename matching. All findings above were re-verified at the function-name level. No code change required.

## Systemic patterns

The four actionable findings (1-4) share one root cause: **assertion depth lags code depth in shared, internal infrastructure.** Each implicated file is either a numerical core (transport solvers), a cross-cutting gate (encoder validators), a reuse hub (encoder utilities), or an intermediate-output producer (BOLDZip texture fit), and each is currently validated only through a downstream end-to-end test with a tolerance. The recurring risk is the same: a regression that perturbs an internal relationship — path equivalence, input rejection, column accounting, lag/sparsity selection — but stays inside the end-to-end tolerance band escapes detection. The remediation pattern is uniform across all four: add direct, path-isolated unit tests (metamorphic equivalence, exact condition-class assertions, shape/support recovery) on small fixed inputs, rather than relying on the round-trip tolerance to catch internal drift.

Finding 5 documents a separate, audit-methodology pattern: because tests are named by symbol rather than by source file, coverage assessment in this repo must use function-name grep, never basename matching.
