# Correctness

## Summary

The package's correctness profile is bimodal: the high-level S4 data model, the
encode dispatch pipeline, and the round-trip transform contracts are largely
sound, but several **leaf computational kernels carry silent, data-dependent
defects** that pass the existing test suite because the tests exercise only
clean, regular inputs (full rectangular masks, even/power-of-two dimensions,
default parameters). The most dangerous findings — a Morton-vs-sequential block
ordering mismatch in the Haar wavelet, a non-symmetric normalized Laplacian, an
off-by-one B-spline knot formula, and a non-square-root eigen fallback — all
produce *plausible-looking but numerically wrong* output on realistic
neuroimaging data while remaining invisible on toy fixtures. A second large class
of findings is defensive: missing positivity/finiteness/`NA` guards that defer
errors into C++ or `RSpectra`/`Matrix` internals, yielding cryptic failures
rather than classed conditions. Round-trip invertibility is repeatedly used as
the de-facto correctness check, which is misleading for transforms whose
forward/inverse share the same flawed permutation or normalization.

---

## Critical

### C1. Haar wavelet block ordering: scalings are z-inner, data is Morton — wrong normalization on real masks
**`critical`** · `R/haar_wavelet.R:210-241`, `src/haar_wavelet_rcpp.cpp:179-202,225,289`

`precompute_haar_scalings` (R and C++) enumerates 2×2×2 blocks in
x-outer/y-middle/z-inner order, but `forward_lift_*`/`inverse_lift_*` consume
voxels in **Morton (Z-order)** block order. The two orderings agree only at block
(0,0,0). When blocks hold equal voxel counts (full mask, even dims) all
normalization factors are identical and the mismatch is silent — which is exactly
what the round-trip tests use. On any real brain mask, boundary blocks hold
different valid-voxel counts, so `sqrt_nvalid[b]` is applied to the wrong
Morton block, corrupting averages, detail coefficients, and invertibility. The
C++ path explicitly discards the Morton mask (`(void)mask_flat_morton;`),
confirming the reconciliation step was never implemented.

**Why it matters:** every spatial Haar transform on real data is wrong; the
forward/inverse pair shares the flaw, so round-trip tests cannot detect it.

**Fix:** enumerate blocks in Morton order in `precompute_haar_scalings*` (build
block coords, sort by Morton code, then fill `sqrt_nvalid`/`sqrt_nvalid_div_8`),
or carry a z-inner→Morton permutation and apply it when consuming `counts`. Apply
identically to R and C++ paths. (Root cause shared with the spatial-semantics
finding at `haar_wavelet_rcpp.cpp:179-202` — they are the same bug viewed at two
levels.)

### C2. Normalized Laplacian scales rows twice instead of rows-then-columns
**`critical`** · `R/hierarchical_helpers.R:398`

`L <- diag(n) - (D_inv_sqrt * W) * D_inv_sqrt` computes element `[i,j]` as
`d[i]·W[i,j]·d[i]` (vector recycling scales rows in both products) instead of the
symmetric `d[i]·W[i,j]·d[j]`. The result is non-symmetric (empirical
`max|L - t(L)| ≈ 0.5`), so `RSpectra::eigs` returns complex eigenvectors,
`stats::dist()` silently drops the imaginary parts, and Ward clustering runs on
corrupted coordinates — producing wrong hierarchical parcellations with no error.

**Fix:** `L <- diag(n) - outer(D_inv_sqrt, D_inv_sqrt) * W` (or the two-step
`sweep`). See also C9: the downstream `eigs()` should be `eigs_sym()`.

### C3. B-spline interior-knot formula has inverted intercept correction → k±1 columns
**`critical`** · `R/bspline_basis.R:47` (with `:54`, `:64-75`)

`n_int <- max(0L, k - degree - (if (include_intercept) 0L else 1L))` subtracts on
the wrong side. With the default `include_intercept=FALSE`, `bs()` returns `k-1`
columns; the silent zero-pad branch then appends a zero column and (with the
default `orthonormalize=TRUE`) `qr.Q()` fabricates an arbitrary non-spline
direction in the last column. With `include_intercept=TRUE`, `bs()` returns `k+1`
and a real spline column is silently trimmed. Empirically reproduced for k=6,
degree=3.

**Why it matters:** the default code path always returns a corrupted basis whose
last component is not a B-spline at all.

**Fix:** `n_int <- max(0L, k - degree - as.integer(include_intercept))`, set
`df_use <- k` unconditionally (finding M-bspline-df), and replace the silent
zero-pad with a classed abort so the QR-fabrication path is unreachable.

### C4. DCT/Haar/B-spline kernels: `k <= 0` / `n_time <= 0` not guarded → `seq.int(0, k-1)` returns wrong-length vector
**`critical`** · `R/dct_basis.R:15-23` (also `n_time=0` at `:15-22`; B-spline guards `R/bspline_basis.R:33-39`)

`build_dct_basis` guards only `k > n_time`. For `k=0`, `seq.int(0L, -1L)` is
`c(0L, -1L)` (length 2, not `integer(0)`), so `outer()` produces an 8×2 matrix of
negative-frequency cosines and the DC-normalization runs on garbage — no error.
The same `seq.int` countdown pathology fires for `n_time=0`. (The negative-`k`
and `n_time<1` variants reaching the same path are lower-severity duplicates.)

**Fix:** add `if (k < 1L) abort(class="fmrilatent_error_invalid_count")` and a
matching `n_time < 1L` guard before any `seq.int`. Route both exported builders
through the existing `.validate_positive_count()` helper so direct callers are
protected too.

### C5. `latent_dct_heatwavelet` builds `NeuroSpace` without spacing/origin → fails on every real fMRI mask
**`critical`** · `R/latent_dct_heatwavelet.R:91`

`spc <- NeuroSpace(c(dim(mask_arr), n_time))` produces unit spacing / zero
origin, but `mask_vol` retains the real `LogicalNeuroVol` geometry. The
`LatentNeuroVec` constructor compares the two spaces and hard-errors. Every other
encoder uses `.space_with_time_from_mask()`; this function is the sole outlier.

**Fix:** `spc <- .space_with_time_from_mask(mask, n_time, "latent_dct_heatwavelet")`.

---

## High

### H1. AWPT last-scale filter emits low-pass smoother → filter bank is not a partition of unity
**`high`** · `R/awpt.R:284-288`

The final scale uses `cur` (`exp(-s_last·λ)`) instead of `prev` (the low-pass
residual), so filters sum to `1 - exp(-s_{S-1}λ) + exp(-s_Sλ)` rather than 1.
Low-frequency atom energies are systematically underweighted (sum as low as
~0.78), distorting reconstruction and roughness penalties.

**Fix:** use `prev` for the last scale; verified to restore exact partition of unity.

### H2. AWPT scale weights cycled across columns instead of per-scale block
**`high`** · `R/awpt.R:694`

Loadings are scale-major (`n_centers` columns per scale), but
`Diagonal(x = rep_len(sw, ncol))` cycles weights `[w1,w2,…,w1,w2,…]` instead of
`rep(sw, each = n_centers)`. The roughness penalty is wrong whenever there are
multiple scales and centers.

**Fix:** `Diagonal(x = rep(sw, each = ncol(loadings) %/% length(sw)))`.

### H3. Eigen fallback in `.symmetric_matrix_factor` computes the wrong matrix square root
**`high`** · `R/encode_template.R:94-96`

When Cholesky fails, the fallback returns `V·diag(sqrt(D))`, which satisfies
`tᵀ·factor = D`, not `mat`. The Cholesky path returns `U` with `Uᵀ U = mat`; all
downstream callers (`to_raw`, `.analysis_loadings_from_transform`,
`.transform_quadratic_form`) rely on that convention, so the PSD fallback path
corrupts loadings, transforms, and quadratic forms.

**Fix:** return `V·diag(sqrt(vals))·Vᵀ` (symmetric sqrt) and use a pseudoinverse
for the rank-deficient inversions that follow.

### H4. `concat()` silently discards offsets from all objects except the first
**`high`** · `R/latent_methods.R:206-321`

`concat()` checks dims/masks/k/loadings but never offsets, and always stores
`x@offset`. Concatenating PCA-encoded scans (whose `center=TRUE` default produces
non-empty per-voxel offsets) from different sessions silently shifts every
reconstructed time point from the later objects by `(x@offset - y@offset)`.

**Fix:** add an offset-equality check; fall back to `NeuroVecSeq` on mismatch.

### H5. `project_hierarchical` returns OLS coefficients incompatible with stored analysis loadings
**`high`** · `R/hierarchical_template.R:234-238`

`encode_hierarchical` stores analysis-space coefficients `X·B·R⁻¹` and loadings
`B·R⁻¹`; `project_hierarchical` returns raw OLS `X·B·G⁻¹`. Both reconstruct via
`coeff %*% t(B)`, but reconstructing a `project_hierarchical` result against the
stored `analysis_loadings` gives the wrong answer — a public API contract break.

**Fix:** route `project_hierarchical` through `.template_projection_payload`
(center=FALSE) so both paths share the analysis space, or document the raw-OLS
convention and expose a matching field on `encode_hierarchical`.

### H6. `gram_factor` slot `"ANY"` stores incompatible factorizations; QR fallback breaks `Matrix::solve()`
**`high`** · `R/all_class.R:234` (consumer `R/hierarchical_template.R:115-122,235`)

`build_hierarchical_template` stores either a `CHMfactor` (Cholesky) or a base-R
`qr` object in the untyped `gram_factor` slot, but `project_hierarchical` calls
`Matrix::solve(G_factor, proj)` unconditionally. `Matrix::solve` cannot dispatch
on a base-R `qr` object, so the QR fallback (reachable on near-singular Gram)
errors or produces wrong results.

**Fix:** constrain the slot via a `setClassUnion`/`setValidity`, and branch the
consumer on factorization type (`qr.coef` for QR).

### H7. `is_explicit_latent` decoder validation omits `levels_keep` → every prediction crashes for compliant decoders
**`high`** · `R/implicit_latent.R:148-157`

The constructor checks only `time_idx`/`roi_mask`, but `predict.ImplicitLatent`
(line 218) always passes `levels_keep`. A decoder declared
`(time_idx, roi_mask)` without `...` passes construction then fails at every call
with `unused argument (levels_keep)`.

**Fix:** add `levels_keep` to the checked set, or call the decoder defensively
via `do.call`.

### H8. `neuroarchive_handoff_contract` stores manifest lists instead of template objects
**`high`** · `R/shared_structure.R:691`

`templates = lapply(templates, validate_template_protocol)` stores the validator's
plain manifest list, not the original templates. Downstream consumers calling
`template_loadings()`/`basis_decoder()` on `contract$templates` fail at dispatch.

**Fix:** `templates = templates` (validation already happens at line 665); store
manifests under a separate key if needed.

### H9. `make_cluster_reduction`: `NA` cluster ids silently dropped from `cluster_ids` but kept in `map`
**`high`** · `R/reduction.R:126-127` (also `as_cluster_reduction:151-152`)

`sort(unique(map))` drops `NA` (default `na.last=NA`), so `cluster_ids` never
contains `NA` while `map` does. Downstream `lift()` (slepian/pca/diffusion) skip
`NA`-mapped voxels via `which(map == cid)`/`match`, silently excluding them from
the basis or corrupting the prolongation matrix.

**Fix:** `if (anyNA(map)) abort(...)` in both constructors.

### H10. Inverted `orthonormalize` condition in `spec_st` bspline → always force-orthonormalizes
**`high`** · `R/encode_methods_st.R:31-34`

`build_bspline_basis` already orthonormalizes internally; the post-hoc guard
`if (!spec$time$orthonormalize)` fires exactly when the user requested
`FALSE`, overriding their choice. `spec_time_bspline(orthonormalize=FALSE)` inside
`spec_st` is silently a no-op.

**Fix:** remove the redundant QR block entirely.

### H11. `wavelet_active` `t(apply())` transposes coeffs when exactly one voxel is masked
**`high`** · `R/wavelet_active.R:29-31,42,48-52`

For a 1-voxel mask, `apply(X_mat, 1, f)` returns a length-`n_time` vector, so
`t(...)` yields `1×n_time` instead of `n_time×1`. `meta$n_time` is then stored as
1 and the decoder produces a fully transposed reconstruction. (Empty-mask case at
`:42` is the same `t(apply())` family.)

**Fix:** `do.call(rbind, lapply(seq_len(nrow(X_mat)), …))` in both encoder and
decoder; capture `n_time <- nrow(X_mat)` before the transform.

### H12. `plot_spatial_atom` ggplot2 branch never displays the plot
**`high`** · `R/spatial_plot.R:29-35`

The ggplot object is a bare statement; the function returns `invisible(arr)` and
the plot is silently discarded (auto-print only happens at top level). The
near-universal ggplot2 path shows nothing.

**Fix:** `print(p)` before `invisible(arr)`, or return the ggplot.

### H13. `k_neighbors` not clamped before `rgsp::graph_knn` → opaque crash on small parcels
**`high`** · `R/graph_bridge.R:21` (caller `R/hierarchical_template.R:412`; sibling `R/heat_wavelet.R:45-50`)

`graph_knn` requires `k < nrow(coords)`. A 2-voxel parcel with default
`k_neighbors=6` errors inside rgsp; the hierarchical caller computes `k_use` but
forwards the unclamped value. The same missing clamp exists in the heat-wavelet
lift path.

**Fix:** `k_use <- min(k_neighbors, nrow(coords) - 1L)`; abort if `< 1`. Clamp
per-cluster in `heat_wavelet` too.

### H14. `inverse_lift_rcpp` scalar `root_coeff` truncates the multi-block root vector
**`high`** · `src/haar_wavelet_rcpp.cpp:289` (caller `R/haar_wavelet.R:408`)

C++ signature `double root_coeff` silently coerces the caller's `root_coeff[tt,]`
vector to element [1]. When the decomposition retains >1 top-level lowpass block
(non-power-of-two / partial masks), all but the first are dropped, giving wrong
output or a hard stop. The R reference handles the full vector.

**Fix:** change the parameter to `NumericVector`, seed `current` from the full
vector, add a non-empty check.

### H15. Validators accept integer-overflowing counts → `NA_integer_` reaches C++
**`high`** · `R/encoder_validators.R:12-18,22-28`

`.validate_positive_count`/`.validate_nonnegative_count` pass `is.finite`,
`< 1`, and `all.equal(x, round(x))` for `x >= 2^31`, then `as.integer()` returns
`NA` with only a warning. The `NA` dimension propagates into DPSS/Haar/HRBF C++
entry points.

**Fix:** add `if (x > .Machine$integer.max) abort(...)` before coercion in both
functions. (The relative-tolerance variant — large non-integers accepted by
`all.equal` — is the same validator; replace the integer test with `x != floor(x)`.)

### H16. Heat-wavelet NULL-reduction default creates one singleton cluster per voxel
**`high`** · `R/heat_wavelet.R:78-81`

`reduction=NULL` builds `map <- seq_len(n_vox)`, so every voxel is its own
cluster, every singleton hits the degenerate path, and the loadings matrix becomes
`n_vox × (n_vox · n_scales)` of all-ones — catastrophic memory, useless basis.
This is the default path.

**Fix:** `map <- rep(1L, n_vox)` for a single global cluster, or remove the NULL
default and require a meaningful reduction.

### H17. `eigs()` used instead of `eigs_sym()` for symmetric Laplacian
**`high`** · `R/hierarchical_helpers.R:406`

After fixing C2, `L` is symmetric PSD, but `RSpectra::eigs` (general solver) may
still return complex-typed results that `stats::dist()` silently realizes.

**Fix:** `RSpectra::eigs_sym(L, k = k_use, which = "SM")`.

### H18. `gather()` unconditionally discards Poisson samples for small components
**`high`** · `R/hrbf.R:235-241`

In `poisson_disk_sample_neuroim2`, the closure overwrites the Poisson-disk result
with a single centroid whenever `nrow(coords) <= 8`, regardless of whether
sampling succeeded. The compat sibling has the correct `nrow(out) == 0` guard.
Every ≤8-voxel component always yields exactly one atom.

**Fix:** `if (nrow(out) == 0L && nrow(coords) <= 8L) out <- centroid`.

### H19. `benchmark_roundtrip` unit conversion off by 10⁹ (seconds treated as nanoseconds)
**`high`** · `R/benchmark_roundtrip.R:66`

`bench::mark` medians are in **seconds**; dividing by `1e6` yields nanoseconds.
The `median_ms` column is ~9 orders of magnitude wrong.

**Fix:** `median_ms = as.numeric(mark_res$median) * 1e3`.

---

## Medium

> The medium tier is dominated by three recurring root causes (see Systemic
> Patterns). Grouped by theme below.

### M-noise. `noise_scale` computed from the wrong (untrimmed / pre-projection) value
**`medium`** · `R/codec_boldzip.R:151-156`, `R/codec_boldzip_decompose.R:387-392`, `R/codec_boldzip_diagnostics.R:137-144`

Three BOLDZip quantization sites derive the quantization step from a distribution
other than the values being quantized: `theta` quantized with `sd(z_raw)` (the
pre-projection carrier); event amplitudes quantized with the full pre-trim
`amplitude` vector instead of `amplitude[ord]`; and per-timepoint reliability
weights recycled column-major against a voxel-major error matrix. All miscalibrate
the step or weighting; impact gated on quantization being enabled.

**Fix:** pass `.boldzip_noise_scale(<the trimmed values being quantized>)` at each
site; validate `reliability_weights` length against `nrow(X)`.

### M-symmetry. AWPT symmetry/normalization guards are vacuous or skipped
**`medium`** · `R/awpt.R:123-130` (vacuous post-symmetrize check), `R/awpt.R:677` (rw-Laplacian → asymmetric Q)

`.awpt_enforce_symmetric` checks symmetry of `0.5(M+Mᵀ)`, which is always
symmetric, so severely asymmetric inputs pass silently. Separately, passing an
`rw`-normalized (`D⁻¹L`, non-symmetric) operator as `anatomical_operator`
produces a non-symmetric roughness `Q = BᵀL_rw B` that is stored unsymmetrized.

**Fix:** check the *input* against its transpose before symmetrizing; symmetrize
`Q` after computing it (or reject non-symmetric operators).

### M-hrbf. HRBF orientation/index/projection inconsistencies
**`medium`** · `R/hrbf.R:703-705` (square-matrix transpose), `:797-948` (sparseMatrix col dim = active not full), `:62-74` (`hrbf_reconstruct_partial` index space), `:379-383/365` (silent `wendland_c4→c6` alias), `:215-233` (inert seed)

A cluster of independent HRBF defects: (1) the orientation guard
`nrow==n_active && ncol==k_actual` transposes an already-correct K×N matrix when
`n_active == k_actual`; (2) the R-fallback `sparseMatrix` uses active-voxel count
as the column dimension while indexing with raw grid indices (errors on
non-contiguous masks; Rcpp fast path unaffected); (3) `hrbf_reconstruct_partial`
validates `voxel_idx` against the full grid but indexes active-voxel columns; (4)
`kernel_type="wendland_c4"` is silently coerced to C6 (different smoothness); (5)
`sample_component` calls `set.seed` then iterates unshuffled coords, so the seed
is inert and centre selection is scan-order biased.

**Fix:** enforce the C++ K×N contract and drop the shape heuristic; remap indices
to a single convention; implement the C4 polynomial and remove the alias; shuffle
coords after `set.seed` (matching the compat variant).

### M-validate. Missing positivity/finiteness guards defer errors into internals
**`medium`** · `R/diffusion_wavelet.R:110,129` (`max_scales=0`→`Matrix(NULL)`), `R/pca_spatial.R:90` (`k<=0`→empty loadings), `R/hierarchical_helpers.R:405-406` (`n<3`→`eigs` crash), `R/slepian_handles.R:30-35` (unvalidated tr/bandwidth → bad cached handle), `R/latent_dct_heatwavelet.R:31,85,91` (unvalidated `n_time`), `R/slepian_temporal.R:59-77` (unvalidated tr/bandwidth → unclassed `if(NA)` error), `R/hierarchical_template.R:406-426` & `R/slepian_spatial.R:21-34` (`k==n_loc` → undocumented RSpectra fallback / crash; missing `requireNamespace` guard at `:34`)

A broad family: scalar parameters reach `RSpectra::eigs`, `Matrix::Matrix`, or
C++ without lower-bound/`NA`/finiteness checks, producing cryptic
internal errors, silent degenerate output (empty loadings, all-zero latent), or
polluted handle caches instead of classed `fmrilatent_*` conditions.

**Fix:** validate scalars at the function boundary (positive, finite, scalar,
non-`NA`); cap `k_use <- min(k, n_loc - 1L)` before RSpectra and route the
all-eigenvectors case to `base::eigen`; add the missing `requireNamespace`
abort.

### M-codepath. Validated/accepted arguments silently ignored
**`medium`** · `R/encode_methods_space.R:173-186` (`materialize` ignored in wavelet_active), `R/slepian_temporal.R:53-74` (`denoise` never read), `R/encode_methods_time.R:23-29` (`meta$k` unclamped vs actual rank), `R/transport_latent.R:406-424` (adjoint-convention warning suppressed)

Several public arguments are resolved/documented but never consulted:
`materialize` is validated then dropped; `denoise` is a no-op so a large `k` is
never capped to the Shannon number; `meta$k` records the unclamped value while the
stored basis has `min(k, n_time)` columns; and the composed adjoint convention
can become non-euclidean without a warning (error surfaces only later at
covariance time).

**Fix:** honour or remove each argument; clamp `k_use` before recording `meta$k`;
warn whenever the composed convention is non-euclidean.

### M-roi. ROI/coordinate handling reimplemented or mis-shaped
**`medium`** · `R/slepian_spatiotemporal.R:52-57` & `R/encode_methods_st.R:84-87` (inline ROI subsetting bypasses `roi_subset_columns`), `R/implicit_latent.R:321-325` (`reconstruct_array` drops spacing/origin), `R/latent_block_vector.R:179-184` (non-list `roi_mask` silently dropped)

ROI subsetting is duplicated inline (mishandling mask-indexed vector ROIs that the
canonical helper would reject) in two files; `reconstruct_array` converts the
`LogicalNeuroVol` mask to a plain array before `.wrap_decoded_volume`, so the
returned `DenseNeuroVec` has default spacing/zero origin; and a non-list
`roi_mask` is silently replaced with an empty list, skipping all spatial
subsetting.

**Fix:** delegate to `roi_subset_columns`; pass the `LogicalNeuroVol` (or
`.space_with_time_from_mask`) through; error or broadcast non-list `roi_mask`.

### M-diffusion. Diffusion-wavelet operator and cache defects
**`medium`** · `R/diffusion_wavelet.R:63-65,122-125` (compressed operator is `QᵀTᵀTQ` not `QᵀTQ`), `:40-43,97-103` (function `T_mat` accepted then always aborted), `:78-79` (`seed=NULL` replaced by `1L`), `R/diffusion_wavelet_handle.R:26-35` (`data` hashed into id but ignored by `lift()`)

The multi-scale compressed operator uses `crossprod(B)` (Gram of T) instead of the
Galerkin restriction `crossprod(Q, T_op(Q))`, breaking hierarchical scales
(`max_scales > 1` only; default 1 masks it). The function-operator branch in
`as_linear_operator` is dead because the guard checks `dim(T_mat)` not `T_op`.
`seed=NULL` cannot opt out of seeding. `data` is hashed into the cache id but
never consumed, causing spurious cache misses/bloat.

**Fix:** Galerkin restriction; check `T_op`/pass an explicit `n`; drop the
`%||% 1L`; exclude `data` from the id hash.

### M-quad. AWPT FISTA objective double-counts the Gram term
**`medium`** · `R/encode_transport_solve.R:58-70`

`.awpt_objective` adds `0.5·tr(ZᵀZ·A)` with `A = gram + λQ` on top of the residual
term, which already contains `0.5·tr(ZᵀZ·gram)`. The gradient is correct (solution
converges), but the stopping criterion and divergence guard evaluate an inflated
objective at the wrong scale.

**Fix:** pass only the penalty part (`λQ`) to the objective.

### M-dead. Misordered/dead validation guards yield base-R errors instead of classed conditions
**`medium`** · `R/encode_template.R:79-83,100-104` (square check after symmetrize), `R/latent_methods.R:94-99` (`if(NA)` on i,j,k), `R/all_generic.R:233` (`coef_metric` default `"raw"` vs peers' `"analysis"`), `R/latent_neurovec_materialize.R:100,126` (`Matrix()` without `sparse=`), `R/latent_indexing.R:405-408` (integer overflow in 4-col linear index)

Several guards are unreachable (square check runs after the non-conformable
symmetrization that throws first) or produce base-R errors (`if(NA)` on missing
i,j,k coords). `coef_metric` defaults to `"raw"` while all peer generics default
to `"analysis"` (the documented Euclidean contract). `Matrix()` without `sparse=`
yields a data-dependent class. The 4-column matrix index path overflows 32-bit
integers for large volumes.

**Fix:** move dimension checks before symmetrization; guard `NA` coords; align the
`coef_metric` default; specify `sparse=`; cast to double before the index
multiplication.

### M-misc. Other medium correctness gaps
**`medium`** · `R/haar_wavelet.R:560-580` (`perform_haar_lift_synthesis` hardcodes `root[tt,1]`), `:699-700` (`relative_to_root_std` mixes root + coarsest detail), `R/encode_operator.R:311` (eager `template_roughness` regardless of `spatial_lambda`), `R/reduction.R:39-45` (`CoarsenedReduction` no constructor/validity), `R/searchlight_utils.R:20-23,52-53` (no bounds/k-match validation), `R/encoder_utilities.R:11-20` (swallowed mask-conversion error), `R/encode_methods_space.R:74` (HRBF encode uses LS while `hrbf_latent` uses dot-product projection — divergent objects), `R/latent_bilateral_surface_vector.R:140-155` (NeuroSurface vs NeuroSurfaceVector round-trip)

Assorted single-site defects: the Haar synthesis side drops multi-block roots
(analysis side guards `top_blocks==1`, synthesis does not); the
`relative_to_root_std` threshold conflates the root with the coarsest detail band;
`encode_awpt` eagerly materializes the roughness matrix even at `spatial_lambda=0`;
`CoarsenedReduction` has no dimension-consistency enforcement; searchlight helpers
don't validate neighborhood indices or `ncol(B)==ncol(L)`; mask-conversion errors
are swallowed, hiding the root cause; and the two HRBF encode entry points use
different projection operators, producing different `LatentNeuroVec` objects for
identical inputs.

---

## Low

The low tier (≈40 findings) is overwhelmingly **dead code, redundant guards, and
diagnostic-quality issues** with no current runtime impact:

- **Dead/vacuous conditions:** `is.null()` on environment/list slots
  (`R/reduction.R:143-148`); `is.matrix(x) && !is(x,"Matrix")`
  (`R/latent_surface_vector.R:73-80`, `R/latent_indexing.R:57-66`); always-true
  Matrix density guards; `!is.numeric(i)` on integer-typed dispatch
  (`R/latent_indexing.R:158`); `drop` in `...` shadowing a formal
  (`R/latent_methods.R:33-35`); `simplify` guard structurally always TRUE
  (`R/searchlight_utils.R:25`); `isGeneric(..., where=topenv(parent.frame()))`
  resolving to base (`R/all_generic.R:74-76`).
- **Condition-object hygiene:** duplicate class strings in `.encoder_cli_warn`
  (`R/encode.R:3-10`, `R/encoder_registry.R:65-68`); `rlang::caller_env()`
  inside a tryCatch handler (`R/latent_handles.R:469-476`).
- **Missing diagnostic/early-validation:** unguarded `k`/`n_time`/`dims`/`seed`
  in plotting, registry, and test helpers; misleading validator reuse
  (`.boldzip_validate_matrix` on adjacency); opaque `getClass()` dependency probe
  (`R/latent_bilateral_surface_vector.R:130-131`); missing `@include all_class.R`
  in `all_generic.R`.
- **Defensive C++ gaps (require malformed direct `.Call`):** unchecked length
  validation in `forward_lift_rcpp`/`inverse_lift_rcpp`
  (`src/haar_wavelet_rcpp.cpp`); duplicate-coordinate and `levels<=0` acceptance
  in `active_pencil_wavelet.cpp`; unbounded `n×n` DPSS allocation
  (`src/slepian_dpss_rcpp.cpp`); exported-but-blocked dense prolate solver;
  unvalidated `kernel_type` falling through to Wendland.
- **Edge-case-only:** `t(apply())`/empty-mask transposition for degenerate
  inputs; `as.matrix()` on a `LoadingsHandle` (latent, no caller triggers it);
  `gram_factor` computed but unread in surface/parcel templates;
  `pair_time` asymmetry for odd `n_time` halves; `sd()` of length-1 inputs.

Each is documented in the raw findings; none blocks correct operation on
well-formed inputs but several mask future regressions (e.g. the zero-pad+QR and
`t(apply())` dead branches).

---

## Systemic Patterns

1. **Round-trip ≠ correctness.** The Haar wavelet (C1, the C++ ordering twin, H14,
   M-misc synthesis), DCT, and active-pencil transforms all pass round-trip tests
   because forward and inverse share the same flawed permutation/normalization.
   Block ordering, spatial localization, and multi-root handling are wrong while
   reconstruction is exact. **Recommendation:** add tests that assert
   *intermediate* invariants (partition of unity, per-block voxel counts on
   irregular masks, coefficient values against a dense reference), not just
   `decode(encode(x)) ≈ x`.

2. **`noise_scale` calibrated from the wrong distribution** (M-noise, three sites)
   — a copy-paste of the explicit-override pattern that defeats the correct
   `%||% .boldzip_noise_scale(values)` default. **Recommendation:** remove the
   explicit overrides and rely on the helper's default; add an assertion that the
   scale source has the same length as the quantized values.

3. **Vector-recycling and `apply()`-shape bugs** (C2 Laplacian, H2 scale weights,
   H11/M wavelet transpose) — all stem from R's silent recycling/dimension-drop
   semantics. **Recommendation:** prefer `outer()`/`rep(…, each=)`/`vapply`/
   `do.call(rbind, …)` over `vec * matrix` and `t(apply())`.

4. **Validation deferred into internals** (H13, H15, the entire M-validate group,
   most of Low) — scalar guards live in some functions (`dpss_time_basis`,
   `make_cluster_reduction`) but are missing in their siblings, so bad inputs
   surface as cryptic `RSpectra`/`Matrix`/C++ errors. **Recommendation:**
   centralize boundary validation in the `.validate_*` helpers and call them
   uniformly; promote the geometry helper `.space_with_time_from_mask` and the
   ROI helper `roi_subset_columns` to the single source of truth (C5, M-roi).

5. **`gram_factor` / metric-factorization convention drift** (H3, H5, H6, M-quad,
   the surface/parcel dead-factor low findings) — multiple code paths recompute or
   store Gram factorizations under inconsistent conventions (ridged vs unridged,
   Cholesky vs QR vs eigen, weighted vs unweighted). **Recommendation:** define
   one factorization contract (`factorᵀ·factor = metric`), enforce it via slot
   typing/validity, and have every projector consume the same stored factor.
