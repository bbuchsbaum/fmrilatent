# Performance

The package is functionally correct on the performance axis — none of the verified findings is a correctness defect — but it carries a consistent set of memory-and-compute inefficiencies that surface specifically at neuroimaging scale (10^5 voxels, 10^2–10^3 components, 10^2–10^3 time points). The dominant theme is premature dense materialization: sparse and lazy (handle-backed) representations are repeatedly forced to dense base matrices for validation, identity comparison, reconstruction, and arithmetic, defeating the package's own lazy/sparse storage design. The second theme is wasted computation around the handle/registry cache: every spatial-wavelet handle constructor runs its expensive `lift()` before consulting the cache that would have made it unnecessary. Individually most findings are medium or low severity; collectively they mean that the lightweight-handle and sparse-storage abstractions do not deliver their intended savings on large brains, and a few paths can drive peak memory into the multi-GB range or trigger OOM.

## Findings (by severity)

### 1. Spatial-wavelet handle constructors run `lift()` before checking the cache
**`high`** · `R/diffusion_wavelet_handle.R:37-48`, `R/heat_wavelet_handle.R:38-48`, `R/slepian_handles.R:100-103`

*Merged: three constructors share an identical root cause.* Each handle constructor computes a deterministic `id` (available before the lift), then calls the expensive `lift()` — a full voxel×components dense materialization, including kNN-graph build, Chebyshev/diffusion expansion, or graph-Laplacian eigendecomposition — and only afterward calls `.latent_register_matrix(..., overwrite = FALSE)`. On a fingerprint-matching cache hit that registration call silently discards the freshly computed result (`latent_handles.R:230-232`), so the entire `lift()` was wasted. The `loadings_mat()` access path correctly checks the cache first; these constructors do not.

**Why it matters:** For large brains (50k+ voxels, hundreds of components) `lift()` is the dominant cost — seconds of wall time, and a full dense allocation. Any pattern that re-constructs a handle for an already-cached spec (loops, multiple `LatentNeuroVec` instances sharing a dictionary) silently recomputes everything.

**Fix:** Consult `.latent_get_matrix(id, type = "loadings")` before calling `lift()`. On a hit, reconstruct the handle from `dim(cached)` and return immediately:
```r
cached <- .latent_get_matrix(id, type = "loadings")
if (!is.null(cached)) {
  return(new("LoadingsHandle", id = id, dim = as.integer(dim(cached)),
             kind = "lifted", spec = spec_payload, label = label))
}
L <- lift(reduction, basis_spec, data = data, k_neighbors = k_neighbors)
# ... construct + register as before
```
The diffusion/heat constructors are rated highest because the wasted computation is largest; slepian is the same defect with smaller per-call cost.

The `diffusion_wavelet_handle.R:37` instance is rated `high` on its own; the heat and slepian variants are `medium` individually but share the identical remedy.

---

### 2. Premature dense materialization of sparse / lazy matrices
**`medium`** · multiple files (themed cluster — see table)

A family of findings where a sparse `Matrix` or handle-backed loadings/basis is forced to a full dense `matrix` for an operation that does not require it. Each is a real allocation of `O(n_vox × k)` (or larger) that is immediately discarded or could be done sparse-aware.

| File:line | What's wrong | Fix |
|---|---|---|
| `R/latent_neurovector.R:213-220` | `all(is.finite(basis))` / `all(is.finite(loadings))` triggers `Matrix`'s `is.finite,sparseMatrix` method, which allocates a dense logical vector of `prod(dim)` (~100 MB for 50k×500) on **every** `LatentNeuroVec()` construction, before the density-based conversion. | `.all_finite()` helper checking only `m@x` for `sparseMatrix` subclasses (structural zeros are always finite). |
| `R/shared_structure.R:38-43` | `any(!is.finite(as.matrix(x)))` densifies the whole matrix in a validator called from several contract paths; no short-circuit. | Sparse-aware: check `x@x` for `sparseMatrix`. |
| `R/encoder_utilities.R:104` | `as.matrix(X) %*% as.matrix(loadings)` densifies the (often sparse) loadings dictionary before the product. | Drop `as.matrix(loadings)`; let `Matrix` dispatch handle `dense %*% sparse`. |
| `R/latent_methods.R:268-278` | `concat()` does `identical(as.matrix(x_loadings), as.matrix(loadings_mat(obj)))` per object — `O(N × 160 MB)` for whole-brain handles, just for an identity test. | Compare handle `@id` + `@spec` before falling back to matrix equality. |
| `R/latent_methods.R:549-564` | `as.array()` holds the full `nTime×nVox` matrix *and* the `X·Y·Z·nTime` array simultaneously (~4.3 GB peak for MNI×300). | Compute each time slice directly from `B`/`L`; never materialize the full matrix. |
| `R/latent_surface_vector.R:128-141` | `reconstruct_matrix` calls `as.matrix(x)` for all time points then row-subsets `time_idx`. | Pass `time_idx` to `basis_mat(x, i = ...)` (as `latent_methods.R` already does). |
| `R/latent_bilateral_surface_vector.R:81-87` | `loadings()` forces `sparse = FALSE` on the rbind of two hemisphere loadings, densifying block-diagonal-sparse data. | `methods::rbind2(loadings(x@left), loadings(x@right))` to preserve sparsity. |
| `R/slepian_spatiotemporal.R:43-44,70` | `as.matrix(L_s)` densifies a block-sparse lift (potentially ~40 GB at 500k×10k) before the decoder needs it. | Keep `L_s` sparse in `coeff`; coerce inside the decoder multiply only. |
| `R/spatial_plot.R:13` | `loadings_mat(loadings)` materializes the full dictionary to use one column. (Note: a `j=idx` call site change alone does *not* avoid materialization — `materialize_loadings_from_spec` is column-blind; true fix is architectural.) | Add column-selective materialization to `materialize_loadings_from_spec`, or accept the limitation. |

**Why it matters:** These are exactly the operations the lazy/sparse design exists to avoid. At whole-brain scale several of them (notably `as.array()`, `concat()`, the slepian-ST decoder, and the `is.finite` guard on construction) can exhaust memory or add hundreds of MB of avoidable peak allocation.

---

### 3. Redundant recomputation (non-cache)
**`medium`**

| File:line | What's wrong | Fix |
|---|---|---|
| `R/hrbf.R:91-92` | `hrbf_latent` calls `hrbf_generate_basis()` directly, then `hrbf_project_matrix()` which calls `hrbf_generate_basis()` again — the Poisson-sampled atom basis is built twice. | Build the basis once, inline the `tcrossprod`. |
| `R/parcel_basis.R:144-166` | `gram_factor` (Cholesky of LᵀL) is stored but never read; `template_project` looks for `analysis_transform` which is never populated, so the Gram is re-factored `O(k³)` per subject. The "cached" doc claim is false. | Compute and store `analysis_transform` in `meta`, or pass `gram_factor` into `template_project`. |

---

### 4. Avoidable dense / `O(n²)` allocations on inactive or guarded paths
**`medium`**

| File:line | What's wrong | Fix |
|---|---|---|
| `R/encode_transport_solve.R:337-341` | When `temporal_lambda <= 0`, `Lt` is a dense `n_time×n_time` zero matrix that is never used; the matrix-free sibling correctly uses `NULL`. | `} else { NULL }`. |
| `R/transport_latent.R:776-779` | `decode_covariance(diag_only = FALSE)` materializes a dense `n_target×n_source` decoder, then `D Σ Dᵀ` is `n_target×n_target` (~80 GB at 100k voxels). No guard. Opt-in only (default `diag_only = TRUE`). | Add a size-based warning before materializing. |
| `R/heat_wavelet.R:47` | `diag(n_loc)` allocates a dense `n_loc×n_loc` identity per cluster (~32 MB at 2000 voxels). *Note:* the reviewer's `Matrix::Diagonal` fix is invalid — rgsp's C++ backend requires `arma::mat`; the identity is the intended full-dictionary extraction. | Add a `n_loc` size guard; densification is inherent to the algorithm. |

---

### 5. Quadratic / nested-loop scaling
**`medium`**

| File:line | What's wrong | Fix |
|---|---|---|
| `R/latent_indexing.R:454-463` | `[` ANY/ANY method does `which(mask_indices == ...)` inside a `t_idx × v_idx` double loop — `O(n_time × n_unique_vox × n_vox_req)` (~2.5 B comparisons for a full-brain query). The numeric `[` method solves the same scatter in `O(n_vox)`. | Precompute `col_map <- match(mask_indices, valid_mask_idx)` once; vectorized scatter per time point. |
| `R/shared_structure.R:610-635` | `render_shared_events` grows `rows`/`cols`/`values` with `c()` inside a per-event loop — classic `O(n²)` vector-growing. | Collect into pre-sized lists, `unlist()` once after the loop. |

---

### 6. Numerical reliability (spectral)
**`medium`** · `R/hierarchical_template.R:425`

`RSpectra::eigs(L, k, which = "SM")` on a graph Laplacian (always PSD with a zero eigenvalue) asks ARPACK to converge near a singularity — the hardest regime — and can fail or return inaccurate eigenvectors for elongated parcels or near-disconnected kNN graphs. The failure is detectable (error / fewer vectors), so this is a reliability rather than silent-correctness issue.

**Fix:** Shift-invert mode `RSpectra::eigs(L, k, which = "LM", sigma = 0)`, the standard method for Laplacian/Fiedler eigenvectors.

---

### 7. C++ Gaussian HRBF atom: no distance cutoff
**`medium`** · `src/hrbf_atoms_rcpp.cpp:54-55` (and duplicated 84-85)

The Wendland branch short-circuits with `if (dist >= 1.0) continue;`; the Gaussian branch evaluates `std::exp(...)` for every voxel unconditionally, gating only afterward via `value_threshold`. At brain scale (N≈200k voxels × K≈5k atoms) this is ~10^9 transcendental calls, the vast majority producing sub-threshold values.

**Fix:** Precompute `max_dist2 = -2.0 * sigma * sigma * std::log(value_threshold)` and `continue` when `dist2 > max_dist2`, matching the Wendland early-exit.

---

### 8. Validity functions materialize bases on every `validObject()`
**`low`** · `R/latent_bilateral_surface_vector.R:268-291`, `R/latent_block_vector.R:319-323` (+ `R/latent_shared_validation.R:16-45`)

*Merged: shared root cause across two S4 validity functions.* Both `setValidity` handlers run a full value-equality check (`as.matrix(basis(...))` + `all.equal`) on every `validObject()`/`new()`, materializing handle-backed or large dense bases. For `BlockLatentNeuroVector` the constructor invokes the check twice. After the constructor's first materialization the result is usually an LRU-cache hit, so practical cost is bounded, but value equality does not belong in a validity function.

**Fix:** Reduce the validity check to dimension agreement (`.latent_basis_dim()`, which never materializes); keep value-equality in the constructor only.

---

### 9. Redundant self-comparison / unconditional densify in validators
**`low`**

| File:line | What's wrong | Fix |
|---|---|---|
| `R/latent_shared_validation.R:33-43` | Loop runs `seq_len(n)`, so `i=1` re-materializes `bases[[1]]` and compares the reference to itself. | `seq.int(2L, n)`. |
| `R/surface_template.R:49-55` | Full dense coercion of `loadings` happens *before* the row-count guard. | Check `nrow` before coercing. |
| `R/surface_template.R:77` | `Matrix::Matrix(as.matrix(roughness), sparse = FALSE)` double-copies a `k×k` matrix the consumer densifies anyway. | Store as-is / drop the intermediate `as.matrix`. |
| `R/pca_spatial.R:79` | `as.matrix(data)` copies even when `data` is already a base matrix (full pass is genuinely needed for centering, so deferral helps only the non-default path). | `if (is.matrix(data)) data else as.matrix(data)`. |
| `R/latent_block_vector.R:162-168` | `loadings()` hardcodes `sparse = FALSE`, densifying sparse block loadings (consistent with the class's dense assumptions elsewhere). | Drop `as.matrix()`/`sparse=FALSE`; rbind sparse blocks. |

---

### 10. Minor compute waste
**`low`**

| File:line | What's wrong | Fix |
|---|---|---|
| `R/hierarchical_helpers.R:195` | `max(levels_lh[[lvl]])` recomputed inside the per-voxel loop (constant across iterations). Loop overhead dominates, but the redundant `max` is real. | Hoist `lh_max <- vapply(levels_lh, max, integer(1))`; ideally vectorize the loop. |
| `R/encode_transport_solve.R:207-229` | Power iteration applies the (expensive matrix-free) quadratic operator twice per step; one application suffices (`lambda_est <- norm_w`). | Drop the second application. |
| `R/codec_boldzip_decompose.R:53-57` | Orthonormality check computes a `k×k` Gram with no `ncol <= nrow` guard; overcomplete bases trigger two `O(k²)` allocations before rejection. | Guard `ncol(basis) > nrow(basis)` before `crossprod`. |
| `R/codec_boldzip_validate.R:153` | Zero-variance guard uses `.Machine$double.eps` (~2.2e-16); `cor()` is unreliable well above that, and matrix inputs are flattened by `as.numeric()`. Diagnostic-only (no data-integrity impact). | Use `sqrt(.Machine$double.eps)` (~1.5e-8). |

---

## Systemic patterns

1. **Cache-after-compute in handle constructors (Finding 1).** The three spatial-wavelet handle constructors all run `lift()` before consulting a cache that would obviate it, relying on `.latent_register_matrix(overwrite = FALSE)` to silently drop the result. The correct early-exit (`.latent_get_matrix` first) already exists on the read path and should be mirrored in every constructor. A single shared helper (e.g. `.with_loadings_cache(id, spec, fn = lift_closure)`) would eliminate the class of bug.

2. **`as.matrix()` as a reflex (Findings 2, 8, 9).** The single most common defect is calling `as.matrix()` / `Matrix(..., sparse = FALSE)` on sparse or handle-backed data for validation, identity, reconstruction, or arithmetic that `Matrix` dispatch or `@x`-slot inspection would handle without densifying. This directly undermines the package's lazy/sparse storage premise. A small set of sparse-aware helpers (`.all_finite()`, dimension-only basis comparison, time-subset-aware reconstruction) would retire most of these.

3. **Value-equality in S4 validity (Finding 8).** Validity functions should verify shape, not content; value comparison forces materialization on a path invoked implicitly throughout dispatch. Move it to constructors.

4. **Inactive-branch allocation (Finding 4).** Dense `O(n²)` placeholders are allocated on logically-inactive code paths (`Lt` zero matrix) where a sibling implementation already uses `NULL`. Aligning the two transport-solver paths removes the waste for free.
