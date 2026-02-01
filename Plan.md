# HRBF Implementation Plan (fmrilatent)

## Goal
Port the hierarchical radial basis (HRBF) representation from `neuroarchive` into `fmrilatent`, producing a fast, deterministic path to create `LatentNeuroVec` objects (basis = time × atoms coefficients, loadings = voxels × atoms) while supporting shared/implicit dictionaries via lightweight handles.

# Hierarchical Laplacian Frame Plan (template-backed, non-orthogonal)

## Goal
Ship a sparse, multiscale Laplacian frame built from nested parcellations, encoded via cached sparse Cholesky of the Gram matrix. The structure is template-dependent (atlas + mask) but data-agnostic, enabling reuse across studies.

## Milestones
1) **Template class and helpers**
   - Add S4 `HierarchicalBasisTemplate` with slots: `mask`, `space`, `parcellations` (list of int labels per level), `loadings` (sparse B), `gram_factor` (Cholesky/QR), `meta` (atlas id, version, k_per_level, build hash). 
   - Accessors: `is_hierarchical_template()`, `template_meta()`, `template_loadings()`, `template_solver()`; keep consistent with `LatentNeuroVec` helper naming.

2) **Builder pipeline (offline)**
   - `build_hierarchical_basis(mask, parcellations, k_per_level, ridge = 1e-8, solver = c("chol","qr"))`.
   - Validate nesting (parent is union of children), reindex labels per level, and assert coverage of in-mask voxels.
   - Use existing `basis_slepian` + `make_cluster_reduction` + `lift` to assemble per-level blocks; `cbind` into sparse B.
   - Compute `G = crossprod(B)`, add ridge if requested, factor with `Matrix::Cholesky(perm = TRUE, LDL = TRUE)`; fallback to sparse QR if Cholesky fails.

3) **Encoding/decoding helpers (runtime)**
   - `encode_hierarchical(X, template)` → `LatentNeuroVec` with `basis = coeff (time×k)`, `loadings = template@loadings`, `meta$family = "hierarchical_laplacian"`, caching `gram_factor` in meta for reuse.
   - `decode_hierarchical(latent)` simply multiplies `loadings %*% t(latent@basis)`.
   - `project_hierarchical(template, X)` returns coefficients only (time×k) for reuse in factories.

4) **Template assets strategy**
   - Provide packaged `.rds` templates under `inst/extdata` (e.g., `mni_hier_schaefer100_400_subcort_v1.rds`) built offline from MNI mask, Schaefer 100/400 cortex, and a subcortical atlas (aseg/HO/CIT168). Keep voxel order documented.
   - Include a reproducible build script (not run on CRAN) under `inst/scripts/build_hier_templates.R` describing atlas sources and seeds.
   - `load_hierarchical_template(name = "mni_v1")` convenience loader used by `latent_factory(family = "hierarchical_mni")`.
   - Size presets (initial):
     * **small (L0–L2)**: whole brain → ~32 networks → ~64 parcels; k_per_level ≈ c(8,3,1).
     * **medium (L0–L3)**: whole brain → ~32 → ~64 → Schaefer-400 leaves; k_per_level ≈ c(8,5,3,1).
     * **large (optional L4)**: add an ~800 tile leaf if we fabricate supervoxels; reuse k_per_level c(8,5,3,1).
   - Semantics: cortical leaves piggy-back on Schaefer-400 (names, Yeo17 tags); coarser levels built by agglomerating leaves while recording parent mapping so meaning is preserved.

5) **Subcortex handling**
   - Build cortical and subcortical blocks separately; concatenate in voxel order to keep B block-diagonal and Gram sparse. Allow per-block `k_per_level`.
   - Permit small parcels/nuclei to use k = 1–3 or identity to avoid unstable eigenmodes.

6) **Testing**
   - Tiny mask round-trip: encode→decode within tolerance; check sparsity (nnz) and Gram bandwidth.
   - Mask/space mismatch errors; degenerate parcel handling (zero-size cluster) surfaces a clear message.
   - Performance smoke: solve time scales with K and remains sparse (symbolic nnz tracked).

7) **Docs**
   - Roxygen for template builder/loader and encoding helpers; vignette section describing the frame vs orthogonal trade-off and the cached-solver workflow.
   - Document atlas/licensing assumptions and regeneration steps for templates.

## Parcellation construction (Schaefer-driven)
- Inputs: surface Schaefer L/R (fsaverage) for geodesic + boundary contact; volumetric Schaefer-400 in MNI for voxel labelmaps; subcortical atlas (aseg/HO/CIT168).
- Similarity on parcel graph (nodes = 400 cortex parcels): w_ij = α*boundary_contact + β*exp(-geodesic/d0) + γ*1{Yeo17_i==Yeo17_j}; defaults α=0.5, β=0.3, γ=0.2, d0≈30mm.
- Clustering: compute Laplacian spectral embedding (2–3 dims), then Ward linkage on the embedding; restrict merges within hemisphere+Yeo17 until target counts met (e.g., 64→32), then optionally relax Yeo17 for top level (~16).
- Levels shipped: L0 whole brain; L1 ~16–20 (Yeo17-like, hemi-split); L2 ~32–64; L3 = Schaefer-400 leaves; optional L4 = supervoxel/voxel leaf per parcel (high-res SKU).
- Subcortex: keep atlas labels as leaves; optional splits for big nuclei; low k_per_level (1–3).

## Milestones
1) **Extract & Adapt Core R Helpers**
   - Import `generate_hrbf_atom`, `poisson_disk_sample_neuroim2`, and `hrbf_basis_from_params` (trim HDF5/pipeline hooks).
   - Replace neuroarchive options with fmrilatent-local flags (e.g., `fmri.hrbf.use_rcpp = TRUE/FALSE`).
   - Ensure inputs are plain mask arrays/`LogicalNeuroVol` and return `Matrix` objects.

2) **LatentNeuroVec Wiring**
   - Add a user-facing constructor/wrapper `hrbf_latent(mask, data, params)` that returns `LatentNeuroVec`.
   - Implement encode (data → coeff) and decode (coeff → data) helpers using the HRBF basis.
   - Support kernel choice (`gaussian`, `wendland_c6`), levels, extra fine levels, radius factor, seed.
   - Store family metadata in `LatentNeuroVec@meta` with helpers (`as_hrbf_latent()`, `hrbf_meta()`, `is_hrbf_latent()`); use this pattern for future basis types.

3) **Performance Path**
   - Port `hrbf_atoms_rcpp` (Eigen/OpenMP) and expose a toggle.
   - Benchmark vs R path on small (32³) and medium (64³) masks; document thresholds to switch.

4) **Testing**
   - Round-trip tests on tiny deterministic masks (≤8³) and medium masks with set seeds.
   - Component access tests: reconstruct selected voxels/timepoints; compare encode/decode equality within tolerance.
   - Performance smoke tests (microbenchmark) to guard regressions.

5) **Docs & Examples**
   - Roxygen for public helpers; short vignette/example script showing HRBF latent construction, encoding, reconstruction, and component inspection.
   - Codify the `meta$family` tagging pattern for other basis families (`as_<family>_latent()`, `<family>_meta()`, `is_<family>_latent()`); avoid new classes unless dispatch is needed.

6) **Dictionary handles (cross-cutting)**
   - Add BasisHandle/LoadingsHandle unions so LatentNeuroVec can store explicit matrices or shared/implicit dictionaries (DCT, lifted heat wavelets).
   - Provide registry helpers and operator-aware accessors `basis_mat()` / `loadings_mat()`; refactor core access paths to use them.
   - Ship first handles: `dct_basis_handle()`, `heat_wavelet_loadings_handle()`, and a convenience constructor `latent_dct_heatwavelet()`.

7) **Slepian / DPSS pathway**
   - Temporal: wrap existing DPSS/Slepian generator in `slepian_temporal_handle()` and add `slepian_temporal_latent()` (handle-backed basis; matrix loadings).
   - Spatial: implement `lift(ClusterReduction, spec_slepian)` to produce sparse voxel×k loadings from local graph Laplacians; add `slepian_spatial_loadings_handle()` for sharing.
   - Spatiotemporal: prototype `slepian_spatiotemporal_latent()` using handles for both temporal and spatial Slepians plus a small core; start with explicit core tensor, consider `ImplicitLatent` if decoder form is cleaner.
   - Tests: handle vs explicit parity; sparsity checks; round-trip on tiny masks; partial-access (linear/matricized) without full reconstruction.

8) **Encoding API & factory**
   - Add small spec constructors (`spec_time_*`, `spec_space_*`, `spec_st`) to make basis choices declarative.
   - Provide a single `encode(x, spec, ...)` generic that builds the appropriate latent (handles by default).
   - Optional newcomer helper `latent_factory(x, family, ...)` that maps simple family strings to specs and calls `encode`.
   - Cover with tests for parity (handle vs matrix), class/type correctness, and a tiny spatiotemporal path.

## Forward-looking (mixable bases)
- **Temporal families:** DPSS (band-limited), B-splines (smooth, low-order), DCT (orthonormal, fast). Provide `*_latent()` wrappers with `meta$family = "<name>_time"`.
- **Spatial families:** HRBF (analytic), Slepian (graph concentration, supports binary or weighted masks e.g., gray-matter prob), maybe DCT on voxel grid for fast orthonormal bases. Tag as `<name>_space`.
- **Mix-and-match:** Allow pairing a temporal family with a spatial family (e.g., DPSS_time + HRBF_space, DCT_time + Slepian_space). Encode compatibility rules in docs; keep `meta$families` as a vector (e.g., `c("dpss_time","hrbf_space")`). Avoid new classes; reuse `LatentNeuroVec` plus meta tags.
## Deliverables
- R implementations (and optional Rcpp) under `R/` and `src/`.
- testthat coverage for HRBF paths.
- Updated `NAMESPACE` and man pages via `devtools::document()`.


# Diffusion / Heat Wavelets (rgsp) Plan

Goal: add a graph-based, multiscale basis that bridges Haar (implicit, blocky) and Slepian (spectral) using diffusion/heat wavelets built via `rgsp`, while keeping the latent API consistent with existing constructors.

Conventions (shared): expose user-friendly constructors (`*_latent()`) that return `LatentNeuroVec` for explicit bases or `ImplicitLatent` (with `meta$family`) for implicit/overcomplete transforms. Avoid requiring users to pass basis/loadings directly.

Milestones
1) Graph bridge helper
   - `voxel_subset_to_gsp(mask, voxel_indices, k_neighbors = 6)` using `neuroim2::indices` → `rgsp::graph_knn`; ensure Laplacian ready for filters.

2) Basis spec
   - `basis_heat_wavelet(scales = c(1,2,4,8), order = 30, threshold = 1e-6)`; purely a spec object tagged `spec_heat_wavelet`.

3) Supervoxel lift (explicit basis path)
   - Method `lift(reduction = "ClusterReduction", basis_spec = "spec_heat_wavelet")`:
     * For each cluster: build local graph, run `rgsp::wavelet_heat_transform` on identity to materialize atoms per scale.
     * Threshold small entries; assemble block-diagonal sparse loadings (global voxel rows, concatenated atoms by scale/cluster).
   - Constructor `heat_wavelet_latent(X, mask, reduction, spec)` returning `LatentNeuroVec` (basis = coeff, loadings = sparse atoms); mirror HRBF pattern.

4) Whole-brain implicit path
   - `rgsp_heat_implicit(X, mask, scales, order = 30, threshold = ...)`:
     * Build global graph; forward transform (time x voxels → coeff tensor); sparsify coefficients.
     * Wrap in `ImplicitLatent` with `meta$family = "heat_wavelet"` and decoder that calls inverse heat transform.

5) Tests & benchmarks
   - Small masks: round-trip within tolerance; thresholding reduces nnz; deterministic with fixed seeds.
   - Performance smoke vs Haar/HRBF on modest masks (e.g., 10k voxels).

6) Docs
   - Add roxygen for bridge, basis spec, constructors; briefly document explicit vs implicit modes and constructor naming convention.

# Spatial/Temporal Basis Integration Plan (Slepian, PCA, etc.)

## Goal
Add DPSS/Slepian plus other bases (local PCA, flat/means) as first-class options using the shared pipeline: reduction → basis → lift → latent. Keep `LatentNeuroVec`/`ImplicitLatent` as the only latent shapes; return sparse loadings whenever possible.

## Milestones
1) **Temporal DPSS (band-limited time basis)**
   - RcppEigen tridiagonal DPSS generator `generate_dpss_basis_rcpp`.
   - `slepian_temporal_latent(X, mask, tr, bandwidth, k = NULL, denoise = TRUE)` with meta tagging.
   - Tests: compare to multitaper output on small T; round-trip encode/decode.

2) **Reduction & Lift API (basis-agnostic)**
   - S4 `GraphReduction` base; `ClusterReduction` (voxel→cluster map) and `CoarsenedReduction` (prolongation `P_matrix`).
   - Basis specs (`basis_slepian`, `basis_pca`, `basis_flat`, extensible) and `lift()` generic to build voxel×components loadings (sparse when possible) regardless of basis choice.
   - Tests: block-diagonal sparsity for clusters; P-based lift correctness on toy graphs; PCA vs Slepian parity on small parcels.

3) **Supervoxel/Cluster Spatial Bases (default spatial path)**
   - `graph_reduce(method="supervoxels")` + `lift(basis_slepian)` for spectral Slepians.
   - Reuse the same path for local PCA/means: `lift(basis_pca)` / `lift(basis_flat)` to get parcel PCA or parcel averages.
   - Controls for sparsity/compression (k_per_voxel, n_supervoxels); meta tags for chosen family (e.g., `"slepian_supervoxel"`, `"pca_supervoxel"`).
   - Tests: reconstruction on small masks; sparsity checks; seed determinism.

4) **Spatiotemporal Combination**
   - Option A (separable): pair DPSS_time with any lifted spatial basis (Slepian/PCA/flat) via LatentNeuroVec meta tags.
   - Option B (clean tensor): `*_spatiotemporal_latent()` returning `ImplicitLatent` with time_basis, space_basis, core C = t(B_t) %*% X %*% L_s.
   - Tests: toy ROI compression and decode correctness.

5) **Optional Graph Coarsening Path (advanced)**
   - `graph_reduce(method="coarsen")` + `lift` for coarse global gradients; mark experimental.
   - Benchmarks on density/memory vs supervoxel path.

6) **Docs & Vignettes**
   - Update roxygen for new constructors; vignette comparing PCA vs DPSS temporal; spectral supervoxel example; note default stance (supervoxel first, coarsen optional).
