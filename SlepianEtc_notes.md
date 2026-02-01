# SlepianEtc Notes

## Part 1 — 2025-11-29

- Why Slepian/DPSS: BOLD is band-limited; DPSS maximizes energy
  concentration in a target bandwidth/region, giving physics-informed
  compression and denoising versus PCA variance-only bases.

- Two modes

  - **Temporal Slepian → LatentNeuroVec**: fixed DPSS dictionary `B`
    (time × k) derived from TR and bandwidth; voxel loadings `L` found
    by projection. Reconstruction: `X ≈ B_dpss %*% t(L)`; automatically
    removes high-frequency physio noise.
  - **Spatiotemporal Graph Slepian → ImplicitLatent**: temporal DPSS
    plus spatial graph eigenmodes (smallest Laplacian eigenvalues)
    respecting voxel adjacency/ROI topology. Reconstruction:
    `X ≈ B_dpss %*% C_core %*% t(U_graph)`; strong compression for ROIs
    but computationally heavier.

- DPSS generation (avoid heavy deps): implement tridiagonal eigen-solver
  with RcppEigen.

- Temporal constructor sketch (`R/slepian.R`)

  ``` r
  slepian_temporal_latent <- function(X, mask, tr, bandwidth = 0.1, k = NULL) {
    n_time <- nrow(X)
    W  <- bandwidth * tr           # normalized half-bandwidth
    NW <- n_time * W                # time-bandwidth product
    if (is.null(k)) k <- max(1, floor(2 * NW) - 1)  # Shannon number

    B <- generate_dpss_basis_rcpp(n_time, NW, k)    # orthogonal DPSS basis
    L <- Matrix::crossprod(X, B)                    # loadings (voxels × k)

    LatentNeuroVec(
      basis    = Matrix::Matrix(B),
      loadings = Matrix::Matrix(L),
      space    = neuroim2::space(mask),
      mask     = mask,
      offset   = colMeans(X),
      meta     = list(family = "slepian_temporal", tr = tr, bandwidth = bandwidth, nw = NW)
    )
  }
  ```

- Graph Slepian sketch (`R/slepian_spatial.R`): build voxel adjacency →
  Laplacian; `RSpectra::eigs_sym(laplacian, k_space, which = "SM")` for
  smooth spatial modes `B_s`; temporal DPSS `B_t`; core
  `C <- t(B_t) %*% X %*% B_s`; decoder multiplies subsets on demand.

- Comparative value: denoising “for free” by keeping first `K = 2NW - 1`
  tapers; no data SVD needed (basis analytic); coefficients may store
  sparsely; PCA remains general, Haar provides locality, Slepian gives
  bandwidth-aware frequency interpretability.

- Roadmap: implement temporal Slepian first; add Rcpp tridiagonal DPSS
  generator; optional `denoise = TRUE` to drop out-of-band tapers; later
  add graph Slepian for ROI compression.

## Part 2 — 2025-11-29

- Core hurdle: full-brain spatial Slepian (≈1e5 voxels) naive eigensolve
  is O(N^3); need sparsity/locality tricks.

- Approach 1 — Cluster-Slepian (recommended): parcel mask (e.g.,
  Schaefer 400 or K-means) into ~500–1000 clusters; for each parcel
  solve tiny Laplacian eigenproblem (k_per_cluster = 1 gives parcel
  mean; 2–5 adds gradients). Assemble block-diagonal sparse `L`;
  temporal basis `B <- X %*% L`; yields LatentNeuroVec with sparse
  loadings and natural “lifting” (voxel signal → smooth parcel modes).

  ``` r
  slepian_cluster_latent <- function(X, mask, clusters, k_per_cluster = 3, tr, bw) {
    cluster_ids <- sort(unique(clusters[clusters > 0]))
    i <- j <- x <- list(); col_off <- 0
    for (cid in cluster_ids) {
      vox <- which(clusters == cid)
      lap <- get_laplacian(get_local_adjacency(mask, vox))      # small (<= few hundred)
      k <- min(k_per_cluster, length(vox))
      vecs <- RSpectra::eigs_sym(lap, k = k, which = "SM")$vectors
      for (cc in seq_len(k)) {
        i[[length(i)+1]] <- vox
        j[[length(j)+1]] <- rep(col_off + cc, length(vox))
        x[[length(x)+1]] <- vecs[, cc]
      }
      col_off <- col_off + k
    }
    L <- Matrix::sparseMatrix(i = unlist(i), j = unlist(j), x = unlist(x),
                              dims = c(sum(mask), col_off))
    B <- as.matrix(X %*% L)
    LatentNeuroVec(basis = Matrix::Matrix(B), loadings = L,
                   space = neuroim2::space(mask), mask = mask)
  }
  ```

- Approach 2 — Nyström: pick ~5% voxels as landmarks, solve eigenproblem
  there (~5k nodes), extend eigenvectors to remaining voxels via
  Nyström/interpolation (HRBF/linear); optionally threshold to sparsify
  `L`; keeps smooth whole-brain gradients without hard parcel
  boundaries.

- Approach 3 — Tensor-product lattice: assume regular 3D grid; take
  Kronecker of 1D Slepians/DCT along x/y/z, then mask out non-brain
  voxels; avoids eigendecomp entirely, extreme speed on regular grids.

- Helper idea: `lift_parcellation(parcellation, mask, k=3)` builds the
  sparse block-diagonal lifting matrix (per-parcel eigenmodes) for reuse
  across runs.

## Part 3 — 2025-11-29

- Clustering vs Graph Coarsening: partition (supervoxels) makes
  Laplacian block-diagonal; coarsening builds multiscale coarse graph
  plus prolongation `P` to lift eigenvectors
  (`v_fine = P %*% v_coarse`).

- Tradeoffs

  - Clustering: extreme sparsity, embarrassingly parallel tiny
    eigenproblems, hard boundaries that cut inter-parcel edges—great for
    storage/compression; kills very global gradients.
  - Coarsening: preserves topology and global gradients; lifting via `P`
    spreads support → denser loadings, heavier implementation
    (multigrid/Kron reduction).

- Fit for fmrilatent: prefer supervoxel/cluster Slepians—maximally
  sparse loadings suit `LatentNeuroVec`; blockiness mitigated by keeping
  k\>1 modes per parcel (mean + gradients).

- “Soft supervoxel” compromise: keep supervoxel solves but allow
  boundary voxels to share weights across adjacent parcels, turning
  block-diagonal into lightly banded loading matrix (smooths seams
  without dense blow-up).

- Supervoxel constructor sketch

  ``` r
  slepian_supervoxel_latent <- function(X, mask, n_supervoxels = 500, k_per_voxel = 3) {
    sv <- neuroim2::cluster_voxels(mask, k = n_supervoxels, method = "slic")
    i <- j <- x <- list(); col <- 0
    for (id in sort(unique(sv[sv > 0]))) {
      vox <- which(sv == id)
      lap <- get_local_laplacian(mask, vox)
      vecs <- RSpectra::eigs_sym(lap, k = min(k_per_voxel, length(vox)), which = "SM")$vectors
      for (cc in seq_len(ncol(vecs))) {
        i[[length(i)+1]] <- vox
        j[[length(j)+1]] <- rep(col + cc, length(vox))
        x[[length(x)+1]] <- vecs[, cc]
      }
      col <- col + ncol(vecs)
    }
    L <- Matrix::sparseMatrix(i = unlist(i), j = unlist(j), x = unlist(x),
                              dims = c(sum(mask), col))
    B <- as.matrix(X %*% L)
    LatentNeuroVec(basis = B, loadings = L, space = neuroim2::space(mask),
                   mask = mask, meta = list(method = "slepian_supervoxel"))
  }
  ```

- Guideline: use clustering for compression/denoising pipelines;
  consider coarsening only when smooth whole-brain gradients are primary
  and added density is acceptable.

## Part 4 — 2025-11-29

- Goal: general Lifting API decoupling reduction (topology
  simplification) from basis math; plug-in architecture for Slepians,
  PCA, flat averaging, etc.

- Pipeline: Reduce (ambient voxels → reduced graph via
  clusters/coarsening/landmarks) → Define Basis (e.g., Laplacian
  eigenvectors, local PCA, identity) → Lift (pull reduced basis back to
  ambient; produce sparse/dense loadings).

- Classes: abstract `GraphReduction` (mask, adjacency optional);
  `ClusterReduction` (map voxel→cluster, cluster_ids);
  `CoarsenedReduction` (prolongation `P_matrix`, coarse coords). Basis
  specs are light structs: `basis_slepian(k,type)`,
  `basis_pca(k,whiten)`, [`basis_flat()`](reference/basis_slepian.md).

- Generic lifter

  ``` r
  setGeneric("lift", function(reduction, basis_spec, data = NULL, ...) standardGeneric("lift"))
  ```

  - `lift(ClusterReduction, spec_slepian)`: loop clusters, build local
    adjacency/laplacian, solve `k` smallest eigs (RSpectra), assemble
    block-diagonal `dgCMatrix`.
  - `lift(CoarsenedReduction, spec_slepian)`: solve eigs on coarse
    laplacian, lift via `P_matrix %*% V_coarse` (denser, band-sparse).

- Workflow: `mask` → `graph_reduce(..., method="supervoxels"/"coarsen")`
  → `spec <- basis_slepian(k=5)` → `L <- lift(reducer, spec)` →
  `LatentNeuroVec(basis = X %*% L, loadings = L, ...)`.

- Elegance: separation of concerns (topology vs basis vs storage),
  extensible (add wavelets/new reducers without touching lifter users),
  defaults return
  [`Matrix::sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html)
  to stay memory-light; implicit latents could return decoder closures
  under same API.

## Part 6 — 2025-11-29

- Plug-in fit to current `LatentNeuroVec`: spatial atoms live in
  `loadings`, temporal atoms in `basis`; Slepian_space can replace HRBF
  with no new class. Pattern:

  ``` r
  B_atoms_vox <- slepian_space_basis(mask, roi_mask, params)   # atoms x vox
  coeff       <- slepian_project_matrix(X, mask, roi_mask, params) # time x atoms
  loadings    <- Matrix::t(B_atoms_vox)
  LatentNeuroVec(basis = coeff, loadings = loadings,
                 space = neuroim2::NeuroSpace(...),
                 mask = mask, meta = list(family = "slepian_space", params = params))
  ```

- Spatial graph Slepian recipe: build voxel graph; take first K
  Laplacian eigenvectors (band-limit); solve Slepian concentration in
  ROI to get K’ orthonormal atoms concentrated in ROI. Orthogonal,
  bandwidth- and ROI-controllable; fewer hyperparams than HRBF.

- Temporal DPSS family: treat 1D Slepians as a temporal family
  (`meta$family = "dpss_time"`). You can combine DPSS_time ×
  Slepian_space while keeping separable
  `basis[,k] * loadings[,k]`—encode pairing in meta (time_family,
  space_family, params). Full joint spatiotemporal Slepians would break
  separability and would need a new class (time_basis, space_basis,
  coeff grid)—defer as future phase.

- API sketches mirroring HRBF:

  - `slepian_space_basis(mask, roi_mask=mask, bandwidth=K, n_slepian=K', laplacian=c("unnormalized","normalized")) -> atoms x vox`
  - `slepian_project_matrix(X, mask, roi_mask, params) -> time x atoms`
  - `slepian_latent(X, mask, roi_mask=mask, params, label="")` returns
    `LatentNeuroVec` tagged `family="slepian_space"`.
  - Helpers: `as_slepian_latent()`, `is_slepian_latent()`,
    `slepian_params()`.

- Practical checklist:

  1.  `make_voxel_graph(mask, connectivity=6/18/26, weight_scheme="binary"|"distance")`.
  2.  `slepian_space_basis(...)` (store eigenvalues/concentration).
  3.  `slepian_project_matrix(...)` (time × atoms).
  4.  Temporal DPSS: `dpss_time_basis(T, nw, K)`; optional
      `as_dpss_time_latent`.
  5.  Mixing helpers to combine/compare HRBF vs Slepian vs DPSS latents
      using meta tags.

## Part 7 — 2025-11-29

- Core mantra: `reduction → basis → lift → latent`; Slepian is one basis
  choice, not a bolt-on. Keep single latent shapes (LatentNeuroVec;
  ImplicitLatent when a core tensor makes sense) and separate topology
  from math; default to sparse loadings.

- Layered architecture

  - Layer 0 primitives: DPSS generator in RcppEigen
    `generate_dpss_basis_rcpp(N, NW, K)` with R wrapper
    `dpss_time_basis(n_time, tr, bandwidth, k=NULL)`.
  - Layer 1 engine: `GraphReduction` subclasses (`ClusterReduction` map
    voxel→cluster; `CoarsenedReduction` with `P_matrix`), `basis_*`
    specs (`basis_slepian`, `basis_pca`, `basis_flat`), and
    [`lift()`](reference/lift.md) generic. Cluster+Slepian →
    block-sparse L; Coarsen+Slepian → denser via P; Cluster+PCA →
    data-driven local basis.

- User workflows

  - Temporal DPSS (Mode A):
    `slepian_temporal_latent(X, mask, tr, bandwidth=0.1, k=NULL, denoise=TRUE)`
    builds DPSS B, L = X’ %\*% B, returns LatentNeuroVec with meta
    {family=“slepian_temporal”, tr, bandwidth, nw}.
  - Spectral supervoxel Slepian (Mode B default spatial):
    `slepian_supervoxel_latent(X, mask, n_supervoxels=500, k_per_voxel=3)`
    → reduction=supervoxels, spec=basis_slepian, L = lift(…), B = X %\*%
    L, LatentNeuroVec meta {family=“slepian_supervoxel”}; columns:
    parcel mean + gradients.
  - Spatiotemporal: combine DPSS_time × Slepian_space. Simple path keeps
    separable LatentNeuroVec (pairing encoded in meta). Clean path uses
    ImplicitLatent with core C = t(B_t) %*% X %*% L_s, time_basis = B_t,
    space_basis = L_s, meta tags {time_family=“slepian_temporal”,
    space_family=“slepian_supervoxel”}.

- Official stance on whole-brain eigensolves: default to
  cluster/supervoxel reductions for scalability and sparsity; offer
  CoarsenedReduction explicitly for global gradients (warn about denser
  loadings); full whole-brain Slepians out-of-scope for now.

- Roadmap

  1.  Phase 1: Temporal DPSS (C++ gen, slepian_temporal_latent,
      tests/vignette vs PCA).
  2.  Phase 2: Reduction + lift API (GraphReduction, graph_reduce, basis
      specs, lift; ClusterReduction + basis_slepian path working).
  3.  Phase 3: Spectral supervoxel Slepian (k_per_voxel/n_supervoxels
      knobs, denoise/compression controls).
  4.  Phase 4: Spatiotemporal ImplicitLatent (B_t, L_s, C) with
      benchmarks.
  5.  Phase 5: Optional graph coarsening path (advanced/experimental).
