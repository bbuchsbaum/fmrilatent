Diffusion Wavelets (Section 2.3 & Figure 2)
===========================================

- Goal: build data-adaptive multiscale bases on a graph/manifold by diffusing a Markov operator `T` and compressing its range as dyadic powers filter high-frequency detail.
- Dyadic powers: examine `T^(2^j)`; as `j` grows, diffusion spreads and numerical rank drops.
- Compression/basis update: at each scale, run a modified, rank-revealing QR on the current operator representation to find a coarse orthonormal scaling basis `phi_{j+1}` within tolerance `epsilon`. Scaling functions are stored “extended” back to the original basis `phi_0`. Wavelets are the orthogonal complement.
- Downsampling the operator: project/square the compressed operator to form `T^(2^{j+1})` in the new basis, iterating until the space collapses (rank 1 or stopping rule).

Faithful Pseudocode (Figure 2-aligned)
--------------------------------------

```python
def DiffusionWavelets(T, epsilon):
    # Inputs: T (N x N transition/similarity matrix), epsilon (compression precision)
    # Outputs: extended scaling bases phis, extended wavelet bases psis (each in terms of phi_0)
    phi_0 = I_N                          # Dirac basis on original nodes
    current_T = T                        # [T]_phi0^phi0
    extended_phis = [phi_0]              # [phi_j]_phi0
    extended_psis = []
    j = 0
    while size(current_T) > 1:           # or other stop (rank collapse / max_J)
        Q, R = ModifiedQR(current_T, epsilon)        # Q ~ [phi_{j+1}]_phi_j ; R ~ [T^(2^j)]_phi_j^phi_{j+1}
        phi_next_global = Q @ extended_phis[j]       # [phi_{j+1}]_phi0
        extended_phis.append(phi_next_global)

        Q_wavelet, _ = ModifiedQR(I - Q @ Q.T, epsilon)  # wavelet complement in phi_j coordinates
        psi_next_global = Q_wavelet @ extended_phis[j]   # [psi_j]_phi0
        extended_psis.append(psi_next_global)

        compressed_T = R @ Q                          # [T^(2^j)]_phi_{j+1}^phi_{j+1}
        current_T = compressed_T @ compressed_T       # square for dyadic step: [T^(2^{j+1})]_phi_{j+1}^phi_{j+1}
        j += 1
    return extended_phis, extended_psis
```

Scalable QR Proposal
--------------------

- Bottleneck: deterministic rank-revealing QR on dense `N x N` matrices is `O(N^3)` and diffusion causes fill-in. Intractable for many-thousand–million node graphs.
- Remedy: randomized rank-revealing QR (Halko–Martinsson–Tropp style) to approximate the column space before factoring.
- Steps (target rank `k`, oversample `p`):
  1. Draw Gaussian `Omega` of shape `n x (k+p)`.
  2. Form sketch `Y = A @ Omega` (fast with sparse/implicit `A`).
  3. Orthonormalize `Y -> Q = qr(Y)`.
  4. Project small matrix `B = Q.T @ A`.
  5. Factor small `B` with deterministic QR/SVD: `U_small, R_small = qr(B)`.
  6. Reconstruct basis `Q_final = Q @ U_small`, giving `A ≈ Q_final @ R_small`.
- Complexity drops toward `O(N * k^2)` or `O(N^2 log k)`, with the heavy lift in the highly parallelizable multiply `A @ Omega`.

Scaling for 100K-node brain connectivity
----------------------------------------

- Short answer: vanilla DWT (and basic randomized QR) will not scale on a 100K-node small-world connectome unless extra constraints are added.
- Fill-in from small-world topology: degree ~100 gives ~10K neighbors at `T^2` and near-full coverage by `T^4`, making `T^(2^j)` dense (~80 GB at 100K x 100K) and QR `O(1e15)` FLOPs.
- Slow spectral decay: connectome spectra decay slowly, so epsilon-driven truncation may still require tens of thousands of basis vectors (e.g., 40K), making QR of a 100K x 40K matrix infeasible.

Constraints to make it workable
- Never form `T^2` explicitly: evaluate products as `Y = T(T Omega)` so cost scales with `nnz(T)` rather than `N^2`.
- Aggressive sparsification: threshold small entries (e.g., `< 1e-5`) after compression to prevent density growth.
- Fix a target rank (e.g., 2K–5K) instead of epsilon-based rank discovery; accept loss of some high-frequency detail to stay tractable.

Sketch of a scalable loop (brain-sized)

```python
def Scalable_Brain_DWT(T_sparse, target_rank=2000, oversample=20, max_scales=...):
    current_T = T_sparse  # 100K x 100K sparse (~10M nnz)
    for _ in range(max_scales):
        Omega = randn(current_T.shape[0], target_rank + oversample)
        Y = current_T.dot(current_T.dot(Omega))          # implicit T^2 @ Omega, sparse-dense
        Q, _ = np.linalg.qr(Y, mode="reduced")           # tall-skinny QR
        B = current_T.dot(Q)                             # sparse-dense multiply
        T_compressed = B.T.dot(B)                        # small dense (target_rank x target_rank)
        T_compressed[np.abs(T_compressed) < 1e-5] = 0    # enforce sparsity if desired
        save_basis(Q)                                    # store basis analog of [phi_{j+1}]_phi_j
        current_T = T_compressed
        if current_T.shape[0] < target_rank:
            break
```

- Effect: the Figure 2 flow would blow up at the first squaring; this constrained randomized path can finish on realistic hardware by capping rank and keeping operations implicit.

Applicability to voxel vs cluster/coarsened graphs
--------------------------------------------------

- Works on raw voxel graphs (100K+ nodes) if using the implicit, capped-rank path above; sparsity and thresholding are essential.
- Directly usable on reduced/clustered/coarsened graphs; smaller N eases the rank cap and may allow epsilon-driven truncation.
- Bases and compressed operators can be composed across resolutions: run DWT on a coarse graph, then lift or restrict via cluster assignments when mixing with other methods (e.g., spectral embeddings, HRBF, GCNs).
- The randomized, implicit multiplies mean the same code path can operate on any graph representation as long as `T` is provided as a sparse/linear operator.

Implementation plan for fmrilatent
----------------------------------

- API/spec: add `basis_diffusion_wavelet()` spec (params: target_rank/oversample, threshold, max_scales, epsilon optional) in a new `R/diffusion_wavelet.R`; mirror patterns in `basis_heat_wavelet` and use S3 class `spec_diffusion_wavelet`.
- Operator interface: accept `T` as `Matrix` or as a function `f(x)` returning `T %*% x`; never materialize `T^2`. Provide helper `implicit_power_mult(T_op, X)` to evaluate `T(T X)`.
- Core routine: implement `randomized_diffusion_step(T_op, target_rank, oversample, threshold)` returning `(Q, T_compressed)` using the implicit two-pass multiply `Y = T(T Omega)`, tall-skinny QR, and small `T_compressed = B^T B`; include optional thresholding of `T_compressed`.
- Integration with reductions: add a `lift` method `lift(ClusterReduction, spec_diffusion_wavelet, ...)` that builds `T_op` from the cluster graph (or raw voxel graph via `voxel_subset_to_gsp` / adjacency), then produces stacked scaling functions `[phi_{j+1}]_phi0` analogous to other basis lifts.
- Latent constructor: add `diffusion_wavelet_latent(X, mask, reduction = NULL, spec = basis_diffusion_wavelet(), ...)` to assemble a `LatentNeuroVec` (basis = `X %*% loadings`, loadings = scaling funcs); store meta with family = "diffusion_wavelet" and spec parameters.
- Interop: ensure loadings can be composed with HRBF or spectral reductions by exposing them as `Matrix` objects and keeping the operator interface linear; support both voxel-space and coarsened graphs by swapping the `T_op` builder.
- Testing: add testthat specs with small deterministic graphs (voxel mask toy) to check rank-capped compression, orthonormality of `Q`, and reconstruction error vs. direct dense computation; include perf smoke (microbenchmark) for `T(T Omega)` sparse-dense multiplies.
