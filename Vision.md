# Vision: Latent fMRI in fmrilatent

We aim to provide lightweight, high-performance latent representations
for fMRI data with clean R interfaces and optional Rcpp acceleration.
HRBF is the first flavor; others (PCA/ICA, wavelets, DPSS, B-splines,
DCT, spatial Slepian) follow a shared contract: explicit-basis families
use `LatentNeuroVec`; implicit/dictionary-free families (e.g., Haar
lifting) use `ImplicitLatent` with `meta$family` tags plus
`as_<family>_latent()`, `<family>_meta()`, `is_<family>_latent()`
helpers. Spatial Slepian should accept binary or weighted domains (e.g.,
gray-matter probability); temporal DPSS/DCT/B-splines can mix with
spatial families under clear compatibility rules. Shared or fixed
dictionaries should be representable via handles rather than duplicating
matrices.

## Template-backed hierarchical frames

- Accept that locality + multiscale beats orthogonality: hierarchical
  Laplacian bases stay sparse; we solve in the dual using cached sparse
  Cholesky of `G = B'B`.
- Introduce a template-only asset (`HierarchicalBasisTemplate`)
  capturing mask, parcellations, sparse loadings `B`, and the solver; it
  is data-agnostic and reusable across studies.
- Packaged templates (e.g., MNI + Schaefer100/400 + subcortex) live in
  `inst/extdata` with documented build scripts; factories
  (`latent_factory(family="hierarchical_mni")`) load and project without
  rebuilding.
- Subcortex handled as a separate block to keep `B` block-diagonal;
  small nuclei may use tiny k or identity modes to avoid unstable
  eigenvectors.

## Principles

- **Speed-first:** Favor sparse matrices, vectorized math, and C++ paths
  where they clearly win.
- **Deterministic by default:** Seeds and sampling should yield
  repeatable bases.
- **Minimal dependencies:** Keep external requirements small; avoid
  heavy pipelines or HDF5 bindings here.
- **Composable:** All basis/encoding functions should return objects
  usable directly with `LatentNeuroVec`, either as explicit matrices or
  as handles that can materialize lazily and be shared.

## Architectural Lens: reduction → basis → lift → latent

- **Reduction (topology):** how we simplify or partition voxels
  (supervoxels/atlas via `ClusterReduction`; optional coarse graphs via
  `CoarsenedReduction`). Defaults favor partitions that yield
  block-sparse loadings.
- **Basis (math):** what functions live on the reduced topology
  (Slepian, PCA, flat/identity, wavelets; extensible); specified via
  lightweight `basis_*` specs.
- **Lift:** generic [`lift()`](reference/lift.md) turns (reduction,
  basis_spec) into a voxel×components loading matrix (sparse when
  possible) or a decoder closure for implicit cases—works identically
  for Slepian, PCA, or flat means.
- **Latent:** keep only two shapes—`LatentNeuroVec` for explicit or
  handle-backed bases/loadings; `ImplicitLatent` when a core tensor plus
  decoders is cleaner (e.g., spatiotemporal factorizations).

## Defaults and Stance

- Temporal denoising/compression: DPSS (band-limited) as the first-class
  temporal family—no data SVD, physics-informed.
- Spatial compression at scale: supervoxel reduction + lift as the
  default spatial path; basis choice within that path is pluggable
  (Slepian for physics/spectral control, local PCA for data-driven, flat
  for means). Coarsened/global modes are optional and documented as
  denser/advanced.
- Sparsity-first: prefer constructions that keep loadings sparse; dense
  paths are opt-in and clearly marked.

## Roadmap Snapshot (augmented)

- Deliver HRBF latent pipeline (R + optional Eigen/OpenMP).
- Add temporal DPSS and the reduction/basis/lift API, then spectral
  supervoxel Slepian as the default spatial Slepian implementation.
- Stand up the hierarchical Laplacian template path: class, builder,
  packaged MNI template, encode/decode helpers, and tests for
  sparsity/round-trip.
- Provide spatiotemporal pairing (DPSS_time × Slepian_space) via
  `LatentNeuroVec` meta tags or an `ImplicitLatent` core when tensor
  structure is clearer.
- Add benchmarks and guidance for choosing R vs C++ backend and for
  selecting supervoxel vs coarsened reductions.
- Add dictionary handles + registry (BasisHandle/LoadingsHandle) so
  fixed/shared dictionaries (DCT, heat-wavelet) can be stored by
  reference.
- Expand to additional basis families once HRBF is stable, reusing the
  same encode/decode contracts.
