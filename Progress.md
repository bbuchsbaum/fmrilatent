# Progress (HRBF + Basis Pipeline in fmrilatent)

## Current Status
- Repo has `LatentNeuroVec` core classes and now accepts Matrix or handle-backed bases/loadings (BasisHandle/LoadingsHandle unions + registry).
- HRBF R path added (basis generation, encode/decode, LatentNeuroVec wrapper).
- HRBF C++ path (Eigen/OpenMP) ported; compiled via Rcpp attributes.
- Tests in `tests/testthat` cover small-mask roundtrip and kernel alias; package builds and tests pass locally.
- Slepian design recorded in `SlepianEtc_notes.md`; plan added to `Plan.md` and vision updated for reduction→basis→lift pipeline.
- DPSS temporal path scaffolded: Rcpp DPSS generator, `dpss_time_basis()`, and `slepian_temporal_latent()` added (no tests yet).
- Abstract reduction/basis/lift scaffolding added (`GraphReduction`, `ClusterReduction`, `CoarsenedReduction`, basis specs, `lift` generic); concrete implementations remain external.
- DPSS tridiagonal backend implemented (LAPACK `dstev`); tests for dense vs tridiag parity and projection correctness now pass.
- Haar lifting ported: pure R + Rcpp path, implicit Morton-based transform with `ImplicitLatent`/`HaarLatent`, partial reconstruction (`levels_keep`, ROI/time), option flag `fmrilatent.haar.use_rcpp`.
- Introduced `ImplicitLatent` container and helpers for basis-free decoders; Haar now uses this track while keeping predict/as.matrix entry points.
- New operator-aware accessors `basis_mat()` / `loadings_mat()` route all decoding to either explicit matrices or handles; registry avoids duplicating shared dictionaries.
- Added first handle-backed families: `dct_basis_handle()` (lazy temporal DCT), `heat_wavelet_loadings_handle()` (shared spatial dictionary), and a convenience constructor `latent_dct_heatwavelet()`.
- Added spatiotemporal Slepian to near-term roadmap: temporal handles (`slepian_temporal_handle`), spatial lift (`lift(..., spec_slepian)` + `slepian_spatial_loadings_handle`), and a prototype `slepian_spatiotemporal_latent()` using handles + core; tests planned for round-trip/sparsity/partial access.
- New encoding API planned: spec constructors (`spec_time_*`, `spec_space_*`, `spec_st`), a unified `encode()` front-end (handle by default), and an optional `latent_factory()` for newcomers; tests will cover handle vs matrix parity and tiny spatiotemporal paths.

## Next Actions
- Benchmark R vs C++ on medium/large masks; document recommended toggle defaults. Initial runs (16³, 24³) show C++ slower; default backend switched to R. Need to retest on larger masks after C++ tuning.
- Extend tests to medium masks and both kernels; add encode/decode tolerances with noise.
- Generate/update man pages with `devtools::document()` once imports are stable.
- Begin Slepian/DPSS phase:
  - Add tests for DPSS temporal (vs multitaper on toy data) and slepian_temporal_latent round-trip.
  - Wire concrete lift methods via external reduction providers (supervoxel/atlas), or stub a minimal internal path for flat/PCA if needed.
  - Prototype `slepian_supervoxel_latent()` once reduction provider is available; include deterministic seeding and sparsity checks.
- Clean up roxygen warnings (unknown imports, example braces) and run full test suite after document rebuild.
- New plan added for diffusion/heat wavelets (rgsp): build voxel→gsp bridge, spec object, supervoxel lift (explicit `LatentNeuroVec`) and whole-brain implicit variant; add tests/benchmarks.
- Unify constructor convention: keep user-facing `*_latent()` wrappers (explicit or implicit) so users never pass basis/loadings directly; document in forthcoming man pages.
- Add tests for handle-backed paths (DCT + heat-wavelet) ensuring reconstruction matches explicit equivalents and registry reuse works.

- Hierarchical Laplacian frame work started: added `HierarchicalBasisTemplate` class, builder/encoder helpers, RSpectra dependency, and offline template build script stub. Added generic hierarchy helpers (`cut_hclust_nested`, nesting validation, parent maps). Plan updated with Schaefer-driven hierarchical construction.
- **Schaefer geodesic integration added**: new functions `build_schaefer_levels()` and `build_schaefer_hierarchical_template()` in `hierarchical_helpers.R`. These bridge neuroatlas (Schaefer surface atlas) and neurosurf (geodesic distances, boundary contact) to build hierarchical templates with geodesic-informed parcel clustering. Requires pending neurosurf functions: `parcel_geodesic_distance_matrix()`, `parcel_boundary_contact()`.

## 2025-02 diffusion wavelets follow-ups
- Current graph build is internal (`build_cluster_graph()` -> `rgsp::graph_knn(..., k_neighbors, weight = "distance", sym = "union")` + row-normalization).
- Potential API additions:
  - Expose k-NN options (sym rule, distance metric, epsilon radius) and allow passing a precomputed adjacency/transition.
  - Hook to override graph construction in spec/constructors for voxel-level or cluster-level inputs.

## Notes / Risks
- Determinism: Poisson sampling seeded; verify after refactors.
- Performance: R path may be slow on ≥64³ masks; need benchmarks before promising defaults. Current C++ path slower on 16³/24³ (R≈0.46s vs C++≈0.76s; R≈3.8s vs C++≈7.4s).
- Slepian spatial pieces depend on robust adjacency builders; performance and memory need profiling once lift is wired.
- Diffusion/heat wavelets depend on `rgsp` availability and graph construction cost; supervoxel strategy chosen to keep atom matrices sparse and tractable.
