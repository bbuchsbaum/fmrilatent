# BOLDZip-SR PRD Notes

**Status**: proposal notes, not implemented.
**Feature name**: BOLDZip-SR.
**Long name**: Split-Reliable Graph Carrier Compression.
**Target artifact**: detailed source notes sufficient to generate a future `PRD.json`.
**Source**: user-provided proposal text on 2026-05-09.
**Repo fit**: `fmrilatent` already represents fMRI data as basis/loadings/offset and includes graph, wavelet, temporal basis, handle, and encode/decode infrastructure. BOLDZip-SR is a higher-level codec proposal that would likely sit above or beside the existing `LatentNeuroVec` and `encode()` systems rather than replacing them.

---

## 1. PRD Identity

Suggested JSON fields:

```json
{
  "id": "boldzip_sr",
  "name": "BOLDZip-SR",
  "title": "Split-Reliable Graph Carrier Compression",
  "status": "proposal",
  "domain": "fMRI compression",
  "package": "fmrilatent",
  "primary_goal": "Compress preprocessed fMRI into carrier signals, high-resolution spatial texture loadings, and sparse reliable innovation events while preserving biologically reproducible signal."
}
```

Short description:

BOLDZip-SR is a subject-adaptive, training-light fMRI codec. It stops treating each parcel time series as the stored object. Instead, it stores a small bank of temporally compressed whole-brain carrier signals, a high-resolution but mostly time-independent sparse spatial texture map describing how fine spatial details ride on those carriers, and a sparse event stream for reliable deviations.

One-line product claim:

At matched or lower bitrate than parcel time series, BOLDZip-SR should reconstruct full-resolution voxel or grayordinate time series while preserving reliable functional connectivity, task betas, event timing, and subject-specific spatial detail better than parcel averages.

---

## 2. Problem Statement

Current parcel-average representations store `P x T` values and discard within-parcel spatial structure. This is a spatial low-pass projection with little or no temporal compression. The reconstruction is blocky and parcel-constant, so subject-specific fine-scale spatial features, boundary shifts, gradients, and focal activations are removed even when they are reliable.

Generic low-rank methods compress time and space but can blur high-spatial-frequency structure. Generic wavelet or DCT codecs do not allocate bits according to fMRI-specific reliability or downstream neuroscience objectives. Per-subject implicit neural representations may compress fMRI but require neural optimization and may be slow to encode or decode.

BOLDZip-SR addresses the gap by separating fMRI signal into:

- carriers: few global or network-level temporal degrees of freedom;
- texture: high-resolution spatial detail stored mostly once as sparse loadings on those carriers;
- events: rare, reliable deviations from the carrier-texture model;
- noise: unreliable variance that should be coarsely quantized or dropped.

---

## 3. Core Hypothesis

Fine spatial detail in fMRI usually does not require an independent full time series per detail atom. Many high-resolution spatial atoms can be represented as static signed mixtures, local loadings, or small-lag versions of a small set of whole-brain carrier signals.

Therefore, instead of storing `D x T` fine-detail time series, the codec can store roughly `D x q` static loading parameters with `q << T`, plus sparse exceptions. Split-half reliability should drive which detail atoms, coefficients, temporal basis coefficients, and events receive bits.

The expected scaling target is:

```text
payload_size ~= K * M + q * D_kept + S
```

where:

- `K`: number of carrier signals;
- `M`: number of temporal coefficients per carrier after temporal compression;
- `D_kept`: number of reliable fine spatial detail atoms kept;
- `q`: average number of carrier loadings per detail atom;
- `S`: number of sparse residual events.

---

## 4. Model Specification

Input:

```text
X in R^(V x T)
```

where:

- `V`: number of voxels or grayordinates;
- `T`: number of time points;
- `X`: preprocessed subject fMRI matrix.

Spatial graph basis:

```text
Phi = [Phi_c, Phi_d]
```

where:

- `Phi_c`: coarse scaling functions, roughly parcel-like;
- `Phi_d`: fine graph-wavelet/detail atoms;
- basis should be invertible or reconstructive enough for codec synthesis;
- basis should preserve locality, graph boundaries, cortical sheet adjacency, and subcortical volumetric adjacency.

Codec model:

```text
X ~= mu + Phi_c C + Phi_d A Z + Phi_d R
```

Terms:

- `mu`: mean or baseline map, needed for raw intensity reconstruction or optional centering reversal;
- `C`: coarse coefficient matrix for parcel-like signal, shape approximately `n_coarse x T`;
- `Z`: carrier time courses, shape `K x T`;
- `A`: sparse spatial texture loading matrix, shape `D x K`;
- `R`: sparse residual event or innovation stream in fine-detail coefficient space;
- `Phi_d A Z`: high-resolution spatial texture riding on carrier dynamics;
- `Phi_d R`: sparse reliable deviations not captured by carrier-texture reconstruction.

Important nuance:

The coarse signal is not the final representation. It is the source or support for learning carrier dynamics. High-resolution signal is reintroduced through sparse texture loadings and sparse events.

---

## 5. Data Product And Payload Components

A future `PRD.json` should describe a BOLDZip-SR encoded object as a structured payload with at least these components:

```json
{
  "payload_components": [
    "header",
    "graph_basis_descriptor",
    "temporal_basis_descriptor",
    "baseline_map",
    "coarse_coefficients_or_carrier_fit_metadata",
    "carrier_temporal_coefficients",
    "texture_loadings",
    "innovation_events",
    "quantization_metadata",
    "reliability_metadata"
  ]
}
```

Required payload details:

- `header`: codec version, dimensions, spatial domain type, temporal length, preprocessing assumptions, endian/encoding details if a real bitstream is implemented.
- `graph_basis_descriptor`: graph construction settings, basis family, levels, atom ordering, coarse/detail split, surface or grayordinate support metadata, any deterministic seeds.
- `temporal_basis_descriptor`: DCT, lapped cosine, wavelet, spline, HRF-aware atom set, basis length, coefficient ordering.
- `baseline_map`: optional `mu`; may be absent for centered data.
- `carrier_temporal_coefficients`: quantized `Theta` for `Z ~= Theta B_t^T`.
- `texture_loadings`: sparse tuples such as `(detail_atom_id, carrier_id, amplitude, lag, quantization_bin, reliability_score)`.
- `innovation_events`: sparse tuples such as `(spatial_atom_or_component_id, t0, duration, amplitude, shape_id, reliability_score)`.
- `quantization_metadata`: step sizes, reliability-to-step mapping, entropy coding tables if applicable.
- `reliability_metadata`: split method, shrinkage assumptions, thresholds, validation scores.

---

## 6. Encoder Pipeline

### Step A: Build graph basis

Inputs:

- subject cortical surface, grayordinate mesh, or voxel mask;
- cortex graph edges following cortical sheet topology;
- subcortex graph edges from local volumetric adjacency;
- optional hybrid cortex plus subcortex support once that infrastructure exists.

Output:

- implicit or explicit multiresolution basis `Phi = [Phi_c, Phi_d]`;
- coarse scaling functions;
- fine detail atoms;
- block iteration support for `Phi_d` so fine coefficients can be processed without materializing full matrices.

Product requirement:

The graph basis is not just a transform choice. It is the coordinate system for the codec. It must be deterministic, locality-preserving, boundary-sensitive, and invertible or near-invertible under the intended reconstruction path.

Candidate implementations in this repo:

- existing graph bridge code for voxel subsets;
- existing heat/diffusion wavelet machinery;
- existing Haar/wavelet-active code for fast wavelet transforms;
- existing hierarchical template machinery for multi-scale parcel templates;
- future CIFTI-aware grayordinate support from `dev/round3_cifti_design.md`.

### Step B: Estimate split-half reliability

Split strategy options:

- odd/even TRs after temporal correction;
- interleaved chunks;
- run halves if there are multiple comparable runs;
- task-design-aware splits to avoid confounding stimulus timing with split identity.

Reliability score option 1:

```text
r_i = max(0, cov(c_i^(1), c_i^(2)) / sqrt(var(c_i^(1)) * var(c_i^(2))))
```

Reliability score option 2:

```text
r_i = sigma_signal_i^2 / (sigma_signal_i^2 + sigma_noise_i^2)
```

Purpose:

The codec should not spend bits on non-reproducible high-frequency noise. Reliability scoring is the anti-overfitting and bit-allocation rule.

PRD requirement:

The implementation must define reliability separately for carrier candidates, temporal packet coefficients, fine detail atoms/loadings, and residual events. A single reliability API should be able to return scores and metadata for each coefficient family.

### Step C: Learn carrier signals

Projection:

```text
Y_c = Phi_c^T X
```

Carrier learning:

```text
Y_c ~= U_K Z
```

Candidate methods:

- reliability-weighted randomized SVD;
- robust PCA;
- weighted low-rank approximation with shrinkage;
- rank range initially `K = 32..128`.

Carrier output:

- `Z in R^(K x T)`;
- each row is a carrier signal;
- carriers should represent global or network-level temporal degrees of freedom.

Temporal compression:

```text
Z ~= Theta B_t^T
```

where:

- `B_t`: temporal packet basis;
- `Theta`: temporal coefficients;
- keep only reliable or rate-distortion-beneficial coefficients.

Candidate temporal bases:

- DCT;
- lapped cosine;
- temporal wavelets;
- splines;
- HRF-aware atoms;
- existing `fmrilatent` temporal DCT, B-spline, and Slepian/DPSS machinery.

### Step D: Fit high-resolution texture loadings

Fine projection:

```text
Y_d = Phi_d^T X
```

For each detail atom `i`:

```text
Y_d[i, t] ~= sum_{k in K_i} a_ik Z_k(t - tau_ik)
```

where:

- `K_i`: small set of top carriers for detail atom `i`;
- `q = |K_i|`, commonly 1, 2, or 3;
- `a_ik`: static loading;
- `tau_ik`: optional small discrete lag for task or vascular latency.

Selection rule:

Keep a loading only if it predicts the held-out split. Store sparse tuples:

```text
(detail_atom_id, carrier_id, loading, lag)
```

Possible implementation detail:

Process `Phi_d` in blocks. For each block, compute `Yd_block`, score carriers via correlations or sparse regression, validate against held-out split, quantize accepted loading tuples, and append to payload.

### Step E: Encode sparse innovations

Residual:

```text
E = Y_d - A Z
```

or with temporal lag expansion:

```text
E = Y_d - A lagged(Z)
```

Event candidates:

- local graph-time connected components;
- high-amplitude cofluctuation events;
- task-locked residuals;
- sharp residual changes that validate across split halves;
- spatially clustered detail atoms with common timing.

Event tuple:

```text
(spatial_atom_or_component_id, t0, duration, amplitude, shape_id)
```

Role:

Events are the innovation stream, not the whole codec. The method should avoid representing all fMRI as a point process.

### Step F: Noise-shaped quantization

Quantization rule:

```text
Delta_i proportional to sigma_noise_i / sqrt(r_i + epsilon)
```

Expected behavior:

- reliable coefficients receive finer quantization;
- unreliable coefficients receive coarse quantization or are zeroed;
- compression acts partly as denoising rather than blind sample damage.

PRD requirement:

Quantization should be tied to reliability metadata and downstream objective weights. The product spec should expose configurable rate targets or quality targets, not only fixed thresholds.

---

## 7. Decoder Pipeline

Decoder should require no fitting or neural optimization.

Steps:

1. Decode metadata, graph basis descriptor, temporal basis descriptor, baseline, carrier temporal coefficients, sparse texture loadings, and sparse events.
2. Reconstruct carriers:

   ```text
   Z_hat = Theta_hat B_t^T
   ```

3. Reconstruct fine detail coefficients:

   ```text
   Y_d_hat = A_hat Z_hat + R_hat
   ```

4. Reconstruct full data:

   ```text
   X_hat = mu_hat + Phi_c C_hat + Phi_d Y_d_hat
   ```

5. Return either:

   - a full matrix `V x T`;
   - a `LatentNeuroVec`-like object;
   - an implicit object that supports block/time/ROI decoding without full materialization.

Implementation priority:

The first version should prefer an implicit decoder object. Full `V x T` materialization is useful for evaluation but should not be required for normal use.

---

## 8. Compression Budget Example

Example values:

```text
V = 90,000
T = 1,200
full_data_samples = 108,000,000
```

Parcel average baseline:

```text
P = 400
parcel_samples = 400 * 1,200 = 480,000
```

BOLDZip-SR target:

```text
K = 96
M = 64
carrier_coefficients = 96 * 64 = 6,144
D_kept = 20,000
q = 1.5
texture_loadings = 20,000 * 1.5 = 30,000
S = 5,000
event_payload = 5,000
total_scalarish_payload = 6,144 + 30,000 + 5,000 = 41,144
```

Approximate scalar-count comparison:

```text
480,000 / 41,144 ~= 11.7
```

Interpretation:

This is about 11.7x fewer scalar-ish values than a 400-parcel representation while reconstructing full-resolution signal instead of parcel-constant signal. Actual byte ratio will depend on index coding, quantization, entropy coding, headers, baseline storage, graph basis metadata, and event representation.

PRD caveat:

The scalar budget is a target hypothesis, not a proven guarantee. The first experiment should report actual bytes, scalar counts, and reconstruction metrics separately.

---

## 9. Objective Function

Do not optimize ordinary voxelwise MSE alone because it rewards preserving scanner noise. The rate-distortion objective should emphasize reliable biological signal.

Suggested loss:

```text
L =
  lambda_1 * ||W_r o (X - X_hat)||_F^2
+ lambda_2 * ||FC(X) - FC(X_hat)||_F^2
+ lambda_3 * ||beta_task(X) - beta_task(X_hat)||^2
+ lambda_4 * d_events(X, X_hat)
+ lambda_5 * bits
```

where:

- `W_r`: split-half reliability weighting;
- `FC`: functional connectivity summary;
- `beta_task`: GLM beta or contrast map for task data;
- `d_events`: event timing, cofluctuation, or burst preservation distance;
- `bits`: actual encoded payload size or a differentiable proxy.

Resting-state priorities:

- full-resolution FC preservation;
- spectra and autocorrelation preservation;
- cofluctuation event timing;
- test-retest fingerprinting;
- downstream prediction performance;
- held-out reliable coefficient reconstruction.

Task priorities:

- GLM beta map preservation;
- contrast preservation;
- stimulus-locked residual preservation;
- latency structure;
- focal activation preservation;
- task decoding or downstream prediction.

---

## 10. Baselines For First Experiment

Compare at matched bitrates or matched storage budgets:

1. 400-parcel average.
2. 1,000-parcel average.
3. PCA/SVD rank-K.
4. Generic 4D wavelet/DCT compression.
5. Graph-wavelet thresholding.
6. INR-style per-subject compressor.
7. BOLDZip-SR.

Primary comparison claim:

BOLDZip-SR should not be judged only by raw MSE against the original noisy matrix. It should beat parcel baselines on full-resolution, biologically weighted reconstruction at substantially lower rate.

---

## 11. Evaluation Metrics

Suggested PRD JSON fields:

```json
{
  "metrics": {
    "compression": [],
    "reconstruction": [],
    "resting_state": [],
    "task": [],
    "runtime": [],
    "robustness": []
  }
}
```

Compression metrics:

- actual bytes on disk;
- scalar payload count;
- compression ratio versus full data;
- compression ratio versus 400-parcel and 1,000-parcel baselines;
- payload breakdown by carriers, texture, events, metadata, baseline.

Reconstruction metrics:

- `corr(X, X_hat)` on held-out reliable coefficients;
- reliability-weighted MSE;
- unweighted MSE reported as secondary;
- spatially stratified error by coarse/detail scale;
- temporal spectra preservation;
- ROI and full-resolution reconstruction error.

Resting-state metrics:

- full-resolution FC preservation;
- parcel FC preservation for common atlases;
- edge cofluctuation/event timing preservation;
- test-retest fingerprinting;
- subject identification;
- downstream behavioral prediction if available.

Task metrics:

- beta-map correlation;
- contrast-map correlation;
- task effect localization;
- event/stimulus-locked residual preservation;
- latency estimates if lags are enabled.

Runtime metrics:

- encode time;
- decode time;
- peak memory;
- block-processing throughput;
- full reconstruction throughput;
- ROI/time-slice decode throughput.

Robustness metrics:

- behavior under low SNR;
- behavior under different split strategies;
- sensitivity to graph basis family;
- sensitivity to `K`, `M`, `q`, reliability thresholds, event thresholds;
- reproducibility across random seeds.

---

## 12. Minimum Viable Product

MVP goal:

Implement an experimental BOLDZip-SR prototype that can encode and decode a preprocessed grayordinate or matrix dataset, then evaluate it against parcel and SVD baselines at matched scalar budgets.

MVP scope:

- matrix input `X` with `V x T` orientation explicitly checked;
- graph/basis object supplied by caller or built with existing fmrilatent graph/wavelet utilities;
- deterministic split-half reliability function;
- coarse projection and randomized SVD carrier learning;
- DCT-based temporal compression for carriers;
- sparse per-detail texture loading selection with `q` top carriers and optional no-lag first version;
- sparse residual event detector with a simple reliable-threshold rule;
- decoder that returns matrix reconstruction and optional implicit object;
- evaluation script producing scalar counts and core metrics.

MVP non-goals:

- production entropy-coded binary file format;
- neural INR baseline implementation inside fmrilatent;
- full CIFTI file I/O;
- multi-subject dictionary learning;
- optimized C++ kernels for all phases;
- perfect reconstruction of raw intensity;
- clinical-grade archival compression claims.

Recommended first simplifications:

- no lag in first texture model (`tau_ik = 0`);
- DCT temporal basis first;
- reliability-weighted hard thresholding before advanced rate-distortion optimization;
- simple event atoms before connected graph-time component coding;
- matrix-based evaluation before new S4 class design.

---

## 13. Acceptance Criteria

Prototype acceptance criteria:

- Given a small synthetic or real matrix `X`, encoder returns a structured BOLDZip-SR object with carrier coefficients, texture loadings, events, and metadata.
- Decoder reconstructs `X_hat` with dimensions identical to `X`.
- The encoded object reports a payload breakdown with carrier, texture, event, metadata, and baseline sizes.
- Reliability scoring is deterministic under fixed seed and split strategy.
- Loadings are only kept if held-out split validation passes the configured threshold.
- Events are only kept if residual event validation passes the configured threshold.
- Baseline comparison script can evaluate 400-parcel, 1,000-parcel, SVD, graph-threshold, and BOLDZip-SR outputs on the same metric set.
- The first benchmark reports actual bytes and scalar counts, not only theoretical formula counts.
- Unit tests cover encode/decode dimensions, empty events, zero retained details, full retained details on toy data, quantization edge cases, and deterministic seeds.

Research acceptance criteria:

- At one or more matched bitrates, BOLDZip-SR preserves held-out reliable fine-detail coefficients better than parcel averages.
- At comparable or lower payload size than 400-parcel time series, BOLDZip-SR reconstructs full-resolution signal with better reliability-weighted reconstruction than parcel expansion.
- For resting-state data, FC preservation and event timing preservation are not worse than parcel baselines at the target bitrate.
- For task data, beta or contrast preservation is not worse than parcel baselines at the target bitrate.
- Encode/decode runtime is practical enough for iterative experimentation on one subject without neural fitting loops.

---

## 14. Suggested API Surface

Possible high-level API:

```r
fit <- boldzip_sr_encode(
  X,
  graph = graph,
  spatial_basis = NULL,
  k_carriers = 96,
  temporal_basis = spec_time_dct(k = 64),
  q_texture = 2,
  split = c("odd_even", "chunks", "runs"),
  reliability = boldzip_reliability(),
  quantization = boldzip_quantization(),
  events = boldzip_events(),
  rate_target = NULL,
  quality_target = NULL,
  seed = 1
)

X_hat <- boldzip_sr_decode(fit)
summary <- boldzip_sr_payload_summary(fit)
metrics <- evaluate_boldzip_sr(X, X_hat, baselines = ...)
```

Possible object classes:

- `BoldZipSR`: encoded object or S3 list;
- `BoldZipSRBasis`: graph basis descriptor and operators;
- `BoldZipSRReliability`: split strategy and reliability scores;
- `BoldZipSRPayload`: encoded carrier, texture, event, quantization metadata;
- `BoldZipSRDecoder`: implicit decoder for block/time/ROI reconstruction.

Possible integration with existing classes:

- map carrier temporal basis and spatial texture loadings into a `LatentNeuroVec`-like representation where practical;
- use `BasisHandle` and `LoadingsHandle` for deferred materialization;
- use existing `encode()` specs for temporal basis choices;
- add experimental functions under a separate file before public export;
- keep full production binary codec out of the initial package API until validated.

---

## 15. Implementation Phases

### Phase 0: Design and simulation harness

Deliverables:

- synthetic data generator with known carriers, sparse spatial texture, events, and noise;
- reliability split utilities;
- metric functions;
- design review of how current `fmrilatent` basis/handle infrastructure maps to BOLDZip-SR.

Exit criteria:

- synthetic generator can recover known structure in small examples;
- objective metrics run without large neuroimaging dependencies.

### Phase 1: Matrix prototype

Deliverables:

- matrix input encoder;
- supplied graph basis or simple wavelet basis;
- carrier SVD;
- DCT carrier compression;
- sparse texture regression with no lag;
- simple residual event encoder;
- decoder to matrix.

Exit criteria:

- small toy tests pass;
- scalar budget and reconstruction metrics are generated.

### Phase 2: Graph/grayordinate basis integration

Deliverables:

- block-wise detail atom processing;
- graph wavelet or lifting basis integration;
- support for grayordinate-like arrays;
- no full materialization of huge `Phi_d`.

Exit criteria:

- one realistic subject/run can be encoded in bounded memory.

### Phase 3: Reliability-aware rate-distortion

Deliverables:

- reliability-weighted quantization;
- threshold/rate search;
- validation split reporting;
- payload quality presets.

Exit criteria:

- compression can target a scalar or byte budget;
- results include reliability-weighted metric improvements over baselines.

### Phase 4: Evaluation suite

Deliverables:

- baseline implementations;
- matched-bitrate evaluation;
- resting/task metric adapters;
- reproducible benchmark report.

Exit criteria:

- results decide whether BOLDZip-SR warrants productizing.

### Phase 5: Package productization

Deliverables:

- stable API;
- tests;
- docs/vignette;
- optional binary serialization;
- performance optimization and possible C++ kernels.

Exit criteria:

- documented experimental feature or exported stable feature depending on results.

---

## 16. Risks And Mitigations

Risk: fine detail does require independent dynamics in some datasets.

Mitigation: allow event stream and variable `q`; report residual energy by detail scale; compare against graph-wavelet thresholding and SVD.

Risk: split-half reliability underestimates task-locked or nonstationary effects.

Mitigation: support task-aware splits, run-wise splits, and design-aware reliability scoring.

Risk: carrier bank learned from coarse coefficients misses high-resolution local dynamics.

Mitigation: optionally include a small set of reliable fine-detail candidates in carrier learning, or add carrier refinement after initial texture fit.

Risk: payload indices and metadata erase scalar-count advantage.

Mitigation: report actual bytes early; delta-code sorted atom IDs; keep first implementation honest about headers and metadata.

Risk: graph basis construction is expensive or unstable.

Mitigation: cache basis descriptors; support deterministic seeds; start with supplied bases; add block operators before production use.

Risk: raw MSE looks worse than generic codecs.

Mitigation: define success around reliability-weighted and neuroscience-weighted metrics while still reporting raw MSE transparently.

Risk: quantization introduces bias into downstream GLMs or FC.

Mitigation: evaluate downstream metrics directly; include task beta and FC preservation in acceptance criteria.

Risk: too many API concepts for users.

Mitigation: expose a high-level `boldzip_sr_encode()` with presets; keep internals inspectable but not required for basic use.

---

## 17. Open Questions

Scientific/model questions:

- How stable is the assumption that fine detail can be modeled as static loadings on global carriers?
- Does lagged texture improve task/vascular latency preservation enough to justify extra payload?
- Should carriers be learned only from coarse functions, or from reliability-filtered full graph coefficients?
- What is the best split strategy for resting-state versus task fMRI?
- How should event shape dictionaries be learned or fixed?
- Does the event stream improve FC/event metrics after carrier-texture fitting, or is it mostly overhead?

Engineering questions:

- Should the first encoded object be an S3 list, S4 class, or an extension of existing latent objects?
- Should graph basis descriptors be serializable before binary payload work begins?
- Which existing `fmrilatent` wavelet or graph basis path is closest to the needed invertible detail/coarse split?
- How should block-wise `Phi_d^T X` be exposed?
- Should quantization be implemented as actual integer bins in the prototype or simulated by coefficient thresholding first?
- Where should CIFTI/grayordinate support land relative to the hybrid support design in `dev/round3_cifti_design.md`?

Product questions:

- Is the initial user a researcher comparing codecs, or an analyst needing a compact exchange/storage format?
- Is success defined by scalar budget, actual byte budget, or downstream analytic preservation?
- Should BOLDZip-SR be an experimental research module or a package-level flagship feature?
- What datasets should be named as required first benchmarks?

Evidence questions:

- Literature claims in the proposal should be verified and cited before PRD finalization:
  - parcellation choice affects FC and individual-difference interpretation;
  - ROI averaging ignores within-ROI heterogeneity;
  - graph/spectral wavelets help represent fMRI on brain geometry;
  - low-rank/temporal-subspace constraints work for accelerated fMRI reconstruction;
  - high-amplitude BOLD events preserve important voxelwise FC information;
  - INR-style fMRI compression exists and has the stated tradeoffs;
  - lossy fMRI compression has prior research support.

---

## 18. PRD JSON Blueprint

A future `PRD.json` could use this top-level shape:

```json
{
  "id": "boldzip_sr",
  "title": "BOLDZip-SR: Split-Reliable Graph Carrier Compression",
  "status": "proposal",
  "summary": "",
  "problem": "",
  "goals": [],
  "non_goals": [],
  "users": [],
  "core_hypothesis": "",
  "model": {
    "input": "X in R^(V x T)",
    "basis": "Phi = [Phi_c, Phi_d]",
    "reconstruction": "X ~= mu + Phi_c C + Phi_d A Z + Phi_d R",
    "components": {}
  },
  "payload": {
    "formula": "K*M + q*D_kept + S",
    "components": [],
    "metadata": []
  },
  "encoder": {
    "steps": []
  },
  "decoder": {
    "steps": []
  },
  "objective": {
    "loss": "",
    "resting_state_terms": [],
    "task_terms": []
  },
  "baselines": [],
  "metrics": {},
  "mvp": {
    "scope": [],
    "non_goals": [],
    "acceptance_criteria": []
  },
  "phases": [],
  "risks": [],
  "open_questions": [],
  "evidence_to_verify": [],
  "repo_integration": {
    "likely_files": [],
    "existing_systems": [],
    "test_strategy": []
  }
}
```

---

## 19. Repo Integration Notes

Existing systems that appear relevant:

- `LatentNeuroVec` model: existing `basis x loadings + offset` representation is conceptually aligned with carrier/texture reconstruction.
- `BasisHandle` and `LoadingsHandle`: useful for avoiding eager materialization of large carrier or spatial basis objects.
- `encode()` specs: existing temporal and spatial basis specs can inform high-level API style.
- DCT/B-spline/Slepian temporal bases: candidate `B_t` implementations.
- heat/diffusion/active/Haar wavelets and graph bridge: candidate graph detail basis infrastructure.
- hierarchical templates and parcel basis: candidate `Phi_c` and multi-scale basis starting points.
- transport/implicit latent machinery: possible route for decoding without full matrix materialization.
- `dev/round3_cifti_design.md`: important if first real experiment uses grayordinate arrays.

Likely implementation files if productized later:

- `R/boldzip_sr.R`: high-level encode/decode API and object constructors.
- `R/boldzip_sr_reliability.R`: split-half scoring and validation.
- `R/boldzip_sr_texture.R`: sparse texture regression.
- `R/boldzip_sr_events.R`: residual event detection and coding.
- `R/boldzip_sr_quantization.R`: reliability-aware quantization.
- `R/boldzip_sr_metrics.R`: evaluation metrics and payload summaries.
- `tests/testthat/test-boldzip-sr*.R`: prototype and edge-case tests.
- `vignettes/boldzip-sr.Rmd`: only after prototype results justify docs.

Testing strategy:

- synthetic ground-truth tests for carrier recovery, texture loading recovery, and event recovery;
- dimension and reconstruction tests;
- deterministic seed tests;
- empty/degenerate payload tests;
- split-half validation behavior tests;
- payload accounting tests;
- small integration benchmark with a mock graph basis.

---

## 20. Minimal Pseudocode For PRD

Encoder:

```text
X1, X2 = split_halves(X)
Phi_c, Phi_d = build_or_load_graph_basis(G)

Y_c = Phi_c^T X
U, Z = randomized_reliable_svd(Y_c, X1, X2, rank = K)

Theta = temporal_packet_transform(Z, basis = B_t)
Theta_hat = quantize_keep_reliable(Theta)
Z_hat_for_fit = inverse_temporal_packet(Theta_hat, B_t)

payload_A = []
for block in Phi_d.blocks():
    Yd_block = block^T X
    loadings = split_reliable_sparse_regression(
        Yd_block,
        Z_hat_for_fit,
        q = q,
        lags = allowed_lags
    )
    payload_A.append(quantize(loadings))

Yd_pred = sparse_A_times_Z(payload_A, Z_hat_for_fit)
Residual = Phi_d^T X - Yd_pred
events = encode_reliable_graph_time_events(Residual, X1, X2)

bitstream_or_object = encode_payload(
    graph_basis = graph_basis_descriptor,
    temporal_basis = temporal_basis_descriptor,
    theta = Theta_hat,
    texture = payload_A,
    events = events,
    reliability = reliability_metadata,
    quantization = quantization_metadata
)
```

Decoder:

```text
payload = decode_payload(bitstream_or_object)
Phi_c, Phi_d = rebuild_graph_basis(payload.graph_basis_descriptor)
B_t = rebuild_temporal_basis(payload.temporal_basis_descriptor)

Z_hat = payload.Theta_hat B_t^T
Yd_hat = payload.A_hat Z_hat + payload.R_hat
X_hat = payload.mu_hat + Phi_c C_hat + Phi_d Yd_hat
```

---

## 21. Bottom-Line PRD Framing

BOLDZip-SR is sharper than parcellation because it does not average away detail. It stores temporal degrees of freedom once as carriers, stores high-resolution spatial detail as cheap sparse texture loadings, and reserves additional bits for rare reliable events. Split-half reliability is the central policy: preserve reproducible biological signal, not raw variance.

The proposal should be treated first as an experimental codec and evaluation program. The decisive question is empirical: at matched actual byte rates, does carrier-plus-texture-plus-event reconstruction preserve reliable full-resolution fMRI information better than parcel averages, low-rank compression, and graph-wavelet thresholding?
