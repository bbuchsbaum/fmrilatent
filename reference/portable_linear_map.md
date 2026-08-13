# Portable linear-map contract

The portable linear map is the seam \`fmrilatent\` uses to consume
subject field operators from neurofunctor (or any other producer)
without pulling in file-format-specific warp code. A portable linear map
represents a single linear operator \\A : \mathbb{R}^{n\_\text{source}}
\to \mathbb{R}^{n\_\text{target}}\\ together with enough metadata to
compose, adjoint, and (optionally) materialize it.

## Details

An object satisfies the contract if it is either

1.  a base `matrix` or
    [`Matrix::Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html), or

2.  a `list` with the following fields:

    - `n_source`, `n_target`:

      Positive integer scalar dimensions of the source and target sample
      vectors.

    - `forward(x, ...)`:

      Function applying the operator to source-space data. Accepts a
      vector of length `n_source` or a matrix with `n_source` rows;
      returns a vector/matrix with `n_target` rows.

    - `adjoint_apply(y, ...)` (alias `adjoint`):

      Function applying the discrete adjoint (or a declared alternative)
      to target-space data.

    - `source_domain_id`, `target_domain_id`:

      Stable character identifiers for the source and target domains.
      Used for composition safety and provenance digests.
      `.compose_linear_maps()` requires matching domain ids across the
      junction unless at least one side is an empty string
      ("unspecified").

    - `source_support`, `target_support`:

      Optional descriptors of the sample layout on each side (for
      example a `LogicalNeuroVol`, a 3D logical mask, surface vertex
      indices, or a vector of sample ids). When present these let
      `fmrilatent` reconstruct domain-aware outputs; when absent it
      falls back to `provenance$target_support` /
      `provenance$target_mask`.

    - `adjoint_convention`:

      Character tag identifying the adjoint convention. Defaults to
      `"euclidean_discrete"` (the discrete transpose under the standard
      Euclidean inner product). Covariance pushforward paths such as
      [`decode_covariance()`](decode_covariance.md) and
      `.project_covariance_diag()` currently require this value;
      operators using another convention must be re-normalized before
      being handed to the pushforward code.

    - `provenance`:

      Optional named list of provenance metadata. `fmrilatent`
      interprets a small set of reserved keys:

      `target_mask`, `target_support`, `target_domain`

      :   Legacy fallback locations read by
          `.resolve_transport_target_support()` when the top-level
          fields are absent. New producers should set the top-level
          fields instead.

      `source_support`

      :   Legacy fallback for the top-level field of the same name.

      `coordinates`

      :   Appended with value `"raw"` when the map has been wrapped with
          `.transform_linear_map_coordinates()`.

      All other keys are advisory and carried through composition
      verbatim. Producers should prefer top-level fields over provenance
      entries for anything fmrilatent actually reads.

    - `materialize(...)` (optional):

      Returns the dense `n_target x n_source` matrix representation.
      Used only when a caller explicitly asks for it; the main encode
      and decode paths are matrix-free and only call
      `forward`/`adjoint_apply`.

    - `contract_version`:

      Integer scalar identifying the contract revision. Auto-filled by
      [`as_portable_linear_map()`](as_portable_linear_map.md) to `1L` in
      this release. Future revisions will bump this tag so producers and
      consumers can negotiate compatibility.

**Producer-form vs canonical form on `TransportLatent`.**
[`transport_latent()`](transport_latent.md) records the raw producer
object at `x$field_operator` (and its legacy alias
`x$observation_operator`) as a back-reference for debugging and
round-trip. The canonical portable-linear-map representation used
internally by all encode/decode math is at `x$transport$field_operator`.
Consumers that need to introspect the operator should always read the
canonical form.

Producers (such as `neurofunctor::compile_observation_operator`) can
return either a bare callback list or a wrapped object. Call
[`as_portable_linear_map()`](as_portable_linear_map.md) to coerce any
accepted input to the canonical contract and
[`validate_portable_linear_map()`](validate_portable_linear_map.md) to
verify that an input already satisfies it without coercing.
