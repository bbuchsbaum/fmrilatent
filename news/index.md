# Changelog

## fmrilatent 0.2.0

- Added layered [`latent_space_id()`](../reference/latent_space_id.md)
  records that distinguish analysis-coordinate identity, decoder-domain
  identity, and decoder support while preserving a compatibility
  identifier for downstream consumers.

- Added immutable
  [`latent_units_record()`](../reference/latent_units_record.md)
  declarations and [`latent_units()`](../reference/latent_units.md)
  accessors. Legacy representations fail closed with explicitly
  undeclared units instead of receiving inferred scientific semantics.

- Made decoder-domain metadata explicit so integration layers can
  validate spatial identity from the latent object rather than
  placeholder transport labels.

- Updated decoded surface matrices to the current `neurosurf`
  full-domain row contract while retaining the ordered support indices.

- Corrected the `rgsp` remote source and made the hierarchical Laplacian
  solver use its symmetric smallest-magnitude path for stable low modes.

- Fixed active-pencil wavelet round-trips for sparse, non-contiguous
  masks by lifting compact active runs rather than padded full-grid
  pencils.

- Fixed HRBF tiny-ROI sparse fallback construction for non-contiguous
  masks.

- Clamped inferred temporal Slepian handle ranks to the valid
  `[1, n_time]` range.

- Preserved one-column matrix shape in
  `decode_coefficients(wrap = "none")` when the caller supplied a
  matrix.

- Derived hierarchical-template and spatial PCA component counts from
  returned eigensolver columns rather than requested ranks.

- Guarded bilateral surface wrapping by checking the
  `BilatNeuroSurfaceVector` class before construction.
