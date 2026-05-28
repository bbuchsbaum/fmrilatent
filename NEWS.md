# fmrilatent 0.1.0

## Development

* Fixed active-pencil wavelet round-trips for sparse, non-contiguous masks by
  lifting compact active runs rather than padded full-grid pencils.
* Fixed HRBF tiny-ROI sparse fallback construction for non-contiguous masks.
* Clamped inferred temporal Slepian handle ranks to the valid `[1, n_time]`
  range.
* Preserved one-column matrix shape in `decode_coefficients(wrap = "none")`
  when the caller supplied a matrix.
* Derived hierarchical-template and spatial PCA component counts from returned
  eigensolver columns rather than requested ranks.
* Guarded bilateral surface wrapping by checking the `BilatNeuroSurfaceVector`
  class before construction.
