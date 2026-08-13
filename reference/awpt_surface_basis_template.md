# Build a surface AWPT basis template

Build a surface AWPT basis template

## Usage

``` r
awpt_surface_basis_template(
  geometry,
  basis_spec = basis_awpt_wavelet(),
  support = NULL,
  loadings = NULL,
  centers = NULL,
  anatomical_operator = NULL,
  conductance = NULL,
  coefficient_roughness = NULL,
  measure = NULL,
  center = FALSE,
  ridge = 1e-08,
  label = "surface_awpt_wavelet"
)
```

## Arguments

- geometry:

  A
  [`neurosurf::SurfaceGeometry`](https://bbuchsbaum.github.io/neurosurf/reference/SurfaceGeometry.html)
  or `neurosurf::SurfaceSet`.

- basis_spec:

  An AWPT basis specification created by
  [`basis_awpt_wavelet()`](basis_awpt_wavelet.md).

- support:

  Surface support as vertex indices or a logical vector over all
  vertices.

- loadings:

  Optional explicit surface decoder loadings.

- centers:

  Optional support-local center indices used for automatic wave-packet
  construction.

- anatomical_operator:

  Optional field-space roughness operator on the supported surface
  domain.

- conductance:

  Optional conductance matrix on the supported surface graph.

- coefficient_roughness:

  Optional coefficient-space roughness matrix.

- measure:

  Optional support-aligned weighting or mass information.

- center:

  Logical; if `TRUE`, center data before projection.

- ridge:

  Small positive ridge added to the Gram diagonal if needed.

- label:

  Optional label stored in metadata.

## Value

A `SurfaceAWPTBasisTemplate`.
