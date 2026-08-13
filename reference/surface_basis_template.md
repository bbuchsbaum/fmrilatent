# Build a shared surface basis template

Build a shared surface basis template

## Usage

``` r
surface_basis_template(
  geometry,
  loadings,
  support = NULL,
  roughness = NULL,
  measure = NULL,
  ridge = 1e-08,
  label = "surface_basis",
  meta = list()
)
```

## Arguments

- geometry:

  A
  [`neurosurf::SurfaceGeometry`](https://bbuchsbaum.github.io/neurosurf/reference/SurfaceGeometry.html)
  or `neurosurf::SurfaceSet`.

- loadings:

  A matrix-like decoder basis with rows aligned to `support` and columns
  aligned to latent coefficients.

- support:

  Surface support as vertex indices or a logical vector over all surface
  nodes. Defaults to the full surface.

- roughness:

  Optional coefficient-space roughness matrix.

- measure:

  Optional support-aligned weighting or mass information.

- ridge:

  Small diagonal ridge added to the Gram matrix if needed.

- label:

  Optional label stored in metadata.

- meta:

  Optional additional metadata.

## Value

A `SurfaceBasisTemplate` object.
