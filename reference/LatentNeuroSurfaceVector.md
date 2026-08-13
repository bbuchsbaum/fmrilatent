# Construct a LatentNeuroSurfaceVector

Construct a LatentNeuroSurfaceVector

## Usage

``` r
LatentNeuroSurfaceVector(
  basis,
  loadings,
  geometry,
  support = NULL,
  offset = NULL,
  label = "",
  meta = list()
)
```

## Arguments

- basis:

  Temporal basis matrix (`n_time x k`) or `BasisHandle`.

- loadings:

  Surface loadings matrix (`n_support x k`) or `LoadingsHandle`.

- geometry:

  A
  [`neurosurf::SurfaceGeometry`](https://bbuchsbaum.github.io/neurosurf/reference/SurfaceGeometry.html)
  or `neurosurf::SurfaceSet`.

- support:

  Surface support as vertex indices or a logical vector over all
  vertices.

- offset:

  Optional numeric vector of length `n_support`.

- label:

  Optional label.

- meta:

  Optional metadata list.

## Value

A `LatentNeuroSurfaceVector`.
