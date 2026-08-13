# LatentNeuroSurfaceVector Class

An explicit latent representation for surface-domain neuroimaging data.
The data are stored as a temporal basis matrix and a surface loadings
matrix over a supported set of vertices on a `neurosurf` geometry.

## Slots

- `basis`:

  A `Matrix` or `BasisHandle` with dimensions (nTime x k).

- `loadings`:

  A `Matrix` or `LoadingsHandle` with dimensions (nVerticesInSupport x
  k).

- `offset`:

  Optional numeric vector of length `nVerticesInSupport`.

- `geometry`:

  Surface domain object, typically a
  [`neurosurf::SurfaceGeometry`](https://bbuchsbaum.github.io/neurosurf/reference/SurfaceGeometry.html)
  or `neurosurf::SurfaceSet`.

- `support`:

  Integer vector of vertex indices within `geometry`.

- `label`:

  Character label.

- `meta`:

  Lightweight metadata list.
