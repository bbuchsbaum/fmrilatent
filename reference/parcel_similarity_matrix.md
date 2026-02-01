# Build a similarity matrix for parcel clustering (Schaefer-like)

Combines geometric and functional priors into a parcel–parcel similarity
matrix W used to cluster Schaefer parcels into coarser levels. Real data
glue: - \`boundary_contact\`: from surface meshes (fsaverage), compute
fraction of shared boundary between parcels. - \`geo_dist\`: geodesic
distance between parcel centroids on the surface (fallback: Euclidean in
MNI). - \`yeo17\`: network labels from the Schaefer/Yeo annotation
(hemisphere-specific).

## Usage

``` r
parcel_similarity_matrix(
  boundary_contact,
  geo_dist,
  yeo17,
  alpha = 0.5,
  beta = 0.3,
  gamma = 0.2,
  d0 = 30
)
```

## Arguments

- boundary_contact:

  Numeric matrix \[n_parc x n_parc\], symmetric, contact fraction
  (0..1).

- geo_dist:

  Numeric matrix \[n_parc x n_parc\], symmetric geodesic distances (mm).

- yeo17:

  Integer or factor vector of length n_parc with Yeo17 network labels.

- alpha:

  Weight for boundary contact (default 0.5).

- beta:

  Weight for geodesic kernel (default 0.3).

- gamma:

  Weight for Yeo17 agreement (default 0.2).

- d0:

  Scale for geodesic exponential (default 30 mm).

## Value

Similarity matrix W (symmetric, zero diagonal).
