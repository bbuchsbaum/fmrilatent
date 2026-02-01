# Cut an hclust into nested label vectors

Utility to turn a hierarchical clustering (on parcels) into a nested set
of parcel label vectors that downstream code can lift into voxel space.
Real data glue: you will pass the hclust built on a parcel-level
similarity graph (e.g., Schaefer-400 parcel graph constructed from
surface geodesics and boundary contacts); then map these parcel labels
back to voxels using the volumetric atlas.

## Usage

``` r
cut_hclust_nested(hc, k_levels)
```

## Arguments

- hc:

  hclust object (e.g., from spectral_ward_hclust()).

- k_levels:

  Integer vector of cluster counts, ordered coarse→fine (e.g., c(16, 32,
  64, 400)).

## Value

List of integer label vectors (same length as hc\$labels), one per
level, nested by construction.
