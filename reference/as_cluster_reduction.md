# Convert a ClusteredNeuroVol to a ClusterReduction

Bridges the
[`neuroim2::ClusteredNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.html)
parcellation representation to fmrilatent's `ClusterReduction` class,
preserving label metadata.

## Usage

``` r
as_cluster_reduction(cvol)
```

## Arguments

- cvol:

  A `ClusteredNeuroVol` object (from neuroim2).

## Value

A `ClusterReduction` object with label metadata in `info`.
