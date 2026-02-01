# Derive parent maps for a nested set of parcellations

Derive parent maps for a nested set of parcellations

## Usage

``` r
parent_maps_from_levels(levels)
```

## Arguments

- levels:

  List of integer label vectors (all same length), coarse→fine or
  fine→coarse; order agnostic.

## Value

List of integer named vectors: parents\[\[lvl\]\] maps child ids (names)
at level lvl to parent ids at lvl-1; parents\[\[1\]\] is integer(0).
