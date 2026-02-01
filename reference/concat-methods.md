# Concatenate LatentNeuroVec Objects

Concatenates two or more `LatentNeuroVec` objects along the temporal
dimension.

## Usage

``` r
# S4 method for class 'LatentNeuroVec,LatentNeuroVec'
concat(x, y, ...)
```

## Arguments

- x:

  First `LatentNeuroVec`.

- y:

  Second `LatentNeuroVec`.

- ...:

  Additional `LatentNeuroVec` objects to concatenate.

## Value

A new `LatentNeuroVec` if all objects are compatible, otherwise a
[`NeuroVecSeq-class`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVecSeq-class.html).

## Examples

``` r
if (FALSE) { # \dontrun{
combined <- concat(lvec1, lvec2, lvec3)
dim(combined)
} # }
```
