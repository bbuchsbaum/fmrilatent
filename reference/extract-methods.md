# Extract Elements from LatentNeuroVec

Extract elements from a `LatentNeuroVec`. Use `[[` to extract a single
volume as a `SparseNeuroVol`, or `[` for 4D subsetting.

## Usage

``` r
# S4 method for class 'LatentNeuroVec,numeric'
x[[i]]

# S4 method for class 'LatentNeuroVec,numeric,numeric,ANY'
x[i, j, k, l, ..., drop = TRUE]

# S4 method for class 'LatentNeuroVec,ANY,ANY,ANY'
x[i, j, k, l, ..., drop = TRUE]
```

## Arguments

- x:

  A [`LatentNeuroVec-class`](LatentNeuroVec-class.md) object.

- i:

  Numeric index for first dimension (x-axis) or volume index for `[[`.

- j:

  Numeric index for second dimension (y-axis).

- k:

  Numeric index for third dimension (z-axis).

- l:

  Numeric index for fourth dimension (time).

- ...:

  Additional arguments (unused).

- drop:

  Logical, whether to drop dimensions (default TRUE).

## Value

For `[[`: A
[`SparseNeuroVol-class`](https://bbuchsbaum.github.io/neuroim2/reference/SparseNeuroVol-class.html)
containing the computed volume. For `[`: An array of extracted values.

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract volumes
vol1 <- lvec[[1]]
vol_mid <- lvec[[dim(lvec)[4] / 2]]
vol_last <- lvec[[dim(lvec)[4]]]
} # }
```
