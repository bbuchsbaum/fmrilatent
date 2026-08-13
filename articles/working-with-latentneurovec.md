# Working With LatentNeuroVec Objects

``` r

library(fmrilatent)
```

A `LatentNeuroVec` behaves like a neuroimaging time series, but it
stores a factorization instead of a dense time-by-voxel matrix. This
vignette shows the basic access pattern: encode data, inspect the
factors, pull out a voxel time series, and reconstruct only when you
need the dense matrix.

## What do we start with?

The examples use a small 4 x 4 x 4 mask and a 24-time-point matrix. Rows
are time points and columns are in-mask voxels.

``` r

mask <- array(TRUE, dim = c(4, 4, 4))
mask_vol <- neuroim2::LogicalNeuroVol(mask, neuroim2::NeuroSpace(dim(mask)))

set.seed(7)
X <- matrix(rnorm(24 * sum(mask)), nrow = 24)
dim(X)
#> [1] 24 64
```

## How is the object encoded?

A temporal DCT basis stores shared time functions in
[`basis()`](../reference/basis-methods.md) and voxel weights in
[`loadings()`](../reference/loadings-methods.md).

``` r

lat <- encode(X, spec_time_dct(k = 8), mask = mask_vol, materialize = "matrix")

data.frame(
  time_points = nrow(basis(lat)),
  components = ncol(basis(lat)),
  voxels = nrow(loadings(lat))
)
#>   time_points components voxels
#> 1          24          8     64
```

The dense reconstruction has the same shape as the input, but it is
lossy because eight DCT components cannot represent every 24-point time
course.

``` r

recon <- as.matrix(lat)

data.frame(
  rows = nrow(recon),
  cols = ncol(recon),
  rmse = sqrt(mean((recon - X)^2))
)
#>   rows cols      rmse
#> 1   24   64 0.7947202
```

## How do you read one voxel?

Use [`series()`](../reference/series-methods.md) when you want a time
course for one voxel or a small group of voxels without thinking about
the full matrix.

``` r

voxel_one <- series(lat, 1L)
head(voxel_one)
#> [1]  0.9893044  0.3916249 -0.4362176 -1.0104257 -1.0466904 -0.6113938
```

Coordinate access is available when you know the voxel location in the
image space.

``` r

same_voxel <- series(lat, 1L, 1L, 1L)
all.equal(as.numeric(voxel_one), as.numeric(same_voxel))
#> [1] TRUE
```

## How do you slice it like an array?

A `LatentNeuroVec` behaves like the 4D `(x, y, z, time)` array it stands
in for, so the usual `neuroim2` indexing verbs work — each one
reconstructs only the slice you ask for.

Use `[[` to pull one time point as a 3D volume:

``` r

vol1 <- lat[[1]]
class(vol1)
#> [1] "SparseNeuroVol"
#> attr(,"package")
#> [1] "neuroim2"
dim(vol1)
#> [1] 4 4 4
```

Use `[` with four indices to extract a spatial-by-time sub-array:

``` r

block <- lat[1:2, 1:2, 1:2, 1:5]
dim(block)
#> [1] 2 2 2 5
```

You never materialize the full dense array to get a corner of it; the
slice is reconstructed on demand from the factors.

## When should handles be used?

The default handle mode stores the basis as a lightweight `BasisHandle`
reference instead of a dense matrix, and only materializes it the first
time you ask for it. Matrix mode materializes up front, which is useful
when you plan to reconstruct or access many series repeatedly.

To see the difference you have to look at the stored slot directly:
[`basis()`](../reference/basis-methods.md) *materializes* a handle on
access, so it always hands back a dense matrix in either mode. The slot
is where the two representations diverge.

``` r

lat_handle <- encode(X, spec_time_dct(k = 8), mask = mask_vol,
                     materialize = "handle")

data.frame(
  eager_slot = class(lat@basis)[1],
  handle_slot = class(lat_handle@basis)[1],
  same_reconstruction = isTRUE(all.equal(as.matrix(lat), as.matrix(lat_handle)))
)
#>   eager_slot handle_slot same_reconstruction
#> 1  dgeMatrix BasisHandle                TRUE
```

The eager object stores a materialized matrix; the handle object stores
a `BasisHandle`. Both reconstruct to the same values — the handle just
defers the cost until you call
[`basis()`](../reference/basis-methods.md),
[`series()`](../reference/series-methods.md), or
[`as.matrix()`](https://rdrr.io/r/base/matrix.html).

Use
[`vignette("compression-diagnostics")`](../articles/compression-diagnostics.md)
when you want to compare error and storage budgets across component
counts.
