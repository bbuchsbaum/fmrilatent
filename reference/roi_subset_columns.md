# Subset reconstruction matrix columns by ROI mask

Extracts columns from a reconstruction matrix that correspond to voxels
inside an ROI mask. This is a shared utility used by decoder functions
to handle ROI subsetting consistently.

## Usage

``` r
roi_subset_columns(rec_mat, mask_arr, roi_mask = NULL)
```

## Arguments

- rec_mat:

  Numeric matrix (time x voxels-in-mask) to subset.

- mask_arr:

  Logical array (the full brain mask).

- roi_mask:

  Logical array (the ROI mask), or `NULL` to return `rec_mat` unchanged.

## Value

The subsetted matrix, or `rec_mat` unchanged when `roi_mask` is `NULL`.

## Examples

``` r
mask <- array(c(TRUE, TRUE, FALSE, TRUE), dim = c(2, 2, 1))
rec  <- matrix(1:9, nrow = 3, ncol = 3)  # 3 time x 3 masked voxels
roi  <- array(c(TRUE, FALSE, FALSE, TRUE), dim = c(2, 2, 1))
roi_subset_columns(rec, mask, roi)  # keeps columns 1 and 3
#>      [,1] [,2]
#> [1,]    1    7
#> [2,]    2    8
#> [3,]    3    9
```
