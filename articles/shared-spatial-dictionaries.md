# Shared Spatial Dictionaries Across Subjects

``` r

library(fmrilatent)
```

Shared spatial dictionaries are useful when multiple subjects should be
represented in the same spatial coordinates. The dictionary is built
once from a template reduction or training sample; each subject gets its
own coefficient time series and offset.

## What is shared?

We start with a nine-voxel mask and three parcels. The training matrix
has one dominant signal in each parcel.

``` r

mask <- array(TRUE, dim = c(3, 3, 1))
mask_vol <- neuroim2::LogicalNeuroVol(mask, neuroim2::NeuroSpace(dim(mask)))
parcel_map <- rep(1:3, each = 3L)
reduction <- make_cluster_reduction(mask_vol, parcel_map)

n_time <- 8L
t <- seq_len(n_time)
signals <- cbind(
  sin(2 * pi * t / n_time),
  cos(2 * pi * t / n_time),
  seq(-1, 1, length.out = n_time)
)
weights <- matrix(c(
  1.0, 0.8, 1.2,
  1.1, 0.9, 1.3,
  0.7, 1.0, 1.4
), nrow = 3, byrow = TRUE)

X_train <- matrix(0, nrow = n_time, ncol = sum(mask))
for (parcel in 1:3) {
  cols <- which(parcel_map == parcel)
  X_train[, cols] <- signals[, parcel] %*% t(weights[parcel, ])
}
```

Build the template once. This example uses one PCA component per parcel
so the render is deterministic and fast.

``` r

tmpl <- parcel_basis_template(
  reduction,
  basis_pca(k = 1),
  data = X_train,
  center = TRUE
)

validate_template_protocol(tmpl)[c("rank", "n_features", "support_cardinality")]
#> $rank
#> [1] 3
#> 
#> $n_features
#> [1] 9
#> 
#> $support_cardinality
#> [1] 9
```

## How do subjects differ?

Two subject matrices can share the same dictionary while having
different coefficient time series.

``` r

set.seed(19)
X_subject1 <- X_train + matrix(rnorm(length(X_train), sd = 0.01), nrow = n_time)
X_subject2 <- 0.8 * X_train + matrix(rnorm(length(X_train), sd = 0.01), nrow = n_time)

lat1 <- encode(X_subject1, spec_space_parcel(tmpl), mask = mask_vol)
lat2 <- encode(X_subject2, spec_space_parcel(tmpl), mask = mask_vol)

data.frame(
  subject = c("subject1", "subject2"),
  coefficient_rows = c(nrow(coef_time(lat1)), nrow(coef_time(lat2))),
  coefficient_cols = c(ncol(coef_time(lat1)), ncol(coef_time(lat2))),
  rmse = round(c(
    sqrt(mean((as.matrix(lat1) - X_subject1)^2)),
    sqrt(mean((as.matrix(lat2) - X_subject2)^2))
  ), 6)
)
#>    subject coefficient_rows coefficient_cols     rmse
#> 1 subject1                8                3 0.008307
#> 2 subject2                8                3 0.008162
```

The spatial loadings are the shared part. The coefficients are the
subject-specific part.

``` r

data.frame(
  same_loading_shape = identical(dim(loadings(lat1)), dim(loadings(lat2))),
  max_loading_difference = max(abs(loadings(lat1) - loadings(lat2))),
  coefficient_difference = max(abs(coef_time(lat1) - coef_time(lat2)))
)
#>   same_loading_shape max_loading_difference coefficient_difference
#> 1               TRUE                      0              0.3927821
```

Use
[`vignette("shared-structure-boldzip")`](../articles/shared-structure-boldzip.md)
when you need the lower-level shared-structure contracts and archive
handoff boundaries.
