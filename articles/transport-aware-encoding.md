# Transport-Aware Encoding

``` r

library(fmrilatent)
```

Transport-aware encoding separates the shared basis from the
subject-specific field operator. `fmrilatent` owns coefficient recovery
and decode/projection; the operator that maps template fields into the
subject target domain is compiled outside the package.

## What basis are we transporting?

This small example uses six voxels split into two parcels. A
one-component PCA dictionary per parcel gives a two-column shared basis
asset.

``` r

mask <- array(TRUE, dim = c(3, 2, 1))
mask_vol <- neuroim2::LogicalNeuroVol(mask, neuroim2::NeuroSpace(dim(mask)))
parcel_map <- rep(1:2, each = 3L)
reduction <- make_cluster_reduction(mask_vol, parcel_map)

X_template <- matrix(c(
  1, 2, 2, 1, 0, 0,
  0, 1, 1, 0, 1, 1,
  1, 0, 0, 1, 2, 2,
  0, 2, 2, 0, 1, 1,
  1, 1, 0, 0, 2, 1,
  0, 0, 1, 1, 1, 2
), nrow = 6, byrow = TRUE)

tmpl <- parcel_basis_template(reduction, basis_pca(k = 1),
                              data = X_template, center = TRUE)
```

## What does the field operator look like?

For a real study, this object would come from registration and sampling
code. Here the template and subject grids match, so the operator is
identity.

``` r

identity_operator <- list(
  n_source = sum(mask),
  n_target = sum(mask),
  source_domain_id = "template:demo",
  target_domain_id = "subject01:bold",
  forward = function(z) as.matrix(z),
  adjoint_apply = function(y) as.matrix(y),
  provenance = list(operator = "identity-demo", target_mask = mask_vol)
)
```

## How are coefficients recovered?

[`encode_transport()`](../reference/encode_transport.md) solves for
coefficient time series in the shared analysis coordinates and carries
run metadata forward.

``` r

run_info <- list(run_id = rep(c(1L, 2L), each = 3L))

lat_tx <- encode_transport(
  X_template,
  basis_asset = tmpl,
  field_operator = identity_operator,
  lambda = 0,
  run_info = run_info
)

data.frame(
  time_points = nrow(coef_time(lat_tx, "analysis")),
  analysis_components = ncol(coef_time(lat_tx, "analysis")),
  finite = all(is.finite(coef_time(lat_tx, "analysis")))
)
#>   time_points analysis_components finite
#> 1           6                   2   TRUE
```

## What is handed downstream?

Downstream model-fitting code should consume analysis coordinates. When
a model returns a coefficient-space effect,
[`decode_coefficients()`](../reference/decode_coefficients.md) projects
it back to native or template space.

``` r

gamma <- coef_time(lat_tx, "analysis")[1, ]
native_map <- decode_coefficients(lat_tx, gamma, space = "native", wrap = "none")

data.frame(
  coefficient_length = length(gamma),
  native_values = length(native_map),
  finite_native = all(is.finite(native_map))
)
#>   coefficient_length native_values finite_native
#> 1                  2             6          TRUE
```

This is the package boundary: external code builds the field operator
and fits the statistical model; `fmrilatent` keeps the shared basis,
coefficients, and decode semantics coherent.
