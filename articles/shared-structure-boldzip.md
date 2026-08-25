# Shared Structure and BOLDZip-SR

``` r

library(fmrilatent)
#> 
#> Attaching package: 'fmrilatent'
#> The following object is masked from 'package:stats':
#> 
#>     offset
library(Matrix)
```

This article covers two research-facing surfaces that build on the core
`LatentNeuroVec` workflow:

- shared structure contracts for reusable spatial dictionaries and
  downstream handoff;
- BOLDZip-SR, an experimental matrix codec that stores carrier signals,
  static high-resolution texture loadings, and sparse residual events.

The examples are intentionally small. They demonstrate object flow and
sanity checks rather than compression claims for production-scale fMRI.

## Shared Basis Assets

A shared basis asset is a spatial dictionary that can be reused across
subjects. The subject changes the coefficient time series and offset;
the template loadings stay fixed.

We start with a tiny nine-voxel mask, three parcels, and a training
matrix with six time points:

``` r

mask <- array(TRUE, dim = c(3, 3, 1))
mask_vol <- neuroim2::LogicalNeuroVol(mask, neuroim2::NeuroSpace(dim(mask)))

parcel_map <- rep(1:3, each = 3L)
reduction <- make_cluster_reduction(mask_vol, parcel_map)

n_time <- 6L
t <- seq_len(n_time)
parcel_signal <- cbind(
  sin(2 * pi * t / n_time),
  cos(2 * pi * t / n_time),
  seq(-1, 1, length.out = n_time)
)
voxel_weights <- c(1.0, 0.8, 1.2, 1.1, 0.9, 1.3, 0.7, 1.0, 1.4)
X <- matrix(0, nrow = n_time, ncol = sum(mask))
for (parcel in 1:3) {
  cols <- which(parcel_map == parcel)
  X[, cols] <- parcel_signal[, parcel] %*% t(voxel_weights[cols])
}
dim(X)
#> [1] 6 9
```

[`parcel_basis_template()`](../reference/parcel_basis_template.md)
builds the reusable dictionary. Here we use a one-component PCA
dictionary per parcel so the example is deterministic and does not
require graph-neighbor tuning.

``` r

tmpl <- parcel_basis_template(
  reduction,
  basis_pca(k = 1),
  data = X,
  center = TRUE
)

protocol <- validate_template_protocol(tmpl)
data.frame(
  rank = protocol$rank,
  features = protocol$n_features,
  support = protocol$support_cardinality,
  has_roughness = protocol$has_roughness
)
#>   rank features support has_roughness
#> 1    3        9       9         FALSE
```

The template protocol is the public surface shared-structure consumers
use: loadings, rank, domain/support, measure, roughness, projection, and
decoder behavior. The same template can encode a subject matrix:

``` r

lat <- encode(X, spec_space_parcel(tmpl), mask = mask_vol)

data.frame(
  coefficient_rows = nrow(coef_time(lat, "analysis")),
  coefficient_cols = ncol(coef_time(lat, "analysis")),
  loading_rows = nrow(loadings(lat)),
  loading_cols = ncol(loadings(lat)),
  rmse = sqrt(mean((as.matrix(lat) - X)^2))
)
#>   coefficient_rows coefficient_cols loading_rows loading_cols         rmse
#> 1                6                3            9            3 1.003728e-16
```

The shared component contract records method-neutral facts about the
spatial dictionary. It does not say whether the basis came from PCA,
graph Slepians, AWPT, or another method.

``` r

component <- shared_component_contract(
  template_loadings(tmpl),
  family = "parcel_pca",
  domain_id = "toy-template",
  support = template_support(tmpl),
  measure = template_measure(tmpl)
)

data.frame(
  id_prefix = substr(component$id, 1, 18),
  family = component$family,
  features = component$n_features,
  components = component$n_components,
  digest_chars = nchar(component$digest)
)
#>            id_prefix     family features components digest_chars
#> 1 component:add7f59c parcel_pca        9          3           32
```

## In-Session References and Subject Deltas

[`shared_reference()`](../reference/shared_reference.md) is
intentionally an in-session cache. It can avoid duplicating a matrix
within an R session, but it does not create archive paths, checksums,
HDF5 datasets, or lazy file readers.

``` r

shared_reference_clear()
loadings_ref <- shared_reference(template_loadings(tmpl), kind = "loadings")
shared_reference_info(loadings_ref)[c("kind", "has_cached_value", "persistence")]
#> $kind
#> [1] "loadings"
#> 
#> $has_cached_value
#> [1] TRUE
#> 
#> $persistence
#> [1] "session"
```

Group-plus-delta loadings represent a subject-specific dictionary as a
shared group matrix plus a small subject delta.

``` r

delta <- Matrix(0, nrow = nrow(template_loadings(tmpl)), ncol = template_rank(tmpl),
                sparse = TRUE)
delta[1, 1] <- 0.05

subject_loadings <- group_delta_loadings(loadings_ref, delta, scale = 1)
materialized <- materialize_group_delta_loadings(subject_loadings)

data.frame(
  rows = nrow(materialized),
  cols = ncol(materialized),
  max_abs_delta = max(abs(as.matrix(materialized - template_loadings(tmpl))))
)
#>   rows cols max_abs_delta
#> 1    9    3          0.05
```

The neuroarchive handoff contract makes the package boundary explicit:
`fmrilatent` owns in-memory representation contracts and decode/project
behavior; archive packages own persistent stores and locators.

``` r

handoff <- neuroarchive_handoff_contract(
  representation = lat,
  components = list(component),
  templates = list(tmpl),
  references = list(loadings_ref)
)

data.frame(
  fmrilatent_items = length(handoff$fmrilatent_owns),
  neuroarchive_items = length(handoff$neuroarchive_owns),
  persistent = handoff$persistent,
  has_archive_locator = !is.null(handoff$archive_locator)
)
#>   fmrilatent_items neuroarchive_items persistent has_archive_locator
#> 1                6                  6      FALSE               FALSE
```

## Shared Temporal Specs and Event Dictionaries

The shared temporal descriptor captures a reusable temporal basis. A DCT
descriptor can be materialized directly:

``` r

time_spec <- shared_temporal_spec("dct", n_time = 8L, rank = 3L)
time_basis <- materialize_shared_temporal_spec(time_spec)

data.frame(
  rows = nrow(time_basis),
  cols = ncol(time_basis),
  max_orthogonality_error = max(abs(crossprod(time_basis) - diag(3)))
)
#>   rows cols max_orthogonality_error
#> 1    8    3            1.110223e-16
```

An event dictionary turns sparse event rows into a time-by-atom event
matrix. This is useful for sparse innovation streams such as BOLDZip-SR
residual events.

``` r

dictionary <- shared_event_dictionary(
  shapes = list(impulse = 1, box2 = c(1, 1)),
  n_time = 6L
)
events <- data.frame(
  atom = c(1L, 2L),
  time = c(2L, 4L),
  amplitude = c(3, -2),
  shape_id = c("impulse", "box2")
)
event_matrix <- render_shared_events(dictionary, events, n_atoms = 2L)

data.frame(
  rows = nrow(event_matrix),
  cols = ncol(event_matrix),
  nonzero = Matrix::nnzero(event_matrix),
  finite = all(is.finite(as.matrix(event_matrix)))
)
#>   rows cols nonzero finite
#> 1    2    6       3   TRUE
```

## BOLDZip-SR

BOLDZip-SR works at the matrix level with rows as voxels or
grayordinates and columns as time points. This orientation is different
from [`encode()`](../reference/encode.md), which uses time x voxels.

The model stores:

- a voxel-wise baseline map;
- temporally compressed carrier signals;
- sparse static texture loadings from fine detail atoms onto carriers;
- sparse residual events for reliable deviations.

``` r

n_voxels <- 8L
n_time <- 24L
tt <- seq_len(n_time)
carrier <- sin(2 * pi * tt / n_time)
texture <- seq(0.25, 2, length.out = n_voxels)
baseline <- seq(-1, 1, length.out = n_voxels)

X_vt <- texture %*% t(carrier)
X_vt <- sweep(X_vt, 1L, baseline, "+")
dim(X_vt)
#> [1]  8 24
```

The example is exactly rank one after centering, so one carrier and a
full DCT temporal budget should reconstruct it to numerical precision.

``` r

fit <- boldzip_sr_encode(
  X_vt,
  k_carriers = 1L,
  temporal_k = n_time,
  q_texture = 1L,
  reliability = boldzip_reliability(min_texture_reliability = 0.1),
  events = boldzip_events(max_events = 0L)
)

X_hat <- boldzip_sr_decode(fit)
metrics <- evaluate_boldzip_sr(X_vt, fit)

data.frame(
  mse = signif(metrics[["mse"]], 3),
  rmse = signif(metrics[["rmse"]], 3),
  payload_scalars = metrics[["payload_scalars"]],
  finite_reconstruction = all(is.finite(X_hat))
)
#>        mse     rmse payload_scalars finite_reconstruction
#> 1 3.16e-30 1.78e-15             650                  TRUE
```

The payload summary reports where the stored scalar budget goes. It is
an accounting diagnostic, not a binary file-size estimate.

``` r

boldzip_sr_payload_summary(fit)
#>          component scalar_count bytes
#> 1   carriers_theta           24   408
#> 2 texture_loadings           42  1464
#> 3           events            0  1352
#> 4      baseline_mu            8   112
#> 5   basis_metadata          576  5720
#> 6     total_object          650 21184
```

The codec can also use a graph-derived coarse/detail spatial basis. The
small path graph below creates orthonormal Laplacian eigenvectors; low
frequencies become carriers and higher frequencies become detail atoms.

``` r

adjacency <- matrix(0, nrow = 6L, ncol = 6L)
for (idx in seq_len(5L)) {
  adjacency[idx, idx + 1L] <- 1
  adjacency[idx + 1L, idx] <- 1
}

graph_basis <- boldzip_graph_spatial_basis(
  adjacency,
  n_coarse = 2L,
  n_detail = 3L
)
gram_error <- max(abs(crossprod(cbind(graph_basis$phi_c, graph_basis$phi_d)) - diag(5L)))

data.frame(
  coarse_atoms = ncol(graph_basis$phi_c),
  detail_atoms = ncol(graph_basis$phi_d),
  max_gram_error = signif(gram_error, 3)
)
#>   coarse_atoms detail_atoms max_gram_error
#> 1            2            3       1.07e-15
```

## What to Use When

Use shared basis contracts when you need a reusable representation
boundary: the same spatial dictionary, explicit projection behavior, and
a clean handoff to downstream model-fitting or archive code.

Use BOLDZip-SR when you want to experiment with a compression-oriented
carrier/texture/event decomposition. In its current form it is an
experimental matrix codec, not a production bitstream format.
