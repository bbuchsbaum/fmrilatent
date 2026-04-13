# Tests for parcel_basis.R: shared parcel basis templates and projection encoding

library(testthat)
library(Matrix)
library(neuroim2)

# --- Test helpers ---

make_test_data <- function(n_time = 50, dims = c(5, 5, 5), n_parcels = 3) {
  mask <- array(TRUE, dim = dims)
  mask_vol <- LogicalNeuroVol(mask, NeuroSpace(dims))
  n_vox <- prod(dims)

  # Create parcel map (roughly equal-sized parcels)
  map <- as.integer(rep_len(seq_len(n_parcels), n_vox))

  # Generate data
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)

  list(X = X, mask = mask, mask_vol = mask_vol, map = map, n_time = n_time,
       n_vox = n_vox, n_parcels = n_parcels)
}

# --- as_cluster_reduction tests ---

test_that("as_cluster_reduction converts ClusteredNeuroVol to ClusterReduction", {
  dims <- c(5, 5, 5)
  mask_vol <- LogicalNeuroVol(array(TRUE, dim = dims), NeuroSpace(dims))
  n_vox <- prod(dims)
  clusters <- as.integer(rep_len(1:4, n_vox))

  cvol <- ClusteredNeuroVol(mask_vol, clusters)
  red <- as_cluster_reduction(cvol)

  expect_s4_class(red, "ClusterReduction")
  expect_equal(length(red@map), n_vox)
  expect_equal(sort(red@cluster_ids), 1:4)
  expect_s4_class(red@mask, "LogicalNeuroVol")
})

test_that("as_cluster_reduction preserves label_map", {
  dims <- c(4, 4, 4)
  mask_vol <- LogicalNeuroVol(array(TRUE, dim = dims), NeuroSpace(dims))
  clusters <- as.integer(rep_len(1:2, prod(dims)))
  lmap <- list("RegionA" = 1L, "RegionB" = 2L)

  cvol <- ClusteredNeuroVol(mask_vol, clusters, label_map = lmap)
  red <- as_cluster_reduction(cvol)

  expect_true(!is.null(red@info$label_map))
  expect_equal(red@info$label_map[["RegionA"]], 1L)
})

test_that("as_cluster_reduction preserves cluster_map when available", {
  dims <- c(4, 4, 4)
  mask_vol <- LogicalNeuroVol(array(TRUE, dim = dims), NeuroSpace(dims))
  clusters <- as.integer(rep_len(1:2, prod(dims)))

  cvol <- ClusteredNeuroVol(mask_vol, clusters)
  red <- as_cluster_reduction(cvol)

  expect_true(!is.null(red@info$cluster_map))
  expect_length(red@info$cluster_map, num_clusters(cvol))
})

test_that("as_cluster_reduction rejects non-ClusteredNeuroVol", {
  expect_error(as_cluster_reduction("not a cvol"), "ClusteredNeuroVol")
})

# --- parcel_basis_template tests with Slepian (geometric, no data needed) ---

test_that("parcel_basis_template builds Slepian template", {
  skip_if_not_installed("rgsp")
  td <- make_test_data(n_time = 30, dims = c(4, 4, 4), n_parcels = 2)
  red <- make_cluster_reduction(td$mask, td$map)

  tmpl <- parcel_basis_template(red, basis_slepian(k = 3))

  expect_s3_class(tmpl, "ParcelBasisTemplate")
  expect_equal(nrow(tmpl$loadings), td$n_vox)
  expect_true(ncol(tmpl$loadings) > 0)
  expect_equal(tmpl$meta$family, "spec_slepian")
  expect_true(tmpl$center)
  expect_true(!is.null(tmpl$gram_factor))
})

test_that("parcel_basis_template builds PCA template from data", {
  td <- make_test_data(n_time = 30, dims = c(4, 4, 4), n_parcels = 2)
  red <- make_cluster_reduction(td$mask, td$map)

  tmpl <- parcel_basis_template(red, basis_pca(k = 3), data = td$X)

  expect_s3_class(tmpl, "ParcelBasisTemplate")
  expect_equal(nrow(tmpl$loadings), td$n_vox)
  # 2 parcels x 3 components each = 6 total atoms
  expect_equal(ncol(tmpl$loadings), 6L)
  expect_equal(tmpl$meta$family, "spec_pca")
})

test_that("parcel_basis_template accepts ClusteredNeuroVol", {
  skip_if_not_installed("rgsp")
  dims <- c(4, 4, 4)
  mask_vol <- LogicalNeuroVol(array(TRUE, dim = dims), NeuroSpace(dims))
  n_vox <- prod(dims)
  clusters <- as.integer(rep_len(1:2, n_vox))
  cvol <- ClusteredNeuroVol(mask_vol, clusters)

  tmpl <- parcel_basis_template(cvol, basis_slepian(k = 2))

  expect_s3_class(tmpl, "ParcelBasisTemplate")
  expect_equal(length(tmpl$reduction@cluster_ids), 2L)
})

test_that("parcel_basis_template rejects invalid parcellation", {
  expect_error(parcel_basis_template("bad"), "ClusterReduction or ClusteredNeuroVol")
})

test_that("parcel_basis_template center=FALSE works", {
  skip_if_not_installed("rgsp")
  td <- make_test_data(n_time = 20, dims = c(4, 4, 4), n_parcels = 2)
  red <- make_cluster_reduction(td$mask, td$map)
  tmpl <- parcel_basis_template(red, basis_slepian(k = 2), center = FALSE)
  expect_false(tmpl$center)
})

test_that("parcel_basis_template passes center=FALSE through to PCA lift", {
  td <- make_test_data(n_time = 20, dims = c(4, 2, 2), n_parcels = 4)
  red <- make_cluster_reduction(td$mask, td$map)

  tmpl <- parcel_basis_template(
    red,
    basis_pca(k = 1),
    data = td$X,
    center = FALSE,
    backend = "svd"
  )
  L_nocenter <- lift(red, basis_pca(k = 1), data = td$X, center = FALSE, backend = "svd")

  expect_equal(as.matrix(tmpl$loadings), as.matrix(L_nocenter))
})

test_that("parcel_basis_template rejects unsupported PCA scale preprocessing", {
  td <- make_test_data(n_time = 20, dims = c(4, 2, 2), n_parcels = 4)
  red <- make_cluster_reduction(td$mask, td$map)

  expect_error(
    parcel_basis_template(red, basis_pca(k = 1), data = td$X, scale = TRUE),
    "scale = TRUE"
  )
})

test_that("parcel_basis_template rejects unsupported PCA whitening", {
  td <- make_test_data(n_time = 20, dims = c(4, 2, 2), n_parcels = 4)
  red <- make_cluster_reduction(td$mask, td$map)

  expect_error(
    parcel_basis_template(red, basis_pca(k = 1, whiten = TRUE), data = td$X),
    "whiten = TRUE"
  )
})

# --- print method ---

test_that("print.ParcelBasisTemplate runs without error", {
  skip_if_not_installed("rgsp")
  td <- make_test_data(n_time = 20, dims = c(4, 4, 4), n_parcels = 2)
  red <- make_cluster_reduction(td$mask, td$map)
  tmpl <- parcel_basis_template(red, basis_slepian(k = 2))
  expect_output(print(tmpl), "ParcelBasisTemplate")
  expect_output(print(tmpl), "Atoms:")
  expect_output(print(tmpl), "Parcels:")
})

# --- spec_space_parcel + encode tests ---

test_that("spec_space_parcel creates valid spec", {
  skip_if_not_installed("rgsp")
  td <- make_test_data(n_time = 20, dims = c(4, 4, 4), n_parcels = 2)
  red <- make_cluster_reduction(td$mask, td$map)
  tmpl <- parcel_basis_template(red, basis_slepian(k = 2))

  spec <- spec_space_parcel(tmpl)
  expect_s3_class(spec, "spec_space_parcel")
  expect_identical(spec$template, tmpl)
})

test_that("spec_space_parcel rejects non-template", {
  expect_error(spec_space_parcel(list()), "ParcelBasisTemplate")
})

test_that("encode with spec_space_parcel produces valid LatentNeuroVec", {
  skip_if_not_installed("rgsp")
  td <- make_test_data(n_time = 30, dims = c(4, 4, 4), n_parcels = 2)
  red <- make_cluster_reduction(td$mask, td$map)
  tmpl <- parcel_basis_template(red, basis_slepian(k = 3))

  lvec <- encode(td$X, spec_space_parcel(tmpl), mask = td$mask_vol)

  expect_s4_class(lvec, "LatentNeuroVec")
  expect_equal(dim(lvec)[4], td$n_time)
  expect_true(length(offset(lvec)) > 0)
  expect_equal(lvec@meta$family, "parcel_basis")
})

test_that("encode with spec_space_parcel preserves atlas labels in metadata", {
  skip_if_not_installed("rgsp")
  dims <- c(4, 4, 4)
  mask_vol <- LogicalNeuroVol(array(TRUE, dim = dims), NeuroSpace(dims))
  n_vox <- prod(dims)
  clusters <- as.integer(rep_len(1:2, n_vox))
  lmap <- list("RegionA" = 1L, "RegionB" = 2L)
  cvol <- ClusteredNeuroVol(mask_vol, clusters, label_map = lmap)
  tmpl <- parcel_basis_template(cvol, basis_slepian(k = 2))
  X <- matrix(rnorm(20 * n_vox), nrow = 20, ncol = n_vox)

  lvec <- encode(X, spec_space_parcel(tmpl), mask = mask_vol)

  expect_true(!is.null(lvec@meta$label_map))
  expect_equal(lvec@meta$label_map[["RegionA"]], 1L)
})

test_that("encode with spec_space_parcel preserves custom mask geometry", {
  skip_if_not_installed("rgsp")
  dims <- c(2, 2, 1)
  spc3 <- NeuroSpace(dims, spacing = c(2, 3, 4), origin = c(10, 20, 30))
  mask_vol <- LogicalNeuroVol(array(TRUE, dim = dims), spc3)
  map <- c(1L, 1L, 2L, 2L)
  X <- matrix(rnorm(5 * 4), nrow = 5, ncol = 4)
  tmpl <- parcel_basis_template(make_cluster_reduction(mask_vol, map), basis_slepian(k = 1))

  lvec <- encode(X, spec_space_parcel(tmpl), mask = mask_vol)

  expect_equal(neuroim2::spacing(neuroim2::space(mask(lvec))), c(2, 3, 4))
  expect_equal(neuroim2::origin(neuroim2::space(mask(lvec))), c(10, 20, 30))
})

test_that("encode with spec_space_parcel center=FALSE produces zero offset", {
  skip_if_not_installed("rgsp")
  td <- make_test_data(n_time = 30, dims = c(4, 4, 4), n_parcels = 2)
  red <- make_cluster_reduction(td$mask, td$map)
  tmpl <- parcel_basis_template(red, basis_slepian(k = 3), center = FALSE)

  lvec <- encode(td$X, spec_space_parcel(tmpl), mask = td$mask_vol)

  expect_equal(length(offset(lvec)), 0)
})

test_that("shared template gives same loadings for different subjects", {
  skip_if_not_installed("rgsp")
  td <- make_test_data(n_time = 30, dims = c(4, 4, 4), n_parcels = 2)
  red <- make_cluster_reduction(td$mask, td$map)
  tmpl <- parcel_basis_template(red, basis_slepian(k = 3))

  X2 <- matrix(rnorm(30 * td$n_vox), nrow = 30, ncol = td$n_vox)

  lvec1 <- encode(td$X, spec_space_parcel(tmpl), mask = td$mask_vol)
  lvec2 <- encode(X2, spec_space_parcel(tmpl), mask = td$mask_vol)

  # Loadings must be identical (same template)
  expect_identical(loadings(lvec1), loadings(lvec2))
  # Basis (scores) should differ
  expect_false(identical(as.matrix(basis(lvec1)), as.matrix(basis(lvec2))))
})

test_that("reconstruction via shared Slepian template is reasonable", {
  skip_if_not_installed("rgsp")
  td <- make_test_data(n_time = 30, dims = c(4, 4, 4), n_parcels = 2)
  red <- make_cluster_reduction(td$mask, td$map)
  # Use enough components to get a decent reconstruction
  tmpl <- parcel_basis_template(red, basis_slepian(k = 5))
  lvec <- encode(td$X, spec_space_parcel(tmpl), mask = td$mask_vol)

  # Reconstruct
  recon <- as.matrix(basis(lvec)) %*% t(as.matrix(loadings(lvec)))
  if (length(offset(lvec)) > 0) {
    recon <- sweep(recon, 2, offset(lvec), "+")
  }

  # With Laplacian basis and centered data, correlation should be positive
  cors <- sapply(seq_len(ncol(recon)), function(v) cor(td$X[, v], recon[, v]))
  expect_true(mean(cors, na.rm = TRUE) > 0)
})

test_that("encode with shared PCA template produces valid LatentNeuroVec", {
  td <- make_test_data(n_time = 40, dims = c(4, 4, 4), n_parcels = 3)
  red <- make_cluster_reduction(td$mask, td$map)
  tmpl <- parcel_basis_template(red, basis_pca(k = 3), data = td$X)

  X2 <- matrix(rnorm(40 * td$n_vox), nrow = 40, ncol = td$n_vox)
  lvec <- encode(X2, spec_space_parcel(tmpl), mask = td$mask_vol)

  expect_s4_class(lvec, "LatentNeuroVec")
  expect_equal(ncol(as.matrix(loadings(lvec))), 9L)  # 3 parcels x 3 components
})

test_that("encode rejects a mask that differs from the template mask", {
  skip_if_not_installed("rgsp")
  dims <- c(4, 4, 4)
  mask_arr <- array(FALSE, dim = dims)
  mask_arr[1:2, , ] <- TRUE
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(dims))
  n_vox <- sum(mask_arr)
  map <- as.integer(rep_len(1:2, n_vox))
  X <- matrix(rnorm(20 * n_vox), nrow = 20, ncol = n_vox)
  red <- make_cluster_reduction(mask_vol, map)
  tmpl <- parcel_basis_template(red, basis_slepian(k = 2))

  bad_mask_arr <- array(FALSE, dim = dims)
  bad_mask_arr[3:4, , ] <- TRUE
  bad_mask <- LogicalNeuroVol(bad_mask_arr, NeuroSpace(dims))

  expect_false(identical(as.array(bad_mask), as.array(mask_vol)))
  expect_equal(sum(as.array(bad_mask)), sum(as.array(mask_vol)))

  expect_error(
    encode(X, spec_space_parcel(tmpl), mask = bad_mask),
    "template mask"
  )
})

test_that("encode rejects data with wrong number of voxels", {
  skip_if_not_installed("rgsp")
  td <- make_test_data(n_time = 20, dims = c(4, 4, 4), n_parcels = 2)
  red <- make_cluster_reduction(td$mask, td$map)
  tmpl <- parcel_basis_template(red, basis_slepian(k = 2))

  X_bad <- matrix(rnorm(20 * 10), nrow = 20, ncol = 10)
  expect_error(
    encode_spec(X_bad, spec_space_parcel(tmpl), mask = td$mask_vol,
                reduction = NULL, materialize = "matrix", label = ""),
    "voxels"
  )
})

# --- latent_factory integration ---

test_that("latent_factory with parcel_space works", {
  skip_if_not_installed("rgsp")
  td <- make_test_data(n_time = 30, dims = c(4, 4, 4), n_parcels = 2)
  red <- make_cluster_reduction(td$mask, td$map)
  tmpl <- parcel_basis_template(red, basis_slepian(k = 3))

  lvec <- latent_factory("parcel_space", td$X, mask = td$mask_vol,
                         template = tmpl, materialize = "matrix")

  expect_s4_class(lvec, "LatentNeuroVec")
})

# --- make_cluster_reduction moved to reduction.R ---

test_that("make_cluster_reduction still works after move from heat_wavelet.R", {
  mask <- array(TRUE, dim = c(3, 3, 3))
  map <- rep_len(1:3, 27)
  red <- make_cluster_reduction(mask, map)

  expect_s4_class(red, "ClusterReduction")
  expect_equal(length(red@cluster_ids), 3L)
  expect_equal(length(red@map), 27L)
})

test_that("make_cluster_reduction rejects maps with wrong length", {
  mask <- array(TRUE, dim = c(2, 2, 1))
  expect_error(make_cluster_reduction(mask, 1:2), "number of voxels")
})
