# Tests for hierarchical_helpers.R
# Comprehensive test suite for hierarchical parcellation helper functions

# =============================================================================
# validate_nested_parcellations tests
# =============================================================================

test_that("validate_nested_parcellations accepts valid nested parcellations", {
  # Two levels: coarse (2 clusters) and fine (4 clusters)
  # Each fine cluster maps to exactly one coarse cluster
  coarse <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)
  fine <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)

  levels <- list(coarse, fine)

  expect_silent(validate_nested_parcellations(levels))
  expect_true(validate_nested_parcellations(levels))
})

test_that("validate_nested_parcellations rejects non-nested parcellations", {
  # Invalid: fine cluster 2 spans two coarse clusters
  coarse <- c(1L, 1L, 2L, 2L, 2L, 2L)
  fine <- c(1L, 2L, 2L, 3L, 3L, 4L)

  levels <- list(coarse, fine)

  expect_error(validate_nested_parcellations(levels), "not nested")
})

test_that("validate_nested_parcellations rejects different length vectors", {
  coarse <- c(1L, 1L, 2L, 2L)
  fine <- c(1L, 1L, 2L, 2L, 3L)

  levels <- list(coarse, fine)

  expect_error(validate_nested_parcellations(levels), "same length")
})

test_that("validate_nested_parcellations accepts single level", {
  single <- c(1L, 1L, 2L, 2L, 3L, 3L)
  levels <- list(single)

  expect_silent(validate_nested_parcellations(levels))
  expect_true(validate_nested_parcellations(levels))
})

test_that("validate_nested_parcellations rejects empty list", {
  expect_error(validate_nested_parcellations(list()), "non-empty list")
})

test_that("validate_nested_parcellations rejects non-list input", {
  expect_error(validate_nested_parcellations(c(1L, 2L, 3L)), "non-empty list")
  expect_error(validate_nested_parcellations(NULL), "non-empty list")
})

test_that("validate_nested_parcellations accepts three-level hierarchy", {
  # Three levels: coarse -> medium -> fine
  coarse <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)
  medium <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)
  fine <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L)

  levels <- list(coarse, medium, fine)

  expect_silent(validate_nested_parcellations(levels))
  expect_true(validate_nested_parcellations(levels))
})

test_that("validate_nested_parcellations handles levels with NAs in nested structure", {
  # Test that NAs are handled correctly during parent mapping
  coarse <- c(1L, 1L, NA_integer_, 2L, 2L)
  fine <- c(1L, 2L, NA_integer_, 3L, 4L)

  levels <- list(coarse, fine)

  # Should work since NAs are filtered out in parent mapping
  expect_silent(validate_nested_parcellations(levels))
})

test_that("validate_nested_parcellations detects non-nesting at middle level", {
  # Three levels where middle level breaks nesting with fine
  coarse <- c(1L, 1L, 1L, 2L, 2L, 2L)
  medium <- c(1L, 1L, 2L, 3L, 3L, 4L)
  fine <- c(1L, 2L, 2L, 3L, 4L, 5L)  # cluster 2 in fine maps to both 1 and 2 in medium

  levels <- list(coarse, medium, fine)

  expect_error(validate_nested_parcellations(levels), "not nested")
})

# =============================================================================
# cut_hclust_nested tests
# =============================================================================

test_that("cut_hclust_nested produces nested label vectors", {
  # Create a simple hierarchical clustering
  set.seed(42)
  d <- dist(matrix(rnorm(20), nrow = 10))
  hc <- hclust(d, method = "ward.D2")

  k_levels <- c(2, 4, 6)

  result <- cut_hclust_nested(hc, k_levels)

  expect_type(result, "list")
  expect_length(result, 3)

  # Each level should have 10 elements (number of objects)
  for (lvl in result) {
    expect_length(lvl, 10)
  }

  # Verify nesting
  expect_silent(validate_nested_parcellations(result))
})

test_that("cut_hclust_nested warns on unsorted k_levels", {
  set.seed(123)
  d <- dist(matrix(rnorm(12), nrow = 6))
  hc <- hclust(d)

  k_levels <- c(4, 2, 3)  # not sorted

  expect_warning(cut_hclust_nested(hc, k_levels), "not sorted")
})

test_that("cut_hclust_nested rejects non-hclust input", {
  expect_error(cut_hclust_nested("not_hclust", c(2, 4)), "must be an hclust")
})

test_that("cut_hclust_nested rejects k_levels < 1", {
  set.seed(456)
  d <- dist(matrix(rnorm(8), nrow = 4))
  hc <- hclust(d)

  expect_error(cut_hclust_nested(hc, c(0, 2)), "must be >= 1")
})

test_that("cut_hclust_nested handles single k value", {
  set.seed(100)
  d <- dist(matrix(rnorm(16), nrow = 8))
  hc <- hclust(d, method = "ward.D2")

  k_levels <- c(3)

  result <- cut_hclust_nested(hc, k_levels)

  expect_type(result, "list")
  expect_length(result, 1)
  expect_length(result[[1]], 8)
  expect_equal(length(unique(result[[1]])), 3)
})

test_that("cut_hclust_nested handles k equal to number of observations", {
  set.seed(200)
  d <- dist(matrix(rnorm(10), nrow = 5))
  hc <- hclust(d)

  # k = n means each observation is its own cluster
  k_levels <- c(2, 5)

  result <- cut_hclust_nested(hc, k_levels)

  expect_length(result, 2)
  expect_equal(length(unique(result[[2]])), 5)
})

test_that("cut_hclust_nested handles duplicate k values after sorting", {
  set.seed(300)
  d <- dist(matrix(rnorm(12), nrow = 6))
  hc <- hclust(d)

  # After sorting and unique, should become c(2, 3, 4)
  k_levels <- c(4, 2, 3, 2)

  expect_warning(result <- cut_hclust_nested(hc, k_levels), "not sorted")

  # Should have 3 unique levels after deduplication
  expect_length(result, 3)
})

test_that("cut_hclust_nested returns integer vectors", {
  set.seed(400)
  d <- dist(matrix(rnorm(20), nrow = 10))
  hc <- hclust(d)

  result <- cut_hclust_nested(hc, c(2, 5))

  for (lvl in result) {
    expect_type(lvl, "integer")
  }
})

test_that("cut_hclust_nested coerces k_levels to integer", {
  set.seed(500)
  d <- dist(matrix(rnorm(12), nrow = 6))
  hc <- hclust(d)

  # Pass numeric k_levels
  result <- cut_hclust_nested(hc, c(2.9, 4.1))

  expect_length(result, 2)
  # 2.9 becomes 2, 4.1 becomes 4
  expect_equal(length(unique(result[[1]])), 2)
  expect_equal(length(unique(result[[2]])), 4)
})

# =============================================================================
# parcel_similarity_matrix tests
# =============================================================================

test_that("parcel_similarity_matrix computes valid similarity", {
  n <- 5
  boundary_contact <- matrix(runif(n * n), n, n)
  boundary_contact <- (boundary_contact + t(boundary_contact)) / 2
  diag(boundary_contact) <- 0

  geo_dist <- matrix(runif(n * n) * 50, n, n)
  geo_dist <- (geo_dist + t(geo_dist)) / 2
  diag(geo_dist) <- 0

  yeo17 <- c(1, 1, 2, 2, 3)

  W <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17,
                                alpha = 0.5, beta = 0.3, gamma = 0.2, d0 = 30)

  expect_true(is.matrix(W))
  expect_equal(dim(W), c(n, n))

  # W should be symmetric
  expect_equal(W, t(W), tolerance = 1e-10)

  # Diagonal should be zero
  expect_equal(diag(W), rep(0, n))

  # All values should be non-negative (similarity)
  expect_true(all(W >= 0))
})

test_that("parcel_similarity_matrix rejects mismatched dimensions", {
  boundary <- matrix(1, 4, 4)
  geo <- matrix(1, 5, 5)
  yeo <- c(1, 2, 3, 4)


  expect_error(parcel_similarity_matrix(boundary, geo, yeo), "same dimensions")
})

test_that("parcel_similarity_matrix rejects mismatched yeo length", {
  boundary <- matrix(1, 4, 4)
  geo <- matrix(1, 4, 4)
  yeo <- c(1, 2, 3)

  expect_error(parcel_similarity_matrix(boundary, geo, yeo), "length must match")
})

test_that("parcel_similarity_matrix rejects non-matrix inputs", {
  n <- 4
  boundary <- matrix(1, n, n)
  geo <- matrix(1, n, n)
  yeo <- c(1, 2, 3, 4)

  expect_error(parcel_similarity_matrix(as.data.frame(boundary), geo, yeo), "must be matrices")
  expect_error(parcel_similarity_matrix(boundary, as.data.frame(geo), yeo), "must be matrices")
})

test_that("parcel_similarity_matrix respects weight parameters", {
  n <- 3
  # Create predictable inputs
  boundary_contact <- matrix(c(0, 1, 0,
                               1, 0, 1,
                               0, 1, 0), n, n, byrow = TRUE)
  geo_dist <- matrix(c(0, 10, 20,
                       10, 0, 10,
                       20, 10, 0), n, n, byrow = TRUE)
  yeo17 <- c(1, 1, 2)

  # Test with only boundary contribution
  W_boundary <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17,
                                          alpha = 1.0, beta = 0.0, gamma = 0.0, d0 = 30)

  # Off-diagonal should reflect boundary contact only
  expect_equal(W_boundary[1, 2], 1.0)
  expect_equal(W_boundary[1, 3], 0.0)

  # Test with only geodesic contribution
  W_geo <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17,
                                     alpha = 0.0, beta = 1.0, gamma = 0.0, d0 = 10)

  # Geodesic kernel: exp(-geo_dist / d0)
  expect_equal(W_geo[1, 2], exp(-10/10), tolerance = 1e-10)
  expect_equal(W_geo[1, 3], exp(-20/10), tolerance = 1e-10)

  # Test with only network contribution
  W_net <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17,
                                     alpha = 0.0, beta = 0.0, gamma = 1.0, d0 = 30)

  # Same network (1,2) should have value 1, different network (1,3) should have 0
  expect_equal(W_net[1, 2], 1.0)
  expect_equal(W_net[1, 3], 0.0)
})

test_that("parcel_similarity_matrix handles factor yeo17 input", {
  n <- 4
  boundary_contact <- matrix(0.5, n, n)
  diag(boundary_contact) <- 0
  geo_dist <- matrix(10, n, n)
  diag(geo_dist) <- 0

  # Pass as factor
  yeo17 <- factor(c("Network1", "Network1", "Network2", "Network2"))

  W <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17,
                                alpha = 0.5, beta = 0.3, gamma = 0.2, d0 = 30)

  expect_true(is.matrix(W))
  expect_equal(dim(W), c(n, n))
  expect_equal(diag(W), rep(0, n))
})

test_that("parcel_similarity_matrix handles single parcel case", {
  n <- 1
  boundary_contact <- matrix(0, n, n)
  geo_dist <- matrix(0, n, n)
  yeo17 <- c(1)

  W <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17,
                                alpha = 0.5, beta = 0.3, gamma = 0.2, d0 = 30)

  expect_equal(dim(W), c(1, 1))
  expect_equal(W[1, 1], 0)
})

test_that("parcel_similarity_matrix handles large d0 (weak geodesic effect)", {
  n <- 3
  boundary_contact <- matrix(0, n, n)
  geo_dist <- matrix(c(0, 100, 200,
                       100, 0, 100,
                       200, 100, 0), n, n, byrow = TRUE)
  yeo17 <- c(1, 2, 3)

  # Very large d0 should make geodesic kernel nearly uniform
  W <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17,
                                alpha = 0.0, beta = 1.0, gamma = 0.0, d0 = 10000)

  # All off-diagonal elements should be close to 1
  expect_true(all(W[upper.tri(W)] > 0.98))
})

# =============================================================================
# spectral_ward_hclust tests
# =============================================================================

test_that("spectral_ward_hclust produces valid hclust", {
  skip_if_not_installed("RSpectra")

  set.seed(789)
  n <- 10
  W <- matrix(runif(n * n), n, n)
  W <- (W + t(W)) / 2
  diag(W) <- 0

  hc <- spectral_ward_hclust(W, k_embed = 3)

  expect_s3_class(hc, "hclust")
  expect_length(hc$order, n)
})

test_that("spectral_ward_hclust respects hemisphere constraint", {
  skip_if_not_installed("RSpectra")

  set.seed(101)
  n <- 8
  W <- matrix(runif(n * n), n, n)
  W <- (W + t(W)) / 2
  diag(W) <- 0

  hemi <- c("L", "L", "L", "L", "R", "R", "R", "R")

  hc <- spectral_ward_hclust(W, k_embed = 2, hemi = hemi)

  expect_s3_class(hc, "hclust")
})

test_that("spectral_ward_hclust rejects non-matrix input", {
  expect_error(spectral_ward_hclust("not_matrix"), "must be a matrix")
})

test_that("spectral_ward_hclust rejects non-square matrix", {
  W <- matrix(1, 4, 5)
  expect_error(spectral_ward_hclust(W), "must be square")
})

test_that("spectral_ward_hclust respects network constraint", {
  skip_if_not_installed("RSpectra")

  set.seed(202)
  n <- 12
  W <- matrix(runif(n * n), n, n)
  W <- (W + t(W)) / 2
  diag(W) <- 0

  network <- c("Net1", "Net1", "Net1", "Net2", "Net2", "Net2",
               "Net3", "Net3", "Net3", "Net4", "Net4", "Net4")

  hc <- spectral_ward_hclust(W, k_embed = 3, network = network)

  expect_s3_class(hc, "hclust")
  expect_length(hc$order, n)
})

test_that("spectral_ward_hclust handles combined hemi and network constraints", {
  skip_if_not_installed("RSpectra")

  set.seed(303)
  n <- 8
  W <- matrix(runif(n * n), n, n)
  W <- (W + t(W)) / 2
  diag(W) <- 0

  hemi <- c("L", "L", "L", "L", "R", "R", "R", "R")
  network <- c("A", "A", "B", "B", "A", "A", "B", "B")

  hc <- spectral_ward_hclust(W, k_embed = 2, hemi = hemi, network = network)

  expect_s3_class(hc, "hclust")
})

test_that("spectral_ward_hclust handles small matrix (k_embed > n-1)", {
  skip_if_not_installed("RSpectra")

  set.seed(404)
  n <- 3
  W <- matrix(runif(n * n), n, n)
  W <- (W + t(W)) / 2
  diag(W) <- 0

  # k_embed = 5 is larger than n-1 = 2, should be clamped
  # RSpectra may fall back to eigen() with a message, so allow warnings
  hc <- suppressWarnings(spectral_ward_hclust(W, k_embed = 5))

  expect_s3_class(hc, "hclust")
  expect_length(hc$order, n)
})

test_that("spectral_ward_hclust handles k_embed = 2", {
  skip_if_not_installed("RSpectra")

  set.seed(505)
  n <- 6
  W <- matrix(runif(n * n), n, n)
  W <- (W + t(W)) / 2
  diag(W) <- 0

  hc <- spectral_ward_hclust(W, k_embed = 2)

  expect_s3_class(hc, "hclust")
})

test_that("spectral_ward_hclust produces consistent output for same seed", {
  skip_if_not_installed("RSpectra")

  n <- 8
  set.seed(606)
  W <- matrix(runif(n * n), n, n)
  W <- (W + t(W)) / 2
  diag(W) <- 0

  set.seed(707)
  hc1 <- spectral_ward_hclust(W, k_embed = 3)

  set.seed(707)
  hc2 <- spectral_ward_hclust(W, k_embed = 3)

  expect_equal(hc1$merge, hc2$merge)
  expect_equal(hc1$height, hc2$height)
})

# =============================================================================
# parent_maps_from_levels tests
# =============================================================================

test_that("parent_maps_from_levels produces valid parent maps", {
  coarse <- c(1L, 1L, 1L, 2L, 2L, 2L)
  fine <- c(1L, 2L, 2L, 3L, 4L, 4L)

  levels <- list(coarse, fine)
  parents <- parent_maps_from_levels(levels)

  expect_type(parents, "list")
  expect_length(parents, 2)

  # First level has no parents
  expect_length(parents[[1]], 0)

  # Second level maps each fine parcel to its coarse parent
  expect_equal(parents[[2]][["1"]], 1L)
  expect_equal(parents[[2]][["2"]], 1L)
  expect_equal(parents[[2]][["3"]], 2L)
  expect_equal(parents[[2]][["4"]], 2L)
})

test_that("parent_maps_from_levels handles single level", {
  single <- c(1L, 1L, 2L, 2L)
  levels <- list(single)

  parents <- parent_maps_from_levels(levels)

  expect_length(parents, 1)
  expect_length(parents[[1]], 0)
})

test_that("parent_maps_from_levels handles three-level hierarchy", {
  coarse <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)
  medium <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)
  fine <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L)

  levels <- list(coarse, medium, fine)
  parents <- parent_maps_from_levels(levels)

  expect_length(parents, 3)
  expect_length(parents[[1]], 0)

  # Level 2: medium -> coarse
  expect_equal(parents[[2]][["1"]], 1L)
  expect_equal(parents[[2]][["2"]], 1L)
  expect_equal(parents[[2]][["3"]], 2L)
  expect_equal(parents[[2]][["4"]], 2L)

  # Level 3: fine -> medium
  expect_equal(parents[[3]][["1"]], 1L)
  expect_equal(parents[[3]][["2"]], 1L)
  expect_equal(parents[[3]][["3"]], 2L)
  expect_equal(parents[[3]][["4"]], 2L)
  expect_equal(parents[[3]][["5"]], 3L)
  expect_equal(parents[[3]][["6"]], 3L)
  expect_equal(parents[[3]][["7"]], 4L)
  expect_equal(parents[[3]][["8"]], 4L)
})

test_that("parent_maps_from_levels returns named integer vectors", {
  coarse <- c(1L, 1L, 2L, 2L)
  fine <- c(1L, 2L, 3L, 4L)

  levels <- list(coarse, fine)
  parents <- parent_maps_from_levels(levels)

  expect_type(parents[[2]], "integer")
  expect_true(!is.null(names(parents[[2]])))
  expect_equal(sort(names(parents[[2]])), c("1", "2", "3", "4"))
})

test_that("parent_maps_from_levels handles non-contiguous cluster IDs", {
  # Cluster IDs don't have to be 1, 2, 3, ...
  coarse <- c(5L, 5L, 10L, 10L)
  fine <- c(100L, 200L, 300L, 400L)

  levels <- list(coarse, fine)
  parents <- parent_maps_from_levels(levels)

  expect_equal(parents[[2]][["100"]], 5L)
  expect_equal(parents[[2]][["200"]], 5L)
  expect_equal(parents[[2]][["300"]], 10L)
  expect_equal(parents[[2]][["400"]], 10L)
})

test_that("parent_maps_from_levels rejects invalid nesting", {
  coarse <- c(1L, 1L, 2L, 2L)
  fine <- c(1L, 2L, 2L, 3L)  # cluster 2 spans both coarse clusters

  levels <- list(coarse, fine)

  expect_error(parent_maps_from_levels(levels), "not nested")
})

# =============================================================================
# Edge cases and integration tests
# =============================================================================

test_that("full pipeline: similarity matrix to hclust to nested levels", {
  skip_if_not_installed("RSpectra")

  set.seed(999)
  n <- 12

  # Create synthetic parcel data
  boundary_contact <- matrix(runif(n * n), n, n)
  boundary_contact <- (boundary_contact + t(boundary_contact)) / 2
  diag(boundary_contact) <- 0

  geo_dist <- matrix(runif(n * n) * 100, n, n)
  geo_dist <- (geo_dist + t(geo_dist)) / 2
  diag(geo_dist) <- 0

  yeo17 <- rep(1:3, each = 4)

  # Compute similarity matrix
  W <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17,
                                alpha = 0.5, beta = 0.3, gamma = 0.2, d0 = 30)

  # Perform spectral clustering
  hc <- spectral_ward_hclust(W, k_embed = 3)

  # Cut into nested levels
  k_levels <- c(2, 4, 8)
  levels <- cut_hclust_nested(hc, k_levels)

  # Validate nesting
  expect_silent(validate_nested_parcellations(levels))

  # Get parent maps
  parents <- parent_maps_from_levels(levels)

  expect_length(parents, 3)
  expect_length(parents[[1]], 0)

  # Verify cluster counts
  expect_equal(length(unique(levels[[1]])), 2)
  expect_equal(length(unique(levels[[2]])), 4)
  expect_equal(length(unique(levels[[3]])), 8)
})

test_that("validate_nested_parcellations handles all-same-cluster level", {
  # All voxels in one cluster at coarse level
  coarse <- rep(1L, 8)
  fine <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)

  levels <- list(coarse, fine)

  expect_silent(validate_nested_parcellations(levels))
})

test_that("parcel_similarity_matrix handles zero boundary contact", {
  n <- 4
  boundary_contact <- matrix(0, n, n)
  geo_dist <- matrix(c(0, 10, 20, 30,
                       10, 0, 10, 20,
                       20, 10, 0, 10,
                       30, 20, 10, 0), n, n, byrow = TRUE)
  yeo17 <- c(1, 1, 2, 2)

  W <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17,
                                alpha = 0.5, beta = 0.3, gamma = 0.2, d0 = 30)

  expect_true(is.matrix(W))
  # No boundary contribution, only geodesic and network
  expect_true(all(W >= 0))
})

test_that("parcel_similarity_matrix handles zero geodesic distance", {
  n <- 3
  boundary_contact <- matrix(c(0, 0.5, 0.3,
                               0.5, 0, 0.4,
                               0.3, 0.4, 0), n, n, byrow = TRUE)
  geo_dist <- matrix(0, n, n)  # All distances are zero
  yeo17 <- c(1, 1, 2)

  W <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17,
                                alpha = 0.5, beta = 0.3, gamma = 0.2, d0 = 30)

  expect_true(is.matrix(W))
  expect_equal(diag(W), rep(0, n))

  # With zero distance, geodesic kernel = exp(0) = 1 for all pairs
  # So off-diagonal should include beta * 1 = 0.3 contribution
  expect_true(all(W[upper.tri(W)] >= 0.3))
})

test_that("cut_hclust_nested handles negative k values after coercion", {
  set.seed(111)
  d <- dist(matrix(rnorm(12), nrow = 6))
  hc <- hclust(d)

  expect_error(cut_hclust_nested(hc, c(-1, 2)), "must be >= 1")
})

test_that("spectral_ward_hclust handles uniform similarity matrix", {
  skip_if_not_installed("RSpectra")

  n <- 5
  W <- matrix(0.5, n, n)
  diag(W) <- 0

  # Should still produce valid hclust even with uniform similarities
  hc <- spectral_ward_hclust(W, k_embed = 2)

  expect_s3_class(hc, "hclust")
  expect_length(hc$order, n)
})

test_that("parcel_similarity_matrix returns numeric matrix", {
  n <- 4
  boundary_contact <- matrix(1, n, n)
  diag(boundary_contact) <- 0
  geo_dist <- matrix(10, n, n)
  diag(geo_dist) <- 0
  yeo17 <- c(1, 1, 2, 2)

  W <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17)

  expect_true(is.numeric(W))
  expect_true(is.matrix(W))
  expect_false(is.integer(W[1, 2]))
})

test_that("parent_maps_from_levels preserves all child cluster IDs",
 {
  coarse <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L)
  fine <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L)

  levels <- list(coarse, fine)
  parents <- parent_maps_from_levels(levels)

  # All 9 fine cluster IDs should be in the parent map
  expect_length(parents[[2]], 9)
  expect_equal(sort(as.integer(names(parents[[2]]))), 1:9)
})

# =============================================================================
# Additional edge case tests for .parent_map_one (internal)
# =============================================================================

test_that(".parent_map_one detects invalid nesting with multiple parents", {
  # Child cluster 2 maps to both parent 1 and parent 2
  child <- c(1L, 2L, 2L, 3L)
  parent <- c(1L, 1L, 2L, 2L)

  expect_error(fmrilatent:::.parent_map_one(child, parent, lvl = 2),
               "not nested")
})

test_that(".parent_map_one handles NAs in parent mapping", {
  # NAs should be filtered out
  child <- c(1L, 1L, 2L, 2L, 3L)
  parent <- c(1L, 1L, NA_integer_, 2L, 2L)

  # Child 2 has positions 3 and 4, with parent values NA and 2
  # After filtering NAs, child 2 maps to parent 2 (valid)
  result <- fmrilatent:::.parent_map_one(child, parent, lvl = 2)
  expect_equal(result[["2"]], 2L)

  # Test case where child truly has only NA parents
  child2 <- c(1L, 1L, 2L, 3L, 3L)
  parent2 <- c(1L, 1L, NA_integer_, 2L, 2L)
  # Child 2 only at position 3 with parent NA - should error
  expect_error(fmrilatent:::.parent_map_one(child2, parent2, lvl = 2))
})

test_that(".parent_map_one creates named integer vector with all child IDs", {
  child <- c(10L, 10L, 20L, 20L, 30L, 30L)
  parent <- c(1L, 1L, 1L, 1L, 2L, 2L)

  result <- fmrilatent:::.parent_map_one(child, parent, lvl = 2)

  expect_type(result, "integer")
  expect_equal(sort(names(result)), c("10", "20", "30"))
  expect_equal(result[["10"]], 1L)
  expect_equal(result[["20"]], 1L)
  expect_equal(result[["30"]], 2L)
})

test_that(".parent_map_one handles single-element clusters", {
  child <- c(1L, 2L, 3L, 4L)
  parent <- c(1L, 1L, 2L, 2L)

  result <- fmrilatent:::.parent_map_one(child, parent, lvl = 2)

  expect_length(result, 4)
  expect_equal(result[["1"]], 1L)
  expect_equal(result[["2"]], 1L)
  expect_equal(result[["3"]], 2L)
  expect_equal(result[["4"]], 2L)
})

# =============================================================================
# Additional parcel_similarity_matrix tests
# =============================================================================

test_that("parcel_similarity_matrix handles non-symmetric inputs by symmetrizing", {
  n <- 3
  # Create non-symmetric boundary matrix
  boundary_contact <- matrix(c(0, 0.5, 0.2,
                               0.6, 0, 0.3,
                               0.4, 0.5, 0), n, n, byrow = TRUE)
  geo_dist <- matrix(c(0, 10, 20,
                       10, 0, 15,
                       20, 15, 0), n, n, byrow = TRUE)
  yeo17 <- c(1, 1, 2)

  W <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17,
                                alpha = 0.5, beta = 0.3, gamma = 0.2, d0 = 30)

  # Result should be symmetric
  expect_equal(W, t(W), tolerance = 1e-10)
})

test_that("parcel_similarity_matrix handles all weights equal to zero", {
  n <- 4
  boundary_contact <- matrix(runif(n * n), n, n)
  boundary_contact <- (boundary_contact + t(boundary_contact)) / 2
  diag(boundary_contact) <- 0
  geo_dist <- matrix(runif(n * n) * 50, n, n)
  geo_dist <- (geo_dist + t(geo_dist)) / 2
  diag(geo_dist) <- 0
  yeo17 <- c(1, 1, 2, 2)

  # All weights zero
  W <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17,
                                alpha = 0.0, beta = 0.0, gamma = 0.0, d0 = 30)

  expect_equal(W, matrix(0, n, n))
})

test_that("parcel_similarity_matrix handles very small d0", {
  n <- 3
  boundary_contact <- matrix(0, n, n)
  geo_dist <- matrix(c(0, 10, 20,
                       10, 0, 10,
                       20, 10, 0), n, n, byrow = TRUE)
  yeo17 <- c(1, 2, 3)

  # Very small d0 should make geodesic kernel decay rapidly
  W <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17,
                                alpha = 0.0, beta = 1.0, gamma = 0.0, d0 = 0.1)

  # All off-diagonal should be very close to zero
  expect_true(all(W[upper.tri(W)] < 1e-30))
})

test_that("parcel_similarity_matrix validates alpha + beta + gamma don't need to sum to 1", {
  n <- 3
  boundary_contact <- matrix(0.5, n, n)
  diag(boundary_contact) <- 0
  geo_dist <- matrix(10, n, n)
  diag(geo_dist) <- 0
  yeo17 <- c(1, 1, 2)

  # Weights sum to 1.5
  W <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17,
                                alpha = 0.5, beta = 0.5, gamma = 0.5, d0 = 30)

  expect_true(is.matrix(W))
  expect_equal(dim(W), c(n, n))
})

test_that("parcel_similarity_matrix handles negative geodesic distances gracefully", {
  n <- 3
  boundary_contact <- matrix(0, n, n)
  # Include negative distances (shouldn't happen in practice, but test robustness)
  geo_dist <- matrix(c(0, -10, 20,
                       -10, 0, 15,
                       20, 15, 0), n, n, byrow = TRUE)
  yeo17 <- c(1, 2, 3)

  # exp(-(-10)/30) = exp(0.333) > 1
  W <- parcel_similarity_matrix(boundary_contact, geo_dist, yeo17,
                                alpha = 0.0, beta = 1.0, gamma = 0.0, d0 = 30)

  expect_true(is.matrix(W))
  # With negative distance, kernel will be > 1
  expect_true(W[1, 2] > 1.0)
})

# =============================================================================
# Additional validate_nested_parcellations tests
# =============================================================================

test_that("validate_nested_parcellations handles parcellations with gaps in cluster IDs", {
  # Cluster IDs skip values: 1, 3, 5 (no 2, 4)
  coarse <- c(1L, 1L, 3L, 3L, 5L, 5L)
  fine <- c(10L, 11L, 30L, 31L, 50L, 51L)

  levels <- list(coarse, fine)

  expect_silent(validate_nested_parcellations(levels))
})

test_that("validate_nested_parcellations rejects when fine level is coarser", {
  # "Coarse" has more clusters than "fine"
  coarse <- c(1L, 2L, 3L, 4L, 5L, 6L)
  fine <- c(1L, 1L, 2L, 2L, 3L, 3L)

  levels <- list(coarse, fine)

  # This should work if nesting is valid (each fine cluster maps to one coarse)
  # But here cluster 1 in fine maps to both 1 and 2 in coarse
  expect_error(validate_nested_parcellations(levels), "not nested")
})

test_that("validate_nested_parcellations handles all zeros", {
  coarse <- rep(0L, 6)
  fine <- rep(0L, 6)

  levels <- list(coarse, fine)

  # All voxels in cluster 0 at both levels
  expect_silent(validate_nested_parcellations(levels))
})

test_that("validate_nested_parcellations handles mixed NA and valid values", {
  coarse <- c(1L, 1L, NA_integer_, NA_integer_, 2L, 2L)
  fine <- c(1L, 2L, NA_integer_, NA_integer_, 3L, 4L)

  levels <- list(coarse, fine)

  # NAs should be filtered in parent mapping, so this should work
  expect_silent(validate_nested_parcellations(levels))
})

# =============================================================================
# Additional cut_hclust_nested tests
# =============================================================================

test_that("cut_hclust_nested handles very large k (k > n)", {
  set.seed(888)
  d <- dist(matrix(rnorm(12), nrow = 6))
  hc <- hclust(d)

  # k=10 is larger than n=6; cutree will error
  k_levels <- c(2, 10)

  # cutree errors when k > n
  expect_error(cut_hclust_nested(hc, k_levels), "elements of 'k' must be between")

  # But k=n should work fine
  k_levels2 <- c(2, 6)
  result <- cut_hclust_nested(hc, k_levels2)
  expect_length(result, 2)
  expect_equal(length(unique(result[[2]])), 6)
})

test_that("cut_hclust_nested preserves cluster labels across levels", {
  set.seed(777)
  d <- dist(matrix(rnorm(20), nrow = 10))
  hc <- hclust(d, method = "ward.D2")

  k_levels <- c(2, 3, 5)
  result <- cut_hclust_nested(hc, k_levels)

  # Each observation should have a label at each level
  for (lvl in result) {
    expect_false(any(is.na(lvl)))
    expect_true(all(lvl > 0))
  }
})

test_that("cut_hclust_nested handles k_levels with length 1", {
  set.seed(666)
  d <- dist(matrix(rnorm(16), nrow = 8))
  hc <- hclust(d)

  result <- cut_hclust_nested(hc, c(4))

  expect_length(result, 1)
  expect_equal(length(unique(result[[1]])), 4)
})

# =============================================================================
# Additional spectral_ward_hclust tests
# =============================================================================

test_that("spectral_ward_hclust handles sparse similarity matrix", {
  skip_if_not_installed("RSpectra")

  set.seed(555)
  n <- 10
  # Create sparse W (mostly zeros)
  W <- matrix(0, n, n)
  # Add a few edges
  W[1, 2] <- W[2, 1] <- 0.8
  W[3, 4] <- W[4, 3] <- 0.7
  W[5, 6] <- W[6, 5] <- 0.6

  hc <- spectral_ward_hclust(W, k_embed = 2)

  expect_s3_class(hc, "hclust")
  expect_length(hc$order, n)
})

test_that("spectral_ward_hclust handles identity-like similarity matrix", {
  skip_if_not_installed("RSpectra")

  n <- 6
  W <- diag(n)
  diag(W) <- 0  # Zero diagonal as required

  # All parcels have zero similarity - degenerate case
  hc <- spectral_ward_hclust(W, k_embed = 2)

  expect_s3_class(hc, "hclust")
})

test_that("spectral_ward_hclust handles k_embed = 1", {
  skip_if_not_installed("RSpectra")

  set.seed(444)
  n <- 8
  W <- matrix(runif(n * n), n, n)
  W <- (W + t(W)) / 2
  diag(W) <- 0

  hc <- spectral_ward_hclust(W, k_embed = 1)

  expect_s3_class(hc, "hclust")
  expect_length(hc$order, n)
})

test_that("spectral_ward_hclust handles very dense similarity matrix", {
  skip_if_not_installed("RSpectra")

  set.seed(333)
  n <- 7
  # All similarities high
  W <- matrix(0.9, n, n)
  diag(W) <- 0

  hc <- spectral_ward_hclust(W, k_embed = 3)

  expect_s3_class(hc, "hclust")
  # hc$merge is a matrix with (n-1) rows and 2 columns
  expect_equal(nrow(hc$merge), n - 1)
})

test_that("spectral_ward_hclust with network constraint produces valid clustering", {
  skip_if_not_installed("RSpectra")

  set.seed(222)
  n <- 9
  W <- matrix(runif(n * n), n, n)
  W <- (W + t(W)) / 2
  diag(W) <- 0

  network <- rep(c("A", "B", "C"), each = 3)

  hc <- spectral_ward_hclust(W, k_embed = 3, network = network)

  expect_s3_class(hc, "hclust")
  # Check that height increases monotonically
  expect_true(all(diff(hc$height) >= 0))
})

test_that("spectral_ward_hclust with both constraints applies combined penalties", {
  skip_if_not_installed("RSpectra")

  set.seed(111)
  n <- 12
  W <- matrix(runif(n * n), n, n)
  W <- (W + t(W)) / 2
  diag(W) <- 0

  hemi <- rep(c("L", "R"), each = 6)
  network <- rep(c("Net1", "Net2", "Net3"), times = 4)

  hc <- spectral_ward_hclust(W, k_embed = 3, hemi = hemi, network = network)

  expect_s3_class(hc, "hclust")
  expect_length(hc$labels, n)
})

# =============================================================================
# Additional parent_maps_from_levels tests
# =============================================================================

test_that("parent_maps_from_levels handles four-level hierarchy", {
  L0 <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)
  L1 <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)
  L2 <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L)
  L3 <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L)  # Same as L2

  levels <- list(L0, L1, L2, L3)
  parents <- parent_maps_from_levels(levels)

  expect_length(parents, 4)
  expect_length(parents[[1]], 0)
  expect_length(parents[[2]], 4)
  expect_length(parents[[3]], 8)
  expect_length(parents[[4]], 8)
})

test_that("parent_maps_from_levels returns empty list for first level", {
  single <- c(1L, 1L, 2L, 2L, 3L, 3L)
  levels <- list(single)

  parents <- parent_maps_from_levels(levels)

  expect_equal(parents[[1]], integer(0))
})

test_that("parent_maps_from_levels handles unequal cluster sizes", {
  # Coarse cluster 1 has 5 voxels, cluster 2 has 3
  coarse <- c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L)
  fine <- c(1L, 1L, 2L, 2L, 3L, 4L, 5L, 5L)

  levels <- list(coarse, fine)
  parents <- parent_maps_from_levels(levels)

  expect_equal(parents[[2]][["1"]], 1L)
  expect_equal(parents[[2]][["2"]], 1L)
  expect_equal(parents[[2]][["3"]], 1L)
  expect_equal(parents[[2]][["4"]], 2L)
  expect_equal(parents[[2]][["5"]], 2L)
})

# =============================================================================
# Integration tests with realistic edge cases
# =============================================================================

test_that("integration: validate_nested -> parent_maps -> validate again", {
  coarse <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L)
  medium <- c(1L, 1L, 2L, 3L, 3L, 4L, 5L, 6L, 6L)
  fine <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L)

  levels <- list(coarse, medium, fine)

  # Validate nesting
  expect_silent(validate_nested_parcellations(levels))

  # Get parent maps
  parents <- parent_maps_from_levels(levels)

  # Check parent relationships
  expect_equal(parents[[2]][["1"]], 1L)
  expect_equal(parents[[2]][["2"]], 1L)
  expect_equal(parents[[3]][["1"]], 1L)
  expect_equal(parents[[3]][["2"]], 1L)
  expect_equal(parents[[3]][["3"]], 2L)
})

test_that("integration: full workflow with degenerate single-cluster level", {
  # Start with all voxels in one cluster
  level0 <- rep(1L, 10)
  level1 <- c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L)
  level2 <- c(1L, 1L, 2L, 2L, 3L, 4L, 4L, 5L, 5L, 6L)

  levels <- list(level0, level1, level2)

  expect_silent(validate_nested_parcellations(levels))

  parents <- parent_maps_from_levels(levels)

  # Level 1 parents: both clusters map to cluster 1
  expect_equal(parents[[2]][["1"]], 1L)
  expect_equal(parents[[2]][["2"]], 1L)
})

# =============================================================================
# Additional coverage tests for spectral_ward_hclust
# =============================================================================

test_that("spectral_ward_hclust Laplacian and embedding are correct", {
  skip_if_not_installed("RSpectra")

  # Test with a block-diagonal W to verify spectral clustering finds blocks
  n <- 8
  W <- matrix(0, n, n)
  # Block 1: nodes 1-4

  W[1:4, 1:4] <- 0.9
  # Block 2: nodes 5-8
  W[5:8, 5:8] <- 0.9
  # Weak inter-block connections
  W[4, 5] <- W[5, 4] <- 0.01
  diag(W) <- 0

  hc <- spectral_ward_hclust(W, k_embed = 2)
  expect_s3_class(hc, "hclust")

  # Cutting at k=2 should separate the two blocks
  labels <- cutree(hc, k = 2)
  # Nodes 1-4 should be in one cluster, 5-8 in another
  expect_equal(length(unique(labels[1:4])), 1L)
  expect_equal(length(unique(labels[5:8])), 1L)
  expect_true(labels[1] != labels[5])
})

test_that("spectral_ward_hclust with hemi penalty pushes cross-hemi merges later", {
  skip_if_not_installed("RSpectra")

  set.seed(808)
  n <- 6
  # Uniform similarity
  W <- matrix(0.5, n, n)
  diag(W) <- 0

  hemi <- c("L", "L", "L", "R", "R", "R")

  hc <- spectral_ward_hclust(W, k_embed = 2, hemi = hemi)
  expect_s3_class(hc, "hclust")

  # At k=2, the two groups should correspond to the two hemispheres
  labels <- cutree(hc, k = 2)
  expect_equal(length(unique(labels[1:3])), 1L)
  expect_equal(length(unique(labels[4:6])), 1L)
})

test_that("spectral_ward_hclust with network penalty only", {
  skip_if_not_installed("RSpectra")

  set.seed(909)
  n <- 9
  W <- matrix(0.3, n, n)
  diag(W) <- 0

  network <- c("A", "A", "A", "B", "B", "B", "C", "C", "C")

  hc <- spectral_ward_hclust(W, k_embed = 3, hemi = NULL, network = network)
  expect_s3_class(hc, "hclust")

  # At k=3, clusters should roughly match networks
  labels <- cutree(hc, k = 3)
  # Each network group should have a single cluster label
  expect_equal(length(unique(labels[1:3])), 1L)
  expect_equal(length(unique(labels[4:6])), 1L)
  expect_equal(length(unique(labels[7:9])), 1L)
})

test_that("spectral_ward_hclust merge matrix has correct dimensions", {
  skip_if_not_installed("RSpectra")

  set.seed(1111)
  n <- 10
  W <- matrix(runif(n * n), n, n)
  W <- (W + t(W)) / 2
  diag(W) <- 0

  hc <- spectral_ward_hclust(W, k_embed = 3)
  # merge matrix should be (n-1) x 2
  expect_equal(dim(hc$merge), c(n - 1, 2))
  # height should be length n-1
  expect_equal(length(hc$height), n - 1)
  # heights should be non-decreasing
  expect_true(all(diff(hc$height) >= -1e-10))
})

# =============================================================================
# Additional coverage tests for .parent_map_one
# =============================================================================

test_that(".parent_map_one returns correct map for simple case", {
  child <- c(1L, 1L, 2L, 2L)
  parent <- c(10L, 10L, 20L, 20L)

  result <- fmrilatent:::.parent_map_one(child, parent, lvl = 2)
  expect_equal(result[["1"]], 10L)
  expect_equal(result[["2"]], 20L)
})

test_that(".parent_map_one includes error message with level number", {
  child <- c(1L, 2L, 2L, 3L)
  parent <- c(1L, 1L, 2L, 2L)

  expect_error(
    fmrilatent:::.parent_map_one(child, parent, lvl = 5),
    "level 5"
  )
})

test_that(".parent_map_one handles all-NA parent for a child (error)", {
  child <- c(1L, 2L, 2L)
  parent <- c(1L, NA_integer_, NA_integer_)

  # Child 2 maps only to NAs -> after filtering, 0 unique parents -> error
  expect_error(
    fmrilatent:::.parent_map_one(child, parent, lvl = 2),
    "not nested"
  )
})

test_that(".parent_map_one with many children mapping to one parent", {
  child <- c(1L, 2L, 3L, 4L, 5L)
  parent <- rep(1L, 5)

  result <- fmrilatent:::.parent_map_one(child, parent, lvl = 2)
  expect_length(result, 5L)
  expect_true(all(result == 1L))
})

# =============================================================================
# Additional integration tests
# =============================================================================

test_that("spectral_ward_hclust + cut_hclust_nested pipeline produces valid hierarchy", {
  skip_if_not_installed("RSpectra")

  set.seed(1234)
  n <- 15

  # Build a similarity matrix with 3 clear groups
  W <- matrix(0.1, n, n)
  W[1:5, 1:5] <- 0.8
  W[6:10, 6:10] <- 0.8
  W[11:15, 11:15] <- 0.8
  diag(W) <- 0

  hc <- spectral_ward_hclust(W, k_embed = 3)
  levels <- cut_hclust_nested(hc, c(3, 5, 10))

  # Validate nesting
  expect_silent(validate_nested_parcellations(levels))

  # Get parent maps
  parents <- parent_maps_from_levels(levels)
  expect_length(parents, 3L)
  expect_length(parents[[1]], 0L)
  # All child IDs at level 2 should map to a parent at level 1
  expect_true(all(parents[[2]] %in% unique(levels[[1]])))
  expect_true(all(parents[[3]] %in% unique(levels[[2]])))
})

test_that("parent_maps_from_levels handles identical adjacent levels", {
  # Two levels that are identical
  labels <- c(1L, 1L, 2L, 2L, 3L, 3L)
  levels <- list(labels, labels)

  parents <- parent_maps_from_levels(levels)
  # Each cluster maps to itself
  expect_equal(parents[[2]][["1"]], 1L)
  expect_equal(parents[[2]][["2"]], 2L)
  expect_equal(parents[[2]][["3"]], 3L)
})

test_that("validate_nested_parcellations handles large number of levels", {
  # 5 levels, each a refinement
  L0 <- rep(1L, 16)
  L1 <- rep(1:2, each = 8)
  L2 <- rep(1:4, each = 4)
  L3 <- rep(1:8, each = 2)
  L4 <- 1:16

  levels <- list(L0, L1, L2, L3, L4)
  expect_silent(validate_nested_parcellations(levels))

  parents <- parent_maps_from_levels(levels)
  expect_length(parents, 5L)
  # At level 5 (L4), each of 16 children maps to one of 8 parents (L3)
  expect_length(parents[[5]], 16L)
})
