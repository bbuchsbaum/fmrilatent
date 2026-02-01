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
