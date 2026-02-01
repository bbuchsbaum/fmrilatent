# Tests for fmrilatent registry lifecycle API

test_that("fmrilatent_registry_clear removes all entries", {
  # Start fresh
  fmrilatent_registry_clear()
  
  # Verify empty
  stats <- fmrilatent_registry_stats()
  expect_equal(stats$total$count, 0)
  
  # Register some matrices
  .latent_register_matrix("test-basis-1", matrix(1:4, 2, 2), type = "basis")
  .latent_register_matrix("test-basis-2", matrix(1:6, 2, 3), type = "basis")
  .latent_register_matrix("test-loadings-1", matrix(1:9, 3, 3), type = "loadings")
  
  # Verify populated
  stats <- fmrilatent_registry_stats()
  expect_equal(stats$basis$count, 2)
  expect_equal(stats$loadings$count, 1)
  expect_equal(stats$total$count, 3)
  
  # Clear all
  removed <- fmrilatent_registry_clear()
  expect_equal(removed, 3)
  
  # Verify empty again
  stats <- fmrilatent_registry_stats()
  expect_equal(stats$total$count, 0)
})

test_that("fmrilatent_registry_clear respects type argument", {
  fmrilatent_registry_clear()
  
  .latent_register_matrix("test-b", matrix(1, 2, 2), type = "basis")
  .latent_register_matrix("test-l", matrix(1, 3, 3), type = "loadings")
  
  # Clear only basis
  removed <- fmrilatent_registry_clear("basis")
  expect_equal(removed, 1)
  
  stats <- fmrilatent_registry_stats()
  expect_equal(stats$basis$count, 0)
  expect_equal(stats$loadings$count, 1)
  
  # Cleanup
  fmrilatent_registry_clear()
})

test_that("fmrilatent_registry_list returns correct IDs", {
  fmrilatent_registry_clear()
  
  .latent_register_matrix("alpha", matrix(1, 2, 2), type = "basis")
  .latent_register_matrix("beta", matrix(1, 2, 2), type = "basis")
  .latent_register_matrix("gamma", matrix(1, 3, 3), type = "loadings")
  
  # List all
  all_ids <- fmrilatent_registry_list()
  expect_length(all_ids, 3)
  expect_true("alpha" %in% all_ids)
  expect_true("beta" %in% all_ids)
  expect_true("gamma" %in% all_ids)
  
  # List by type
  basis_ids <- fmrilatent_registry_list("basis")
  expect_length(basis_ids, 2)
  expect_true(all(names(basis_ids) == "basis"))
  
  loadings_ids <- fmrilatent_registry_list("loadings")
  expect_length(loadings_ids, 1)
  expect_equal(loadings_ids[[1]], "gamma")
  
  fmrilatent_registry_clear()
})

test_that("fmrilatent_registry_stats reports memory usage", {
  fmrilatent_registry_clear()
  
  # Register a matrix
  m <- matrix(rnorm(1000), 100, 10)
  .latent_register_matrix("large-matrix", m, type = "basis")
  
  stats <- fmrilatent_registry_stats()
  
  expect_equal(stats$basis$count, 1)
  expect_gt(stats$basis$bytes, 0)
  expect_equal(stats$loadings$count, 0)
  expect_equal(stats$loadings$bytes, 0)
  
  fmrilatent_registry_clear()
})

test_that("clearing registry forces handle re-materialization", {
  skip_if_not_installed("neuroim2")
  
  fmrilatent_registry_clear()
  
  # Create a handle
  handle <- dct_basis_handle(n_time = 10L, k = 3L)
  
  # First materialization
  mat1 <- basis_mat(handle)
  expect_true(.latent_has_matrix(handle@id, "basis"))
  
  # Clear registry
  fmrilatent_registry_clear()
  expect_false(.latent_has_matrix(handle@id, "basis"))
  
  # Re-materialize
  mat2 <- basis_mat(handle)
  expect_true(.latent_has_matrix(handle@id, "basis"))
  
  # Results should be identical
  expect_equal(mat1, mat2)
  
  fmrilatent_registry_clear()
})

test_that("cache hit returns same object without recomputation", {
  fmrilatent_registry_clear()
  
  # Register a matrix
  m <- matrix(1:6, 2, 3)
  .latent_register_matrix("cached-matrix", m, type = "basis")
  
  # Get it twice
  retrieved1 <- .latent_get_matrix("cached-matrix", "basis")
  retrieved2 <- .latent_get_matrix("cached-matrix", "basis")
  
  # Should be identical (same object)
  expect_identical(retrieved1, retrieved2)
  
  fmrilatent_registry_clear()
})

test_that("cache miss returns NULL", {
  fmrilatent_registry_clear()
  
  result <- .latent_get_matrix("nonexistent-id", "basis")
  expect_null(result)
  
  result <- .latent_get_matrix("nonexistent-id", "loadings")
  expect_null(result)
})
