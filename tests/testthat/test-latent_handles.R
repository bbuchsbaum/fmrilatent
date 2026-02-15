# Tests for latent_handles.R - Registry lifecycle API

# -----------------------------------------------------------------------------
# Helper to ensure clean registry state between tests
# -----------------------------------------------------------------------------

setup_clean_registry <- function() {
  fmrilatent_registry_clear()
}

# -----------------------------------------------------------------------------
# Tests for fmrilatent_registry_clear
# -----------------------------------------------------------------------------

test_that("fmrilatent_registry_clear clears all registries by default", {
  setup_clean_registry()

  # Add some items to the registry
  fmrilatent:::.latent_register_matrix("test_basis_1", matrix(1:4, 2), "basis")
  fmrilatent:::.latent_register_matrix("test_loadings_1", matrix(1:6, 3), "loadings")

  # Verify they exist
  expect_true(fmrilatent:::.latent_has_matrix("test_basis_1", "basis"))
  expect_true(fmrilatent:::.latent_has_matrix("test_loadings_1", "loadings"))

  # Clear all
  removed <- fmrilatent_registry_clear()

  # Verify removed
  expect_false(fmrilatent:::.latent_has_matrix("test_basis_1", "basis"))
  expect_false(fmrilatent:::.latent_has_matrix("test_loadings_1", "loadings"))

  # Returned count should be 2

  expect_equal(removed, 2L)
})

test_that("fmrilatent_registry_clear('basis') only clears basis registry", {
  setup_clean_registry()

  # Add items to both registries
  fmrilatent:::.latent_register_matrix("test_basis", matrix(1:4, 2), "basis")
  fmrilatent:::.latent_register_matrix("test_loadings", matrix(1:6, 3), "loadings")

  # Clear only basis
  removed <- fmrilatent_registry_clear("basis")

  # Basis should be cleared
  expect_false(fmrilatent:::.latent_has_matrix("test_basis", "basis"))

  # Loadings should still exist
  expect_true(fmrilatent:::.latent_has_matrix("test_loadings", "loadings"))

  expect_equal(removed, 1L)

  # Cleanup
  fmrilatent_registry_clear()
})

test_that("fmrilatent_registry_clear('loadings') only clears loadings registry", {
  setup_clean_registry()

  # Add items to both registries
  fmrilatent:::.latent_register_matrix("test_basis", matrix(1:4, 2), "basis")
  fmrilatent:::.latent_register_matrix("test_loadings", matrix(1:6, 3), "loadings")

  # Clear only loadings
  removed <- fmrilatent_registry_clear("loadings")

  # Basis should still exist
  expect_true(fmrilatent:::.latent_has_matrix("test_basis", "basis"))

  # Loadings should be cleared
  expect_false(fmrilatent:::.latent_has_matrix("test_loadings", "loadings"))

  expect_equal(removed, 1L)

  # Cleanup
  fmrilatent_registry_clear()
})

test_that("fmrilatent_registry_clear returns 0 when registry is empty", {
  setup_clean_registry()

  removed <- fmrilatent_registry_clear()

  expect_equal(removed, 0L)
})

# -----------------------------------------------------------------------------
# Tests for fmrilatent_registry_list
# -----------------------------------------------------------------------------

test_that("fmrilatent_registry_list returns empty vector when registry is empty", {
  setup_clean_registry()

  ids <- fmrilatent_registry_list()

  expect_length(ids, 0L)
  expect_type(ids, "character")
})

test_that("fmrilatent_registry_list returns all IDs by default", {
  setup_clean_registry()

  # Add items
  fmrilatent:::.latent_register_matrix("basis_1", matrix(1:4, 2), "basis")
  fmrilatent:::.latent_register_matrix("basis_2", matrix(1:4, 2), "basis")
  fmrilatent:::.latent_register_matrix("loadings_1", matrix(1:6, 3), "loadings")

  ids <- fmrilatent_registry_list()

  expect_length(ids, 3L)
  expect_true("basis_1" %in% ids)
  expect_true("basis_2" %in% ids)
  expect_true("loadings_1" %in% ids)

  # Cleanup
  fmrilatent_registry_clear()
})

test_that("fmrilatent_registry_list('basis') returns only basis IDs", {
  setup_clean_registry()

  fmrilatent:::.latent_register_matrix("my_basis", matrix(1:4, 2), "basis")
  fmrilatent:::.latent_register_matrix("my_loadings", matrix(1:6, 3), "loadings")

  ids <- fmrilatent_registry_list("basis")

  expect_length(ids, 1L)
  expect_equal(ids[[1]], "my_basis")
  expect_equal(names(ids), "basis")

  # Cleanup
  fmrilatent_registry_clear()
})

test_that("fmrilatent_registry_list('loadings') returns only loadings IDs", {
  setup_clean_registry()

  fmrilatent:::.latent_register_matrix("my_basis", matrix(1:4, 2), "basis")
  fmrilatent:::.latent_register_matrix("my_loadings", matrix(1:6, 3), "loadings")

  ids <- fmrilatent_registry_list("loadings")

  expect_length(ids, 1L)
  expect_equal(ids[[1]], "my_loadings")
  expect_equal(names(ids), "loadings")

  # Cleanup
  fmrilatent_registry_clear()
})

test_that("fmrilatent_registry_list names indicate registry type", {
  setup_clean_registry()

  fmrilatent:::.latent_register_matrix("b1", matrix(1:4, 2), "basis")
  fmrilatent:::.latent_register_matrix("l1", matrix(1:6, 3), "loadings")

  ids <- fmrilatent_registry_list("all")

  # All IDs should have names indicating their type
  expect_true(all(names(ids) %in% c("basis", "loadings")))

  # Cleanup
  fmrilatent_registry_clear()
})

# -----------------------------------------------------------------------------
# Tests for fmrilatent_registry_stats
# -----------------------------------------------------------------------------

test_that("fmrilatent_registry_stats returns zero counts for empty registry", {
  setup_clean_registry()

  stats <- fmrilatent_registry_stats()

  expect_equal(stats$basis$count, 0L)
  expect_equal(stats$loadings$count, 0L)
  expect_equal(stats$total$count, 0L)

  expect_equal(stats$basis$bytes, 0)
  expect_equal(stats$loadings$bytes, 0)
  expect_equal(stats$total$bytes, 0)
})

test_that("fmrilatent_registry_stats returns correct counts", {
  setup_clean_registry()

  fmrilatent:::.latent_register_matrix("b1", matrix(1:4, 2), "basis")
  fmrilatent:::.latent_register_matrix("b2", matrix(1:4, 2), "basis")
  fmrilatent:::.latent_register_matrix("l1", matrix(1:6, 3), "loadings")

  stats <- fmrilatent_registry_stats()

  expect_equal(stats$basis$count, 2L)
  expect_equal(stats$loadings$count, 1L)
  expect_equal(stats$total$count, 3L)

  # Cleanup
  fmrilatent_registry_clear()
})

test_that("fmrilatent_registry_stats reports non-zero bytes for cached objects", {
  setup_clean_registry()

  # Create a larger matrix to ensure measurable bytes
  big_matrix <- matrix(rnorm(1000), 100, 10)
  fmrilatent:::.latent_register_matrix("big_one", big_matrix, "basis")

  stats <- fmrilatent_registry_stats()

  expect_gt(stats$basis$bytes, 0)
  expect_gt(stats$total$bytes, 0)

  # Cleanup
  fmrilatent_registry_clear()
})

test_that("fmrilatent_registry_stats('basis') only reports basis stats", {
  setup_clean_registry()

  fmrilatent:::.latent_register_matrix("b1", matrix(1:4, 2), "basis")
  fmrilatent:::.latent_register_matrix("l1", matrix(1:6, 3), "loadings")

  stats <- fmrilatent_registry_stats("basis")

  expect_true("basis" %in% names(stats))
  expect_false("loadings" %in% names(stats))
  expect_false("total" %in% names(stats))

  # Cleanup
  fmrilatent_registry_clear()
})

test_that("fmrilatent_registry_stats('loadings') only reports loadings stats", {
  setup_clean_registry()

  fmrilatent:::.latent_register_matrix("b1", matrix(1:4, 2), "basis")
  fmrilatent:::.latent_register_matrix("l1", matrix(1:6, 3), "loadings")

  stats <- fmrilatent_registry_stats("loadings")

  expect_false("basis" %in% names(stats))
  expect_true("loadings" %in% names(stats))
  expect_false("total" %in% names(stats))

  # Cleanup
  fmrilatent_registry_clear()
})

# -----------------------------------------------------------------------------
# Tests for cache hit/miss behavior
# -----------------------------------------------------------------------------

test_that("registry caches matrix on first access and returns same object", {
  setup_clean_registry()

  test_matrix <- matrix(rnorm(100), 10, 10)
  test_id <- "cache_test_1"

  # Register matrix
  fmrilatent:::.latent_register_matrix(test_id, test_matrix, "basis")

  # First retrieval
  cached1 <- fmrilatent:::.latent_get_matrix(test_id, "basis")

  # Second retrieval should return identical object
  cached2 <- fmrilatent:::.latent_get_matrix(test_id, "basis")

  expect_identical(cached1, cached2)
  expect_equal(cached1, test_matrix)

  # Cleanup
  fmrilatent_registry_clear()
})

test_that("clear forces re-registration (no stale cache)", {
  setup_clean_registry()

  test_id <- "refresh_test"
  matrix_v1 <- matrix(1:4, 2)
  matrix_v2 <- matrix(5:8, 2)

  # Register first version
  fmrilatent:::.latent_register_matrix(test_id, matrix_v1, "basis")
  cached_v1 <- fmrilatent:::.latent_get_matrix(test_id, "basis")
  expect_equal(cached_v1, matrix_v1)

  # Clear and register second version
  fmrilatent_registry_clear()
  fmrilatent:::.latent_register_matrix(test_id, matrix_v2, "basis")
  cached_v2 <- fmrilatent:::.latent_get_matrix(test_id, "basis")

  # Should get version 2, not version 1
  expect_equal(cached_v2, matrix_v2)

  # Cleanup
  fmrilatent_registry_clear()
})

test_that("cache miss returns NULL", {
  setup_clean_registry()

  result <- fmrilatent:::.latent_get_matrix("nonexistent_id", "basis")

  expect_null(result)
})

test_that(".latent_has_matrix returns TRUE for cached, FALSE for uncached", {
  setup_clean_registry()

  fmrilatent:::.latent_register_matrix("exists", matrix(1:4, 2), "basis")

  expect_true(fmrilatent:::.latent_has_matrix("exists", "basis"))
  expect_false(fmrilatent:::.latent_has_matrix("not_exists", "basis"))

  # Cleanup
  fmrilatent_registry_clear()
})

test_that("overwrite=FALSE does not replace existing entry", {
  setup_clean_registry()

  test_id <- "no_overwrite"
  matrix_v1 <- matrix(1:4, 2)
  matrix_v2 <- matrix(5:8, 2)

  # Register first version
  fmrilatent:::.latent_register_matrix(test_id, matrix_v1, "basis")

  # Attempt to register second version without overwrite
  expect_warning(
    result <- fmrilatent:::.latent_register_matrix(test_id, matrix_v2, "basis", overwrite = FALSE),
    "already registered"
  )

  expect_false(result)

  # Should still have version 1
  cached <- fmrilatent:::.latent_get_matrix(test_id, "basis")
  expect_equal(cached, matrix_v1)

  # Cleanup
  fmrilatent_registry_clear()
})

test_that("overwrite=TRUE replaces existing entry", {
  setup_clean_registry()

  test_id <- "do_overwrite"
  matrix_v1 <- matrix(1:4, 2)
  matrix_v2 <- matrix(5:8, 2)

  # Register first version
  fmrilatent:::.latent_register_matrix(test_id, matrix_v1, "basis")

  # Register second version with overwrite
  result <- fmrilatent:::.latent_register_matrix(test_id, matrix_v2, "basis", overwrite = TRUE)

  expect_true(result)

  # Should now have version 2
  cached <- fmrilatent:::.latent_get_matrix(test_id, "basis")
  expect_equal(cached, matrix_v2)

  # Cleanup
  fmrilatent_registry_clear()
})

# -----------------------------------------------------------------------------
# Tests for internal registry helper functions
# -----------------------------------------------------------------------------

test_that(".latent_register_matrix rejects invalid id", {
  expect_error(
    fmrilatent:::.latent_register_matrix(123, matrix(1:4, 2), "basis"),
    "is.character"
  )

  expect_error(
    fmrilatent:::.latent_register_matrix(c("a", "b"), matrix(1:4, 2), "basis"),
    "length"
  )
})

test_that(".latent_get_matrix rejects invalid id", {
  expect_error(
    fmrilatent:::.latent_get_matrix(123, "basis"),
    "is.character"
  )

  expect_error(
    fmrilatent:::.latent_get_matrix(c("a", "b"), "basis"),
    "length"
  )
})

test_that(".latent_get_registry_env returns correct environments", {
  basis_env <- fmrilatent:::.latent_get_registry_env("basis")
  loadings_env <- fmrilatent:::.latent_get_registry_env("loadings")

  expect_true(is.environment(basis_env))
  expect_true(is.environment(loadings_env))

  # They should be different environments
  expect_false(identical(basis_env, loadings_env))
})

test_that(".latent_get_registry_env rejects invalid type", {
  expect_error(
    fmrilatent:::.latent_get_registry_env("invalid"),
    "'arg' should be one of"
  )
})

# -----------------------------------------------------------------------------
# Tests for registry enable/disable API
# -----------------------------------------------------------------------------

test_that("fmrilatent_registry_enable enables the registry", {
  # Save original state
  original_state <- getOption("fmrilatent.registry.enabled")
  on.exit(options(fmrilatent.registry.enabled = original_state), add = TRUE)

  # Disable first
  fmrilatent_registry_disable()
  expect_false(fmrilatent_registry_enabled())

  # Now enable
  result <- fmrilatent_registry_enable()
  expect_true(result)
  expect_true(fmrilatent_registry_enabled())
})

test_that("fmrilatent_registry_disable disables the registry", {
  # Save original state
  original_state <- getOption("fmrilatent.registry.enabled")
  on.exit(options(fmrilatent.registry.enabled = original_state), add = TRUE)

  # Enable first
  fmrilatent_registry_enable()
  expect_true(fmrilatent_registry_enabled())

  # Now disable
  result <- fmrilatent_registry_disable()
  expect_true(result)
  expect_false(fmrilatent_registry_enabled())
})

test_that("fmrilatent_registry_enabled returns correct state", {
  # Save original state
  original_state <- getOption("fmrilatent.registry.enabled")
  on.exit(options(fmrilatent.registry.enabled = original_state), add = TRUE)

  # Test enabled state
  options(fmrilatent.registry.enabled = TRUE)
  expect_true(fmrilatent_registry_enabled())

  # Test disabled state
  options(fmrilatent.registry.enabled = FALSE)
  expect_false(fmrilatent_registry_enabled())

  # Test NULL/missing state (defaults to TRUE)
  options(fmrilatent.registry.enabled = NULL)
  expect_true(fmrilatent_registry_enabled())
})

test_that(".latent_register_matrix returns FALSE when registry is disabled", {
  # Save original state
  original_state <- getOption("fmrilatent.registry.enabled")
  on.exit(options(fmrilatent.registry.enabled = original_state), add = TRUE)

  setup_clean_registry()
  fmrilatent_registry_disable()

  result <- fmrilatent:::.latent_register_matrix("test_id", matrix(1:4, 2), "basis")

  expect_false(result)

  # Re-enable to verify it wasn't registered
  fmrilatent_registry_enable()
  expect_false(fmrilatent:::.latent_has_matrix("test_id", "basis"))
})

test_that(".latent_get_matrix returns NULL when registry is disabled", {
  # Save original state
  original_state <- getOption("fmrilatent.registry.enabled")
  on.exit(options(fmrilatent.registry.enabled = original_state), add = TRUE)

  setup_clean_registry()

  # Register with registry enabled
  fmrilatent_registry_enable()
  fmrilatent:::.latent_register_matrix("test_id", matrix(1:4, 2), "basis")
  expect_true(fmrilatent:::.latent_has_matrix("test_id", "basis"))

  # Disable registry
  fmrilatent_registry_disable()

  # Should return NULL even though item exists in environment
  result <- fmrilatent:::.latent_get_matrix("test_id", "basis")
  expect_null(result)

  # Cleanup
  fmrilatent_registry_enable()
  fmrilatent_registry_clear()
})

test_that(".latent_get_matrix returns NULL when id doesn't exist", {
  # Save original state
  original_state <- getOption("fmrilatent.registry.enabled")
  on.exit(options(fmrilatent.registry.enabled = original_state), add = TRUE)

  setup_clean_registry()
  fmrilatent_registry_enable()

  result <- fmrilatent:::.latent_get_matrix("nonexistent_matrix_id", "basis")
  expect_null(result)

  result2 <- fmrilatent:::.latent_get_matrix("another_missing_id", "loadings")
  expect_null(result2)
})

# -----------------------------------------------------------------------------
# Tests for mask_to_array
# -----------------------------------------------------------------------------

test_that("mask_to_array converts array to array", {
  arr <- array(c(TRUE, FALSE, TRUE, FALSE), dim = c(2, 2))
  result <- mask_to_array(arr)
  expect_identical(result, arr)
})

test_that("mask_to_array converts matrix to array", {
  mat <- matrix(c(TRUE, FALSE, TRUE, FALSE), 2, 2)
  result <- mask_to_array(mat)
  expect_true(is.array(result))
  expect_equal(dim(result), c(2, 2))
})

test_that("mask_to_array errors with informative message when conversion fails", {
  # Use an environment object which genuinely fails as.array conversion
  bad_obj <- new.env()

  expect_error(
    mask_to_array(bad_obj, location = "test_function"),
    "In test_function: mask must be array-like or LogicalNeuroVol"
  )

  expect_error(
    mask_to_array(bad_obj, location = "test_function"),
    "Underlying error:"
  )
})

test_that("mask_to_array errors when as.array returns NULL", {
  # Define a mock class that returns NULL from as.array (must be global)
  assign("as.array.mock_null_test", function(x, ...) NULL, envir = .GlobalEnv)
  on.exit(rm("as.array.mock_null_test", envir = .GlobalEnv), add = TRUE)

  mock_obj <- structure(list(), class = "mock_null_test")

  expect_error(
    mask_to_array(mock_obj, location = "my_test_func"),
    "In my_test_func: mask must be array-like or LogicalNeuroVol"
  )

  expect_error(
    mask_to_array(mock_obj, location = "my_test_func"),
    "conversion returned NULL"
  )
})

test_that("mask_to_array uses location parameter in error messages", {
  # Use a function object which genuinely fails as.array conversion
  bad_obj <- function() {}

  expect_error(
    mask_to_array(bad_obj, location = "custom_location_name"),
    "In custom_location_name:"
  )
})

# -----------------------------------------------------------------------------
# Tests for .latent_basis_dim
# -----------------------------------------------------------------------------

test_that(".latent_basis_dim returns correct dimensions for matrix", {
  mat <- matrix(1:20, 10, 2)
  result <- fmrilatent:::.latent_basis_dim(mat)
  expect_equal(result, c(10, 2))
})

test_that(".latent_basis_dim returns correct dimensions for 2D array", {
  arr <- array(1:20, dim = c(10, 2))
  result <- fmrilatent:::.latent_basis_dim(arr)
  expect_equal(result, c(10, 2))
})

test_that(".latent_basis_dim returns correct dimensions for BasisHandle", {
  # Create a BasisHandle with known dimensions
  handle <- new("BasisHandle",
                id = "test_basis",
                dim = as.integer(c(100, 5)),
                kind = "dct",
                spec = list(),
                label = "test")

  result <- fmrilatent:::.latent_basis_dim(handle)
  expect_equal(result, c(100, 5))
})

test_that(".latent_basis_dim errors on unsupported type", {
  # Test with a vector (1D array)
  vec <- 1:10
  expect_error(
    fmrilatent:::.latent_basis_dim(vec),
    "Unsupported basis slot type"
  )

  # Test with a list
  lst <- list(a = 1, b = 2)
  expect_error(
    fmrilatent:::.latent_basis_dim(lst),
    "Unsupported basis slot type"
  )

  # Test with a 3D array
  arr3d <- array(1:24, dim = c(2, 3, 4))
  expect_error(
    fmrilatent:::.latent_basis_dim(arr3d),
    "Unsupported basis slot type"
  )
})

# -----------------------------------------------------------------------------
# Tests for .latent_loadings_dim
# -----------------------------------------------------------------------------

test_that(".latent_loadings_dim returns correct dimensions for matrix", {
  mat <- matrix(1:30, 15, 2)
  result <- fmrilatent:::.latent_loadings_dim(mat)
  expect_equal(result, c(15, 2))
})

test_that(".latent_loadings_dim returns correct dimensions for 2D array", {
  arr <- array(1:30, dim = c(15, 2))
  result <- fmrilatent:::.latent_loadings_dim(arr)
  expect_equal(result, c(15, 2))
})

test_that(".latent_loadings_dim returns correct dimensions for LoadingsHandle", {
  # Create a LoadingsHandle with known dimensions
  handle <- new("LoadingsHandle",
                id = "test_loadings",
                dim = as.integer(c(1000, 5)),
                kind = "explicit",
                spec = list(),
                label = "test")

  result <- fmrilatent:::.latent_loadings_dim(handle)
  expect_equal(result, c(1000, 5))
})

test_that(".latent_loadings_dim errors on unsupported type", {
  # Test with a vector (1D array)
  vec <- 1:10
  expect_error(
    fmrilatent:::.latent_loadings_dim(vec),
    "Unsupported loadings slot type"
  )

  # Test with a list
  lst <- list(a = 1, b = 2)
  expect_error(
    fmrilatent:::.latent_loadings_dim(lst),
    "Unsupported loadings slot type"
  )

  # Test with a 3D array
  arr3d <- array(1:24, dim = c(2, 3, 4))
  expect_error(
    fmrilatent:::.latent_loadings_dim(arr3d),
    "Unsupported loadings slot type"
  )
})
