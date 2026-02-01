# Tests for implicit_latent.R
# Covers: implicit_latent constructor, predict.ImplicitLatent, is_implicit_latent, implicit_meta

library(testthat)

# =============================================================================
# Helper functions for creating test fixtures
# =============================================================================

#' Create a minimal ImplicitLatent object for testing
make_test_implicit_latent <- function(family = "test_family",
                                       coeff = list(data = 1:10),
                                       extra_meta = list()) {
  mask <- array(TRUE, dim = c(2, 2, 2))
  meta <- c(list(family = family), extra_meta)
  decoder <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL) {
    # Return a simple matrix with time x voxels
    n_time <- if (is.null(time_idx)) 5L else length(time_idx)
    n_vox <- sum(mask)
    mat <- matrix(seq_len(n_time * n_vox), nrow = n_time, ncol = n_vox)
    if (!is.null(time_idx)) {
      # Just return subset rows
      mat <- mat[seq_len(length(time_idx)), , drop = FALSE]
    }
    if (!is.null(roi_mask)) {
      # Subset columns based on roi
      keep <- which(as.logical(roi_mask)[which(mask)])
      if (length(keep) > 0) {
        mat <- mat[, keep, drop = FALSE]
      }
    }
    mat
  }
  implicit_latent(coeff, decoder, meta, mask)
}

#' Create an ImplicitLatent with a decoder that tracks its arguments
make_tracking_implicit_latent <- function() {
  mask <- array(TRUE, dim = c(2, 2, 2))
  meta <- list(family = "tracker")
  call_log <- list()

  decoder <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL) {
    # Store the call for later inspection
    call_info <- list(
      time_idx = time_idx,
      roi_mask = roi_mask,
      levels_keep = levels_keep
    )
    # Store in parent environment
    call_log <<- c(call_log, list(call_info))
    matrix(1:8, nrow = 1, ncol = 8)
  }

  obj <- implicit_latent(list(), decoder, meta, mask)
  # Attach call log to object for inspection
  obj$call_log <- call_log
  obj$get_calls <- function() call_log
  obj
}

# =============================================================================
# Tests for implicit_latent() constructor
# =============================================================================

test_that("implicit_latent creates object with correct class", {
  obj <- make_test_implicit_latent()
  expect_s3_class(obj, "ImplicitLatent")
  expect_true(inherits(obj, "ImplicitLatent"))
})
test_that("implicit_latent stores all required components", {
  coeff <- list(a = 1, b = 2)
  mask <- array(TRUE, dim = c(2, 2, 2))
  meta <- list(family = "test")
  decoder <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL) {
    matrix(0, 1, 1)
  }

  obj <- implicit_latent(coeff, decoder, meta, mask)

  expect_identical(obj$coeff, coeff)
  expect_identical(obj$meta, meta)

  expect_identical(obj$mask, mask)
  expect_true(is.function(obj$decoder))
})

test_that("implicit_latent requires meta$family", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  decoder <- function(...) matrix(0, 1, 1)

  # Missing family entirely
  expect_error(
    implicit_latent(list(), decoder, list(), mask),
    "meta\\$family required"
  )

  # family = NULL
  expect_error(
    implicit_latent(list(), decoder, list(family = NULL), mask),
    "meta\\$family required"
  )
})

test_that("implicit_latent accepts various coeff types", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  meta <- list(family = "test")
  decoder <- function(...) matrix(0, 1, 1)

  # List coefficient
  obj1 <- implicit_latent(list(x = 1:10), decoder, meta, mask)
  expect_type(obj1$coeff, "list")

  # Matrix coefficient
  obj2 <- implicit_latent(matrix(1:4, 2, 2), decoder, meta, mask)
  expect_true(is.matrix(obj2$coeff))

  # Empty list
  obj3 <- implicit_latent(list(), decoder, meta, mask)
  expect_equal(length(obj3$coeff), 0)
})

test_that("implicit_latent accepts various mask types", {
  meta <- list(family = "test")
  decoder <- function(...) matrix(0, 1, 1)

  # Logical 3D array
  mask1 <- array(TRUE, dim = c(2, 2, 2))
  obj1 <- implicit_latent(list(), decoder, meta, mask1)
  expect_true(is.array(obj1$mask))

  # Integer/numeric array that is truthy
  mask2 <- array(1L, dim = c(2, 2, 2))
  obj2 <- implicit_latent(list(), decoder, meta, mask2)
  expect_true(is.array(obj2$mask))

  # Mixed TRUE/FALSE mask
  mask3 <- array(FALSE, dim = c(3, 3, 3))
  mask3[1, 1, 1] <- TRUE
  mask3[2, 2, 2] <- TRUE
  obj3 <- implicit_latent(list(), decoder, meta, mask3)
  expect_equal(sum(obj3$mask), 2)
})

test_that("implicit_latent preserves meta with extra fields", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  meta <- list(
    family = "test",
    levels = 3,
    z_seed = 42L,
    custom_field = "custom_value"
  )
  decoder <- function(...) matrix(0, 1, 1)

  obj <- implicit_latent(list(), decoder, meta, mask)

  expect_equal(obj$meta$family, "test")
  expect_equal(obj$meta$levels, 3)
  expect_equal(obj$meta$z_seed, 42L)
  expect_equal(obj$meta$custom_field, "custom_value")
})

# =============================================================================
# Tests for predict.ImplicitLatent()
# =============================================================================

test_that("predict.ImplicitLatent calls decoder with no arguments", {
  obj <- make_test_implicit_latent()
  result <- predict(obj)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 5)  # default from our test decoder
  expect_equal(ncol(result), 8)  # 2x2x2 mask = 8 voxels
})

test_that("predict.ImplicitLatent passes time_idx to decoder", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  meta <- list(family = "test")

  # Decoder that uses time_idx
  decoder <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL) {
    n <- if (is.null(time_idx)) 10L else length(time_idx)
    matrix(seq_len(n * 8), nrow = n, ncol = 8)
  }

  obj <- implicit_latent(list(), decoder, meta, mask)

  # Full reconstruction
  full <- predict(obj)
  expect_equal(nrow(full), 10)

  # Subset time
  partial <- predict(obj, time_idx = c(1, 3, 5))
  expect_equal(nrow(partial), 3)
})

test_that("predict.ImplicitLatent passes roi_mask to decoder", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  meta <- list(family = "test")

  # Decoder that respects roi_mask
  decoder <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL) {
    n_time <- 5L
    n_vox <- if (is.null(roi_mask)) 8L else sum(as.logical(roi_mask))
    matrix(1, nrow = n_time, ncol = n_vox)
  }

  obj <- implicit_latent(list(), decoder, meta, mask)

  # Full ROI
  full <- predict(obj)
  expect_equal(ncol(full), 8)

  # Subset ROI (4 voxels)
  roi <- array(FALSE, dim = c(2, 2, 2))
  roi[1, , ] <- TRUE
  partial <- predict(obj, roi_mask = roi)
  expect_equal(ncol(partial), 4)
})

test_that("predict.ImplicitLatent passes levels_keep to decoder", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  meta <- list(family = "test")

  received_levels <- NULL
  decoder <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL) {
    received_levels <<- levels_keep
    matrix(1, nrow = 5, ncol = 8)
  }

  obj <- implicit_latent(list(), decoder, meta, mask)

  # Without levels_keep
  predict(obj)
  expect_null(received_levels)

  # With levels_keep
  predict(obj, levels_keep = c(1, 2))
  expect_equal(received_levels, c(1, 2))

  # Empty levels_keep
  predict(obj, levels_keep = integer(0))
  expect_equal(received_levels, integer(0))
})

test_that("predict.ImplicitLatent combines multiple arguments", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  meta <- list(family = "test")

  call_args <- NULL
  decoder <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL) {
    call_args <<- list(
      time_idx = time_idx,
      roi_mask = roi_mask,
      levels_keep = levels_keep
    )
    n_time <- if (is.null(time_idx)) 10L else length(time_idx)
    n_vox <- if (is.null(roi_mask)) 8L else sum(as.logical(roi_mask))
    matrix(1, nrow = n_time, ncol = n_vox)
  }

  obj <- implicit_latent(list(), decoder, meta, mask)

  roi <- array(FALSE, dim = c(2, 2, 2))
  roi[1, 1, 1] <- TRUE
  roi[2, 2, 2] <- TRUE

  result <- predict(obj, time_idx = 1:3, roi_mask = roi, levels_keep = c(1))

  expect_equal(call_args$time_idx, 1:3)
  expect_equal(sum(call_args$roi_mask), 2)
  expect_equal(call_args$levels_keep, c(1))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 2)
})

test_that("predict.ImplicitLatent ignores extra ... arguments gracefully", {
  obj <- make_test_implicit_latent()
  # Should not error with extra arguments
  result <- predict(obj, extra_arg = "ignored", another = 123)
  expect_true(is.matrix(result))
})

# =============================================================================
# Tests for is_implicit_latent()
# =============================================================================

test_that("is_implicit_latent returns TRUE for ImplicitLatent objects", {
  obj <- make_test_implicit_latent()
  expect_true(is_implicit_latent(obj))
})

test_that("is_implicit_latent returns TRUE for subclasses of ImplicitLatent", {
  obj <- make_test_implicit_latent()
  class(obj) <- c("CustomLatent", "ImplicitLatent")
  expect_true(is_implicit_latent(obj))

  # HaarLatent style subclass
  class(obj) <- c("HaarLatent", "ImplicitLatent")
  expect_true(is_implicit_latent(obj))

  # Multiple inheritance
  class(obj) <- c("A", "B", "ImplicitLatent", "C")
  expect_true(is_implicit_latent(obj))
})

test_that("is_implicit_latent returns FALSE for non-ImplicitLatent objects", {
  expect_false(is_implicit_latent(NULL))
  expect_false(is_implicit_latent(list()))
  expect_false(is_implicit_latent(list(coeff = 1, decoder = function() {}, meta = list())))
  expect_false(is_implicit_latent(1:10))
  expect_false(is_implicit_latent("string"))
  expect_false(is_implicit_latent(matrix(1:4, 2, 2)))
  expect_false(is_implicit_latent(data.frame(a = 1:3)))
  expect_false(is_implicit_latent(function() {}))
  expect_false(is_implicit_latent(TRUE))
  expect_false(is_implicit_latent(NA))
})

test_that("is_implicit_latent returns logical of length 1", {
  obj <- make_test_implicit_latent()
  result <- is_implicit_latent(obj)

  expect_type(result, "logical")
  expect_length(result, 1)
  expect_false(is.na(result))
})

test_that("is_implicit_latent works with objects having similar structure but wrong class", {
  # Object with all the right slots but not the class
  fake <- list(
    coeff = list(),
    decoder = function(...) matrix(0, 1, 1),
    meta = list(family = "fake"),
    mask = array(TRUE, c(2, 2, 2))
  )
  expect_false(is_implicit_latent(fake))

  # Wrong class name (similar but different)
  class(fake) <- "ImplicitLatentX"
  expect_false(is_implicit_latent(fake))

  class(fake) <- "implicit_latent"
  expect_false(is_implicit_latent(fake))
})

# =============================================================================
# Tests for implicit_meta()
# =============================================================================

test_that("implicit_meta extracts meta from ImplicitLatent", {
  obj <- make_test_implicit_latent(family = "test_family", extra_meta = list(levels = 3))
  meta <- implicit_meta(obj)

  expect_type(meta, "list")
  expect_equal(meta$family, "test_family")
  expect_equal(meta$levels, 3)
})

test_that("implicit_meta returns NULL for non-ImplicitLatent objects", {
  expect_null(implicit_meta(NULL))
  expect_null(implicit_meta(list()))
  expect_null(implicit_meta(1:10))
  expect_null(implicit_meta("string"))
  expect_null(implicit_meta(matrix(1:4, 2, 2)))
  expect_null(implicit_meta(data.frame(a = 1)))
})

test_that("implicit_meta returns NULL for object without meta slot", {
  obj <- structure(list(coeff = list(), decoder = function() {}),
                   class = "ImplicitLatent")
  # Missing meta slot
  result <- implicit_meta(obj)
  expect_null(result)
})

test_that("implicit_meta returns NULL for object with NULL meta", {
  obj <- structure(list(coeff = list(), decoder = function() {}, meta = NULL),
                   class = "ImplicitLatent")
  result <- implicit_meta(obj)
  expect_null(result)
})

test_that("implicit_meta preserves all meta fields", {
  extra <- list(
    levels = 5,
    z_seed = 123L,
    mask_dims = c(10, 10, 10),
    custom = list(a = 1, b = 2)
  )
  obj <- make_test_implicit_latent(family = "complex", extra_meta = extra)
  meta <- implicit_meta(obj)

  expect_equal(meta$family, "complex")
  expect_equal(meta$levels, 5)
  expect_equal(meta$z_seed, 123L)
  expect_equal(meta$mask_dims, c(10, 10, 10))
  expect_equal(meta$custom, list(a = 1, b = 2))
})

test_that("implicit_meta works for subclassed ImplicitLatent", {
  obj <- make_test_implicit_latent(family = "haar")
  class(obj) <- c("HaarLatent", "ImplicitLatent")

  meta <- implicit_meta(obj)
  expect_equal(meta$family, "haar")
})

# =============================================================================
# Tests for edge cases and robustness
# =============================================================================

test_that("implicit_latent handles empty mask correctly", {
  # Note: The function doesn't validate the mask, so an empty mask is allowed
  # but may cause issues in the decoder
  mask <- array(FALSE, dim = c(2, 2, 2))
  meta <- list(family = "test")
  decoder <- function(...) matrix(numeric(0), nrow = 1, ncol = 0)

  obj <- implicit_latent(list(), decoder, meta, mask)
  expect_s3_class(obj, "ImplicitLatent")
  expect_equal(sum(obj$mask), 0)
})

test_that("implicit_latent handles large mask dimensions", {
  # Just verify construction works, not actual computation
  mask <- array(TRUE, dim = c(10, 10, 10))
  meta <- list(family = "test")
  decoder <- function(...) matrix(1, nrow = 1, ncol = 1000)

  obj <- implicit_latent(list(), decoder, meta, mask)
  expect_s3_class(obj, "ImplicitLatent")
  expect_equal(dim(obj$mask), c(10, 10, 10))
})

test_that("predict.ImplicitLatent handles decoder that returns vector",
  {
  mask <- array(TRUE, dim = c(2, 2, 2))
  meta <- list(family = "test")
  # Decoder returns a vector instead of matrix
  decoder <- function(...) 1:8

  obj <- implicit_latent(list(), decoder, meta, mask)
  # The predict method just passes through decoder output
  result <- predict(obj)
  expect_equal(result, 1:8)
})

test_that("predict.ImplicitLatent handles decoder that returns NULL", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  meta <- list(family = "test")
  decoder <- function(...) NULL

  obj <- implicit_latent(list(), decoder, meta, mask)
  result <- predict(obj)
  expect_null(result)
})

test_that("implicit_latent decoder can access coeff from closure", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  meta <- list(family = "test")

  # Create decoder that uses coeff
  coeff_data <- list(scale = 2, offset = 10)
  decoder <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL) {
    # This demonstrates that decoder can be a closure
    matrix(coeff_data$scale * 1:8 + coeff_data$offset, nrow = 1)
  }

  obj <- implicit_latent(coeff_data, decoder, meta, mask)
  result <- predict(obj)
  expected <- matrix(2 * 1:8 + 10, nrow = 1)
  expect_equal(result, expected)
})

# =============================================================================
# Tests for consistency with other fmrilatent functions
# =============================================================================

test_that("implicit_latent integrates with haar_latent pattern", {
  skip_if_not_installed("digest")
  # Skip if haar_wavelet_forward is not available (e.g. in isolated coverage runs)
  skip_if_not(exists("haar_wavelet_forward", mode = "function"),
              "haar_wavelet_forward not available")

  # This tests the pattern used by haar_latent
  mask <- array(TRUE, dim = c(2, 2, 2))
  X <- matrix(rnorm(5 * 8), nrow = 5)

  # Simulate what haar_latent does
  fw <- haar_wavelet_forward(X, mask, levels = 1, z_seed = 42L)
  meta <- fw$meta
  meta$family <- "haar"

  decoder <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL) {
    haar_wavelet_inverse(list(coeff = fw$coeff, meta = meta), mask,
                         levels = meta$levels, z_seed = meta$z_seed,
                         roi_mask = roi_mask, time_idx = time_idx,
                         levels_keep = levels_keep)
  }

  obj <- implicit_latent(fw$coeff, decoder, meta, mask)

  expect_true(is_implicit_latent(obj))
  expect_equal(implicit_meta(obj)$family, "haar")

  # Verify roundtrip
  reco <- predict(obj)
  expect_equal(reco, X, tolerance = 1e-6)
})

test_that("implicit_latent integrates with slepian_st pattern", {
  skip_if_not_installed("RSpectra")

  # This tests the pattern used by slepian_spatiotemporal_latent
  mask <- array(TRUE, dim = c(2, 2, 1))
  n_time <- 5
  n_vox <- sum(mask)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time)

  # Minimal separable construction
  B_t <- matrix(rnorm(n_time * 2), nrow = n_time)  # time basis
  L_s <- matrix(rnorm(n_vox * 2), nrow = n_vox)    # spatial loadings
  core <- crossprod(B_t, X) %*% L_s

  decoder <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL) {
    t_sel <- if (is.null(time_idx)) seq_len(n_time) else as.integer(time_idx)
    B_sel <- B_t[t_sel, , drop = FALSE]
    rec <- B_sel %*% core %*% t(L_s)
    if (!is.null(roi_mask)) {
      idx <- which(as.logical(roi_mask))
      rec <- rec[, idx, drop = FALSE]
    }
    rec
  }

  meta <- list(family = "st_separable", k_time = 2, k_space = 2)

  obj <- implicit_latent(
    coeff = list(core = core, B_t = B_t, L_s = L_s),
    decoder = decoder,
    meta = meta,
    mask = mask
  )

  expect_true(is_implicit_latent(obj))
  expect_equal(implicit_meta(obj)$family, "st_separable")

  # Verify dimensions
  reco <- predict(obj)
  expect_equal(dim(reco), dim(X))

  # Verify time subsetting
  reco_sub <- predict(obj, time_idx = c(1, 3))
  expect_equal(nrow(reco_sub), 2)
})

# =============================================================================
# Tests for %||% operator behavior (internal to the module)
# =============================================================================

test_that("implicit_meta handles missing meta gracefully via null coalescing", {
  # The implicit_meta function uses %||% internally
  # Test that it handles NULL meta$something correctly
  obj <- make_test_implicit_latent()

  # Note: In R, setting list element to NULL removes it
  # Instead test that accessing a non-existent field returns NULL
  meta <- implicit_meta(obj)
  expect_true("family" %in% names(meta))
  # Non-existent field returns NULL

  expect_null(meta$nonexistent_field)

  # Test that meta is returned even when it only has family
  obj2 <- make_test_implicit_latent(family = "minimal", extra_meta = list())
  meta2 <- implicit_meta(obj2)
  expect_equal(meta2$family, "minimal")
})
