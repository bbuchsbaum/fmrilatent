library(testthat)

# Tests for backward compatibility aliases and helpers for Haar implicit latent
# These functions are defined in R/zzz_haar_aliases.R

# Helper function to create test data
make_test_haar_latent <- function() {

  set.seed(42)
  mask <- array(TRUE, dim = c(2, 2, 2))
  X <- matrix(rnorm(6 * sum(mask)), nrow = 6)
  haar_latent(X, mask, levels = 1, z_seed = 42L)
}

make_test_implicit_latent_haar <- function() {
  # Create an ImplicitLatent with family = "haar" that is NOT already HaarLatent
  set.seed(42)
  mask <- array(TRUE, dim = c(2, 2, 2))
  X <- matrix(rnorm(6 * sum(mask)), nrow = 6)
  fw <- haar_wavelet_forward(X, mask, levels = 1, z_seed = 42L)
  meta <- fw$meta
  meta$family <- "haar"
  decoder <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL) {
    haar_wavelet_inverse(list(coeff = fw$coeff, meta = meta), mask,
                         levels = meta$levels, z_seed = meta$z_seed,
                         roi_mask = roi_mask, time_idx = time_idx,
                         levels_keep = levels_keep)
  }
  # Return raw ImplicitLatent without HaarLatent class
  implicit_latent(fw$coeff, decoder, meta, mask)
}

make_test_implicit_latent_non_haar <- function() {
  # Create an ImplicitLatent with family != "haar"
  meta <- list(family = "other")
  decoder <- function(...) matrix(0, 1, 1)
  implicit_latent(list(), decoder, meta, array(TRUE, c(2, 2, 2)))
}

# =============================================================================
# Tests for is_haar_latent()
# =============================================================================

test_that("is_haar_latent returns TRUE for HaarLatent objects", {
  hl <- make_test_haar_latent()
  expect_true(is_haar_latent(hl))
})

test_that("is_haar_latent returns TRUE for ImplicitLatent with family='haar'", {
  il <- make_test_implicit_latent_haar()
  # Ensure it's not already HaarLatent

  expect_false(inherits(il, "HaarLatent"))
  expect_true(inherits(il, "ImplicitLatent"))
  # Should still be recognized as haar latent

  expect_true(is_haar_latent(il))
})

test_that("is_haar_latent returns FALSE for ImplicitLatent with non-haar family", {
  il <- make_test_implicit_latent_non_haar()
  expect_false(is_haar_latent(il))
})

test_that("is_haar_latent returns FALSE for non-ImplicitLatent objects", {
  expect_false(is_haar_latent(NULL))
  expect_false(is_haar_latent(list()))
  expect_false(is_haar_latent(1:10))
  expect_false(is_haar_latent("not a latent"))
  expect_false(is_haar_latent(matrix(1:4, 2, 2)))
  expect_false(is_haar_latent(data.frame(a = 1:3)))
})

test_that("is_haar_latent returns FALSE for object without meta", {
  obj <- structure(list(coeff = NULL, decoder = function(...) NULL),
                   class = "ImplicitLatent")
  expect_false(is_haar_latent(obj))
})

# =============================================================================
# Tests for haar_meta()
# =============================================================================

test_that("haar_meta extracts meta from HaarLatent objects", {
  hl <- make_test_haar_latent()
  meta <- haar_meta(hl)
  expect_type(meta, "list")
  expect_equal(meta$family, "haar")
  expect_true("levels" %in% names(meta))
  expect_true("z_seed" %in% names(meta))
  expect_true("mask_dims" %in% names(meta))
})

test_that("haar_meta extracts meta from ImplicitLatent with family='haar'", {
  il <- make_test_implicit_latent_haar()
  meta <- haar_meta(il)
  expect_type(meta, "list")
  expect_equal(meta$family, "haar")
})

test_that("haar_meta returns NULL for ImplicitLatent with non-haar family", {
  il <- make_test_implicit_latent_non_haar()
  expect_null(haar_meta(il))
})

test_that("haar_meta returns NULL for non-ImplicitLatent objects", {
  expect_null(haar_meta(NULL))
  expect_null(haar_meta(list()))
  expect_null(haar_meta(1:10))
  expect_null(haar_meta("not a latent"))
  expect_null(haar_meta(matrix(1:4, 2, 2)))
})

test_that("haar_meta returns NULL for ImplicitLatent without meta", {
  obj <- structure(list(coeff = NULL, decoder = function(...) NULL),
                   class = "ImplicitLatent")
  expect_null(haar_meta(obj))
})

test_that("haar_meta returns NULL for HaarLatent-classed object without haar family meta", {
  # Edge case: object has HaarLatent class but meta$family is not "haar"
  obj <- structure(list(meta = list(family = "other")),
                   class = c("HaarLatent", "ImplicitLatent"))
  # The function should return NULL because family != "haar"
  expect_null(haar_meta(obj))
})

# =============================================================================
# Tests for as_haar_latent()
# =============================================================================

test_that("as_haar_latent returns HaarLatent unchanged", {
  hl <- make_test_haar_latent()
  result <- as_haar_latent(hl)
  expect_identical(result, hl)
  expect_true(inherits(result, "HaarLatent"))
})

test_that("as_haar_latent converts ImplicitLatent with family='haar' to HaarLatent", {
  il <- make_test_implicit_latent_haar()
  # Confirm it's not already HaarLatent
  expect_false(inherits(il, "HaarLatent"))

  result <- as_haar_latent(il)
  expect_true(inherits(result, "HaarLatent"))
  expect_true(inherits(result, "ImplicitLatent"))

  # Verify the conversion preserves data
  expect_equal(result$coeff, il$coeff)
  expect_equal(result$meta, il$meta)
})

test_that("as_haar_latent adds HaarLatent class only once", {
  il <- make_test_implicit_latent_haar()
  result <- as_haar_latent(il)
  # Convert again
  result2 <- as_haar_latent(result)
  # Should not duplicate HaarLatent in class
  expect_equal(sum(class(result2) == "HaarLatent"), 1)
})

test_that("as_haar_latent throws error for ImplicitLatent with non-haar family", {
  il <- make_test_implicit_latent_non_haar()
  expect_error(as_haar_latent(il), "not Haar implicit latent")
})

test_that("as_haar_latent throws error for non-ImplicitLatent objects", {
  expect_error(as_haar_latent(NULL), "not Haar implicit latent")
  expect_error(as_haar_latent(list()), "not Haar implicit latent")
  expect_error(as_haar_latent(1:10), "not Haar implicit latent")
  expect_error(as_haar_latent("not a latent"), "not Haar implicit latent")
  expect_error(as_haar_latent(matrix(1:4, 2, 2)), "not Haar implicit latent")
})

test_that("as_haar_latent throws error for ImplicitLatent without meta", {
  obj <- structure(list(coeff = NULL, decoder = function(...) NULL),
                   class = "ImplicitLatent")
  expect_error(as_haar_latent(obj), "not Haar implicit latent")
})

# =============================================================================
# Tests for consistency between alias functions
# =============================================================================

test_that("alias functions are consistent with each other", {
  hl <- make_test_haar_latent()

  # is_haar_latent should agree with non-NULL haar_meta
  expect_equal(is_haar_latent(hl), !is.null(haar_meta(hl)))

  # as_haar_latent should return object where is_haar_latent is TRUE
  il <- make_test_implicit_latent_haar()
  converted <- as_haar_latent(il)
  expect_true(is_haar_latent(converted))
  expect_false(is.null(haar_meta(converted)))
})

test_that("alias functions work correctly with predict", {
  # Verify that converted objects still work with predict
  il <- make_test_implicit_latent_haar()
  converted <- as_haar_latent(il)

  # Both should produce same reconstruction
  pred_il <- predict(il)
  pred_conv <- predict(converted)
  expect_equal(pred_il, pred_conv, tolerance = 1e-10)
})

# =============================================================================
# Tests for type consistency
# =============================================================================

test_that("is_haar_latent returns logical type", {
  hl <- make_test_haar_latent()
  result <- is_haar_latent(hl)
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that("haar_meta returns list or NULL", {
  hl <- make_test_haar_latent()
  meta <- haar_meta(hl)
  expect_true(is.list(meta) || is.null(meta))

  meta_null <- haar_meta(list())
  expect_null(meta_null)
})

test_that("as_haar_latent returns object with HaarLatent class", {
  il <- make_test_implicit_latent_haar()
  result <- as_haar_latent(il)
  expect_s3_class(result, "HaarLatent")
  expect_s3_class(result, "ImplicitLatent")
})

# =============================================================================
# Tests for edge cases
# =============================================================================

test_that("functions handle objects with partial structures", {
  # Object with meta but missing family
  obj1 <- structure(list(meta = list(levels = 1)),
                    class = "ImplicitLatent")
  expect_false(is_haar_latent(obj1))
  expect_null(haar_meta(obj1))
  expect_error(as_haar_latent(obj1))

  # Object with meta$family = NULL
  obj2 <- structure(list(meta = list(family = NULL)),
                    class = "ImplicitLatent")
  expect_false(is_haar_latent(obj2))
  expect_null(haar_meta(obj2))
  expect_error(as_haar_latent(obj2))
})

test_that("functions handle objects with extra classes", {
  il <- make_test_implicit_latent_haar()
  class(il) <- c("ExtraClass", class(il))

  # Should still work

  expect_true(is_haar_latent(il))
  expect_false(is.null(haar_meta(il)))

  converted <- as_haar_latent(il)
  expect_true(inherits(converted, "HaarLatent"))
  expect_true(inherits(converted, "ExtraClass"))
})
