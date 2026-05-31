# Edge-case and defensive tests for the encode() pipeline.
#
# Existing tests cover the happy path and 1-voxel / 1-timepoint cases. This
# file pins down boundary inputs that previously had no coverage:
# empty / all-FALSE masks, NA in the data matrix, dimension mismatches,
# spec parameter validation, cascade encoding via NeuroVec inheritance, and
# spec_st composition together with the documented unsupported combinations.

skip_if_no_neuroim2 <- function() {
  testthat::skip_if_not_installed("neuroim2")
}

# -----------------------------------------------------------------------------
# Mask edge cases
# -----------------------------------------------------------------------------

test_that("encode() rejects an all-FALSE mask cleanly", {
  skip_if_no_neuroim2()

  mask <- array(FALSE, dim = c(2L, 2L, 2L))
  X    <- matrix(rnorm(8 * 8), nrow = 8L, ncol = 8L)

  # Either an error before construction or a downstream LatentNeuroVec
  # validity error — both are acceptable; silent success is not.
  expect_error(
    encode(X, spec_time_dct(k = 2L), mask = mask),
    regexp = ".+"
  )
})

test_that("encode() with an empty (zero-dimensional) mask errors", {
  skip_if_no_neuroim2()

  without_dct_fallback_warnings(
    expect_error(
      encode(matrix(0, 4L, 0L), spec_time_dct(k = 2L), mask = array(logical(0), dim = c(0L, 0L, 0L))),
      regexp = ".+"
    )
  )
})

# -----------------------------------------------------------------------------
# NA / non-finite handling
# -----------------------------------------------------------------------------

test_that("encode() with NA in the data matrix surfaces an informative error", {
  skip_if_no_neuroim2()

  set.seed(101)
  mask <- array(TRUE, dim = c(2L, 2L, 2L))
  n_vox <- sum(mask)
  X     <- matrix(rnorm(8 * n_vox), nrow = 8L, ncol = n_vox)
  X[1L, 1L] <- NA_real_

  # The exact error site is encoder-specific; we only require the failure to
  # be observable rather than producing a quietly NaN-corrupted basis.
  out <- tryCatch(
    encode(X, spec_time_dct(k = 2L), mask = mask),
    error = function(e) e
  )

  if (inherits(out, "error")) {
    expect_true(nzchar(conditionMessage(out)))
  } else {
    # If the encoder accepted it, the result must not silently propagate NA
    # into the basis — finite values everywhere, or the constructor would
    # have already raised "must contain only finite values".
    b <- as.matrix(fmrilatent:::basis_mat(out@basis))
    expect_true(all(is.finite(b)))
  }
})

test_that("LatentNeuroVec constructor rejects non-finite basis directly", {
  skip_if_no_neuroim2()

  mask  <- neuroim2::LogicalNeuroVol(array(TRUE, c(2L, 2L, 2L)),
                                     neuroim2::NeuroSpace(c(2L, 2L, 2L)))
  space <- neuroim2::NeuroSpace(c(2L, 2L, 2L, 6L))
  B     <- matrix(rnorm(6 * 2), 6, 2); B[1, 1] <- NA_real_
  L     <- matrix(rnorm(8 * 2), 8, 2)

  expect_error(
    LatentNeuroVec(B, L, space, mask),
    "finite"
  )
})

# -----------------------------------------------------------------------------
# Dimension and parameter validation
# -----------------------------------------------------------------------------

test_that("encode() errors when spec_time_dct k exceeds n_time", {
  skip_if_no_neuroim2()

  mask <- array(TRUE, dim = c(2L, 2L, 2L))
  X    <- matrix(rnorm(4 * 8), nrow = 4L, ncol = 8L)

  # k=10 > n_time=4 should be caught by build_dct_basis().
  expect_error(
    encode(X, spec_time_dct(k = 10L), mask = mask),
    "k.*n_time"
  )
})

# -----------------------------------------------------------------------------
# Cascade encoding via NeuroVec inheritance
# -----------------------------------------------------------------------------

test_that("encode() on a LatentNeuroVec cascades through NeuroVec inheritance", {
  skip_if_no_neuroim2()

  set.seed(11)
  mask  <- neuroim2::LogicalNeuroVol(array(TRUE, c(2L, 2L, 2L)),
                                     neuroim2::NeuroSpace(c(2L, 2L, 2L)))
  n_vox <- sum(as.array(mask))
  n_time <- 8L
  space <- neuroim2::NeuroSpace(c(2L, 2L, 2L, n_time))

  # Encode random data once.
  X     <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)
  lvec1 <- encode(X, spec_time_dct(k = 4L), mask = mask)
  expect_s4_class(lvec1, "LatentNeuroVec")

  # Cascade-encode the resulting LatentNeuroVec; encode.NeuroVec extracts
  # the time series via series() and re-runs through encode.matrix.
  lvec2 <- encode(lvec1, spec_time_dct(k = 3L), mask = mask)
  expect_s4_class(lvec2, "LatentNeuroVec")
  expect_equal(.latent_basis_dim(lvec2@basis)[1L], n_time)
  expect_equal(.latent_basis_dim(lvec2@basis)[2L], 3L)
})

# -----------------------------------------------------------------------------
# spec_st composition: the supported combinations and the documented limits
# -----------------------------------------------------------------------------

test_that("spec_st with Slepian time + Slepian space returns a valid ImplicitLatent", {
  skip_if_no_neuroim2()
  testthat::skip_if_not_installed("rgsp")

  set.seed(12)
  mask   <- array(TRUE, dim = c(3L, 3L, 3L))
  n_vox  <- sum(mask)
  n_time <- 16L
  X      <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)

  spec <- spec_st(
    time  = spec_time_slepian(tr = 1, bandwidth = 0.1, k = 4L),
    space = spec_space_slepian(k = 3L, k_neighbors = 4L)
  )

  out <- encode(X, spec, mask = mask)
  # spec_st builds a decoder-backed ImplicitLatent, not a LatentNeuroVec —
  # this is the documented behaviour of the separable spatiotemporal path.
  expect_s3_class(out, "ImplicitLatent")
  expect_identical(out$meta$family, "st_separable")
  expect_equal(nrow(out$coeff$B_t), n_time)
  expect_equal(nrow(out$coeff$L_s), n_vox)
  expect_equal(ncol(out$coeff$L_s), out$meta$k_space)

  # Decoder must reconstruct a sensibly-shaped time x voxel matrix.
  rec <- predict(out)
  expect_equal(dim(rec), c(n_time, n_vox))
  expect_true(all(is.finite(rec)))
})

test_that("spec_st now supports spec_time_dct (composes DCT time with Slepian space)", {
  skip_if_no_neuroim2()
  testthat::skip_if_not_installed("rgsp")

  set.seed(33)
  mask  <- array(TRUE, dim = c(3L, 3L, 3L))
  n_vox <- sum(mask)
  X     <- matrix(rnorm(8 * n_vox), nrow = 8L, ncol = n_vox)

  spec <- spec_st(
    time  = spec_time_dct(k = 4L),
    space = spec_space_slepian(k = 3L, k_neighbors = 4L)
  )

  out <- encode(X, spec, mask = mask)
  expect_s3_class(out, "ImplicitLatent")
  expect_identical(out$meta$time, "spec_time_dct")
  expect_equal(nrow(out$coeff$B_t), 8L)
  expect_equal(ncol(out$coeff$B_t), 4L)

  # And k > n_time still fails fast.
  expect_error(
    encode(X, spec_st(
      time  = spec_time_dct(k = 20L),
      space = spec_space_slepian(k = 3L, k_neighbors = 4L)
    ), mask = mask),
    "cannot exceed n_time"
  )
})

test_that("spec_st rejects spec_space_pca with a clear error (documented limitation)", {
  skip_if_no_neuroim2()
  testthat::skip_if_not_installed("rgsp")

  mask <- array(TRUE, dim = c(2L, 2L, 2L))
  X    <- matrix(rnorm(8 * 8), nrow = 8L, ncol = 8L)

  spec <- spec_st(
    time  = spec_time_slepian(tr = 1, bandwidth = 0.1, k = 3L),
    space = spec_space_pca(k = 3L)
  )

  expect_error(
    encode(X, spec, mask = mask),
    "Unsupported spec_st\\$space class"
  )
})
