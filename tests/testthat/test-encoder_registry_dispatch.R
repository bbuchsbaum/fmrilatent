# Integration tests for the encoder registry / S3 dispatch contract.
#
# The registry is for discovery only; actual dispatch goes through the S3
# generic encode_spec(). External packages must do BOTH (register + define
# encode_spec.spec_*) for encode() to route to their encoder. These tests
# pin the contract so the seam between the two systems can't drift silently.

skip_if_no_neuroim2 <- function() {
  testthat::skip_if_not_installed("neuroim2")
}

# Snapshot + restore the registry so leftover entries never leak across tests.
with_clean_registry <- function(code) {
  reg_env <- environment(register_encoder)$.encoder_registry_env
  saved <- as.list(reg_env, all.names = TRUE)
  on.exit({
    rm(list = ls(reg_env, all.names = TRUE), envir = reg_env)
    list2env(saved, envir = reg_env)
  }, add = TRUE)
  force(code)
}

# -----------------------------------------------------------------------------
# Defining a custom encoder family that follows the documented contract
# -----------------------------------------------------------------------------

# Spec constructor: returns a list with a distinctive class.
spec_test_zeros <- function(k = 3L) {
  structure(
    list(k = as.integer(k)),
    class = c("spec_test_zeros", "list")
  )
}

# S3 method that builds a real LatentNeuroVec from the data.
encode_spec.spec_test_zeros <- function(x, spec, mask, reduction = NULL,
                                        materialize = "handle", label = "", ...) {
  k <- spec$k
  n_time <- nrow(x)
  mask_arr <- if (inherits(mask, "LogicalNeuroVol")) as.array(mask) else mask
  n_vox <- sum(mask_arr)

  basis    <- Matrix::Matrix(matrix(0, n_time, k), sparse = FALSE)
  loadings <- Matrix::Matrix(matrix(0, n_vox,  k), sparse = FALSE)

  space <- if (inherits(mask, "LogicalNeuroVol")) {
    neuroim2::NeuroSpace(c(dim(mask_arr), n_time))
  } else {
    neuroim2::NeuroSpace(c(dim(mask_arr), n_time))
  }
  mask_vol <- if (inherits(mask, "LogicalNeuroVol")) mask else
    neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dim(mask_arr)))

  LatentNeuroVec(
    basis    = basis,
    loadings = loadings,
    space    = space,
    mask     = mask_vol,
    label    = label,
    meta     = list(family = "test_zeros")
  )
}

# Make the S3 method visible to UseMethod() during the test session.
# (In a real package this is done via NAMESPACE; in tests we register manually.)
.S3method_safe <- function(generic, class, fn) {
  registerS3method(generic, class, fn, envir = topenv(parent.frame()))
}
.S3method_safe("encode_spec", "spec_test_zeros", encode_spec.spec_test_zeros)

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------

test_that("registry is discovery-only: register_encoder alone does not enable dispatch", {
  with_clean_registry({
    # Define a spec that has NO encode_spec method.
    spec_unmethoded <- function() {
      structure(list(), class = c("spec_unmethoded_test_only", "list"))
    }
    register_encoder("unmethoded_test_only", spec_unmethoded,
                     description = "Registered but no S3 method", package = "test")

    # It appears in list_encoders().
    expect_true("unmethoded_test_only" %in% list_encoders()$family)

    # But encode() routes via S3 dispatch and falls through to the default.
    mask <- array(TRUE, c(2L, 2L, 2L))
    X    <- matrix(rnorm(8 * 8), nrow = 8L, ncol = 8L)
    expect_error(
      encode(X, spec_unmethoded(), mask = mask),
      "Unknown spec class"
    )
  })
})

test_that("AWPT is an explicit transport API exception, not a registry family", {
  with_clean_registry({
    expect_false("awpt" %in% list_encoders()$family)
    mask <- array(TRUE, c(2L, 2L, 1L))
    X <- matrix(rnorm(4 * 4), nrow = 4L, ncol = 4L)

    expect_error(
      encode(X, basis_awpt_wavelet(), mask = mask),
      "requires a shared basis_asset and subject field_operator",
      class = "fmrilatent_error_unsupported_spec"
    )
  })
})

test_that("the documented contract works: register + S3 method => encode() dispatches", {
  skip_if_no_neuroim2()

  with_clean_registry({
    register_encoder(
      "test_zeros", spec_test_zeros,
      description = "Test encoder returning all-zero LatentNeuroVec",
      package     = "fmrilatent.tests"
    )

    expect_true("test_zeros" %in% list_encoders()$family)
    expect_identical(get_encoder("test_zeros")$spec_fn, spec_test_zeros)

    set.seed(11)
    n_time <- 6L
    mask   <- array(TRUE, dim = c(2L, 2L, 2L))
    n_vox  <- sum(mask)
    X      <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

    out <- encode(X, spec_test_zeros(k = 4L), mask = mask, label = "registry-roundtrip")
    expect_s4_class(out, "LatentNeuroVec")
    expect_identical(out@label, "registry-roundtrip")
    expect_identical(out@meta$family, "test_zeros")
    expect_equal(.latent_basis_dim(out@basis),    c(n_time, 4L))
    expect_equal(.latent_loadings_dim(out@loadings), c(n_vox, 4L))
  })
})

test_that("get_encoder() reports the available list when the family is missing", {
  with_clean_registry({
    register_encoder("alpha_known", identity, "alpha", "test")
    register_encoder("beta_known",  identity, "beta",  "test")

    expect_error(
      get_encoder("nonexistent_family"),
      "Available encoders:.*alpha_known.*beta_known"
    )
  })
})
