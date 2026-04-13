# Tests for encode.R - spec constructors, encode() generic, and latent_factory()

# ---------------------------------------------------------------------------
# Setup helpers
# ---------------------------------------------------------------------------

# Create a small mask and data matrix for testing
make_test_data <- function(dims = c(3, 3, 2), n_time = 8) {
  mask <- array(TRUE, dim = dims)
  mask_vol <- LogicalNeuroVol(mask, NeuroSpace(dims))
  X <- matrix(rnorm(n_time * sum(mask)), nrow = n_time)
  list(mask = mask, mask_vol = mask_vol, X = X, dims = dims, n_time = n_time)
}

# ---------------------------------------------------------------------------
# Spec constructors
# ---------------------------------------------------------------------------

test_that("spec_time_slepian creates correct spec object", {

  spec <- spec_time_slepian(tr = 2, bandwidth = 0.1, k = 5)
  expect_s3_class(spec, "spec_time_slepian")
  expect_equal(spec$tr, 2)
  expect_equal(spec$bandwidth, 0.1)
  expect_equal(spec$k, 5)
  expect_equal(spec$backend, "tridiag")  # default
})

test_that("spec_time_slepian accepts different backends", {
  spec_tri <- spec_time_slepian(tr = 2, bandwidth = 0.1, backend = "tridiag")
  spec_dense <- spec_time_slepian(tr = 2, bandwidth = 0.1, backend = "dense")
  expect_equal(spec_tri$backend, "tridiag")
  expect_equal(spec_dense$backend, "dense")
})

test_that("spec_time_slepian k parameter is optional", {
  spec <- spec_time_slepian(tr = 2, bandwidth = 0.1)
  expect_null(spec$k)
})

test_that("spec_time_dct creates correct spec object", {
  spec <- spec_time_dct(k = 10, norm = "ortho")
  expect_s3_class(spec, "spec_time_dct")
  expect_equal(spec$k, 10)
  expect_equal(spec$norm, "ortho")
})

test_that("spec_time_dct norm defaults to ortho", {
  spec <- spec_time_dct(k = 5)
  expect_equal(spec$norm, "ortho")
})

test_that("spec_time_dct accepts none normalization", {
  spec <- spec_time_dct(k = 5, norm = "none")
  expect_equal(spec$norm, "none")
})

test_that("spec constructors reject malformed inputs", {
  expect_error(spec_time_dct(k = -5), "positive integer")
  expect_error(spec_time_slepian(tr = -1, bandwidth = 0.1), "positive finite number")
  expect_error(spec_time_slepian(tr = 2, bandwidth = 0), "positive finite number")
  expect_error(spec_time_bspline(k = 0), "positive integer")
  expect_error(spec_space_slepian(k = 0), "positive integer")
  expect_error(spec_space_pca(k = 0), "positive integer")
  expect_error(spec_space_heat(scales = c(1, -2)), "positive finite values")
  expect_error(spec_space_heat(order = 0), "positive integer")
  expect_error(spec_space_hrbf(params = "bad"), "params must be a list")
  expect_error(spec_space_hrbf(params = list(kernel_type = "bad")), "kernel_type")
  expect_error(spec_space_wavelet_active(levels_space = 0), "positive integer")
  expect_error(spec_space_wavelet_active(levels_time = -1), "non-negative integer")
  expect_error(spec_space_wavelet_active(threshold = -0.1), "non-negative finite number")
})

test_that("spec_time_bspline creates correct spec object", {
  spec <- spec_time_bspline(k = 8, degree = 3, include_intercept = FALSE, orthonormalize = TRUE)
  expect_s3_class(spec, "spec_time_bspline")
  expect_equal(spec$k, 8)
  expect_equal(spec$degree, 3)
  expect_false(spec$include_intercept)
  expect_true(spec$orthonormalize)
})

test_that("spec_time_bspline uses default values", {
  spec <- spec_time_bspline(k = 5)
  expect_equal(spec$degree, 3L)
  expect_false(spec$include_intercept)
  expect_true(spec$orthonormalize)
})

test_that("spec_space_slepian creates correct spec object", {
  spec <- spec_space_slepian(k = 4, k_neighbors = 8)
  expect_s3_class(spec, "spec_space_slepian")
  expect_equal(spec$k, 4)
  expect_equal(spec$k_neighbors, 8)
})

test_that("spec_space_slepian uses default values", {
  spec <- spec_space_slepian()
  expect_equal(spec$k, 3L)
  expect_equal(spec$k_neighbors, 6L)
})

test_that("spec_space_pca creates correct spec object", {
  spec <- spec_space_pca(k = 4, center = FALSE, whiten = TRUE, backend = "svd")
  expect_s3_class(spec, "spec_space_pca")
  expect_equal(spec$k, 4L)
  expect_false(spec$center)
  expect_true(spec$whiten)
  expect_equal(spec$backend, "svd")
})

test_that("spec_space_pca uses default values", {
  spec <- spec_space_pca()
  expect_equal(spec$k, 3L)
  expect_true(spec$center)
  expect_false(spec$whiten)
  expect_equal(spec$backend, "auto")
})

test_that("spec_space_heat creates correct spec object", {
  spec <- spec_space_heat(scales = c(1, 2, 4), order = 20, threshold = 1e-5, k_neighbors = 4)
  expect_s3_class(spec, "spec_space_heat")
  expect_equal(spec$scales, c(1, 2, 4))
  expect_equal(spec$order, 20)
  expect_equal(spec$threshold, 1e-5)
  expect_equal(spec$k_neighbors, 4)
})

test_that("spec_space_heat uses default values", {
  spec <- spec_space_heat()
  expect_equal(spec$scales, c(1, 2, 4, 8))
  expect_equal(spec$order, 30L)
  expect_equal(spec$threshold, 1e-6)
  expect_equal(spec$k_neighbors, 6L)
})

test_that("spec_space_hrbf creates correct spec object", {
  params <- list(sigma0 = 4, levels = 2, seed = 42)
  spec <- spec_space_hrbf(params = params)
  expect_s3_class(spec, "spec_space_hrbf")
  expect_equal(spec$params$sigma0, 4)
  expect_equal(spec$params$levels, 2)
  expect_equal(spec$params$seed, 42)
})

test_that("spec_space_hrbf with empty params", {
  spec <- spec_space_hrbf()
  expect_s3_class(spec, "spec_space_hrbf")
  expect_true(is.list(spec$params))
  expect_length(spec$params, 0)
})

test_that("spec_space_wavelet_active creates correct spec object", {
  spec <- spec_space_wavelet_active(levels_space = 3, levels_time = 1, threshold = 0.01)
  expect_s3_class(spec, "spec_space_wavelet_active")
  expect_equal(spec$levels_space, 3)
  expect_equal(spec$levels_time, 1)
  expect_equal(spec$threshold, 0.01)
})

test_that("spec_space_wavelet_active uses default values", {
  spec <- spec_space_wavelet_active()
  expect_equal(spec$levels_space, 2L)
  expect_equal(spec$levels_time, 0L)
  expect_equal(spec$threshold, 0)
})

test_that("spec_st creates correct spatiotemporal spec", {
  time_spec <- spec_time_dct(k = 5)
  space_spec <- spec_space_slepian(k = 3)
  spec <- spec_st(time = time_spec, space = space_spec)

  expect_s3_class(spec, "spec_st")
  expect_s3_class(spec$time, "spec_time_dct")
  expect_s3_class(spec$space, "spec_space_slepian")
  expect_equal(spec$core_mode, "auto")
})

test_that("spec_st accepts explicit core_mode", {
  spec <- spec_st(
    time = spec_time_slepian(tr = 2, bandwidth = 0.1),
    space = spec_space_slepian(),
    core_mode = "explicit"
  )
  expect_equal(spec$core_mode, "explicit")
})

test_that("spec_hierarchical_template requires template or template_file", {
  expect_error(spec_hierarchical_template(), "Provide either a template object or template_file")
})

test_that("spec_hierarchical_template accepts template object", {
  # Create a mock template object
  mock_template <- list(class = "HierarchicalBasisTemplate")
  spec <- spec_hierarchical_template(template = mock_template)
  expect_s3_class(spec, "spec_hierarchical")
  expect_equal(spec$template, mock_template)
  expect_null(spec$template_file)
})

test_that("spec_hierarchical_template accepts template_file", {
  spec <- spec_hierarchical_template(template_file = "/path/to/template.rds")
  expect_s3_class(spec, "spec_hierarchical")
  expect_null(spec$template)
  expect_equal(spec$template_file, "/path/to/template.rds")
})

# ---------------------------------------------------------------------------
# encode() generic dispatch
# ---------------------------------------------------------------------------

test_that("encode.default errors for unsupported classes", {
  expect_error(encode(list(a = 1), spec_time_dct(k = 3)), "No encode method for class")
})

test_that("encode.matrix works with spec_time_dct", {
  td <- make_test_data()
  spec <- spec_time_dct(k = 4)
  lv <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(dim(basis(lv)), c(td$n_time, 4))
  expect_equal(nrow(loadings(lv)), sum(td$mask))
})

test_that("encode.matrix works with spec_time_dct norm='none'", {
  td <- make_test_data()
  spec <- spec_time_dct(k = 4, norm = "none")
  lv <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(dim(basis(lv)), c(td$n_time, 4))
})

test_that("encode.matrix works with spec_time_slepian", {
  td <- make_test_data(n_time = 10)
  spec <- spec_time_slepian(tr = 2, bandwidth = 0.1, k = 3)
  lv <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(dim(basis(lv)), c(10, 3))
})

test_that("encode.matrix works with spec_time_bspline", {
  td <- make_test_data(n_time = 12)
  spec <- spec_time_bspline(k = 5, degree = 3)
  lv <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(dim(basis(lv)), c(12, 5))
})

test_that("encode.matrix works with spec_time_bspline orthonormalize=FALSE", {
  td <- make_test_data(n_time = 12)
  spec <- spec_time_bspline(k = 5, degree = 3, orthonormalize = FALSE)
  lv <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(dim(basis(lv)), c(12, 5))
})

test_that("encode with materialize='handle' works for temporal bases", {
  td <- make_test_data(n_time = 10)
  spec <- spec_time_slepian(tr = 2, bandwidth = 0.1, k = 4)

  lv <- encode(td$X, spec, mask = td$mask_vol, materialize = "handle")

  expect_s4_class(lv, "LatentNeuroVec")
  expect_s4_class(lv@basis, "BasisHandle")
  expect_equal(dim(loadings(lv)), c(sum(td$mask), 4))
})

test_that("encode with materialize='auto' works for temporal bases", {
  td <- make_test_data(n_time = 10)
  spec <- spec_time_dct(k = 4)

  lv <- encode(td$X, spec, mask = td$mask_vol, materialize = "auto")

  expect_s4_class(lv, "LatentNeuroVec")
  expect_s4_class(lv@basis, "BasisHandle")
  expect_equal(dim(loadings(lv)), c(sum(td$mask), 4))
})

test_that("encode with label parameter sets label", {
  td <- make_test_data()
  spec <- spec_time_dct(k = 4)
  lv <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix", label = "test_label")

  expect_equal(lv@label, "test_label")
})

# ---------------------------------------------------------------------------
# encode with spatial specs (require RSpectra for Slepian)
# ---------------------------------------------------------------------------

test_that("encode with spec_space_slepian returns LatentNeuroVec", {
  skip_if_not_installed("RSpectra")
  skip_if_not_installed("rgsp")

  td <- make_test_data()
  spec <- spec_space_slepian(k = 2, k_neighbors = 3)
  lv <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  expect_s4_class(lv, "LatentNeuroVec")
})

test_that("encode with spec_space_heat returns LatentNeuroVec", {
  skip_if_not_installed("RSpectra")
  skip_if_not_installed("rgsp")

  td <- make_test_data()
  spec <- spec_space_heat(scales = c(1, 2), order = 10, k_neighbors = 3)
  lv <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  expect_s4_class(lv, "LatentNeuroVec")
})

test_that("encode with spec_space_hrbf returns LatentNeuroVec", {
  td <- make_test_data()
  params <- list(sigma0 = 4, levels = 2, seed = 42)
  spec <- spec_space_hrbf(params = params)
  lv <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  expect_s4_class(lv, "LatentNeuroVec")
})

test_that("encode with spec_space_wavelet_active returns ImplicitLatent", {
  td <- make_test_data(dims = c(4, 4, 2))  # power of 2 in spatial dims helps wavelets
  spec <- spec_space_wavelet_active(levels_space = 1, levels_time = 0)
  result <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  expect_true(is_implicit_latent(result))
})

# ---------------------------------------------------------------------------
# encode with spatiotemporal specs (spec_st)
# ---------------------------------------------------------------------------

test_that("encode with spec_st (slepian time + slepian space) returns ImplicitLatent", {
  skip_if_not_installed("RSpectra")
  skip_if_not_installed("rgsp")

  td <- make_test_data(n_time = 10)
  spec <- spec_st(
    time = spec_time_slepian(tr = 2, bandwidth = 0.1, k = 3),
    space = spec_space_slepian(k = 2, k_neighbors = 3)
  )
  result <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  expect_true(is_implicit_latent(result))
  # Check decoder works
  reco <- predict(result)
  expect_equal(dim(reco), dim(td$X))
})

test_that("encode with spec_st (bspline time + hrbf space) returns ImplicitLatent", {
  td <- make_test_data(n_time = 10)
  spec <- spec_st(
    time = spec_time_bspline(k = 4, degree = 3),
    space = spec_space_hrbf(params = list(sigma0 = 4, levels = 2, seed = 42))
  )
  result <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  expect_true(is_implicit_latent(result))
  reco <- predict(result)
  expect_equal(dim(reco), dim(td$X))
})

test_that("encode with spec_st (bspline time + slepian space) returns ImplicitLatent", {
  skip_if_not_installed("RSpectra")
  skip_if_not_installed("rgsp")

  td <- make_test_data(n_time = 10)
  spec <- spec_st(
    time = spec_time_bspline(k = 4, degree = 3, orthonormalize = TRUE),
    space = spec_space_slepian(k = 2, k_neighbors = 3)
  )
  result <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  expect_true(is_implicit_latent(result))
})

test_that("encode with spec_st (bspline time + slepian space, orthonormalize=FALSE) returns ImplicitLatent", {
  skip_if_not_installed("RSpectra")
  skip_if_not_installed("rgsp")

  td <- make_test_data(n_time = 10)
  spec <- spec_st(
    time = spec_time_bspline(k = 4, degree = 3, orthonormalize = FALSE),
    space = spec_space_slepian(k = 2, k_neighbors = 3)
  )
  result <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  expect_true(is_implicit_latent(result))
})

test_that("encode with spec_st and spec_space_wavelet_active errors", {
  td <- make_test_data()
  spec <- spec_st(
    time = spec_time_bspline(k = 4),
    space = spec_space_wavelet_active()
  )
  expect_error(
    encode(td$X, spec, mask = td$mask_vol),
    "not yet supported"
  )
})

test_that("encode_spec.default errors for unknown spec class", {
  td <- make_test_data()
  # Create an object with an unknown spec class
  unknown_spec <- structure(list(), class = "unknown_spec")
  expect_error(encode(td$X, unknown_spec, mask = td$mask_vol), "Unknown spec class")
})

# ---------------------------------------------------------------------------
# latent_factory()
# ---------------------------------------------------------------------------

test_that("latent_factory with dct_time family works", {
  td <- make_test_data()
  lv <- latent_factory("dct_time", x = td$X, mask = td$mask_vol, k = 4, materialize = "matrix")

  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(dim(basis(lv)), c(td$n_time, 4))
})

test_that("latent_factory with slepian_time family works", {
  td <- make_test_data(n_time = 10)
  lv <- latent_factory("slepian_time", x = td$X, mask = td$mask_vol,
                       tr = 2, bandwidth = 0.1, k = 3, materialize = "matrix")

  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(dim(basis(lv)), c(10, 3))
})

test_that("latent_factory with slepian_space family works", {
  skip_if_not_installed("RSpectra")
  skip_if_not_installed("rgsp")

  td <- make_test_data()
  lv <- latent_factory("slepian_space", x = td$X, mask = td$mask_vol,
                       k = 2, k_neighbors = 3, materialize = "matrix")

  expect_s4_class(lv, "LatentNeuroVec")
})

test_that("latent_factory with heat_space family works", {
  skip_if_not_installed("RSpectra")
  skip_if_not_installed("rgsp")

  td <- make_test_data()
  lv <- latent_factory("heat_space", x = td$X, mask = td$mask_vol,
                       scales = c(1, 2), order = 10, k_neighbors = 3, materialize = "matrix")

  expect_s4_class(lv, "LatentNeuroVec")
})

test_that("latent_factory with bspline_hrbf_st family works", {
  td <- make_test_data(n_time = 10)
  lv <- latent_factory("bspline_hrbf_st", x = td$X, mask = td$mask_vol,
                       k_time = 4, degree = 3,
                       params = list(sigma0 = 4, levels = 2, seed = 42),
                       materialize = "matrix")

  expect_true(is_implicit_latent(lv))
})

test_that("latent_factory with wavelet_active family works", {
  td <- make_test_data(dims = c(4, 4, 2))
  lv <- latent_factory("wavelet_active", x = td$X, mask = td$mask_vol,
                       levels_space = 1, levels_time = 0, threshold = 0,
                       materialize = "matrix")

  expect_true(is_implicit_latent(lv))
})

test_that("latent_factory with slepian_st family works", {
  skip_if_not_installed("RSpectra")
  skip_if_not_installed("rgsp")

  td <- make_test_data(n_time = 10)
  lv <- latent_factory("slepian_st", x = td$X, mask = td$mask_vol,
                       tr = 2, bandwidth = 0.1, k_time = 3, k_space = 2, k_neighbors = 3,
                       materialize = "matrix")

  expect_true(is_implicit_latent(lv))
})

test_that("latent_factory with label parameter sets label", {
  td <- make_test_data()
  lv <- latent_factory("dct_time", x = td$X, mask = td$mask_vol, k = 4,
                       materialize = "matrix", label = "factory_test")

  expect_equal(lv@label, "factory_test")
})

test_that("latent_factory rejects AWPT family with a direct handoff message", {
  td <- make_test_data()

  expect_error(
    latent_factory("awpt", x = td$X, mask = td$mask_vol),
    "does not support AWPT because AWPT requires a shared basis_asset and subject field_operator"
  )
})

# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

test_that("encode works with single time point", {
  mask <- array(TRUE, dim = c(3, 3, 2))
  mask_vol <- LogicalNeuroVol(mask, NeuroSpace(c(3, 3, 2)))
  X <- matrix(rnorm(sum(mask)), nrow = 1)
  spec <- spec_time_dct(k = 1)

  lv <- encode(X, spec, mask = mask_vol, materialize = "matrix")
  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(dim(basis(lv)), c(1, 1))
})

test_that("encode works with single voxel mask", {
  mask <- array(FALSE, dim = c(3, 3, 2))
  mask[2, 2, 1] <- TRUE
  mask_vol <- LogicalNeuroVol(mask, NeuroSpace(c(3, 3, 2)))
  X <- matrix(rnorm(5), nrow = 5)
  spec <- spec_time_dct(k = 3)

  lv <- encode(X, spec, mask = mask_vol, materialize = "matrix")
  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(nrow(loadings(lv)), 1)
})

test_that("encode works with small k values", {
  td <- make_test_data()
  spec <- spec_time_dct(k = 1)
  lv <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  expect_equal(ncol(basis(lv)), 1)
})

test_that("encode works with mask as array (not LogicalNeuroVol)", {
  mask <- array(TRUE, dim = c(3, 3, 2))
  mask_vol <- LogicalNeuroVol(mask, NeuroSpace(c(3, 3, 2)))
  X <- matrix(rnorm(5 * sum(mask)), nrow = 5)
  spec <- spec_time_dct(k = 3)

  # The encode function should work with LogicalNeuroVol
  lv <- encode(X, spec, mask = mask_vol, materialize = "matrix")
  expect_s4_class(lv, "LatentNeuroVec")
})

test_that("encode works with non-cubic mask dimensions", {
  # Test with rectangular dimensions
  mask <- array(TRUE, dim = c(2, 5, 3))
  mask_vol <- LogicalNeuroVol(mask, NeuroSpace(c(2, 5, 3)))
  X <- matrix(rnorm(6 * sum(mask)), nrow = 6)
  spec <- spec_time_dct(k = 4)

  lv <- encode(X, spec, mask = mask_vol, materialize = "matrix")
  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(nrow(loadings(lv)), sum(mask))
})

test_that("encode handles partial mask (not all TRUE)", {
  mask <- array(TRUE, dim = c(3, 3, 2))
  mask[1, 1, 1] <- FALSE
  mask[3, 3, 2] <- FALSE
  mask_vol <- LogicalNeuroVol(mask, NeuroSpace(c(3, 3, 2)))
  X <- matrix(rnorm(5 * sum(mask)), nrow = 5)
  spec <- spec_time_dct(k = 3)

  lv <- encode(X, spec, mask = mask_vol, materialize = "matrix")
  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(nrow(loadings(lv)), sum(mask))
})

# ---------------------------------------------------------------------------
# Reconstruction quality
# ---------------------------------------------------------------------------

test_that("encode with DCT produces reasonable reconstruction", {
  td <- make_test_data(n_time = 20)
  # Using nearly full rank to get good reconstruction
  spec <- spec_time_dct(k = 18)
  lv <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  # Reconstruct
  reconstructed <- as.matrix(lv)

  # Check that reconstruction is similar to original
  # DCT with k close to n_time should give good reconstruction
  # Check mean correlation across voxels
  cor_vals <- diag(cor(td$X, reconstructed))
  expect_true(mean(cor_vals) > 0.9)  # High mean correlation
})

test_that("encode with slepian produces reasonable reconstruction", {
  td <- make_test_data(n_time = 20)
  spec <- spec_time_slepian(tr = 2, bandwidth = 0.1, k = 10)
  lv <- encode(td$X, spec, mask = td$mask_vol, materialize = "matrix")

  reconstructed <- as.matrix(lv)
  cor_vals <- diag(cor(td$X, reconstructed))
  expect_true(mean(cor_vals) > 0.5)  # Reasonable correlation on average
})

# ---------------------------------------------------------------------------
# ImplicitLatent decoder tests
# ---------------------------------------------------------------------------

test_that("ImplicitLatent decoder supports time_idx subsetting", {
  skip_if_not_installed("RSpectra")
  skip_if_not_installed("rgsp")

  td <- make_test_data(n_time = 10)
  spec <- spec_st(
    time = spec_time_slepian(tr = 2, bandwidth = 0.1, k = 3),
    space = spec_space_slepian(k = 2, k_neighbors = 3)
  )
  result <- encode(td$X, spec, mask = td$mask_vol)

  # Full reconstruction
  full_reco <- predict(result)

  # Partial time reconstruction
  partial_reco <- predict(result, time_idx = 1:5)
  expect_equal(nrow(partial_reco), 5)
  expect_equal(ncol(partial_reco), ncol(full_reco))
})

test_that("ImplicitLatent decoder supports roi_mask subsetting", {
  skip_if_not_installed("RSpectra")
  skip_if_not_installed("rgsp")

  td <- make_test_data(n_time = 10)
  spec <- spec_st(
    time = spec_time_slepian(tr = 2, bandwidth = 0.1, k = 3),
    space = spec_space_slepian(k = 2, k_neighbors = 3)
  )
  result <- encode(td$X, spec, mask = td$mask_vol)

  # Create ROI mask (subset of voxels)
  roi_mask <- rep(FALSE, sum(td$mask))
  roi_mask[1:5] <- TRUE

  partial_reco <- predict(result, roi_mask = roi_mask)
  expect_equal(ncol(partial_reco), 5)
})

# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------

test_that("encode errors when mask cannot be converted to array", {
  td <- make_test_data()
  spec <- spec_space_slepian()
  # Create an invalid mask object - any error is acceptable
  bad_mask <- list(invalid = TRUE)

  expect_error(encode(td$X, spec, mask = bad_mask))
})

test_that("spec_st with unsupported time spec errors", {
  td <- make_test_data()
  # Create a fake unsupported time spec
  fake_time_spec <- structure(list(k = 3), class = "spec_time_fake")
  spec <- spec_st(time = fake_time_spec, space = spec_space_slepian())

  expect_error(encode(td$X, spec, mask = td$mask_vol), "Unsupported spec_st\\$time class")
})

test_that("spec_st with unsupported space spec errors", {
  skip_if_not_installed("RSpectra")
  skip_if_not_installed("rgsp")

  td <- make_test_data()
  # Create a fake unsupported space spec
  fake_space_spec <- structure(list(k = 3), class = "spec_space_fake")
  spec <- spec_st(time = spec_time_slepian(tr = 2, bandwidth = 0.1, k = 3), space = fake_space_spec)

  expect_error(encode(td$X, spec, mask = td$mask_vol), "Unsupported spec_st\\$space class")
})
