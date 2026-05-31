test_that("wavelet_active encoder returns a known latent type", {
  set.seed(1)
  dims <- c(8, 8, 4)
  mask_arr <- array(TRUE, dims)
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dims))
  n_time <- 16
  n_vox <- prod(dims)
  X <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  spec <- spec_space_wavelet_active(levels_space = 1L, levels_time = 0L, threshold = 0)
  lat <- encode(X, spec, mask = mask)
  # wavelet_active_latent() may return either an explicit LatentNeuroVec or an
  # ImplicitLatent depending on the transform; both are accepted at the seam.
  expect_true(is(lat, "LatentNeuroVec") || inherits(lat, "ImplicitLatent"))
})

test_that("wavelet_active seam aborts if the inner latent has the wrong class", {
  # Stub wavelet_active_latent() in the package namespace so the seam sees a
  # non-latent return and must trip the class assertion.
  ns <- asNamespace("fmrilatent")
  orig <- get("wavelet_active_latent", envir = ns)
  on.exit({
    unlockBinding("wavelet_active_latent", ns)
    assign("wavelet_active_latent", orig, envir = ns)
    lockBinding("wavelet_active_latent", ns)
  }, add = TRUE)
  unlockBinding("wavelet_active_latent", ns)
  assign("wavelet_active_latent", function(...) list(not = "a latent"), envir = ns)
  lockBinding("wavelet_active_latent", ns)

  dims <- c(4, 4, 2)
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dims))
  X <- matrix(rnorm(8 * prod(dims)), 8, prod(dims))
  spec <- spec_space_wavelet_active(levels_space = 1L)

  expect_error(
    encode(X, spec, mask = mask),
    class = "fmrilatent_error_invalid_type"
  )
})

test_that("diffusion-global latent emits no dense warning (expect_dense)", {
  set.seed(1)
  dims <- c(6, 6, 4)
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dims))
  n_time <- 20
  n_vox <- prod(dims)
  X <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  red <- make_cluster_reduction(mask, seq_len(n_vox))

  expect_no_warning(
    expect_no_message(
      diffusion_wavelet_latent(X, mask, reduction = red)
    )
  )
})

test_that("genuinely-unexpected dense loadings still emit the dense warning", {
  dims <- c(3, 3, 2)
  n_vox <- prod(dims)
  k <- 2L
  n_time <- 5L
  space <- neuroim2::NeuroSpace(c(dims, n_time))
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dims))
  basis <- Matrix::Matrix(matrix(rnorm(n_time * k), n_time, k), sparse = FALSE)
  loadings <- matrix(runif(n_vox * k) + 1, n_vox, k)  # fully dense base matrix

  expect_warning(
    LatentNeuroVec(basis = basis, loadings = loadings, space = space, mask = mask),
    "dense",
    class = "fmrilatent_warning_dense_storage"
  )
  expect_no_warning(
    expect_no_message(
      LatentNeuroVec(basis = basis, loadings = loadings, space = space, mask = mask,
                     expect_dense = TRUE)
    )
  )
})
