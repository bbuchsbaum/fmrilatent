# Sparse, non-contiguous mask coverage for heat and diffusion wavelet latents,
# analogous to the AWPT and HRBF regressions (bd-01KSQQNPEBC8SE0XK0Q41DXM2B).
#
# These bases are rgsp-backed and do not form a tight frame, so reconstruction
# (basis = X %*% loadings; reco = X %*% L L^T) is an approximation, not an exact
# round-trip. We therefore assert the invariants that DO hold on a scattered
# mask: valid construction, correct loadings/voxel ordering, internal
# reconstruction consistency, and handle materialization.

library(testthat)

# A 4x4x3 grid with 7 scattered, non-contiguous active voxels.
make_scattered_mask <- function() {
  mask_arr <- array(FALSE, dim = c(4, 4, 3))
  mask_arr[cbind(
    c(1L, 3L, 4L, 2L, 4L, 1L, 3L),
    c(1L, 4L, 2L, 3L, 1L, 2L, 4L),
    c(1L, 2L, 3L, 1L, 2L, 3L, 1L)
  )] <- TRUE
  mask_arr
}

# Shared assertions for a *_wavelet_latent constructor on a scattered mask.
expect_sparse_mask_latent <- function(latent_fn, spec, k_neighbors = 4L, seed = 101L) {
  mask_arr <- make_scattered_mask()
  expect_false(all(mask_arr))           # genuinely sparse
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(4, 4, 3)))
  nv <- sum(mask_arr)
  nt <- 5L
  set.seed(seed)
  X <- matrix(rnorm(nt * nv), nrow = nt)

  lv <- latent_fn(X, mask, spec = spec, k_neighbors = k_neighbors)

  # Construction: shape follows the scattered voxel count, not the full grid.
  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(nrow(loadings(lv)), nv)
  expect_equal(nrow(basis(lv)), nt)

  # Reconstruction is finite and correctly shaped (n_time x n_active_voxels).
  reco <- as.matrix(lv)
  expect_equal(dim(reco), c(nt, nv))
  expect_true(all(is.finite(reco)))

  # Internal consistency: the materialized factors reproduce as.matrix() exactly,
  # i.e. the scattered mask did not corrupt the basis/loadings product.
  manual <- as.matrix(as.matrix(basis(lv)) %*% t(as.matrix(loadings(lv))))
  expect_equal(reco, manual, tolerance = 1e-10)

  # Analysis transform identity: basis == X %*% loadings.
  expect_equal(
    as.matrix(basis(lv)),
    as.matrix(X %*% as.matrix(loadings(lv))),
    tolerance = 1e-10
  )

  # Voxel-ordering: series() indexed by full-grid linear voxel index must map to
  # the reconstruction column for that voxel's position among the active voxels.
  active_lin <- which(mask_arr)
  for (k in c(1L, 3L, nv)) {
    expect_equal(
      as.numeric(series(lv, active_lin[k])),
      reco[, k],
      tolerance = 1e-10
    )
  }

  invisible(lv)
}

# Shared assertions for a *_loadings_handle constructor on a scattered mask.
expect_sparse_mask_handle <- function(handle_fn, spec, id, k_neighbors = 4L) {
  mask_arr <- make_scattered_mask()
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(4, 4, 3)))
  nv <- sum(mask_arr)
  red <- make_cluster_reduction(mask, seq_len(nv))

  fmrilatent_registry_clear()
  handle <- handle_fn(red, spec, k_neighbors = k_neighbors, id = id)

  expect_s4_class(handle, "LoadingsHandle")
  expect_equal(handle@dim[1], nv)
  expect_true(fmrilatent:::.latent_has_matrix(handle@id, type = "loadings"))

  mat <- loadings_mat(handle)
  expect_equal(nrow(mat), nv)
  expect_equal(ncol(mat), handle@dim[2])
  expect_true(all(is.finite(as.matrix(mat))))

  fmrilatent_registry_clear()
  invisible(handle)
}

test_that("heat_wavelet_latent constructs and reconstructs on a sparse non-contiguous mask", {
  skip_if_not_installed("rgsp")
  expect_sparse_mask_latent(
    heat_wavelet_latent,
    spec = basis_heat_wavelet(scales = c(1, 2, 4), order = 20, threshold = 0)
  )
})

test_that("diffusion_wavelet_latent constructs and reconstructs on a sparse non-contiguous mask", {
  skip_if_not_installed("rgsp")
  expect_sparse_mask_latent(
    diffusion_wavelet_latent,
    spec = basis_diffusion_wavelet(target_rank = 3L, oversample = 2L,
                                   threshold = 0, max_scales = 2L)
  )
})

test_that("heat_wavelet_loadings_handle materializes on a sparse non-contiguous mask", {
  skip_if_not_installed("rgsp")
  expect_sparse_mask_handle(
    heat_wavelet_loadings_handle,
    spec = basis_heat_wavelet(scales = c(1, 2), order = 16, threshold = 0),
    id = "hw-sparse-mask-test"
  )
})

test_that("diffusion_wavelet_loadings_handle materializes on a sparse non-contiguous mask", {
  skip_if_not_installed("rgsp")
  expect_sparse_mask_handle(
    diffusion_wavelet_loadings_handle,
    spec = basis_diffusion_wavelet(target_rank = 3L, oversample = 2L,
                                   threshold = 0, max_scales = 2L),
    id = "dw-sparse-mask-test"
  )
})
