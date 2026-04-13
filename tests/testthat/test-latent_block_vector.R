library(testthat)
library(Matrix)
library(neuroim2)

.skip_if_no_neurosurf_block <- function() {
  skip_if_not_installed("neurosurf")
}

.make_shared_basis_block <- function() {
  Matrix(matrix(c(
    1, 0,
    0, 1,
    1, 1
  ), nrow = 3, byrow = TRUE), sparse = FALSE)
}

.make_surface_block_block <- function(hemi = c("left", "right")) {
  hemi <- match.arg(hemi)
  .skip_if_no_neurosurf_block()
  geom <- neurosurf::example_surface_geometry()
  geom@hemi <- hemi
  LatentNeuroSurfaceVector(
    basis = .make_shared_basis_block(),
    loadings = Matrix(matrix(c(
      1, 0,
      0, 1,
      1, 1
    ), nrow = 3, byrow = TRUE), sparse = FALSE),
    geometry = geom,
    support = 1:3,
    offset = c(0, 0.5, -0.5),
    meta = list(family = paste0("surface_", hemi))
  )
}

.make_bilateral_block_block <- function() {
  BlockLatentNeuroVector(
    list(
      cortex = BilatLatentNeuroSurfaceVector(
        .make_surface_block_block("left"),
        .make_surface_block_block("right"),
        meta = list(family = "bilat_cortex")
      ),
      subcortex = {
        mask_arr <- array(c(TRUE, FALSE, TRUE, FALSE), dim = c(2, 2, 1))
        mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
        LatentNeuroVec(
          basis = .make_shared_basis_block(),
          loadings = Matrix(matrix(c(
            1, 0,
            0, 1
          ), nrow = 2, byrow = TRUE), sparse = FALSE),
          space = NeuroSpace(c(2, 2, 1, 3)),
          mask = mask_vol,
          offset = c(0.25, -0.25),
          meta = list(family = "volume_subcortex")
        )
      }
    ),
    meta = list(family = "hybrid_block_test")
  )
}

test_that("BlockLatentNeuroVector reconstructs hybrid bilateral-surface plus volume outputs", {
  hybrid <- .make_bilateral_block_block()
  mat <- reconstruct_matrix(hybrid)

  expect_true(is_explicit_latent(hybrid))
  expect_equal(latent_meta(hybrid)$family, "hybrid_block_test")
  expect_named(latent_domain(hybrid), c("cortex", "subcortex"))
  expect_named(latent_support(hybrid), c("cortex", "subcortex"))
  expect_equal(ncol(mat), 8L)
  expect_equal(
    mat,
    cbind(
      reconstruct_matrix(hybrid@blocks$cortex),
      reconstruct_matrix(hybrid@blocks$subcortex)
    ),
    tolerance = 1e-8
  )

  wrapped <- wrap_decoded(hybrid, mat)
  expect_s3_class(wrapped, "BlockDecodedLatent")
  expect_named(wrapped, c("cortex", "subcortex"))
  expect_s4_class(wrapped$cortex, "BilatNeuroSurfaceVector")
  expect_true(is.array(wrapped$subcortex))
})

test_that("BlockLatentNeuroVector supports decoder pushforward and block ROI splitting", {
  hybrid <- .make_bilateral_block_block()
  L <- as.matrix(loadings(hybrid))
  gamma <- c(2, -1)
  Sigma <- matrix(c(2, 0.25,
                    0.25, 1), nrow = 2, byrow = TRUE)

  expect_equal(
    decode_coefficients(hybrid, gamma, space = "native"),
    as.vector(L %*% matrix(gamma, ncol = 1)),
    tolerance = 1e-8
  )
  expect_equal(
    decode_covariance(hybrid, Sigma, space = "native", diag_only = FALSE),
    L %*% Sigma %*% t(L),
    tolerance = 1e-8
  )

  roi <- list(
    cortex = list(left = c(TRUE, FALSE, TRUE), right = c(FALSE, TRUE, FALSE)),
    subcortex = array(c(TRUE, FALSE, FALSE, FALSE), dim = c(2, 2, 1))
  )
  cortex_mat <- reconstruct_matrix(hybrid@blocks$cortex, roi_mask = roi$cortex)
  subcortex_mat <- reconstruct_matrix(hybrid@blocks$subcortex, roi_mask = roi$subcortex)
  expect_equal(
    reconstruct_matrix(hybrid, roi_mask = roi),
    cbind(cortex_mat, subcortex_mat),
    tolerance = 1e-8
  )
})

test_that("BlockLatentNeuroVector rejects blocks with mismatched shared bases", {
  .skip_if_no_neurosurf_block()
  left <- .make_surface_block_block("left")
  right <- .make_surface_block_block("right")
  bad_volume <- {
    mask_arr <- array(c(TRUE, FALSE, TRUE, FALSE), dim = c(2, 2, 1))
    mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
    LatentNeuroVec(
      basis = Matrix(matrix(c(
        1, 0,
        0, 1,
        2, 1
      ), nrow = 3, byrow = TRUE), sparse = FALSE),
      loadings = Matrix(matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE), sparse = FALSE),
      space = NeuroSpace(c(2, 2, 1, 3)),
      mask = mask_vol
    )
  }

  expect_error(
    BlockLatentNeuroVector(list(
      cortex = BilatLatentNeuroSurfaceVector(left, right),
      subcortex = bad_volume
    )),
    "basis matrix must match the shared block basis"
  )
})
