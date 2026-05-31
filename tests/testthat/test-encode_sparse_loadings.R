# Regression tests for sparse-loadings handling in the spatial encode hot paths
# (bd-01KSQQNTXWG2JBBXXQAWRE51GG).
#
# The Slepian and heat spatial encoders previously densified the loadings
# dictionary (as.matrix(loadings_mat(...))) before forming the temporal basis.
# The lift dictionaries are genuinely sparse (dgCMatrix), so the multiply now
# stays in Matrix space and only the small time x k product is densified. These
# tests pin both that the sparse path is exercised and that the result matches
# the dense reference exactly.

library(testthat)

encode_space_sparse_check <- function(spec, mask_dim = c(3, 3, 2), n_time = 6L, seed = 17L) {
  mask_arr <- array(TRUE, dim = mask_dim)
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(mask_dim))
  nv <- sum(mask_arr)
  set.seed(seed)
  X <- matrix(rnorm(n_time * nv), nrow = n_time)

  lv <- without_dense_basis_warning(
    encode(X, spec, mask = mask, materialize = "matrix")
  )

  # The loadings dictionary must actually be sparse, otherwise this test would
  # silently stop guarding the sparse-aware multiply.
  Lmat <- loadings_mat(loadings(lv))
  expect_true(methods::is(Lmat, "sparseMatrix"))

  # Sparse-aware basis == dense reference, to numerical precision.
  reference <- as.matrix(X) %*% as.matrix(Lmat)
  expect_equal(as.matrix(basis(lv)), reference, tolerance = 1e-10)

  invisible(lv)
}

test_that("spec_space_slepian keeps sparse loadings and matches the dense basis", {
  skip_if_not_installed("rgsp")
  encode_space_sparse_check(spec_space_slepian(k = 3L, k_neighbors = 4L))
})

test_that("spec_space_heat keeps sparse loadings and matches the dense basis", {
  skip_if_not_installed("rgsp")
  encode_space_sparse_check(
    spec_space_heat(scales = c(1, 2), order = 16, threshold = 0, k_neighbors = 4L)
  )
})
