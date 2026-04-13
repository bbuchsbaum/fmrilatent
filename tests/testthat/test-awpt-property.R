library(testthat)
library(Matrix)
library(neuroim2)

.make_awpt_mask_property <- function(n) {
  LogicalNeuroVol(array(TRUE, dim = c(n, 1, 1)), NeuroSpace(c(n, 1, 1)))
}

.make_awpt_template_property <- function(loadings, roughness) {
  mask_vol <- .make_awpt_mask_property(nrow(loadings))
  reduction <- make_cluster_reduction(mask_vol, rep(1L, nrow(loadings)))
  structure(
    list(
      loadings = Matrix(loadings, sparse = FALSE),
      gram_factor = Matrix::Cholesky(Matrix::crossprod(loadings), perm = TRUE),
      roughness = Matrix(roughness, sparse = FALSE),
      reduction = reduction,
      basis_spec = basis_awpt_wavelet(scales = rep(1, ncol(loadings))),
      atoms = data.frame(
        col_id = seq_len(ncol(loadings)),
        cluster_id = 1L,
        scale = rep(1, ncol(loadings)),
        scale_index = seq_len(ncol(loadings)),
        roughness_weight = diag(as.matrix(roughness)),
        stringsAsFactors = FALSE
      ),
      center = FALSE,
      meta = list(family = "awpt_wavelet", method = "AWPT", k = ncol(loadings))
    ),
    class = "AWPTBasisTemplate"
  )
}

.random_spd <- function(n, ridge = 0.2) {
  A <- matrix(rnorm(n * n), nrow = n)
  crossprod(A) + diag(ridge, n)
}

test_that("randomized AWPT small systems remain finite and covariance pushforward stays PSD", {
  set.seed(42)

  for (iter in seq_len(5)) {
    loadings <- diag(3) + matrix(rnorm(9, sd = 0.05), nrow = 3)
    roughness <- .random_spd(3, ridge = 0.1)
    tmpl <- .make_awpt_template_property(loadings, roughness)
    op_mat <- diag(3) + matrix(rnorm(9, sd = 0.05), nrow = 3)
    op <- .linear_map_from_matrix(
      op_mat,
      source_domain_id = "template_field",
      target_domain_id = "native_field",
      provenance = list(target_mask = .make_awpt_mask_property(3))
    )
    Y <- matrix(rnorm(12), nrow = 4)

    fit <- encode_awpt(
      x = Y,
      basis_asset = tmpl,
      observation_operator = op,
      center = FALSE,
      spatial_lambda = 0.4,
      temporal_lambda = 0.2,
      temporal_order = 1L,
      run_info = list(n_runs = 2L, run_lengths = c(2L, 2L))
    )

    coef_mat <- as.matrix(coef_time(fit))
    recon <- reconstruct_matrix(fit)
    Sigma <- .random_spd(3, ridge = 0.05)
    cov_native <- decode_covariance(fit, Sigma, space = "native", diag_only = FALSE)

    expect_true(all(is.finite(coef_mat)))
    expect_true(all(is.finite(recon)))
    expect_equal(dim(coef_mat), c(4, 3))
    expect_equal(dim(recon), c(4, 3))
    expect_equal(cov_native, t(cov_native), tolerance = 1e-8)
    expect_gte(min(eigen(cov_native, symmetric = TRUE, only.values = TRUE)$values), -1e-7)
  }
})

test_that("randomized conductance aggregation and operator construction remain finite", {
  set.seed(99)

  conductances <- replicate(4, {
    A <- matrix(rnorm(9), nrow = 3)
    crossprod(A) + diag(0.5, 3)
  }, simplify = FALSE)

  out <- awpt_operator_from_subject_conductances(
    conductances,
    mean_method = "log_euclidean",
    normalize = "sym",
    shrinkage = 0.1,
    enforce_psd = TRUE
  )

  expect_true(is.list(out))
  expect_true(all(is.finite(as.matrix(out$conductance_mean))))
  expect_true(all(is.finite(as.matrix(out$operator))))
  expect_equal(as.matrix(out$conductance_mean), t(as.matrix(out$conductance_mean)), tolerance = 1e-8)
  expect_equal(as.matrix(out$operator), t(as.matrix(out$operator)), tolerance = 1e-8)
})
