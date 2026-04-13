library(testthat)
library(Matrix)
library(neuroim2)

.make_awpt_mask_analytic <- function(n = 2L) {
  LogicalNeuroVol(array(TRUE, dim = c(n, 1, 1)), NeuroSpace(c(n, 1, 1)))
}

.make_awpt_template_analytic <- function(loadings, roughness) {
  mask_vol <- .make_awpt_mask_analytic(nrow(loadings))
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
      meta = list(
        family = "awpt_wavelet",
        method = "AWPT",
        scales = rep(1, ncol(loadings)),
        k = ncol(loadings)
      )
    ),
    class = "AWPTBasisTemplate"
  )
}

.make_identity_operator_analytic <- function(n) {
  .linear_map_from_matrix(
    diag(n),
    source_domain_id = "template_field",
    target_domain_id = "native_field",
    provenance = list(target_mask = .make_awpt_mask_analytic(n))
  )
}

test_that("encode_awpt matches the closed-form quadratic solution without temporal smoothing", {
  tmpl <- .make_awpt_template_analytic(
    loadings = diag(2),
    roughness = diag(c(0.5, 2), nrow = 2)
  )
  op <- .make_identity_operator_analytic(2)
  Y <- matrix(
    c(1, 2,
      3, 4),
    nrow = 2, byrow = TRUE
  )
  spatial_lambda <- 0.7

  fit <- encode_awpt(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    spatial_lambda = spatial_lambda,
    temporal_lambda = 0
  )

  Q <- as.matrix(template_roughness(tmpl))
  expected <- Y %*% solve(diag(2) + spatial_lambda * Q)

  expect_equal(as.matrix(coef_time(fit)), expected, tolerance = 1e-8)
})

test_that("encode_awpt matches the closed-form temporal smoothing solution in the scalar case", {
  tmpl <- .make_awpt_template_analytic(
    loadings = matrix(1, nrow = 1, ncol = 1),
    roughness = matrix(0, nrow = 1, ncol = 1)
  )
  op <- .make_identity_operator_analytic(1)
  Y <- matrix(c(0, 10, -10, 0), ncol = 1)
  temporal_lambda <- 3

  fit <- encode_awpt(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    spatial_lambda = 0,
    temporal_lambda = temporal_lambda,
    temporal_order = 1L,
    run_info = list(n_runs = 1L, run_lengths = 4L)
  )

  D1 <- diff(diag(nrow(Y)))
  Lt <- crossprod(D1)
  expected <- solve(diag(nrow(Y)) + temporal_lambda * Lt, Y)

  expect_equal(as.matrix(coef_time(fit)), expected, tolerance = 1e-8)
})
