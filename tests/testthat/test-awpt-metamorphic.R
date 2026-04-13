library(testthat)
library(Matrix)
library(neuroim2)

.make_awpt_mask_meta <- function(n = 2L) {
  LogicalNeuroVol(array(TRUE, dim = c(n, 1, 1)), NeuroSpace(c(n, 1, 1)))
}

.make_awpt_template_meta <- function(loadings, roughness_diag) {
  mask_vol <- .make_awpt_mask_meta(nrow(loadings))
  reduction <- make_cluster_reduction(mask_vol, rep(1L, nrow(loadings)))
  structure(
    list(
      loadings = Matrix(loadings, sparse = FALSE),
      gram_factor = Matrix::Cholesky(Matrix::crossprod(loadings), perm = TRUE),
      roughness = Matrix::Diagonal(x = roughness_diag),
      reduction = reduction,
      basis_spec = basis_awpt_wavelet(scales = rep(1, ncol(loadings))),
      atoms = data.frame(
        col_id = seq_len(ncol(loadings)),
        cluster_id = 1L,
        scale = rep(1, ncol(loadings)),
        scale_index = seq_len(ncol(loadings)),
        roughness_weight = roughness_diag,
        stringsAsFactors = FALSE
      ),
      center = FALSE,
      meta = list(family = "awpt_wavelet", method = "AWPT", k = ncol(loadings))
    ),
    class = "AWPTBasisTemplate"
  )
}

.make_identity_operator_meta <- function(n) {
  .linear_map_from_matrix(
    diag(n),
    source_domain_id = "template_field",
    target_domain_id = "native_field",
    provenance = list(target_mask = .make_awpt_mask_meta(n))
  )
}

test_that("decode_coefficients and decode_covariance satisfy linearity and homogeneity", {
  tmpl <- .make_awpt_template_meta(diag(2), roughness_diag = c(0, 0))
  op <- .make_identity_operator_meta(2)
  Y <- matrix(
    c(1, 0,
      0, 1),
    nrow = 2, byrow = TRUE
  )
  fit <- encode_awpt(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE
  )

  g1 <- c(1, -2)
  g2 <- c(0.5, 3)
  a <- 1.7
  b <- -0.3
  Sigma <- matrix(c(2, 0.5,
                    0.5, 1), nrow = 2, byrow = TRUE)

  lhs <- decode_coefficients(fit, a * g1 + b * g2, space = "native")
  rhs <- a * decode_coefficients(fit, g1, space = "native") +
    b * decode_coefficients(fit, g2, space = "native")
  expect_equal(lhs, rhs, tolerance = 1e-8)

  cov_full <- decode_covariance(fit, 2 * Sigma, space = "native", diag_only = FALSE)
  cov_scaled <- 2 * decode_covariance(fit, Sigma, space = "native", diag_only = FALSE)
  expect_equal(cov_full, cov_scaled, tolerance = 1e-8)
})

test_that("group-sparse shrinkage is monotone in sparse_lambda", {
  tmpl <- .make_awpt_template_meta(diag(3), roughness_diag = c(0, 0, 0))
  op <- .make_identity_operator_meta(3)
  Y <- rbind(
    c(1, 0.4, 0.05),
    c(1, 0.4, 0.05)
  )

  fit_low <- encode_awpt(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    sparse_lambda = 0.2,
    sparse_mode = "group_l2",
    max_iter = 500L,
    tol = 1e-10
  )
  fit_high <- encode_awpt(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    sparse_lambda = 0.8,
    sparse_mode = "group_l2",
    max_iter = 500L,
    tol = 1e-10
  )

  active_cols <- function(z) sum(sqrt(colSums(as.matrix(z)^2)) > 1e-6)

  expect_lte(active_cols(coef_time(fit_high)), active_cols(coef_time(fit_low)))
})

test_that("spatial roughness decreases monotonically as spatial_lambda increases", {
  tmpl <- .make_awpt_template_meta(diag(2), roughness_diag = c(0, 10))
  op <- .make_identity_operator_meta(2)
  Y <- matrix(c(1, 1), nrow = 1)

  fit_low <- encode_awpt(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    spatial_lambda = 0.2
  )
  fit_high <- encode_awpt(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    spatial_lambda = 2
  )

  Q <- as.matrix(template_roughness(tmpl))
  roughness_value <- function(z) {
    z <- as.matrix(z)
    sum((z %*% Q) * z)
  }

  expect_lte(roughness_value(coef_time(fit_high)),
             roughness_value(coef_time(fit_low)) + 1e-8)
})
