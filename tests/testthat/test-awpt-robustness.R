library(testthat)
library(Matrix)
library(neuroim2)

.make_awpt_mask_robust <- function(n = 2L) {
  LogicalNeuroVol(array(TRUE, dim = c(n, 1, 1)), NeuroSpace(c(n, 1, 1)))
}

.make_awpt_template_robust <- function(loadings, roughness) {
  mask_vol <- .make_awpt_mask_robust(nrow(loadings))
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

.make_identity_operator_robust <- function(n) {
  .linear_map_from_matrix(
    diag(n),
    source_domain_id = "template_field",
    target_domain_id = "native_field",
    provenance = list(target_mask = .make_awpt_mask_robust(n))
  )
}

test_that("awpt_mean_conductance returns symmetric PSD output for random SPD inputs", {
  set.seed(1)
  mats <- replicate(5, {
    A <- matrix(rnorm(9), nrow = 3)
    crossprod(A) + diag(0.2, 3)
  }, simplify = FALSE)

  mean_arith <- awpt_mean_conductance(mats, method = "arithmetic", enforce_psd = TRUE)
  mean_log <- awpt_mean_conductance(mats, method = "log_euclidean", enforce_psd = TRUE)

  for (mat in list(mean_arith, mean_log)) {
    mat <- as.matrix(mat)
    expect_true(all(is.finite(mat)))
    expect_equal(mat, t(mat), tolerance = 1e-8)
    expect_gte(min(eigen(mat, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
  }
})

test_that("conductance averaging with PSD projection regularizes indefinite means", {
  c1 <- matrix(c(0, 3,
                 3, 0), nrow = 2, byrow = TRUE)
  c2 <- matrix(c(0, 1,
                 1, 0), nrow = 2, byrow = TRUE)

  mean_psd <- awpt_mean_conductance(
    list(c1, c2),
    method = "arithmetic",
    enforce_psd = TRUE
  )

  eigvals <- eigen(as.matrix(mean_psd), symmetric = TRUE, only.values = TRUE)$values
  expect_gte(min(eigvals), -1e-8)
})

test_that("awpt_operator_from_conductance preserves core Laplacian identities", {
  W <- matrix(c(0, 2, 1,
                2, 0, 3,
                1, 3, 0), nrow = 3, byrow = TRUE)

  L_none <- as.matrix(awpt_operator_from_conductance(W, normalize = "none"))
  L_sym <- as.matrix(awpt_operator_from_conductance(W, normalize = "sym"))

  expect_equal(L_none, t(L_none), tolerance = 1e-8)
  expect_equal(rowSums(L_none), rep(0, 3), tolerance = 1e-8)
  expect_equal(L_sym, t(L_sym), tolerance = 1e-8)
})

test_that("conductance helpers reject malformed inputs and conflicting roughness specifications", {
  expect_error(
    awpt_mean_conductance(list(diag(2), diag(3))),
    "identical dimensions"
  )
  expect_error(
    awpt_mean_conductance(list(matrix(1:6, nrow = 2))),
    "square"
  )

  reduction <- make_cluster_reduction(.make_awpt_mask_robust(2), c(1L, 2L))
  expect_error(
    awpt_basis_template(
      parcellation = reduction,
      basis_spec = basis_awpt_wavelet(scales = 1),
      loadings = diag(2),
      anatomical_operator = diag(2),
      conductance = matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
    ),
    "Supply at most one"
  )
})

test_that("encode_awpt remains finite with nearly collinear decoders", {
  loadings <- matrix(
    c(1, 1,
      1, 1 + 1e-4),
    nrow = 2, byrow = TRUE
  )
  tmpl <- .make_awpt_template_robust(
    loadings = loadings,
    roughness = diag(c(0.1, 0.2), nrow = 2)
  )
  op <- .make_identity_operator_robust(2)
  Y <- matrix(
    c(1, 2,
      3, 4),
    nrow = 2, byrow = TRUE
  )

  fit <- encode_awpt(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    spatial_lambda = 0.5,
    temporal_lambda = 0.2,
    run_info = list(n_runs = 1L, run_lengths = 2L)
  )

  expect_true(all(is.finite(as.matrix(coef_time(fit)))))
  expect_true(all(is.finite(reconstruct_matrix(fit))))
})
