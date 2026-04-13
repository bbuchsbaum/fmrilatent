library(testthat)
library(Matrix)
library(neuroim2)

make_awpt_mask <- function(n = 2L) {
  LogicalNeuroVol(array(TRUE, dim = c(n, 1, 1)), NeuroSpace(c(n, 1, 1)))
}

make_awpt_template_manual <- function(loadings, roughness_diag, scales = c(4, 1), center = FALSE) {
  mask_vol <- make_awpt_mask(nrow(loadings))
  reduction <- make_cluster_reduction(mask_vol, rep(1L, nrow(loadings)))
  atoms <- data.frame(
    col_id = seq_len(ncol(loadings)),
    cluster_id = 1L,
    scale = scales,
    scale_index = seq_along(scales),
    roughness_weight = roughness_diag,
    stringsAsFactors = FALSE
  )
  structure(
    list(
      loadings = Matrix(loadings, sparse = FALSE),
      gram_factor = Matrix::Cholesky(Matrix::crossprod(loadings), perm = TRUE),
      roughness = Matrix::Diagonal(x = roughness_diag),
      reduction = reduction,
      basis_spec = basis_awpt_wavelet(scales = scales[seq_along(roughness_diag)]),
      atoms = atoms,
      center = center,
      meta = list(
        family = "awpt_wavelet",
        method = "AWPT",
        scales = scales[seq_along(roughness_diag)],
        k = ncol(loadings)
      )
    ),
    class = "AWPTBasisTemplate"
  )
}

make_awpt_operator <- function(n = 2L) {
  .linear_map_from_matrix(
    diag(n),
    source_domain_id = "template_field",
    target_domain_id = "native_field",
    provenance = list(target_mask = make_awpt_mask(n))
  )
}

test_that("awpt_mean_conductance averages subject conductances", {
  c1 <- matrix(c(2, 0.5,
                 0.5, 1), nrow = 2, byrow = TRUE)
  c2 <- matrix(c(4, 1,
                 1, 3), nrow = 2, byrow = TRUE)

  mean_arith <- awpt_mean_conductance(list(c1, c2), method = "arithmetic", enforce_psd = TRUE)
  expect_equal(as.matrix(mean_arith), (c1 + c2) / 2, tolerance = 1e-8)
})

test_that("awpt_operator_from_conductance builds expected Laplacian", {
  conductance <- matrix(c(0, 3,
                          3, 0), nrow = 2, byrow = TRUE)
  expected <- matrix(c(3, -3,
                       -3, 3), nrow = 2, byrow = TRUE)
  L <- awpt_operator_from_conductance(conductance, normalize = "none")

  expect_equal(as.matrix(L), expected, tolerance = 1e-8)
})

test_that("awpt_operator_from_subject_conductances returns mean and operator", {
  c1 <- matrix(c(0, 2,
                 2, 0), nrow = 2, byrow = TRUE)
  c2 <- matrix(c(0, 4,
                 4, 0), nrow = 2, byrow = TRUE)

  out <- awpt_operator_from_subject_conductances(
    list(c1, c2),
    mean_method = "arithmetic",
    normalize = "none",
    enforce_psd = FALSE
  )

  expect_true(is.list(out))
  expect_equal(as.matrix(out$conductance_mean), matrix(c(0, 3, 3, 0), nrow = 2, byrow = TRUE))
  expect_equal(as.matrix(out$operator), matrix(c(3, -3, -3, 3), nrow = 2, byrow = TRUE))
})

test_that("awpt_basis_template accepts an anatomical operator", {
  reduction <- make_cluster_reduction(make_awpt_mask(2), c(1L, 2L))
  L_field <- matrix(c(2, -1,
                      -1, 2), nrow = 2, byrow = TRUE)
  tmpl <- awpt_basis_template(
    parcellation = reduction,
    basis_spec = basis_awpt_wavelet(scales = 1),
    loadings = diag(2),
    anatomical_operator = L_field
  )

  expect_true(is_awpt_template(tmpl))
  expect_equal(as.matrix(template_loadings(tmpl)), diag(2), tolerance = 1e-8)
  expect_equal(as.matrix(template_roughness(tmpl)), L_field, tolerance = 1e-8)
  expect_equal(tmpl$meta$roughness_source, "anatomical_operator")
})

test_that("awpt_basis_template accepts conductance and constructs a Laplacian roughness", {
  reduction <- make_cluster_reduction(make_awpt_mask(2), c(1L, 2L))
  conductance <- matrix(c(0, 3,
                          3, 0), nrow = 2, byrow = TRUE)
  expected_L <- matrix(c(3, -3,
                         -3, 3), nrow = 2, byrow = TRUE)
  tmpl <- awpt_basis_template(
    parcellation = reduction,
    basis_spec = basis_awpt_wavelet(scales = 1),
    loadings = diag(2),
    conductance = conductance
  )

  expect_equal(as.matrix(template_roughness(tmpl)), expected_L, tolerance = 1e-8)
  expect_equal(tmpl$meta$roughness_source, "conductance_laplacian")
})

test_that("awpt_basis_template uses operator derived from subject conductances", {
  reduction <- make_cluster_reduction(make_awpt_mask(2), c(1L, 2L))
  derived <- awpt_operator_from_subject_conductances(
    list(
      matrix(c(0, 2, 2, 0), nrow = 2, byrow = TRUE),
      matrix(c(0, 4, 4, 0), nrow = 2, byrow = TRUE)
    ),
    mean_method = "arithmetic",
    normalize = "none"
  )
  tmpl <- awpt_basis_template(
    parcellation = reduction,
    basis_spec = basis_awpt_wavelet(scales = 1),
    loadings = diag(2),
    anatomical_operator = derived$operator
  )

  expect_equal(as.matrix(template_roughness(tmpl)), as.matrix(derived$operator), tolerance = 1e-8)
})

test_that("AWPT template exposes decoder and roughness contracts", {
  tmpl <- make_awpt_template_manual(diag(2), roughness_diag = c(0.25, 1))
  dec <- basis_decoder(tmpl)

  expect_true(is_awpt_template(tmpl))
  expect_true(is_template(tmpl))
  expect_equal(template_rank(tmpl), 2L)
  expect_equal(dim(as.matrix(template_roughness(tmpl))), c(2, 2))
  expect_equal(diag(as.matrix(template_roughness(tmpl))), c(0.25, 1))
  expect_equal(dec$n_source, 2L)
  expect_equal(dec$n_target, 2L)
})

test_that("encode_awpt applies spatial roughness shrinkage", {
  tmpl <- make_awpt_template_manual(diag(2), roughness_diag = c(0, 10))
  op <- make_awpt_operator(2)
  Y <- matrix(c(1, 1), nrow = 1)

  fit_free <- encode_awpt(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    spatial_lambda = 0
  )
  fit_pen <- encode_awpt(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    spatial_lambda = 1
  )

  coef_free <- as.numeric(coef_time(fit_free))
  coef_pen <- as.numeric(coef_time(fit_pen))

  expect_equal(coef_free, c(1, 1), tolerance = 1e-8)
  expect_true(abs(coef_pen[[2]]) < abs(coef_free[[2]]))
  expect_true(abs(coef_pen[[1]]) >= abs(coef_pen[[2]]))
  expect_equal(latent_meta(fit_pen)$method, "AWPT")
})

test_that("encode_awpt temporal smoothing respects run boundaries", {
  tmpl <- make_awpt_template_manual(matrix(1, nrow = 1, ncol = 1), roughness_diag = 0)
  op <- make_awpt_operator(1)
  Y <- matrix(c(0, 10, -10, 0), ncol = 1)

  fit_single_run <- encode_awpt(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    spatial_lambda = 0,
    temporal_lambda = 10,
    temporal_order = 1L,
    run_info = list(n_runs = 1L, run_lengths = 4L)
  )
  fit_split_runs <- encode_awpt(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    spatial_lambda = 0,
    temporal_lambda = 10,
    temporal_order = 1L,
    run_info = list(n_runs = 2L, run_lengths = c(2L, 2L))
  )

  z_single <- as.numeric(coef_time(fit_single_run))
  z_split <- as.numeric(coef_time(fit_split_runs))

  expect_true(abs(z_split[[2]] - z_split[[3]]) > abs(z_single[[2]] - z_single[[3]]))
})

test_that("encode_awpt supports group-sparse atom shrinkage", {
  tmpl <- make_awpt_template_manual(diag(2), roughness_diag = c(0, 0))
  op <- make_awpt_operator(2)
  Y <- rbind(c(1, 0.2), c(1, 0.2))

  fit_sparse <- encode_awpt(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    spatial_lambda = 0,
    temporal_lambda = 0,
    sparse_lambda = 0.4,
    sparse_mode = "group_l2",
    max_iter = 500L,
    tol = 1e-10
  )

  coef_sparse <- coef_time(fit_sparse)
  expected_col1 <- c(1, 1) * max(0, 1 - 0.4 / sqrt(2))
  expected_col2 <- c(0, 0)

  expect_equal(as.numeric(coef_sparse[, 1]), expected_col1, tolerance = 1e-4)
  expect_equal(as.numeric(coef_sparse[, 2]), expected_col2, tolerance = 1e-6)
  expect_equal(latent_meta(fit_sparse)$sparse_mode, "group_l2")
})
