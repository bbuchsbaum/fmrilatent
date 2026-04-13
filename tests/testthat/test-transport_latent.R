library(testthat)
library(Matrix)
library(neuroim2)

make_transport_mask <- function() {
  LogicalNeuroVol(array(TRUE, dim = c(2, 2, 1)), NeuroSpace(c(2, 2, 1)))
}

make_transport_template <- function() {
  mask_vol <- make_transport_mask()
  loadings <- Matrix(matrix(
    c(1, 0,
      0, 1,
      1, 1,
      0, 1),
    nrow = 4, byrow = TRUE
  ), sparse = FALSE)
  gram_factor <- Matrix::Cholesky(Matrix::crossprod(loadings), perm = TRUE)
  reduction <- make_cluster_reduction(mask_vol, c(1L, 1L, 2L, 2L))
  structure(
    list(
      loadings = loadings,
      gram_factor = gram_factor,
      reduction = reduction,
      basis_spec = basis_slepian(k = 2),
      center = FALSE,
      meta = list(
        family = "spec_slepian",
        k = 2L,
        ridge = 1e-8,
        label_map = list(A = 1L, B = 2L),
        cluster_map = list(`1` = 1:2, `2` = 3:4)
      )
    ),
    class = "ParcelBasisTemplate"
  )
}

make_transport_operator <- function() {
  A <- matrix(
    c(1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 0.5, 0.5,
      0, 0, 0.5, -0.5),
    nrow = 4, byrow = TRUE
  )
  .linear_map_from_matrix(
    A,
    source_domain_id = "template_field",
    target_domain_id = "native_field",
    provenance = list(
      operator = "toy_transport",
      target_mask = make_transport_mask()
    )
  )
}

make_transport_operator_no_mask <- function() {
  A <- matrix(
    c(1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 0.5, 0.5,
      0, 0, 0.5, -0.5),
    nrow = 4, byrow = TRUE
  )
  .linear_map_from_matrix(
    A,
    source_domain_id = "template_field",
    target_domain_id = "native_field",
    provenance = list(operator = "toy_transport")
  )
}

make_transport_operator_no_materialize <- function() {
  A <- matrix(
    c(1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 0.5, 0.5,
      0, 0, 0.5, -0.5),
    nrow = 4, byrow = TRUE
  )
  list(
    mode = "callbacks",
    n_source = ncol(A),
    n_target = nrow(A),
    source_domain_id = "template_field",
    target_domain_id = "native_field",
    provenance = list(
      operator = "toy_transport_callbacks",
      target_mask = make_transport_mask()
    ),
    forward = function(data, ...) A %*% as.matrix(data),
    adjoint_apply = function(data, ...) t(A) %*% as.matrix(data),
    materialize = function(...) stop("materialize forbidden", call. = FALSE)
  )
}

test_that("parcel template exposes basis decoder protocol", {
  tmpl <- make_transport_template()
  dec <- basis_decoder(tmpl)
  payload <- .template_coordinate_payload(
    raw_loadings = template_loadings(tmpl),
    measure = template_measure(tmpl),
    default_measure = "unit"
  )

  expect_equal(template_rank(tmpl), 2L)
  expect_true(!is.null(template_domain(tmpl)))
  expect_equal(dec$n_source, 2L)
  expect_equal(dec$n_target, 4L)

  gamma <- matrix(c(1, 2), ncol = 1)
  expect_equal(as.vector(dec$forward(gamma)), as.vector(payload$analysis_loadings %*% gamma))
})

test_that("encode_transport recovers coefficients and reconstructs data", {
  tmpl <- make_transport_template()
  op <- make_transport_operator()
  B <- .materialize_linear_map(basis_decoder(tmpl))
  A <- .materialize_linear_map(op)
  D <- A %*% B

  Z <- matrix(
    c(1, 2,
      0, 1,
      -1, 3),
    nrow = 3, byrow = TRUE
  )
  Y <- Z %*% t(D)

  tl <- encode_transport(
    x = Y,
    basis_asset = tmpl,
    field_operator = op,
    center = FALSE,
    lambda = 0
  )

  transform <- .template_asset_analysis_transform(tmpl)
  Z_raw <- t(transform$to_raw(t(Z)))

  expect_true(is_transport_latent(tl))
  expect_equal(unname(as.matrix(coef_time(tl, "analysis"))), unname(Z), tolerance = 1e-8)
  expect_equal(unname(as.matrix(coef_time(tl, "raw"))), unname(Z_raw), tolerance = 1e-8)
  expect_equal(unname(reconstruct_matrix(tl)), unname(Y), tolerance = 1e-8)

  dec_native <- decoder(tl, space = "native", coordinates = "analysis")
  dec_template <- decoder(tl, space = "template", coordinates = "analysis")
  expect_equal(dec_native$n_source, 2L)
  expect_equal(dec_native$n_target, 4L)
  expect_equal(dec_template$n_source, 2L)
  expect_equal(dec_template$n_target, 4L)
  expect_equal(latent_meta(tl)$target_mask_source, "operator_provenance")
})

test_that("encode_transport accepts field_operator as the preferred seam name", {
  tmpl <- make_transport_template()
  op <- make_transport_operator()
  B <- .materialize_linear_map(basis_decoder(tmpl))
  A <- .materialize_linear_map(op)
  D <- A %*% B
  Y <- matrix(c(1, 0), nrow = 1) %*% t(D)

  fit <- encode_transport(
    x = Y,
    basis_asset = tmpl,
    field_operator = op,
    center = FALSE,
    lambda = 0
  )

  expect_true(is_transport_latent(fit))
  expect_identical(fit$field_operator, op)
  expect_identical(fit$observation_operator, op)
})

test_that("encode_transport requires target-domain mask when operator omits one", {
  tmpl <- make_transport_template()
  op <- make_transport_operator_no_mask()
  B <- .materialize_linear_map(basis_decoder(tmpl))
  A <- .materialize_linear_map(op)
  D <- A %*% B
  Y <- matrix(c(1, 0), nrow = 1) %*% t(D)

  expect_error(
    encode_transport(
      x = Y,
      basis_asset = tmpl,
      observation_operator = op,
      center = FALSE,
      lambda = 0
    ),
    "requires an explicit target mask"
  )

  fit <- encode_transport(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    mask = make_transport_mask(),
    center = FALSE,
    lambda = 0
  )

  expect_true(is_transport_latent(fit))
  expect_equal(latent_meta(fit)$target_mask_source, "explicit")
})

test_that("encode_transport quadratic path does not materialize the decoder", {
  tmpl <- make_transport_template()
  op <- make_transport_operator_no_materialize()
  B <- .materialize_linear_map(basis_decoder(tmpl))
  A <- matrix(
    c(1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 0.5, 0.5,
      0, 0, 0.5, -0.5),
    nrow = 4, byrow = TRUE
  )
  D <- A %*% B
  Z <- matrix(
    c(1, 2,
      0, 1,
      -1, 3),
    nrow = 3, byrow = TRUE
  )
  Y <- Z %*% t(D)

  fit <- encode_transport(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    lambda = 0
  )

  expect_equal(unname(as.matrix(coef_time(fit))), unname(Z), tolerance = 1e-8)
})

test_that("encode_transport sparse path does not materialize the decoder", {
  tmpl <- make_transport_template()
  op <- make_transport_operator_no_materialize()
  Y <- rbind(c(1, 0.2, 0.1, 0.1), c(1, 0.2, 0.1, 0.1))

  fit <- encode_transport(
    x = Y,
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    lambda = 0,
    sparse_lambda = 0.2,
    sparse_mode = "group_l2",
    max_iter = 500L,
    tol = 1e-10
  )

  expect_true(is_transport_latent(fit))
  expect_true(all(is.finite(as.matrix(coef_time(fit)))))
})

test_that("encode_transport rejects target masks with wrong cardinality", {
  tmpl <- make_transport_template()
  op <- make_transport_operator_no_mask()
  bad_mask <- LogicalNeuroVol(
    array(c(TRUE, TRUE, TRUE, FALSE), dim = c(2, 2, 1)),
    NeuroSpace(c(2, 2, 1))
  )
  Y <- matrix(0, nrow = 1, ncol = 4)

  expect_error(
    encode_transport(
      x = Y,
      basis_asset = tmpl,
      observation_operator = op,
      mask = bad_mask,
      center = FALSE,
      lambda = 0
    ),
    "target mask cardinality 3 does not match field operator target dimension 4"
  )
})

test_that("encode_transport rejects malformed field operators early", {
  tmpl <- make_transport_template()
  Y <- matrix(0, nrow = 1, ncol = 4)

  expect_error(
    encode_transport(
      x = Y,
      basis_asset = tmpl,
      observation_operator = list(
        n_source = 4,
        n_target = 4,
        adjoint_apply = function(data, ...) as.matrix(data),
        provenance = list(target_mask = make_transport_mask())
      ),
      center = FALSE,
      lambda = 0
    ),
    "subject field operator is missing required fields: forward"
  )

  expect_error(
    encode_transport(
      x = Y,
      basis_asset = tmpl,
      observation_operator = list(
        n_source = 0,
        n_target = 4,
        forward = function(data, ...) as.matrix(data),
        adjoint_apply = function(data, ...) as.matrix(data),
        provenance = list(target_mask = make_transport_mask())
      ),
      center = FALSE,
      lambda = 0
    ),
    "positive scalar n_source"
  )
})

test_that("encode_transport detects invalid forward output dimensions", {
  tmpl <- make_transport_template()
  Y <- matrix(c(1, 0, 0, 0), nrow = 1)
  bad_op <- list(
    mode = "callbacks",
    n_source = 4,
    n_target = 4,
    source_domain_id = "template_field",
    target_domain_id = "native_field",
    provenance = list(target_mask = make_transport_mask()),
    forward = function(data, ...) matrix(0, nrow = 3, ncol = ncol(as.matrix(data))),
    adjoint_apply = function(data, ...) matrix(1, nrow = 4, ncol = ncol(as.matrix(data)))
  )

  expect_error(
    encode_transport(
      x = Y,
      basis_asset = tmpl,
      observation_operator = bad_op,
      center = FALSE,
      lambda = 0
    ),
    "transport forward gram returned 3 rows; expected 4"
  )
})

test_that("transport latent projects effects and covariance", {
  tmpl <- make_transport_template()
  op <- make_transport_operator()
  B <- .materialize_linear_map(basis_decoder(tmpl))
  A <- .materialize_linear_map(op)
  D <- A %*% B
  Z <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  tl <- encode_transport(
    x = Z %*% t(D),
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    lambda = 0,
    run_info = list(n_runs = 1L, run_lengths = 2L)
  )

  gamma <- c(2, -1)
  Sigma <- matrix(c(2, 0.5,
                    0.5, 1), nrow = 2, byrow = TRUE)

  expect_equal(decode_coefficients(tl, gamma, space = "template"),
               as.vector(B %*% matrix(gamma, ncol = 1)),
               tolerance = 1e-8)
  expect_equal(decode_coefficients(tl, gamma, space = "native"),
               as.vector(D %*% matrix(gamma, ncol = 1)),
               tolerance = 1e-8)

  expect_equal(project_effect(tl, gamma, space = "template"),
               as.vector(B %*% matrix(gamma, ncol = 1)),
               tolerance = 1e-8)
  expect_equal(project_effect(tl, gamma, space = "native"),
               as.vector(D %*% matrix(gamma, ncol = 1)),
               tolerance = 1e-8)

  expect_equal(project_vcov(tl, Sigma, space = "template", diag_only = TRUE),
               diag(B %*% Sigma %*% t(B)),
               tolerance = 1e-8)
  expect_equal(project_vcov(tl, Sigma, space = "native", diag_only = TRUE),
               diag(D %*% Sigma %*% t(D)),
               tolerance = 1e-8)
  expect_equal(decode_covariance(tl, Sigma, space = "template", diag_only = TRUE),
               diag(B %*% Sigma %*% t(B)),
               tolerance = 1e-8)
  expect_equal(decode_covariance(tl, Sigma, space = "native", diag_only = TRUE),
               diag(D %*% Sigma %*% t(D)),
               tolerance = 1e-8)

  expect_equal(latent_meta(tl)$run_info$run_lengths, 2L)
})

test_that("decode_covariance diag_only path does not materialize the decoder", {
  tmpl <- make_transport_template()
  op <- make_transport_operator_no_materialize()
  B <- .materialize_linear_map(basis_decoder(tmpl))
  A <- matrix(
    c(1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 0.5, 0.5,
      0, 0, 0.5, -0.5),
    nrow = 4, byrow = TRUE
  )
  D <- A %*% B
  Z <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  tl <- encode_transport(
    x = Z %*% t(D),
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    lambda = 0
  )
  Sigma <- matrix(c(2, 0.5,
                    0.5, 1), nrow = 2, byrow = TRUE)

  expect_equal(
    decode_covariance(tl, Sigma, space = "native", diag_only = TRUE),
    diag(D %*% Sigma %*% t(D)),
    tolerance = 1e-8
  )
})

test_that("transport latent respects raw versus analysis coordinate decoding", {
  tmpl <- make_transport_template()
  op <- make_transport_operator()
  B <- .materialize_linear_map(basis_decoder(tmpl))
  A <- .materialize_linear_map(op)
  D <- A %*% B
  M <- diag(c(2, 0.5), nrow = 2)
  transform <- list(
    type = "linear",
    dim = 2L,
    to_analysis = function(x) M %*% as.matrix(x),
    to_raw = function(x) solve(M, as.matrix(x)),
    matrix = M
  )

  tl <- transport_latent(
    coeff_raw = matrix(c(1, 4), nrow = 1),
    coeff_analysis = matrix(c(2, 2), nrow = 1),
    basis_asset = tmpl,
    observation_operator = op,
    mask = make_transport_mask(),
    analysis_transform = transform,
    meta = list(label = "raw-analysis-test")
  )

  gamma_raw <- c(1, 4)
  gamma_analysis <- as.vector(M %*% matrix(gamma_raw, ncol = 1))
  Sigma_raw <- matrix(c(2, 0.25,
                        0.25, 3), nrow = 2, byrow = TRUE)
  Sigma_analysis <- M %*% Sigma_raw %*% t(M)

  expect_equal(coef_time(tl, "raw"), matrix(gamma_raw, nrow = 1))
  expect_equal(coef_time(tl, "analysis"), matrix(gamma_analysis, nrow = 1))
  expect_equal(coef_metric(tl, "analysis"), diag(2))
  expect_equal(coef_metric(tl, "raw"), crossprod(M))

  expect_equal(
    decode_coefficients(tl, gamma_raw, space = "native", coordinates = "raw"),
    as.vector(D %*% matrix(gamma_analysis, ncol = 1)),
    tolerance = 1e-8
  )
  expect_equal(
    decode_coefficients(tl, gamma_raw, space = "native", coordinates = "raw"),
    decode_coefficients(tl, gamma_analysis, space = "native", coordinates = "analysis"),
    tolerance = 1e-8
  )

  expected_cov <- D %*% Sigma_analysis %*% t(D)
  expect_equal(
    decode_covariance(tl, Sigma_raw, space = "native", coordinates = "raw", diag_only = FALSE),
    expected_cov,
    tolerance = 1e-8
  )
  expect_equal(
    reconstruct_matrix(tl),
    matrix(as.vector(D %*% matrix(gamma_analysis, ncol = 1)), nrow = 1),
    tolerance = 1e-8
  )
})

test_that("downstream consumers can rely on the analysis-space handoff contract", {
  downstream_contract <- function(latent) {
    response <- as.matrix(coef_time(latent, "analysis"))
    metric <- coef_metric(latent, "analysis")
    native_decoder <- decoder(latent, space = "native", coordinates = "analysis")
    template_decoder <- decoder(latent, space = "template", coordinates = "analysis")

    list(
      response = response,
      metric = metric,
      basis_asset = basis_asset(latent),
      native_decoder = native_decoder,
      template_decoder = template_decoder,
      project_effect = function(gamma, space = c("native", "template")) {
        decode_coefficients(latent, gamma, space = match.arg(space), coordinates = "analysis")
      },
      project_vcov = function(Sigma, space = c("native", "template"), diag_only = TRUE) {
        decode_covariance(
          latent,
          Sigma,
          space = match.arg(space),
          coordinates = "analysis",
          diag_only = diag_only
        )
      }
    )
  }

  tmpl <- make_transport_template()
  op <- make_transport_operator()
  B <- .materialize_linear_map(basis_decoder(tmpl))
  A <- .materialize_linear_map(op)
  D <- A %*% B
  M <- diag(c(2, 0.5), nrow = 2)
  gamma_raw <- c(1, 4)
  gamma_analysis <- as.vector(M %*% matrix(gamma_raw, ncol = 1))
  Sigma_raw <- matrix(c(2, 0.25,
                        0.25, 3), nrow = 2, byrow = TRUE)
  Sigma_analysis <- M %*% Sigma_raw %*% t(M)

  tl <- transport_latent(
    coeff_raw = matrix(gamma_raw, nrow = 1),
    coeff_analysis = matrix(gamma_analysis, nrow = 1),
    basis_asset = tmpl,
    observation_operator = op,
    mask = make_transport_mask(),
    analysis_transform = list(
      type = "linear",
      dim = 2L,
      to_analysis = function(x) M %*% as.matrix(x),
      to_raw = function(x) solve(M, as.matrix(x)),
      matrix = M
    ),
    meta = list(label = "downstream-contract")
  )

  handoff <- downstream_contract(tl)

  expect_equal(handoff$response, matrix(gamma_analysis, nrow = 1))
  expect_equal(handoff$metric, diag(2))
  expect_s3_class(handoff$basis_asset, "ParcelBasisTemplate")
  expect_equal(handoff$native_decoder$n_source, 2L)
  expect_equal(handoff$template_decoder$n_source, 2L)
  expect_equal(handoff$native_decoder$n_target, 4L)
  expect_equal(handoff$template_decoder$n_target, 4L)
  expect_equal(
    handoff$project_effect(gamma_analysis, "native"),
    as.vector(D %*% matrix(gamma_analysis, ncol = 1)),
    tolerance = 1e-8
  )
  expect_equal(
    handoff$project_effect(gamma_analysis, "template"),
    as.vector(B %*% matrix(gamma_analysis, ncol = 1)),
    tolerance = 1e-8
  )
  expect_equal(
    handoff$project_vcov(Sigma_analysis, "native", diag_only = TRUE),
    diag(D %*% Sigma_analysis %*% t(D)),
    tolerance = 1e-8
  )
  expect_equal(
    handoff$project_vcov(Sigma_analysis, "template", diag_only = TRUE),
    diag(B %*% Sigma_analysis %*% t(B)),
    tolerance = 1e-8
  )
})

test_that("decode_covariance full pushforward is symmetric PSD", {
  tmpl <- make_transport_template()
  op <- make_transport_operator()
  B <- .materialize_linear_map(basis_decoder(tmpl))
  A <- .materialize_linear_map(op)
  D <- A %*% B
  Z <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  tl <- encode_transport(
    x = Z %*% t(D),
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    lambda = 0
  )
  Sigma <- matrix(c(2, 0.5,
                    0.5, 1), nrow = 2, byrow = TRUE)

  cov_native <- decode_covariance(tl, Sigma, space = "native", diag_only = FALSE)
  cov_template <- decode_covariance(tl, Sigma, space = "template", diag_only = FALSE)

  expect_equal(cov_native, t(cov_native), tolerance = 1e-8)
  expect_equal(cov_template, t(cov_template), tolerance = 1e-8)
  expect_gte(min(eigen(cov_native, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
  expect_gte(min(eigen(cov_template, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
})

test_that("LatentNeuroVec satisfies minimal coefficient protocol", {
  mask_vol <- make_transport_mask()
  lv <- LatentNeuroVec(
    basis = Matrix(matrix(c(1, 2, 3, 4), nrow = 2), sparse = FALSE),
    loadings = Matrix(matrix(c(1, 0, 0, 1, 1, 1, 0, 0), nrow = 4, byrow = TRUE), sparse = FALSE),
    space = NeuroSpace(c(2, 2, 1, 2)),
    mask = mask_vol
  )

  expect_equal(coef_time(lv, "analysis"), as.matrix(basis(lv)))
  expect_equal(coef_metric(lv, "analysis"), diag(2))
  expect_equal(analysis_transform(lv)$type, "identity")
  expect_null(basis_asset(lv))
  gamma <- c(2, -1)
  Sigma <- matrix(c(2, 0.5, 0.5, 1), nrow = 2, byrow = TRUE)
  L <- as.matrix(loadings(lv))
  expect_equal(decode_coefficients(lv, gamma, space = "native"),
               as.vector(L %*% matrix(gamma, ncol = 1)))
  expect_equal(decode_covariance(lv, Sigma, space = "native", diag_only = TRUE),
               diag(L %*% Sigma %*% t(L)))
})

test_that("LatentNeuroVec full covariance pushforward is symmetric PSD", {
  mask_vol <- make_transport_mask()
  lv <- LatentNeuroVec(
    basis = Matrix(matrix(c(1, 2, 3, 4), nrow = 2), sparse = FALSE),
    loadings = Matrix(matrix(c(1, 0, 0, 1, 1, 1, 0, 0), nrow = 4, byrow = TRUE), sparse = FALSE),
    space = NeuroSpace(c(2, 2, 1, 2)),
    mask = mask_vol
  )
  Sigma <- matrix(c(2, 0.5,
                    0.5, 1), nrow = 2, byrow = TRUE)

  cov_native <- suppressWarnings(decode_covariance(lv, Sigma, space = "native", diag_only = FALSE))

  expect_equal(cov_native, t(cov_native), tolerance = 1e-8)
  expect_gte(min(eigen(cov_native, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
})

test_that("decoder-centric APIs fail clearly on coefficient dimension mismatch", {
  tmpl <- make_transport_template()
  op <- make_transport_operator()
  B <- .materialize_linear_map(basis_decoder(tmpl))
  A <- .materialize_linear_map(op)
  D <- A %*% B
  Z <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  tl <- encode_transport(
    x = Z %*% t(D),
    basis_asset = tmpl,
    observation_operator = op,
    center = FALSE,
    lambda = 0
  )
  lv <- LatentNeuroVec(
    basis = Matrix(matrix(c(1, 2, 3, 4), nrow = 2), sparse = FALSE),
    loadings = Matrix(matrix(c(1, 0, 0, 1, 1, 1, 0, 0), nrow = 4, byrow = TRUE), sparse = FALSE),
    space = NeuroSpace(c(2, 2, 1, 2)),
    mask = make_transport_mask()
  )

  expect_error(
    decode_coefficients(tl, c(1, 2, 3), space = "native"),
    "must have 2 rows"
  )
  expect_error(
    decode_covariance(tl, diag(3), space = "native"),
    "must have dimensions 2x2"
  )
  expect_error(
    decode_coefficients(lv, c(1, 2, 3), space = "native"),
    "must have 2 rows"
  )
  expect_error(
    decode_covariance(lv, diag(3), space = "native"),
    "must have dimensions 2x2"
  )
})
