make_transport_solver_map <- function(D) {
  D <- as.matrix(D)
  list(
    n_source = ncol(D),
    n_target = nrow(D),
    forward = function(data, ...) D %*% as.matrix(data),
    adjoint_apply = function(data, ...) crossprod(D, as.matrix(data))
  )
}

transport_dense_oracle <- function(X, D, spatial_lambda = 0, spatial_penalty = NULL) {
  X <- as.matrix(X)
  D <- as.matrix(D)
  k <- ncol(D)
  Q <- .normalize_penalty_matrix(spatial_penalty, k, context = "spatial_penalty")
  A <- crossprod(D) + spatial_lambda * Q
  (X %*% D) %*% solve(A)
}

test_that("matrix-free transport coefficients match dense oracle solve", {
  D <- matrix(
    c(
      1.0, 0.2, -0.1,
      0.0, 0.9, 0.3,
      0.4, -0.2, 0.8,
      0.1, 0.5, -0.4,
      -0.3, 0.1, 0.7
    ),
    nrow = 5L,
    byrow = TRUE
  )
  X <- matrix(seq(-0.6, 1.1, length.out = 20L), nrow = 4L)
  Q <- crossprod(matrix(c(1, 0.2, -0.1, 0, 0.8, 0.3, 0, 0, 0.5), nrow = 3L))
  spatial_lambda <- 0.17

  expected <- transport_dense_oracle(X, D, spatial_lambda, Q)
  got <- .solve_transport_coefficients_matrix_free(
    X,
    make_transport_solver_map(D),
    spatial_lambda = spatial_lambda,
    spatial_penalty = Q,
    max_iter = 200L,
    tol = 1e-12
  )

  expect_equal(got, expected, tolerance = 1e-8)
})

test_that("transport Lipschitz estimate applies quadratic system once per iteration", {
  calls <- new.env(parent = emptyenv())
  calls$forward <- 0L
  calls$adjoint <- 0L
  D <- diag(3)
  map <- list(
    n_source = 3L,
    n_target = 3L,
    forward = function(data, ...) {
      calls$forward <- calls$forward + 1L
      D %*% as.matrix(data)
    },
    adjoint_apply = function(data, ...) {
      calls$adjoint <- calls$adjoint + 1L
      crossprod(D, as.matrix(data))
    }
  )

  estimate <- fmrilatent:::.estimate_transport_lipschitz(map, n_time = 2L, k = 3L, n_iter = 4L)

  expect_true(is.finite(estimate))
  expect_equal(calls$forward, 4L)
  expect_equal(calls$adjoint, 4L)
})

test_that("sparse transport path matches dense path when sparse_lambda is zero", {
  raw_basis <- matrix(
    c(
      1, 0, 0,
      1, 1, 0,
      0, 1, 1,
      0, 0, 1,
      1, -1, 1
    ),
    nrow = 5L,
    byrow = TRUE
  )
  D <- qr.Q(qr(raw_basis))[, 1:3, drop = FALSE]
  X <- matrix(c(0.5, -0.2, 1.0, 0.3, 0.7,
                -0.4, 0.8, 0.1, 1.2, -0.6),
              nrow = 2L, byrow = TRUE)

  expected <- transport_dense_oracle(X, D)
  got <- .solve_transport_coefficients_sparse(
    X,
    D,
    sparse_lambda = 0,
    sparse_mode = "lasso",
    max_iter = 50L,
    tol = 1e-12
  )

  expect_equal(got, expected, tolerance = 1e-10)
})

test_that("sparse transport path skips inactive temporal penalty allocation", {
  skip_if_not(exists("local_mocked_bindings", envir = asNamespace("testthat")),
              "testthat local_mocked_bindings unavailable")
  D <- qr.Q(qr(matrix(seq_len(20), nrow = 5L)))[, 1:3, drop = FALSE]
  X <- matrix(seq(-0.5, 0.7, length.out = 20L), nrow = 4L)
  eigen_dims <- character()

  local_mocked_bindings(
    .largest_symmetric_eigenvalue = function(mat) {
      eigen_dims <<- c(eigen_dims, paste(dim(mat), collapse = "x"))
      1
    },
    .package = "fmrilatent"
  )

  got <- .solve_transport_coefficients_sparse(
    X,
    D,
    temporal_lambda = 0,
    sparse_lambda = 0,
    sparse_mode = "none",
    max_iter = 1L,
    tol = 0
  )

  expect_equal(dim(got), c(nrow(X), ncol(D)))
  expect_equal(eigen_dims, "3x3")
})

test_that("sparse AWPT objective receives only the spatial penalty term", {
  skip_if_not(exists("local_mocked_bindings", envir = asNamespace("testthat")),
              "testthat local_mocked_bindings unavailable")
  D <- qr.Q(qr(matrix(seq_len(20), nrow = 5L)))[, 1:3, drop = FALSE]
  X <- matrix(seq(-0.5, 0.7, length.out = 20L), nrow = 4L)
  Q <- crossprod(matrix(c(1, 0.2, 0, 0, 0.7, 0.1, 0, 0, 0.4), nrow = 3L))
  spatial_lambda <- 0.4
  seen_A <- list()

  local_mocked_bindings(
    .awpt_objective = function(Z, X, D_mat, A, ...) {
      seen_A[[length(seen_A) + 1L]] <<- A
      0
    },
    .package = "fmrilatent"
  )

  .solve_transport_coefficients_sparse(
    X,
    D,
    spatial_lambda = spatial_lambda,
    spatial_penalty = Q,
    sparse_lambda = 0.1,
    sparse_mode = "lasso",
    max_iter = 1L,
    tol = 0
  )

  expect_gt(length(seen_A), 0L)
  expect_equal(seen_A[[1L]], spatial_lambda * Q, tolerance = 1e-12)
})

test_that("transport CG matches direct solve on ill-conditioned Gram", {
  singular_values <- c(1, 1e-2, 1e-4)
  D <- diag(singular_values, nrow = 3L)
  gram <- crossprod(D)
  true_coef <- matrix(c(1, -2, 0.5,
                        0.2, 0.4, -0.7),
                      nrow = 2L, byrow = TRUE)
  rhs <- true_coef %*% gram

  got <- .cg_transport_quadratic(
    rhs,
    make_transport_solver_map(D),
    max_iter = 100L,
    tol = 1e-14
  )
  expected <- rhs %*% solve(gram)

  expect_true(all(is.finite(got)))
  expect_equal(got, expected, tolerance = 1e-6)
})

test_that("transport Lipschitz estimate reaches known spectral norm", {
  set.seed(711)
  D <- diag(c(3, 2, 0.5), nrow = 3L)
  true_lipschitz <- max(eigen(crossprod(D), symmetric = TRUE, only.values = TRUE)$values)

  got <- .estimate_transport_lipschitz(
    make_transport_solver_map(D),
    n_time = 4L,
    k = 3L,
    n_iter = 60L
  )

  expect_true(is.finite(got))
  expect_gte(got, true_lipschitz - 1e-8)
})
