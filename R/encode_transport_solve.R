# Transport-coefficient solvers and matrix-free quadratic CG used by the
# encode() pipeline (encode_operator, encode_transport, sparse AWPT).
#
# Internal helpers consumed by these solvers live elsewhere:
#   - .normalize_penalty_matrix, .block_temporal_penalty,
#     .largest_symmetric_eigenvalue, .awpt_objective, .prox_sparse_awpt,
#     .robust_gram_solve   ->  R/encode.R
#   - .normalize_linear_map, .transport_vector_or_matrix_input
#                          ->  R/transport_latent.R
# All resolution happens at call time, so the relative load order within the
# package only needs to ensure encode.R is loaded (it is, via Collate).

#' @include encode.R
NULL

.soft_threshold <- function(x, thresh) {
  sign(x) * pmax(abs(x) - thresh, 0)
}

.prox_sparse_awpt <- function(Z, step, sparse_lambda = 0,
                              sparse_mode = c("none", "group_l2", "lasso")) {
  sparse_mode <- match.arg(sparse_mode)
  if (sparse_lambda <= 0 || sparse_mode == "none") {
    return(Z)
  }
  if (sparse_mode == "lasso") {
    return(.soft_threshold(Z, step * sparse_lambda))
  }

  out <- Z
  for (j in seq_len(ncol(Z))) {
    norm_j <- sqrt(sum(Z[, j]^2))
    if (norm_j == 0) {
      out[, j] <- 0
    } else {
      shrink <- max(0, 1 - (step * sparse_lambda) / norm_j)
      out[, j] <- Z[, j] * shrink
    }
  }
  out
}

.awpt_objective <- function(Z, X, D_mat, A, Lt = NULL,
                            temporal_lambda = 0,
                            sparse_lambda = 0,
                            sparse_mode = c("none", "group_l2", "lasso")) {
  sparse_mode <- match.arg(sparse_mode)
  resid <- X - Z %*% t(D_mat)
  obj <- 0.5 * sum(resid^2)
  obj <- obj + 0.5 * sum((Z %*% A) * Z)
  if (temporal_lambda > 0 && !is.null(Lt)) {
    obj <- obj + 0.5 * temporal_lambda * sum((Lt %*% Z) * Z)
  }
  if (sparse_lambda > 0 && sparse_mode != "none") {
    if (sparse_mode == "lasso") {
      obj <- obj + sparse_lambda * sum(abs(Z))
    } else {
      obj <- obj + sparse_lambda * sum(sqrt(colSums(Z^2)))
    }
  }
  obj
}

.frobenius_inner <- function(x, y) {
  sum(x * y)
}

.apply_transport_forward <- function(map, data, context = "transport forward") {
  prep <- .transport_vector_or_matrix_input(data, map$n_source, context)
  out <- map$forward(prep$data)
  out <- as.matrix(out)
  if (nrow(out) != map$n_target) {
    stop(context, " returned ", nrow(out), " rows; expected ", map$n_target, ".",
         call. = FALSE)
  }
  out
}

.apply_transport_adjoint <- function(map, data, context = "transport adjoint") {
  prep <- .transport_vector_or_matrix_input(data, map$n_target, context)
  out <- map$adjoint_apply(prep$data)
  out <- as.matrix(out)
  if (nrow(out) != map$n_source) {
    stop(context, " returned ", nrow(out), " rows; expected ", map$n_source, ".",
         call. = FALSE)
  }
  out
}

.transport_rhs <- function(X, decoder_map) {
  t(.apply_transport_adjoint(
    decoder_map,
    t(X),
    context = "transport adjoint right-hand side"
  ))
}

.transport_apply_gram <- function(Z, decoder_map) {
  t(.apply_transport_adjoint(
    decoder_map,
    .apply_transport_forward(decoder_map, t(Z), context = "transport forward gram"),
    context = "transport adjoint gram"
  ))
}

.transport_apply_quadratic_system <- function(Z, decoder_map, spatial_lambda = 0,
                                              spatial_penalty = NULL,
                                              temporal_lambda = 0,
                                              Lt = NULL) {
  out <- .transport_apply_gram(Z, decoder_map)
  if (spatial_lambda > 0) {
    out <- out + spatial_lambda * (Z %*% spatial_penalty)
  }
  if (temporal_lambda > 0 && !is.null(Lt)) {
    out <- out + temporal_lambda * (Lt %*% Z)
  }
  out
}

.transport_reconstruct_matrix <- function(Z, decoder_map) {
  t(.apply_transport_forward(
    decoder_map,
    t(Z),
    context = "transport forward reconstruction"
  ))
}

.awpt_objective_matrix_free <- function(Z, X, decoder_map, spatial_penalty = NULL,
                                        spatial_lambda = 0,
                                        temporal_lambda = 0,
                                        Lt = NULL,
                                        sparse_lambda = 0,
                                        sparse_mode = c("none", "group_l2", "lasso")) {
  sparse_mode <- match.arg(sparse_mode)
  resid <- X - .transport_reconstruct_matrix(Z, decoder_map)
  obj <- 0.5 * sum(resid^2)
  if (spatial_lambda > 0) {
    obj <- obj + 0.5 * spatial_lambda * sum((Z %*% spatial_penalty) * Z)
  }
  if (temporal_lambda > 0 && !is.null(Lt)) {
    obj <- obj + 0.5 * temporal_lambda * sum((Lt %*% Z) * Z)
  }
  if (sparse_lambda > 0 && sparse_mode != "none") {
    if (sparse_mode == "lasso") {
      obj <- obj + sparse_lambda * sum(abs(Z))
    } else {
      obj <- obj + sparse_lambda * sum(sqrt(colSums(Z^2)))
    }
  }
  obj
}

.estimate_transport_lipschitz <- function(decoder_map, n_time, k,
                                          spatial_lambda = 0,
                                          spatial_penalty = NULL,
                                          temporal_lambda = 0,
                                          Lt = NULL,
                                          n_iter = 25L) {
  Q <- .normalize_penalty_matrix(spatial_penalty, k, context = "spatial_penalty")
  V <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)
  norm_v <- sqrt(.frobenius_inner(V, V))
  if (!is.finite(norm_v) || norm_v <= 0) {
    V <- matrix(1, nrow = n_time, ncol = k)
    norm_v <- sqrt(.frobenius_inner(V, V))
  }
  V <- V / norm_v
  lambda_est <- 1

  for (iter in seq_len(as.integer(n_iter))) {
    W <- .transport_apply_quadratic_system(
      V,
      decoder_map = decoder_map,
      spatial_lambda = spatial_lambda,
      spatial_penalty = Q,
      temporal_lambda = temporal_lambda,
      Lt = Lt
    )
    norm_w <- sqrt(.frobenius_inner(W, W))
    if (!is.finite(norm_w) || norm_w <= 0) {
      return(1)
    }
    V <- W / norm_w
    AV <- .transport_apply_quadratic_system(
      V,
      decoder_map = decoder_map,
      spatial_lambda = spatial_lambda,
      spatial_penalty = Q,
      temporal_lambda = temporal_lambda,
      Lt = Lt
    )
    lambda_est <- .frobenius_inner(V, AV)
  }

  if (!is.finite(lambda_est) || lambda_est <= 0) 1 else lambda_est
}

.cg_transport_quadratic <- function(rhs, decoder_map, spatial_lambda = 0,
                                    spatial_penalty = NULL,
                                    temporal_lambda = 0,
                                    temporal_order = 1L,
                                    run_info = NULL,
                                    max_iter = 500L,
                                    tol = 1e-8) {
  rhs <- as.matrix(rhs)
  n_time <- nrow(rhs)
  k <- ncol(rhs)
  Q <- .normalize_penalty_matrix(spatial_penalty, k, context = "spatial_penalty")
  Lt <- if (temporal_lambda > 0) {
    .block_temporal_penalty(n_time, temporal_order = temporal_order, run_info = run_info)
  } else {
    NULL
  }

  Z <- matrix(0, nrow = n_time, ncol = k)
  R <- rhs
  P <- R
  rr <- .frobenius_inner(R, R)
  if (!is.finite(rr) || rr <= 0) {
    return(Z)
  }
  rr0 <- rr

  for (iter in seq_len(as.integer(max_iter))) {
    AP <- .transport_apply_quadratic_system(
      P,
      decoder_map = decoder_map,
      spatial_lambda = spatial_lambda,
      spatial_penalty = Q,
      temporal_lambda = temporal_lambda,
      Lt = Lt
    )
    denom <- .frobenius_inner(P, AP)
    if (!is.finite(denom) || denom <= 0) {
      stop("Matrix-free quadratic solver encountered a non-positive curvature direction.",
           call. = FALSE)
    }
    alpha <- rr / denom
    Z <- Z + alpha * P
    R <- R - alpha * AP
    rr_new <- .frobenius_inner(R, R)
    if (!is.finite(rr_new)) {
      stop("Matrix-free quadratic solver diverged.", call. = FALSE)
    }
    if (sqrt(rr_new) <= tol * max(1, sqrt(rr0))) {
      break
    }
    beta <- rr_new / rr
    P <- R + beta * P
    rr <- rr_new
  }

  Z
}

.solve_transport_coefficients_matrix_free <- function(X, decoder_map,
                                                      spatial_lambda = 0,
                                                      spatial_penalty = NULL,
                                                      temporal_lambda = 0,
                                                      temporal_order = 1L,
                                                      run_info = NULL,
                                                      max_iter = 500L,
                                                      tol = 1e-8) {
  X <- as.matrix(X)
  decoder_map <- .normalize_linear_map(decoder_map, context = "transport decoder")
  rhs <- .transport_rhs(X, decoder_map)
  .cg_transport_quadratic(
    rhs = rhs,
    decoder_map = decoder_map,
    spatial_lambda = spatial_lambda,
    spatial_penalty = spatial_penalty,
    temporal_lambda = temporal_lambda,
    temporal_order = temporal_order,
    run_info = run_info,
    max_iter = max_iter,
    tol = tol
  )
}

.solve_transport_coefficients_sparse <- function(X, D_mat,
                                                 spatial_lambda = 0,
                                                 spatial_penalty = NULL,
                                                 temporal_lambda = 0,
                                                 temporal_order = 1L,
                                                 run_info = NULL,
                                                 sparse_lambda = 0,
                                                 sparse_mode = c("none", "group_l2", "lasso"),
                                                 max_iter = 200L,
                                                 tol = 1e-6) {
  sparse_mode <- match.arg(sparse_mode)
  X <- as.matrix(X)
  D_mat <- as.matrix(D_mat)
  n_time <- nrow(X)
  k <- ncol(D_mat)
  gram <- crossprod(D_mat)
  Q <- .normalize_penalty_matrix(spatial_penalty, k, context = "spatial_penalty")
  A <- gram + spatial_lambda * Q
  Lt <- if (temporal_lambda > 0) {
    .block_temporal_penalty(n_time, temporal_order = temporal_order, run_info = run_info)
  } else {
    matrix(0, nrow = n_time, ncol = n_time)
  }
  rhs <- X %*% D_mat
  lipschitz <- .largest_symmetric_eigenvalue(A) +
    temporal_lambda * .largest_symmetric_eigenvalue(Lt)
  if (!is.finite(lipschitz) || lipschitz <= 0) {
    lipschitz <- 1
  }
  step <- 1 / lipschitz

  Z <- matrix(0, nrow = n_time, ncol = k)
  Y <- Z
  t_k <- 1
  prev_obj <- .awpt_objective(Z, X, D_mat, A, Lt, temporal_lambda,
                              sparse_lambda, sparse_mode = sparse_mode)

  for (iter in seq_len(as.integer(max_iter))) {
    grad <- Y %*% A - rhs
    if (temporal_lambda > 0) {
      grad <- grad + temporal_lambda * Lt %*% Y
    }
    Z_next <- .prox_sparse_awpt(
      Y - step * grad,
      step = step,
      sparse_lambda = sparse_lambda,
      sparse_mode = sparse_mode
    )
    t_next <- 0.5 * (1 + sqrt(1 + 4 * t_k^2))
    Y <- Z_next + ((t_k - 1) / t_next) * (Z_next - Z)

    obj_next <- .awpt_objective(Z_next, X, D_mat, A, Lt, temporal_lambda,
                                sparse_lambda, sparse_mode = sparse_mode)
    if (!is.finite(obj_next)) {
      stop("Sparse AWPT solver diverged.", call. = FALSE)
    }
    rel_change <- abs(obj_next - prev_obj) / max(1, abs(prev_obj))
    Z <- Z_next
    t_k <- t_next
    prev_obj <- obj_next
    if (rel_change < tol) {
      break
    }
  }

  Z
}

.solve_transport_coefficients_sparse_matrix_free <- function(X, decoder_map,
                                                             spatial_lambda = 0,
                                                             spatial_penalty = NULL,
                                                             temporal_lambda = 0,
                                                             temporal_order = 1L,
                                                             run_info = NULL,
                                                             sparse_lambda = 0,
                                                             sparse_mode = c("none", "group_l2", "lasso"),
                                                             max_iter = 200L,
                                                             tol = 1e-6) {
  sparse_mode <- match.arg(sparse_mode)
  X <- as.matrix(X)
  decoder_map <- .normalize_linear_map(decoder_map, context = "transport decoder")
  n_time <- nrow(X)
  k <- decoder_map$n_source
  Q <- .normalize_penalty_matrix(spatial_penalty, k, context = "spatial_penalty")
  Lt <- if (temporal_lambda > 0) {
    .block_temporal_penalty(n_time, temporal_order = temporal_order, run_info = run_info)
  } else {
    NULL
  }
  rhs <- .transport_rhs(X, decoder_map)
  lipschitz <- .estimate_transport_lipschitz(
    decoder_map = decoder_map,
    n_time = n_time,
    k = k,
    spatial_lambda = spatial_lambda,
    spatial_penalty = Q,
    temporal_lambda = temporal_lambda,
    Lt = Lt
  )
  step <- 1 / max(lipschitz, 1e-8)

  Z <- matrix(0, nrow = n_time, ncol = k)
  Y <- Z
  t_k <- 1
  prev_obj <- .awpt_objective_matrix_free(
    Z,
    X = X,
    decoder_map = decoder_map,
    spatial_penalty = Q,
    spatial_lambda = spatial_lambda,
    temporal_lambda = temporal_lambda,
    Lt = Lt,
    sparse_lambda = sparse_lambda,
    sparse_mode = sparse_mode
  )

  for (iter in seq_len(as.integer(max_iter))) {
    grad <- .transport_apply_quadratic_system(
      Y,
      decoder_map = decoder_map,
      spatial_lambda = spatial_lambda,
      spatial_penalty = Q,
      temporal_lambda = temporal_lambda,
      Lt = Lt
    ) - rhs
    Z_next <- .prox_sparse_awpt(
      Y - step * grad,
      step = step,
      sparse_lambda = sparse_lambda,
      sparse_mode = sparse_mode
    )
    t_next <- 0.5 * (1 + sqrt(1 + 4 * t_k^2))
    Y <- Z_next + ((t_k - 1) / t_next) * (Z_next - Z)

    obj_next <- .awpt_objective_matrix_free(
      Z_next,
      X = X,
      decoder_map = decoder_map,
      spatial_penalty = Q,
      spatial_lambda = spatial_lambda,
      temporal_lambda = temporal_lambda,
      Lt = Lt,
      sparse_lambda = sparse_lambda,
      sparse_mode = sparse_mode
    )
    if (!is.finite(obj_next)) {
      stop("Sparse matrix-free AWPT solver diverged.", call. = FALSE)
    }
    rel_change <- abs(obj_next - prev_obj) / max(1, abs(prev_obj))
    Z <- Z_next
    t_k <- t_next
    prev_obj <- obj_next
    if (rel_change < tol) {
      break
    }
  }

  Z
}

.solve_transport_coefficients <- function(X, D_mat,
                                          spatial_lambda = 0,
                                          spatial_penalty = NULL,
                                          temporal_lambda = 0,
                                          temporal_order = 1L,
                                          run_info = NULL,
                                          sparse_lambda = 0,
                                          sparse_mode = c("none", "group_l2", "lasso"),
                                          max_iter = 200L,
                                          tol = 1e-6) {
  sparse_mode <- match.arg(sparse_mode)
  X <- as.matrix(X)
  D_mat <- as.matrix(D_mat)
  n_time <- nrow(X)
  k <- ncol(D_mat)
  gram <- crossprod(D_mat)
  Q <- .normalize_penalty_matrix(spatial_penalty, k, context = "spatial_penalty")
  A <- gram + spatial_lambda * Q
  rhs <- X %*% D_mat

  if (sparse_lambda > 0 && sparse_mode != "none") {
    return(.solve_transport_coefficients_sparse(
      X = X,
      D_mat = D_mat,
      spatial_lambda = spatial_lambda,
      spatial_penalty = spatial_penalty,
      temporal_lambda = temporal_lambda,
      temporal_order = temporal_order,
      run_info = run_info,
      sparse_lambda = sparse_lambda,
      sparse_mode = sparse_mode,
      max_iter = max_iter,
      tol = tol
    ))
  }

  if (temporal_lambda <= 0) {
    sol <- rhs %*% .robust_gram_solve(A, diag(k), ridge = max(1e-8, spatial_lambda))
    colnames(sol) <- NULL
    rownames(sol) <- NULL
    return(sol)
  }

  Lt <- .block_temporal_penalty(n_time, temporal_order = temporal_order, run_info = run_info)
  # Guard: the Kronecker product materializes a (n_time*k) x (n_time*k) dense
  # matrix.  For typical fMRI dimensions this is prohibitively large; prefer
  # the matrix-free CG path (.solve_transport_coefficients_matrix_free).
  kron_dim <- as.double(n_time) * as.double(k)
  if (kron_dim > 10000) {
    stop("Temporal-penalty + materialized-basis path would create a ",
         n_time * k, "x", n_time * k,
         " dense matrix. Use the matrix-free solver instead.",
         call. = FALSE)
  }
  system_mat <- kronecker(diag(k), temporal_lambda * Lt) + kronecker(t(A), diag(n_time))
  rhs_vec <- as.vector(rhs)
  z_vec <- .robust_gram_solve(system_mat, matrix(rhs_vec, ncol = 1L),
                              ridge = max(1e-8, spatial_lambda + temporal_lambda))
  matrix(z_vec, nrow = n_time, ncol = k)
}
