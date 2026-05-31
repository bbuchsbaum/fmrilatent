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

# Additive sparsity penalty term shared by the dense and matrix-free AWPT
# objectives. Returns 0 when sparsity is inactive. sparse_mode is assumed
# already validated by the caller (objectives match.arg it up front).
.awpt_penalty_term <- function(Z, sparse_lambda = 0,
                               sparse_mode = c("none", "group_l2", "lasso")) {
  if (sparse_lambda <= 0 || sparse_mode == "none") {
    return(0)
  }
  if (sparse_mode == "lasso") {
    sparse_lambda * sum(abs(Z))
  } else {
    sparse_lambda * sum(sqrt(colSums(Z^2)))
  }
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
  obj + .awpt_penalty_term(Z, sparse_lambda, sparse_mode)
}

# Shared FISTA scaffold for the sparse AWPT solvers. The caller supplies the
# starting point, step, gradient closure grad_fn(Y), objective closure
# obj_fn(Z), and prox closure prox_fn(W, step). The iteration (momentum
# update, divergence guard, relative-objective stopping rule) is identical
# across the dense and matrix-free solvers.
.fista <- function(Z0, step, max_iter, tol, grad_fn, obj_fn, prox_fn,
                   diverge_msg = "FISTA solver diverged.") {
  Z <- Z0
  Y <- Z
  t_k <- 1
  prev_obj <- obj_fn(Z)

  for (iter in seq_len(as.integer(max_iter))) {
    grad <- grad_fn(Y)
    Z_next <- prox_fn(Y - step * grad, step)
    t_next <- 0.5 * (1 + sqrt(1 + 4 * t_k^2))
    Y <- Z_next + ((t_k - 1) / t_next) * (Z_next - Z)

    obj_next <- obj_fn(Z_next)
    if (!is.finite(obj_next)) {
      .encoder_cli_abort(diverge_msg,
                         class = "fmrilatent_error_value", call = rlang::caller_env())
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

.frobenius_inner <- function(x, y) {
  sum(x * y)
}

.apply_transport_forward <- function(map, data, context = "transport forward") {
  prep <- .transport_vector_or_matrix_input(data, map$n_source, context)
  out <- map$forward(prep$data)
  out <- as.matrix(out)
  if (nrow(out) != map$n_target) {
    .encoder_cli_abort(
      paste0(context, " returned ", nrow(out), " rows; expected ", map$n_target, "."),
      class = "fmrilatent_error_dim", call = rlang::caller_env())
  }
  out
}

.apply_transport_adjoint <- function(map, data, context = "transport adjoint") {
  prep <- .transport_vector_or_matrix_input(data, map$n_target, context)
  out <- map$adjoint_apply(prep$data)
  out <- as.matrix(out)
  if (nrow(out) != map$n_source) {
    .encoder_cli_abort(
      paste0(context, " returned ", nrow(out), " rows; expected ", map$n_source, "."),
      class = "fmrilatent_error_dim", call = rlang::caller_env())
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
  obj + .awpt_penalty_term(Z, sparse_lambda, sparse_mode)
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
    lambda_est <- norm_w
    V <- W / norm_w
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
      .encoder_cli_abort(
        "Matrix-free quadratic solver encountered a non-positive curvature direction.",
        class = "fmrilatent_error_value", call = rlang::caller_env())
    }
    alpha <- rr / denom
    Z <- Z + alpha * P
    R <- R - alpha * AP
    rr_new <- .frobenius_inner(R, R)
    if (!is.finite(rr_new)) {
      .encoder_cli_abort("Matrix-free quadratic solver diverged.",
                         class = "fmrilatent_error_value", call = rlang::caller_env())
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
  objective_penalty <- spatial_lambda * Q
  Lt <- if (temporal_lambda > 0) {
    .block_temporal_penalty(n_time, temporal_order = temporal_order, run_info = run_info)
  } else {
    NULL
  }
  rhs <- X %*% D_mat
  lipschitz <- .largest_symmetric_eigenvalue(A)
  if (temporal_lambda > 0) {
    lipschitz <- lipschitz + temporal_lambda * .largest_symmetric_eigenvalue(Lt)
  }
  if (!is.finite(lipschitz) || lipschitz <= 0) {
    lipschitz <- 1
  }
  step <- 1 / lipschitz

  grad_fn <- function(Y) {
    grad <- Y %*% A - rhs
    if (temporal_lambda > 0) {
      grad <- grad + temporal_lambda * Lt %*% Y
    }
    grad
  }
  obj_fn <- function(Z) {
    .awpt_objective(Z, X, D_mat, objective_penalty, Lt, temporal_lambda,
                    sparse_lambda, sparse_mode = sparse_mode)
  }
  prox_fn <- function(W, step) {
    .prox_sparse_awpt(
      W,
      step = step,
      sparse_lambda = sparse_lambda,
      sparse_mode = sparse_mode
    )
  }

  .fista(
    Z0 = matrix(0, nrow = n_time, ncol = k),
    step = step,
    max_iter = max_iter,
    tol = tol,
    grad_fn = grad_fn,
    obj_fn = obj_fn,
    prox_fn = prox_fn,
    diverge_msg = "Sparse AWPT solver diverged."
  )
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

  grad_fn <- function(Y) {
    .transport_apply_quadratic_system(
      Y,
      decoder_map = decoder_map,
      spatial_lambda = spatial_lambda,
      spatial_penalty = Q,
      temporal_lambda = temporal_lambda,
      Lt = Lt
    ) - rhs
  }
  obj_fn <- function(Z) {
    .awpt_objective_matrix_free(
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
  }
  prox_fn <- function(W, step) {
    .prox_sparse_awpt(
      W,
      step = step,
      sparse_lambda = sparse_lambda,
      sparse_mode = sparse_mode
    )
  }

  .fista(
    Z0 = matrix(0, nrow = n_time, ncol = k),
    step = step,
    max_iter = max_iter,
    tol = tol,
    grad_fn = grad_fn,
    obj_fn = obj_fn,
    prox_fn = prox_fn,
    diverge_msg = "Sparse matrix-free AWPT solver diverged."
  )
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
    .encoder_cli_abort(
      paste0("Temporal-penalty + materialized-basis path would create a ",
             n_time * k, "x", n_time * k,
             " dense matrix. Use the matrix-free solver instead."),
      class = "fmrilatent_error_unsupported_operation", call = rlang::caller_env())
  }
  system_mat <- kronecker(diag(k), temporal_lambda * Lt) + kronecker(t(A), diag(n_time))
  rhs_vec <- as.vector(rhs)
  z_vec <- .robust_gram_solve(system_mat, matrix(rhs_vec, ncol = 1L),
                              ridge = max(1e-8, spatial_lambda + temporal_lambda))
  matrix(z_vec, nrow = n_time, ncol = k)
}
