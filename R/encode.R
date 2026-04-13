# High-level encoding API and specs

# Robust Gram matrix solve with ridge regularization
# Adds a small ridge to the diagonal if solve fails
.robust_gram_solve <- function(gram, rhs, ridge = 1e-8) {
  gram_mat <- as.matrix(gram)
  result <- tryCatch(
    solve(gram_mat, rhs),
    error = function(e) NULL
  )
  if (is.null(result)) {
    # Add ridge regularization and retry
    diag(gram_mat) <- diag(gram_mat) + ridge
    result <- tryCatch(
      solve(gram_mat, rhs),
      error = function(e) {
        # Final fallback: pseudoinverse via SVD
        svd_g <- svd(gram_mat)
        tol <- max(dim(gram_mat)) * max(svd_g$d) * .Machine$double.eps
        pos <- svd_g$d > tol
        if (any(pos)) {
          svd_g$v[, pos, drop = FALSE] %*% (t(svd_g$u[, pos, drop = FALSE]) / svd_g$d[pos]) %*% rhs
        } else {
          matrix(0, nrow = ncol(gram_mat), ncol = ncol(rhs))
        }
      }
    )
  }
  result
}

.difference_penalty_matrix <- function(n, order = 1L) {
  n <- as.integer(n)
  order <- as.integer(order)
  if (order < 1L || n <= 1L) {
    return(matrix(0, nrow = n, ncol = n))
  }
  D <- diag(n)
  for (idx in seq_len(order)) {
    D <- D[-1L, , drop = FALSE] - D[-nrow(D), , drop = FALSE]
  }
  crossprod(D)
}

.template_custom_field <- function(template, field) {
  if (isS4(template) || !is.list(template)) {
    return(NULL)
  }
  template[[field]] %||% NULL
}

.template_measure_resolve <- function(measure, n_support,
                                      default = c("unit", "null"),
                                      context = "template_measure") {
  default <- match.arg(default)
  if (is.null(measure)) {
    if (default == "unit") {
      return(rep(1, n_support))
    }
    return(NULL)
  }

  if (is.atomic(measure) && is.null(dim(measure))) {
    measure <- as.numeric(measure)
    if (length(measure) == 1L) {
      return(rep(measure, n_support))
    }
    if (length(measure) != n_support) {
      stop(context, " vector must have length ", n_support, ".", call. = FALSE)
    }
    return(measure)
  }

  measure <- as.matrix(measure)
  if (!identical(dim(measure), c(n_support, n_support))) {
    stop(context, " matrix must have dimensions ",
         n_support, "x", n_support, ".", call. = FALSE)
  }
  measure
}

.template_weighted_crossprod <- function(loadings, measure = NULL) {
  loadings <- as.matrix(loadings)
  if (is.null(measure)) {
    return(crossprod(loadings))
  }
  if (is.atomic(measure) && is.null(dim(measure))) {
    return(crossprod(loadings, loadings * as.numeric(measure)))
  }
  crossprod(loadings, as.matrix(measure) %*% loadings)
}

.template_weighted_right_projection <- function(data, loadings, measure = NULL) {
  data <- as.matrix(data)
  loadings <- as.matrix(loadings)
  if (is.null(measure)) {
    return(data %*% loadings)
  }
  if (is.atomic(measure) && is.null(dim(measure))) {
    return((data * rep(as.numeric(measure), each = nrow(data))) %*% loadings)
  }
  data %*% as.matrix(measure) %*% loadings
}

.symmetric_matrix_factor <- function(mat, tol = 1e-10, context = "matrix factor") {
  mat <- 0.5 * (as.matrix(mat) + t(as.matrix(mat)))
  if (nrow(mat) != ncol(mat)) {
    stop(context, " requires a square matrix.", call. = FALSE)
  }

  if (isTRUE(all.equal(mat, diag(nrow(mat)), tolerance = tol))) {
    return(diag(nrow(mat)))
  }

  chol_try <- try(chol(mat), silent = TRUE)
  if (!inherits(chol_try, "try-error")) {
    return(chol_try)
  }

  eig <- eigen(mat, symmetric = TRUE)
  vals <- pmax(Re(eig$values), tol)
  eig$vectors %*% diag(sqrt(vals), nrow = length(vals))
}

.analysis_transform_from_metric <- function(raw_metric, tol = 1e-10) {
  raw_metric <- 0.5 * (as.matrix(raw_metric) + t(as.matrix(raw_metric)))
  k <- nrow(raw_metric)
  if (!identical(dim(raw_metric), c(k, k))) {
    stop("raw_metric must be square.", call. = FALSE)
  }
  if (isTRUE(all.equal(raw_metric, diag(k), tolerance = tol))) {
    return(.transport_identity_transform(k))
  }

  factor <- .symmetric_matrix_factor(raw_metric, tol = tol,
                                     context = "analysis transform")
  list(
    type = "metric_factor",
    dim = as.integer(k),
    raw_metric = raw_metric,
    matrix = factor,
    to_analysis = function(data) factor %*% as.matrix(data),
    to_raw = function(data) solve(factor, as.matrix(data))
  )
}

.analysis_loadings_from_transform <- function(raw_loadings, transform) {
  raw_loadings <- as.matrix(raw_loadings)
  mat <- transform$matrix %||% diag(ncol(raw_loadings))
  raw_loadings %*% solve(as.matrix(mat))
}

.template_coordinate_payload <- function(raw_loadings, measure = NULL,
                                         analysis_transform = NULL,
                                         default_measure = c("unit", "null"),
                                         tol = 1e-10) {
  raw_loadings <- as.matrix(raw_loadings)
  measure_use <- .template_measure_resolve(
    measure,
    n_support = nrow(raw_loadings),
    default = match.arg(default_measure),
    context = "template measure"
  )

  transform <- analysis_transform
  if (is.null(transform)) {
    raw_metric <- .template_weighted_crossprod(raw_loadings, measure_use)
    transform <- .analysis_transform_from_metric(raw_metric, tol = tol)
  } else {
    raw_metric <- .transport_raw_metric(transform, ncol(raw_loadings))
  }

  list(
    raw_loadings = raw_loadings,
    analysis_loadings = .analysis_loadings_from_transform(raw_loadings, transform),
    measure = measure_use,
    analysis_transform = transform,
    raw_metric = raw_metric
  )
}

.template_projection_payload <- function(data, raw_loadings, measure = NULL,
                                         center = FALSE, offset = NULL,
                                         analysis_transform = NULL,
                                         default_measure = c("unit", "null"),
                                         tol = 1e-10) {
  X <- as.matrix(data)
  offset_out <- numeric(0)
  X_proj <- X
  if (isTRUE(center)) {
    offset_out <- if (!is.null(offset)) {
      as.numeric(offset)
    } else {
      colMeans(X)
    }
    X_proj <- sweep(X, 2L, offset_out, "-")
  }

  payload <- .template_coordinate_payload(
    raw_loadings = raw_loadings,
    measure = measure,
    analysis_transform = analysis_transform,
    default_measure = default_measure,
    tol = tol
  )

  coeff_analysis <- .template_weighted_right_projection(
    X_proj,
    payload$analysis_loadings,
    payload$measure
  )
  coeff_raw <- t(payload$analysis_transform$to_raw(t(coeff_analysis)))

  list(
    coefficients = Matrix::Matrix(coeff_analysis, sparse = FALSE),
    raw_coefficients = Matrix::Matrix(coeff_raw, sparse = FALSE),
    offset = offset_out,
    analysis_transform = payload$analysis_transform,
    analysis_loadings = Matrix::Matrix(payload$analysis_loadings, sparse = FALSE),
    raw_metric = payload$raw_metric,
    measure = payload$measure
  )
}

.template_default_measure_mode <- function(template) {
  if (inherits(template, "SurfaceBasisTemplate") ||
      inherits(template, "SurfaceAWPTBasisTemplate")) {
    "null"
  } else {
    "unit"
  }
}

.template_asset_analysis_transform <- function(template) {
  transform <- .template_custom_field(template, "analysis_transform") %||%
    tryCatch(template_meta(template)$analysis_transform %||% NULL,
             error = function(e) NULL)
  if (!is.null(transform)) {
    return(transform)
  }

  raw_loadings <- tryCatch(template_loadings(template), error = function(e) NULL)
  if (is.null(raw_loadings)) {
    return(NULL)
  }

  measure <- tryCatch(template_measure(template), error = function(e) NULL)
  .template_coordinate_payload(
    raw_loadings = raw_loadings,
    measure = measure,
    default_measure = .template_default_measure_mode(template)
  )$analysis_transform
}

.transform_quadratic_form <- function(mat, transform, coordinates = c("analysis", "raw")) {
  coordinates <- match.arg(coordinates)
  if (is.null(mat)) {
    return(NULL)
  }
  mat <- as.matrix(mat)
  if (coordinates == "raw") {
    return(mat)
  }

  factor <- transform$matrix %||% diag(nrow(mat))
  solve(t(as.matrix(factor)), mat) %*% solve(as.matrix(factor))
}

.run_lengths_from_info <- function(n_time, run_info = NULL) {
  if (is.null(run_info)) {
    return(as.integer(n_time))
  }
  if (is.numeric(run_info) && is.null(dim(run_info))) {
    vals <- as.integer(run_info)
  } else if (is.list(run_info) && !is.null(run_info$run_lengths)) {
    vals <- as.integer(run_info$run_lengths)
  } else {
    vals <- as.integer(n_time)
  }
  if (sum(vals) != n_time) {
    stop("run lengths must sum to the number of time points.", call. = FALSE)
  }
  vals
}

.block_temporal_penalty <- function(n_time, temporal_order = 1L, run_info = NULL) {
  run_lengths <- .run_lengths_from_info(n_time, run_info)
  blocks <- lapply(run_lengths, function(len) .difference_penalty_matrix(len, order = temporal_order))
  if (length(blocks) == 1L) {
    return(blocks[[1L]])
  }
  as.matrix(Matrix::bdiag(blocks))
}

.normalize_penalty_matrix <- function(penalty, k, context = "penalty") {
  if (is.null(penalty)) {
    return(matrix(0, nrow = k, ncol = k))
  }
  if (is.vector(penalty) && is.null(dim(penalty))) {
    if (length(penalty) != k) {
      stop(context, " vector must have length ", k, ".", call. = FALSE)
    }
    return(diag(as.numeric(penalty), nrow = k))
  }
  penalty <- as.matrix(penalty)
  if (!identical(dim(penalty), c(k, k))) {
    stop(context, " must have dimensions ", k, "x", k, ".", call. = FALSE)
  }
  penalty
}

.largest_symmetric_eigenvalue <- function(mat) {
  mat <- as.matrix(mat)
  if (nrow(mat) == 0L) return(0)
  vals <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
  max(Re(vals), 0)
}

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

# --- Spec constructors --------------------------------------------------------

.validate_positive_count <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 1 || !isTRUE(all.equal(x, round(x)))) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  as.integer(round(x))
}

.validate_nonnegative_count <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0 || !isTRUE(all.equal(x, round(x)))) {
    stop(name, " must be a non-negative integer.", call. = FALSE)
  }
  as.integer(round(x))
}

.validate_positive_scalar <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    stop(name, " must be a positive finite number.", call. = FALSE)
  }
  as.numeric(x)
}

.validate_nonnegative_scalar <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
    stop(name, " must be a non-negative finite number.", call. = FALSE)
  }
  as.numeric(x)
}

.validate_flag_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  isTRUE(x)
}

.validate_hrbf_params <- function(params) {
  if (!is.list(params)) {
    stop("params must be a list.", call. = FALSE)
  }
  params_clean <- params
  if (!is.null(params_clean$sigma0)) {
    params_clean$sigma0 <- .validate_positive_scalar(params_clean$sigma0, "params$sigma0")
  }
  if (!is.null(params_clean$levels)) {
    params_clean$levels <- .validate_nonnegative_count(params_clean$levels, "params$levels")
  }
  if (!is.null(params_clean$radius_factor)) {
    params_clean$radius_factor <- .validate_positive_scalar(params_clean$radius_factor, "params$radius_factor")
  }
  if (!is.null(params_clean$num_extra_fine_levels)) {
    params_clean$num_extra_fine_levels <- .validate_nonnegative_count(
      params_clean$num_extra_fine_levels,
      "params$num_extra_fine_levels"
    )
  }
  if (!is.null(params_clean$seed)) {
    params_clean$seed <- .validate_nonnegative_count(params_clean$seed, "params$seed")
  }
  if (!is.null(params_clean$kernel_type) &&
      !params_clean$kernel_type %in% c("gaussian", "wendland_c4", "wendland_c6")) {
    stop("params$kernel_type must be one of: gaussian, wendland_c4, wendland_c6.", call. = FALSE)
  }
  if (!is.null(params_clean$kernel_type_fine_levels) &&
      !params_clean$kernel_type_fine_levels %in% c("gaussian", "wendland_c4", "wendland_c6")) {
    stop(
      "params$kernel_type_fine_levels must be one of: gaussian, wendland_c4, wendland_c6.",
      call. = FALSE
    )
  }
  params_clean
}

#' Temporal Slepian/DPSS spec
#'
#' @param tr Repetition time (seconds).
#' @param bandwidth Half-bandwidth in Hz (default 0.1).
#' @param k Optional number of components (default floor(2*NW)-1).
#' @param backend Backend to use ("tridiag" or "dense").
#' @return A `spec_time_slepian` object for `encode()` / `spec_st()`.
#' @export
spec_time_slepian <- function(tr, bandwidth = 0.1, k = NULL, backend = c("tridiag", "dense")) {
  backend <- match.arg(backend)
  tr <- .validate_positive_scalar(tr, "tr")
  bandwidth <- .validate_positive_scalar(bandwidth, "bandwidth")
  if (!is.null(k)) {
    k <- .validate_positive_count(k, "k")
  }
  structure(list(tr = tr, bandwidth = bandwidth, k = k, backend = backend), class = "spec_time_slepian")
}

#' Temporal DCT spec
#'
#' @param k Components.
#' @param norm Normalization ("ortho" or "none").
#' @return A `spec_time_dct` object.
#' @export
spec_time_dct <- function(k, norm = c("ortho", "none")) {
  norm <- match.arg(norm)
  k <- .validate_positive_count(k, "k")
  structure(list(k = k, norm = norm), class = "spec_time_dct")
}

#' Temporal B-spline spec
#'
#' @param k Components (df).
#' @param degree Spline degree (default 3).
#' @param include_intercept Logical include intercept.
#' @param orthonormalize Logical orthonormalize columns (default TRUE).
#' @return A `spec_time_bspline` object.
#' @export
spec_time_bspline <- function(k, degree = 3L, include_intercept = FALSE, orthonormalize = TRUE) {
  k <- .validate_positive_count(k, "k")
  degree <- .validate_positive_count(degree, "degree")
  include_intercept <- .validate_flag_scalar(include_intercept, "include_intercept")
  orthonormalize <- .validate_flag_scalar(orthonormalize, "orthonormalize")
  structure(list(k = k, degree = degree, include_intercept = include_intercept,
                 orthonormalize = orthonormalize),
            class = "spec_time_bspline")
}

#' Spatial Slepian spec
#'
#' @param k Components per cluster.
#' @param k_neighbors k-NN graph parameter.
#' @return A `spec_space_slepian` object.
#' @export
spec_space_slepian <- function(k = 3L, k_neighbors = 6L) {
  k <- .validate_positive_count(k, "k")
  k_neighbors <- .validate_positive_count(k_neighbors, "k_neighbors")
  structure(list(k = k, k_neighbors = k_neighbors), class = "spec_space_slepian")
}

#' Spatial PCA spec (cluster-local)
#'
#' Computes PCA eigenvectors within each cluster/parcel specified by a
#' `ClusterReduction` and returns a block-sparse spatial dictionary.
#'
#' @param k Components per cluster.
#' @param center Logical; center voxels before PCA (default TRUE). When TRUE,
#'   voxel means are stored in `LatentNeuroVec@offset`.
#' @param whiten Logical; if TRUE, return whitened scores (unit-variance) and
#'   rescaled loadings such that reconstruction is unchanged.
#' @param backend SVD backend: "auto" (default), "svds" (RSpectra), or "svd" (base).
#' @return A `spec_space_pca` object.
#' @export
spec_space_pca <- function(k = 3L, center = TRUE, whiten = FALSE,
                           backend = c("auto", "svds", "svd")) {
  backend <- match.arg(backend)
  k <- .validate_positive_count(k, "k")
  center <- .validate_flag_scalar(center, "center")
  whiten <- .validate_flag_scalar(whiten, "whiten")
  structure(
    list(
      k = k,
      center = center,
      whiten = whiten,
      backend = backend
    ),
    class = "spec_space_pca"
  )
}

#' Spatial heat-wavelet spec (graph diffusion)
#'
#' @param scales Heat scales.
#' @param order Polynomial order.
#' @param threshold Threshold for small coefficients.
#' @param k_neighbors k-NN graph parameter.
#' @return A `spec_space_heat` object.
#' @export
spec_space_heat <- function(scales = c(1, 2, 4, 8), order = 30L, threshold = 1e-6, k_neighbors = 6L) {
  if (!is.numeric(scales) || length(scales) == 0L || any(!is.finite(scales)) || any(scales <= 0)) {
    stop("scales must be a non-empty numeric vector of positive finite values.", call. = FALSE)
  }
  order <- .validate_positive_count(order, "order")
  threshold <- .validate_nonnegative_scalar(threshold, "threshold")
  k_neighbors <- .validate_positive_count(k_neighbors, "k_neighbors")
  structure(list(scales = as.numeric(scales), order = order, threshold = threshold, k_neighbors = k_neighbors),
            class = "spec_space_heat")
}

#' Spatial HRBF spec
#'
#' @param params HRBF parameter list (sigma0, levels, radius_factor, kernel_type, seed).
#' @return A `spec_space_hrbf` object.
#' @export
spec_space_hrbf <- function(params = list()) {
  structure(list(params = .validate_hrbf_params(params)), class = "spec_space_hrbf")
}

#' Spatial wavelet (active pencil) spec
#'
#' @param levels_space Spatial lifting levels.
#' @param levels_time Optional time lifting levels.
#' @param threshold Threshold after space transform.
#' @return A `spec_space_wavelet_active` object.
#' @export
spec_space_wavelet_active <- function(levels_space = 2L, levels_time = 0L, threshold = 0) {
  levels_space <- .validate_positive_count(levels_space, "levels_space")
  levels_time <- .validate_nonnegative_count(levels_time, "levels_time")
  threshold <- .validate_nonnegative_scalar(threshold, "threshold")
  structure(list(levels_space = levels_space, levels_time = levels_time, threshold = threshold),
            class = "spec_space_wavelet_active")
}
#' Spatiotemporal spec (separable)
#'
#' @param time Temporal spec (`spec_time_*`).
#' @param space Spatial spec (`spec_space_*`).
#' @param core_mode Reserved (currently "auto").
#' @return A `spec_st` object.
#' @export
spec_st <- function(time, space, core_mode = c("auto", "explicit")) {
  core_mode <- match.arg(core_mode)
  structure(list(time = time, space = space, core_mode = core_mode), class = "spec_st")
}

# --- Encode generic -----------------------------------------------------------

#' Encode data into a latent representation using a spec
#'
#' @param x Matrix (time x voxels in mask order).
#' @param spec Spec object created by `spec_time_*`, `spec_space_*`, or `spec_st`.
#' @param mask LogicalNeuroVol or logical array (required for spatial pieces).
#' @param reduction Optional GraphReduction (for spatial specs).
#' @param materialize "handle", "matrix", or "auto" (default "handle").
#' @param label Optional label.
#' @param ... Additional arguments passed to methods.
#' @return A `LatentNeuroVec` (explicit bases) or `ImplicitLatent` (separable cases).
#' @export
encode <- function(x, spec, mask, reduction = NULL,
                   materialize = c("handle", "auto", "matrix"),
                   label = "", ...) {
  UseMethod("encode")
}

#' @export
encode.default <- function(x, spec, mask, reduction = NULL,
                           materialize = c("handle", "auto", "matrix"),
                           label = "", ...) {
  stop("No encode method for class: ", paste(class(x), collapse = ","), call. = FALSE)
}

#' @export
encode.matrix <- function(x, spec, mask, reduction = NULL,
                          materialize = c("handle", "auto", "matrix"),
                          label = "", ...) {
  materialize <- match.arg(materialize)
  encode_spec(
    x, spec,
    mask = mask,
    reduction = reduction,
    materialize = materialize,
    label = label,
    ...
  )
}

#' @export
encode.NeuroVec <- function(x, spec, mask, reduction = NULL,
                            materialize = c("handle", "auto", "matrix"),
                            label = "", ...) {
  materialize <- match.arg(materialize)
  X <- t(neuroim2::series(x, mask != 0))  # series returns voxels x time, transpose to time x voxels
  encode(X, spec, mask = mask, reduction = reduction,
         materialize = materialize, label = label, ...)
}

# --- Factory helper -----------------------------------------------------------

#' Simple factory to build a spec and encode in one call
#'
#' @param family One of: "dct_time", "slepian_time", "slepian_space", "heat_space", "slepian_st".
#' @param x Data matrix (time x voxels).
#' @param mask Mask (required for spatial families).
#' @param reduction Optional GraphReduction for spatial specs.
#' @param ... Passed to spec constructors and encode().
#' @param materialize "handle", "matrix", or "auto" (default "handle").
#' @param label Optional label for the resulting object.
#' @return A LatentNeuroVec or ImplicitLatent object.
#' @export
latent_factory <- function(family, x, mask, reduction = NULL, ..., materialize = "handle", label = "") {
  if (identical(family, "awpt")) {
    stop("latent_factory() does not support AWPT because AWPT requires a shared ",
         "basis_asset and subject field_operator. Use encode_awpt() or ",
         "encode_operator() instead.", call. = FALSE)
  }
  family <- match.arg(family, c("dct_time", "slepian_time", "slepian_space",
    "pca_space", "parcel_space", "heat_space", "slepian_st",
    "bspline_hrbf_st", "wavelet_active", "hierarchical"))
  spec <- switch(
    family,
    dct_time = spec_time_dct(...),
    slepian_time = spec_time_slepian(...),
    slepian_space = spec_space_slepian(...),
    pca_space = spec_space_pca(...),
    heat_space = spec_space_heat(...),
    bspline_hrbf_st = {
      args <- list(...)
      time_spec <- args$time %||% spec_time_bspline(k = args$k_time %||% args$k %||% 5L,
                                                    degree = args$degree %||% 3L,
                                                    include_intercept = args$include_intercept %||% FALSE,
                                                    orthonormalize = TRUE)
      space_spec <- args$space %||% spec_space_hrbf(params = args$params %||% list())
      spec_st(time = time_spec, space = space_spec)
    },
    parcel_space = spec_space_parcel(...),
    wavelet_active = spec_space_wavelet_active(...),
    hierarchical = spec_hierarchical_template(...),
    slepian_st = {
      args <- list(...)
      time_spec <- args$time %||% spec_time_slepian(
        tr = args$tr, bandwidth = args$bandwidth %||% 0.1,
        k = args$k_time %||% NULL)
      space_spec <- args$space %||% spec_space_slepian(k = args$k_space %||% 3L, k_neighbors = args$k_neighbors %||% 6L)
      spec_st(time = time_spec, space = space_spec)
    }
  )
  encode(x, spec, mask = mask, reduction = reduction, materialize = materialize, label = label)
}
#' Dispatch encoding based on spec type
#'
#' @param x Data matrix.
#' @param spec Spec object.
#' @param ... Additional arguments passed to methods.
#' @return Encoded representation.
#' @export
encode_spec <- function(x, spec, ...) UseMethod("encode_spec", spec)

#' @exportS3Method
encode_spec.default <- function(x, spec, ...) {
  stop("Unknown spec class: ", paste(class(spec), collapse = ","), call. = FALSE)
}

#' @exportS3Method
encode_spec.spec_awpt_wavelet <- function(x, spec, mask, reduction, materialize, label, ...) {
  stop("AWPT encoding uses encode_awpt() or encode_operator(), not the standard ",
       "encode() pipeline. See ?encode_awpt for details.", call. = FALSE)
}

#' Encode data using a shared basis asset and subject field operator
#'
#' @param x Numeric matrix (time x target samples) or a \code{NeuroVec}.
#' @param template Shared basis asset providing \code{basis_decoder()}.
#' @param field_operator Subject field operator mapping template field
#'   samples to observed native samples.
#' @param observation_operator Legacy alias for \code{field_operator}.
#' @details
#' The external field-operator contract is intentionally narrow. `fmrilatent`
#' consumes an operator-like object with:
#' \describe{
#'   \item{\code{n_source}, \code{n_target}}{Source and target dimensions.}
#'   \item{\code{source_domain_id}, \code{target_domain_id}}{Stable domain identifiers.}
#'   \item{\code{forward(z)}}{Applies the operator to template field samples.}
#'   \item{\code{adjoint_apply(y)}}{Applies the adjoint map.}
#'   \item{\code{provenance$target_mask}}{Optional target-domain mask when
#'   the caller does not pass \code{mask} explicitly.}
#' }
#'
#' On the main quadratic and sparse transport/AWPT paths, coefficient recovery
#' is matrix-free: `fmrilatent` uses the operator's forward and adjoint methods
#' rather than materializing the full subject decoder.
#' @param mask Optional volumetric target mask for the field-operator
#'   target domain.
#' @param domain Optional non-volumetric target domain, for example a surface
#'   geometry.
#' @param support Optional target support aligned to \code{domain}, for example
#'   vertex indices on a surface.
#'   At least one of \code{mask} or \code{support} must be available either
#'   explicitly or via \code{field_operator$provenance}.
#' @param lambda Ridge penalty strength.
#' @param spatial_lambda Strength of the spatial coefficient penalty.
#' @param spatial_penalty Optional coefficient-space roughness matrix or diagonal weights.
#' @param temporal_lambda Strength of temporal smoothing.
#' @param temporal_order Difference order used for temporal smoothing.
#' @param sparse_lambda Strength of optional sparse coefficient shrinkage.
#' @param sparse_mode Sparse penalty mode. Use \code{"group_l2"} for atom-wise group shrinkage.
#' @param max_iter Maximum iterations for sparse AWPT optimization.
#' @param tol Relative convergence tolerance for sparse AWPT optimization.
#' @param center Logical; if \code{TRUE}, center target samples before solving.
#' @param run_info Optional run metadata carried on the resulting latent object.
#' @param label Optional label stored in metadata.
#' @param ... Reserved for future extensions.
#' @return A \code{TransportLatent} object.
#' @export
encode_operator <- function(x, template, field_operator = NULL, observation_operator = NULL, mask = NULL,
                            domain = NULL, support = NULL,
                            lambda = 0, center = TRUE, run_info = NULL,
                            spatial_lambda = lambda,
                            spatial_penalty = NULL,
                            temporal_lambda = 0,
                            temporal_order = 1L,
                            sparse_lambda = 0,
                            sparse_mode = c("none", "group_l2", "lasso"),
                            max_iter = 200L,
                            tol = 1e-6,
                            label = "", ...) {
  sparse_mode <- match.arg(sparse_mode)
  field_operator <- .resolve_field_operator(
    field_operator = field_operator,
    observation_operator = observation_operator,
    context = "encode_operator() field_operator"
  )
  obs_map <- .normalize_field_operator_map(field_operator)
  target_info <- .resolve_transport_target_support(
    mask = mask,
    domain = domain,
    support = support,
    field_operator = obs_map,
    location = "encode_operator"
  )
  mask_use <- target_info$mask
  support_use <- target_info$support
  domain_use <- target_info$domain
  if (inherits(x, "NeuroVec")) {
    if (is.null(mask_use)) {
      stop("NeuroVec inputs currently require a volumetric target mask for encode_operator().",
           call. = FALSE)
    }
    X <- t(neuroim2::series(x, mask_use != 0))
  } else {
    X <- as.matrix(x)
  }

  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("lambda must be a single non-negative finite number.", call. = FALSE)
  }
  if (!is.numeric(spatial_lambda) || length(spatial_lambda) != 1L ||
      !is.finite(spatial_lambda) || spatial_lambda < 0) {
    stop("spatial_lambda must be a single non-negative finite number.", call. = FALSE)
  }
  if (!is.numeric(temporal_lambda) || length(temporal_lambda) != 1L ||
      !is.finite(temporal_lambda) || temporal_lambda < 0) {
    stop("temporal_lambda must be a single non-negative finite number.", call. = FALSE)
  }
  if (!is.numeric(sparse_lambda) || length(sparse_lambda) != 1L ||
      !is.finite(sparse_lambda) || sparse_lambda < 0) {
    stop("sparse_lambda must be a single non-negative finite number.", call. = FALSE)
  }

  basis_map <- .normalize_basis_decoder_map(template)
  D_map <- .compose_linear_maps(basis_map, obs_map, context = "subject decoder")
  analysis_transform_use <- .template_asset_analysis_transform(template) %||%
    .transport_identity_transform(D_map$n_source)

  if (ncol(X) != D_map$n_target) {
    stop("x has ", ncol(X), " target samples but the composed decoder has ",
         D_map$n_target, ".", call. = FALSE)
  }

  offset <- numeric(0)
  X_proj <- X
  if (isTRUE(center)) {
    offset <- colMeans(X)
    X_proj <- sweep(X, 2L, offset, "-")
  }

  penalty_use <- spatial_penalty %||% template_roughness(template, coordinates = "analysis")
  if (is.null(penalty_use) && spatial_lambda > 0) {
    penalty_use <- diag(D_map$n_source)
  }
  coeff_analysis <- if (sparse_lambda > 0 && sparse_mode != "none") {
    .solve_transport_coefficients_sparse_matrix_free(
      X = X_proj,
      decoder_map = D_map,
      spatial_lambda = spatial_lambda,
      spatial_penalty = penalty_use,
      temporal_lambda = temporal_lambda,
      temporal_order = temporal_order,
      run_info = run_info,
      sparse_lambda = sparse_lambda,
      sparse_mode = sparse_mode,
      max_iter = max_iter,
      tol = tol
    )
  } else {
    .solve_transport_coefficients_matrix_free(
      X = X_proj,
      decoder_map = D_map,
      spatial_lambda = spatial_lambda,
      spatial_penalty = penalty_use,
      temporal_lambda = temporal_lambda,
      temporal_order = temporal_order,
      run_info = run_info,
      max_iter = max_iter,
      tol = tol
    )
  }
  coeff_raw <- t(analysis_transform_use$to_raw(t(coeff_analysis)))

  transport_latent(
    coeff_raw = coeff_raw,
    coeff_analysis = coeff_analysis,
    basis_asset = template,
    field_operator = field_operator,
    mask = mask_use,
    domain = domain_use,
    support = support_use,
    analysis_transform = analysis_transform_use,
    offset = offset,
    run_info = run_info,
    meta = list(
      family = "transport",
      label = label,
      lambda = lambda,
      spatial_lambda = spatial_lambda,
      temporal_lambda = temporal_lambda,
      temporal_order = temporal_order,
      sparse_lambda = sparse_lambda,
      sparse_mode = sparse_mode,
      method = if (is_awpt_template(template)) "AWPT" else NULL,
      centered = isTRUE(center),
      target_mask_source = target_info$source
    )
  )
}

#' Encode data using transport-backed latent semantics
#'
#' @param x Numeric matrix (time x target samples) or a \code{NeuroVec}.
#' @param basis_asset Shared basis asset.
#' @param field_operator Subject field operator. See
#'   \code{\link{encode_operator}()} for the required contract.
#' @param observation_operator Legacy alias for \code{field_operator}.
#' @param mask Optional volumetric target mask for the field-operator target domain.
#' @param domain Optional non-volumetric target domain.
#' @param support Optional target support aligned to \code{domain}.
#' @param lambda Ridge penalty strength.
#' @param spatial_lambda Strength of the spatial coefficient penalty.
#' @param spatial_penalty Optional coefficient-space roughness matrix or diagonal weights.
#' @param temporal_lambda Strength of temporal smoothing.
#' @param temporal_order Difference order used for temporal smoothing.
#' @param sparse_lambda Strength of optional sparse coefficient shrinkage.
#' @param sparse_mode Sparse penalty mode. Use \code{"group_l2"} for atom-wise group shrinkage.
#' @param max_iter Maximum iterations for sparse AWPT optimization.
#' @param tol Relative convergence tolerance for sparse AWPT optimization.
#' @param center Logical; if \code{TRUE}, center target samples before solving.
#' @param run_info Optional run metadata carried on the resulting latent object.
#' @param label Optional label stored in metadata.
#' @param ... Reserved for future extensions.
#' @return A \code{TransportLatent} object.
#' @export
encode_transport <- function(x, basis_asset, field_operator = NULL, observation_operator = NULL, mask = NULL,
                             domain = NULL, support = NULL,
                             lambda = 0, center = TRUE, run_info = NULL,
                             spatial_lambda = lambda,
                             spatial_penalty = NULL,
                             temporal_lambda = 0,
                             temporal_order = 1L,
                             sparse_lambda = 0,
                             sparse_mode = c("none", "group_l2", "lasso"),
                             max_iter = 200L,
                             tol = 1e-6,
                             label = "", ...) {
  encode_operator(
    x = x,
    template = basis_asset,
    field_operator = field_operator,
    observation_operator = observation_operator,
    mask = mask,
    domain = domain,
    support = support,
    lambda = lambda,
    spatial_lambda = spatial_lambda,
    spatial_penalty = spatial_penalty,
    temporal_lambda = temporal_lambda,
    temporal_order = temporal_order,
    sparse_lambda = sparse_lambda,
    sparse_mode = sparse_mode,
    max_iter = max_iter,
    tol = tol,
    center = center,
    run_info = run_info,
    label = label,
    ...
  )
}

#' Encode data using the AWPT model
#'
#' @param x Numeric matrix (time x target samples) or a \code{NeuroVec}.
#' @param basis_asset An \code{AWPTBasisTemplate}.
#' @param field_operator Subject field operator. See
#'   \code{\link{encode_operator}()} for the required contract.
#' @param observation_operator Legacy alias for \code{field_operator}.
#' @param mask Optional volumetric target mask for the field-operator target domain.
#' @param domain Optional non-volumetric target domain.
#' @param support Optional target support aligned to \code{domain}.
#' @param spatial_lambda Strength of the anatomical roughness penalty.
#' @param temporal_lambda Strength of temporal smoothing.
#' @param temporal_order Temporal difference order for smoothing.
#' @param sparse_lambda Strength of optional sparse atom selection.
#' @param sparse_mode Sparse penalty mode. Use \code{"group_l2"} for atom-wise group shrinkage.
#' @param max_iter Maximum iterations for sparse AWPT optimization.
#' @param tol Relative convergence tolerance for sparse AWPT optimization.
#' @param center Logical; if \code{TRUE}, center target samples before solving.
#' @param run_info Optional run metadata; \code{run_lengths} control temporal blocks.
#' @param label Optional label stored in metadata.
#' @param ... Reserved for future extensions.
#' @return A \code{TransportLatent} object with AWPT metadata.
#' @export
encode_awpt <- function(x, basis_asset, field_operator = NULL, observation_operator = NULL, mask = NULL,
                        domain = NULL, support = NULL,
                        spatial_lambda = 0, temporal_lambda = 0,
                        temporal_order = 1L,
                        sparse_lambda = 0,
                        sparse_mode = c("none", "group_l2", "lasso"),
                        max_iter = 200L,
                        tol = 1e-6,
                        center = TRUE,
                        run_info = NULL, label = "", ...) {
  sparse_mode <- match.arg(sparse_mode)
  if (!is_awpt_template(basis_asset)) {
    stop("basis_asset must be an AWPT basis template.", call. = FALSE)
  }
  encode_operator(
    x = x,
    template = basis_asset,
    field_operator = field_operator,
    observation_operator = observation_operator,
    mask = mask,
    domain = domain,
    support = support,
    lambda = spatial_lambda,
    spatial_lambda = spatial_lambda,
    spatial_penalty = template_roughness(basis_asset, coordinates = "analysis"),
    temporal_lambda = temporal_lambda,
    temporal_order = temporal_order,
    sparse_lambda = sparse_lambda,
    sparse_mode = sparse_mode,
    max_iter = max_iter,
    tol = tol,
    center = center,
    run_info = run_info,
    label = label,
    ...
  )
}

# Hierarchical template spec -------------------------------------------------

#' Create hierarchical template spec
#' @param template HierarchicalBasisTemplate object
#' @param template_file Path to saved template file
#' @return Spec object of class spec_hierarchical
#' @export
spec_hierarchical_template <- function(template = NULL, template_file = NULL) {
  if (is.null(template) && is.null(template_file)) {
    stop("Provide either a template object or template_file")
  }
  structure(list(template = template, template_file = template_file), class = "spec_hierarchical")
}

#' @exportS3Method
encode_spec.spec_hierarchical <- function(x, spec, mask, materialize, label, ...) {
  tmpl <- spec$template
  if (is.null(tmpl)) {
    if (is.null(spec$template_file)) stop("template_file is required when template is NULL")
    tmpl <- load_hierarchical_template(spec$template_file)
  }
  .assert_template_mask_match(mask, template_mask(tmpl), "encode_spec.spec_hierarchical")
  encode_hierarchical(x, tmpl, label = label, mask = mask)
}

#' @exportS3Method
encode_spec.spec_time_slepian <- function(x, spec, mask, materialize, label, ...) {
  n_time <- nrow(x)
  k_use <- spec$k %||% max(1L, floor(2 * n_time * spec$bandwidth * spec$tr) - 1L)
  if (materialize == "matrix") {
    basis <- dpss_time_basis(n_time, tr = spec$tr, bandwidth = spec$bandwidth, k = k_use, backend = spec$backend)
  } else {
    basis <- slepian_temporal_handle(n_time = n_time, tr = spec$tr,
      bandwidth = spec$bandwidth, k = k_use, backend = spec$backend)
  }
  loadings <- Matrix::Matrix(crossprod(x, basis_mat(basis)), sparse = FALSE)
  spc <- .space_with_time_from_mask(mask, n_time, "encode_spec.spec_time_slepian")
  meta <- list(
    family = "time_slepian",
    k = k_use,
    tr = spec$tr,
    bandwidth = spec$bandwidth,
    backend = spec$backend
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}

#' @exportS3Method
encode_spec.spec_time_dct <- function(x, spec, mask, materialize, label, ...) {
  n_time <- nrow(x)
  k_use <- spec$k
  if (materialize == "matrix") {
    basis <- build_dct_basis(n_time, k = k_use, norm = spec$norm)
  } else {
    basis <- dct_basis_handle(n_time = n_time, k = k_use, norm = spec$norm)
  }
  loadings <- Matrix::Matrix(crossprod(x, basis_mat(basis)), sparse = FALSE)
  spc <- .space_with_time_from_mask(mask, n_time, "encode_spec.spec_time_dct")
  meta <- list(
    family = "time_dct",
    k = k_use,
    norm = spec$norm
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}

#' @exportS3Method
encode_spec.spec_time_bspline <- function(x, spec, mask, materialize, label, ...) {
  n_time <- nrow(x)
  k_use <- spec$k
  if (materialize == "matrix") {
    basis <- build_bspline_basis(
      n_time = n_time,
      k = k_use,
      degree = spec$degree,
      include_intercept = spec$include_intercept,
      orthonormalize = spec$orthonormalize
    )
  } else {
    basis <- bspline_basis_handle(
      n_time = n_time,
      k = k_use,
      degree = spec$degree,
      include_intercept = spec$include_intercept,
      id = NULL,
      label = NULL
    )
    basis@spec$orthonormalize <- spec$orthonormalize
  }
  loadings <- Matrix::Matrix(crossprod(x, basis_mat(basis)), sparse = FALSE)
  spc <- .space_with_time_from_mask(mask, n_time, "encode_spec.spec_time_bspline")
  meta <- list(
    family = "time_bspline",
    k = k_use,
    degree = spec$degree,
    include_intercept = spec$include_intercept,
    orthonormalize = spec$orthonormalize
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}

#' @exportS3Method
encode_spec.spec_space_slepian <- function(x, spec, mask, reduction, materialize, label, ...) {
  mask_arr <- .mask_to_array(mask, "encode_spec.spec_space_slepian")
  if (is.null(reduction)) reduction <- make_cluster_reduction(mask, seq_len(sum(mask_arr)))
  if (materialize == "matrix") {
    loadings <- lift(reduction, basis_slepian(k = spec$k), k_neighbors = spec$k_neighbors)
  } else {
    loadings <- slepian_spatial_loadings_handle(
      reduction,
      basis_spec = basis_slepian(k = spec$k),
      data = NULL,
      k_neighbors = spec$k_neighbors,
      label = "slepian-spatial"
    )
  }
  basis <- as.matrix(x) %*% as.matrix(loadings_mat(loadings))
  spc <- .space_with_time_from_mask(mask, nrow(x), "encode_spec.spec_space_slepian")
  meta <- list(
    family = "space_slepian",
    k = spec$k,
    k_neighbors = spec$k_neighbors
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}

#' @exportS3Method
encode_spec.spec_space_heat <- function(x, spec, mask, reduction, materialize, label, ...) {
  mask_arr <- .mask_to_array(mask, "encode_spec.spec_space_heat")
  if (is.null(reduction)) reduction <- make_cluster_reduction(mask, seq_len(sum(mask_arr)))
  spec_hw <- basis_heat_wavelet(scales = spec$scales, order = spec$order, threshold = spec$threshold)
  if (materialize == "matrix") {
    loadings <- lift(reduction, spec_hw, k_neighbors = spec$k_neighbors)
  } else {
    loadings <- heat_wavelet_loadings_handle(
      reduction,
      basis_spec = spec_hw,
      data = NULL,
      k_neighbors = spec$k_neighbors,
      label = "heat-wavelet"
    )
  }
  basis <- as.matrix(x) %*% as.matrix(loadings_mat(loadings))
  spc <- .space_with_time_from_mask(mask, nrow(x), "encode_spec.spec_space_heat")
  meta <- list(
    family = "space_heat",
    scales = spec$scales,
    order = spec$order,
    threshold = spec$threshold,
    k_neighbors = spec$k_neighbors
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}

#' @exportS3Method
encode_spec.spec_space_hrbf <- function(x, spec, mask, reduction, materialize, label, ...) {
  mask_arr <- .mask_to_array(mask, "encode_spec.spec_space_hrbf")
  params <- spec$params %||% list()
  B_atoms <- hrbf_generate_basis(params, mask) # atoms x vox
  loadings <- Matrix::t(Matrix::Matrix(B_atoms, sparse = TRUE)) # vox x atoms
  gram <- as.matrix(B_atoms %*% Matrix::t(B_atoms))
  rhs <- t(as.matrix(x) %*% Matrix::t(B_atoms))
  coeff <- t(.robust_gram_solve(gram, rhs))
  basis <- Matrix::Matrix(coeff, sparse = FALSE) # time x atoms
  spc <- .space_with_time_from_mask(mask, nrow(x), "encode_spec.spec_space_hrbf")
  meta <- list(
    family = "space_hrbf",
    params = params
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}

#' @exportS3Method
encode_spec.spec_space_pca <- function(x, spec, mask, reduction, materialize, label, ...) {
  mask_arr <- .mask_to_array(mask, "encode_spec.spec_space_pca")
  n_time <- nrow(x)
  n_vox <- sum(mask_arr)

  if (is.null(reduction)) {
    reduction <- make_cluster_reduction(mask, rep.int(1L, n_vox))
  }

  offset <- numeric(0)
  if (isTRUE(spec$center)) {
    offset <- colMeans(x)
  }

  loadings <- lift(
    reduction,
    basis_pca(k = spec$k, whiten = isTRUE(spec$whiten)),
    data = x,
    center = isTRUE(spec$center),
    offset = if (length(offset) > 0) offset else NULL,
    backend = spec$backend %||% "auto",
    ...
  )

  basis <- x %*% loadings
  if (length(offset) > 0) {
    mu_scores <- as.matrix(crossprod(offset, loadings))
    basis <- basis - matrix(1, nrow = n_time, ncol = 1) %*% mu_scores
  }

  if (isTRUE(spec$whiten)) {
    d <- attr(loadings, "fmrilatent.singular_values")
    if (is.null(d) || length(d) != ncol(loadings)) {
      stop("PCA whitening requested, but singular values were not returned by lift().", call. = FALSE)
    }
    if (any(!is.finite(d)) || any(d <= 0)) {
      stop("PCA whitening requires strictly positive finite singular values.", call. = FALSE)
    }
    scale_fac <- sqrt(max(1, n_time - 1))
    basis <- sweep(as.matrix(basis), 2, d, "/") * scale_fac
    loadings <- loadings %*% Matrix::Diagonal(x = d / scale_fac)
  }

  spc <- .space_with_time_from_mask(mask, n_time, "encode_spec.spec_space_pca")
  meta <- list(
    family = "pca_spatial",
    k = spec$k,
    center = isTRUE(spec$center),
    whiten = isTRUE(spec$whiten),
    backend = spec$backend %||% "auto"
  )

  LatentNeuroVec(
    basis = Matrix::Matrix(basis, sparse = FALSE),
    loadings = loadings,
    space = spc,
    mask = mask,
    offset = offset,
    label = label,
    meta = meta
  )
}

#' @exportS3Method
encode_spec.spec_space_wavelet_active <- function(x, spec, mask, reduction, materialize, label, ...) {
  mask_arr <- .mask_to_array(mask, "encode_spec.spec_space_wavelet_active")
  wavelet_active_latent(
    X = x,
    mask = mask,
    levels_space = spec$levels_space,
    levels_time = spec$levels_time,
    threshold = spec$threshold
  )
}

#' @exportS3Method
encode_spec.spec_st <- function(x, spec, mask, reduction, materialize, label, ...) {
  mask_arr <- .mask_to_array(mask, "encode_spec.spec_st")
  n_time <- nrow(x)

  if (inherits(spec$time, "spec_time_slepian")) {
    k_time <- spec$time$k %||% max(1L, floor(2 * n_time * spec$time$bandwidth * spec$time$tr) - 1L)
    B_t <- basis_mat(slepian_temporal_handle(n_time = n_time, tr = spec$time$tr,
                                             bandwidth = spec$time$bandwidth, k = k_time,
                                             backend = spec$time$backend))
  } else if (inherits(spec$time, "spec_time_bspline")) {
    k_time <- spec$time$k
    B_t <- as.matrix(build_bspline_basis(
      n_time = n_time,
      k = k_time,
      degree = spec$time$degree,
      include_intercept = spec$time$include_intercept,
      orthonormalize = spec$time$orthonormalize
    ))
    if (!spec$time$orthonormalize) {
      q <- qr(B_t)
      B_t <- qr.Q(q)
    }
  } else {
    stop("Unsupported spec_st$time class: ", class(spec$time)[1])
  }

  if (is.null(reduction)) reduction <- make_cluster_reduction(mask, seq_len(sum(mask_arr)))
  if (inherits(spec$space, "spec_space_slepian")) {
    L_s <- as.matrix(lift(reduction, basis_slepian(k = spec$space$k), k_neighbors = spec$space$k_neighbors))
    B_atoms <- t(L_s)
  } else if (inherits(spec$space, "spec_space_hrbf")) {
    params <- spec$space$params %||% list()
    B_atoms <- as.matrix(hrbf_generate_basis(params, mask))
    L_s <- t(B_atoms)
  } else if (inherits(spec$space, "spec_space_wavelet_active")) {
    stop("spec_st with space = spec_space_wavelet_active not yet supported; ",
      "use spec_space_wavelet_active directly via encode().", call. = FALSE)
  } else {
    stop("Unsupported spec_st$space class: ", class(spec$space)[1])
  }

  gram <- B_atoms %*% t(B_atoms)
  rhs_st <- t(as.matrix(x) %*% t(B_atoms))
  C_atoms <- t(.robust_gram_solve(as.matrix(gram), rhs_st))

  core <- crossprod(B_t, C_atoms)

  decoder <- function(time_idx = NULL, roi_mask = NULL, ...) {
    t_sel <- if (is.null(time_idx)) seq_len(n_time) else as.integer(time_idx)
    B_sel <- B_t[t_sel, , drop = FALSE]
    rec <- B_sel %*% core %*% t(L_s)
    if (!is.null(roi_mask)) {
      global_idx <- which(as.logical(mask_arr))
      roi_global <- which(as.logical(roi_mask))
      col_keep <- which(global_idx %in% roi_global)
      rec <- rec[, col_keep, drop = FALSE]
    }
    rec
  }

  meta <- list(
    family = "st_separable",
    time = class(spec$time)[1],
    space = class(spec$space)[1],
    k_time = ncol(B_t),
    k_space = ncol(L_s),
    label = label
  )

  implicit_latent(
    coeff = list(core = core, B_t = B_t, L_s = L_s),
    decoder = decoder,
    meta = meta,
    mask = mask_arr
  )
}
