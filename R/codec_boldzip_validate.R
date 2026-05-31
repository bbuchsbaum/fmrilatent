# BOLDZip-SR codec: validators and linear-algebra primitives
#
# Split out of R/codec_boldzip.R (see that file's header for the full map).
# Functions in this file:
#   .boldzip_check_scalar_integer
#   .boldzip_check_scalar_number
#   .boldzip_validate_matrix
#   .boldzip_validate_basis_matrix
#   .boldzip_validate_orthonormal_columns
#   .boldzip_validate_cross_orthogonal
#   .boldzip_validate_spatial_basis_input
#   .boldzip_project
#   .boldzip_synthesize
#   .boldzip_orthonormalize_columns
#   .boldzip_canonicalize_eigenvectors
#   .boldzip_corr

.boldzip_cli_abort <- function(..., class = "fmrilatent_error_boldzip", call = NULL) {
  .encoder_cli_abort(paste0(...), class = class, call = call)
}

.boldzip_check_scalar_integer <- function(x, name, min = 1L) {
  if (length(x) != 1L || is.na(x) || x != as.integer(x) || x < min) {
    .boldzip_cli_abort(name, " must be a scalar integer >= ", min, ".")
  }
  as.integer(x)
}

.boldzip_check_scalar_number <- function(x, name, min = -Inf, max = Inf) {
  if (length(x) != 1L || !is.finite(x) || x < min || x > max) {
    .boldzip_cli_abort(name, " must be a finite scalar in [", min, ", ", max, "].")
  }
  as.numeric(x)
}

.boldzip_validate_matrix <- function(x, name) {
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    .boldzip_cli_abort(name, " must be a numeric matrix.")
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    .boldzip_cli_abort(name, " must contain only finite values.")
  }
  if (nrow(x) < 1L || ncol(x) < 1L) {
    .boldzip_cli_abort(name, " must have at least one row and one column.")
  }
  x
}

.boldzip_validate_basis_matrix <- function(x, name, n_voxels) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    .boldzip_cli_abort(name, " must be a numeric matrix.")
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    .boldzip_cli_abort(name, " must contain only finite values.")
  }
  if (nrow(x) < 1L || ncol(x) < 1L) {
    .boldzip_cli_abort(name, " must have at least one row and one column.")
  }
  if (nrow(x) != n_voxels) {
    .boldzip_cli_abort(name, " must have ", n_voxels, " rows.")
  }
  x
}

.boldzip_validate_orthonormal_columns <- function(x, name, tol = 1e-8) {
  if (is.null(x)) {
    return(invisible(TRUE))
  }
  gram <- crossprod(x)
  if (!isTRUE(all.equal(gram, diag(ncol(x)), tolerance = tol))) {
    .boldzip_cli_abort(name, " columns must be orthonormal for BOLDZip-SR spatial coding.")
  }
  invisible(TRUE)
}

.boldzip_validate_cross_orthogonal <- function(phi_c, phi_d, tol = 1e-8) {
  if (is.null(phi_c) || is.null(phi_d)) {
    return(invisible(TRUE))
  }
  cross <- crossprod(phi_c, phi_d)
  if (!isTRUE(all.equal(cross, matrix(0, nrow = ncol(phi_c), ncol = ncol(phi_d)),
                        tolerance = tol))) {
    .boldzip_cli_abort("phi_c and phi_d must be mutually orthogonal for BOLDZip-SR spatial coding.")
  }
  invisible(TRUE)
}

.boldzip_validate_spatial_basis_input <- function(x, name) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    .boldzip_cli_abort(name, " must be a numeric matrix.")
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    .boldzip_cli_abort(name, " must contain only finite values.")
  }
  if (nrow(x) < 1L || ncol(x) < 1L) {
    .boldzip_cli_abort(name, " must have at least one row and one column.")
  }
  x
}

.boldzip_project <- function(phi, x) {
  if (is.null(phi)) {
    return(x)
  }
  crossprod(phi, x)
}

.boldzip_synthesize <- function(phi, coef) {
  if (is.null(phi)) {
    return(coef)
  }
  phi %*% coef
}

.boldzip_orthonormalize_columns <- function(x, tol = 1e-10) {
  x <- as.matrix(x)
  if (ncol(x) == 0L) {
    return(x)
  }
  q <- qr(x, tol = tol)
  rank <- q$rank
  if (rank < 1L) {
    .boldzip_cli_abort("basis matrix must have at least one non-zero independent column.")
  }
  if (rank < ncol(x)) {
    warning(
      "orthonormalize: dropped ", ncol(x) - rank,
      " linearly dependent column(s); returning ", rank, " columns.",
      call. = FALSE
    )
  }
  qr.Q(q)[, seq_len(rank), drop = FALSE]
}

.boldzip_canonicalize_eigenvectors <- function(vectors) {
  for (idx in seq_len(ncol(vectors))) {
    pivot <- which.max(abs(vectors[, idx]))
    if (length(pivot) && vectors[pivot, idx] < 0) {
      vectors[, idx] <- -vectors[, idx]
    }
  }
  vectors
}

.boldzip_corr <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  min_sd <- sqrt(.Machine$double.eps)
  if (.boldzip_noise_scale(x) <= min_sd || .boldzip_noise_scale(y) <= min_sd) {
    return(NA_real_)
  }
  stats::cor(x, y)
}
