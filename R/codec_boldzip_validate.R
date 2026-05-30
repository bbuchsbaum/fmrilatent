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

.boldzip_check_scalar_integer <- function(x, name, min = 1L) {
  if (length(x) != 1L || is.na(x) || x != as.integer(x) || x < min) {
    stop(name, " must be a scalar integer >= ", min, ".", call. = FALSE)
  }
  as.integer(x)
}

.boldzip_check_scalar_number <- function(x, name, min = -Inf, max = Inf) {
  if (length(x) != 1L || !is.finite(x) || x < min || x > max) {
    stop(name, " must be a finite scalar in [", min, ", ", max, "].",
         call. = FALSE)
  }
  as.numeric(x)
}

.boldzip_validate_matrix <- function(x, name) {
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    stop(name, " must be a numeric matrix.", call. = FALSE)
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  if (nrow(x) < 1L || ncol(x) < 2L) {
    stop(name, " must have at least one row and two columns.", call. = FALSE)
  }
  x
}

.boldzip_validate_basis_matrix <- function(x, name, n_voxels) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    stop(name, " must be a numeric matrix.", call. = FALSE)
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop(name, " must have at least one row and one column.", call. = FALSE)
  }
  if (nrow(x) != n_voxels) {
    stop(name, " must have ", n_voxels, " rows.", call. = FALSE)
  }
  x
}

.boldzip_validate_orthonormal_columns <- function(x, name, tol = 1e-8) {
  if (is.null(x)) {
    return(invisible(TRUE))
  }
  gram <- crossprod(x)
  if (!isTRUE(all.equal(gram, diag(ncol(x)), tolerance = tol))) {
    stop(name, " columns must be orthonormal for BOLDZip-SR spatial coding.",
         call. = FALSE)
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
    stop("phi_c and phi_d must be mutually orthogonal for BOLDZip-SR spatial coding.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.boldzip_validate_spatial_basis_input <- function(x, name) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    stop(name, " must be a numeric matrix.", call. = FALSE)
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop(name, " must have at least one row and one column.", call. = FALSE)
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
    stop("basis matrix must have at least one non-zero independent column.",
         call. = FALSE)
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
  if (stats::sd(x) <= .Machine$double.eps || stats::sd(y) <= .Machine$double.eps) {
    return(NA_real_)
  }
  stats::cor(x, y)
}
