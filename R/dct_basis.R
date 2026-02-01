# Discrete Cosine Transform (DCT-II) temporal basis and handle

#' Build an orthonormal DCT-II basis matrix
#'
#' @param n_time Number of time points.
#' @param k      Number of components (columns); must satisfy k <= n_time.
#' @param norm   Normalization: "ortho" (default) or "none".
#' @return Dense Matrix (n_time x k).
#' @keywords internal
build_dct_basis <- function(n_time, k = n_time, norm = c("ortho", "none")) {
  norm   <- match.arg(norm)
  n_time <- as.integer(n_time)
  k      <- as.integer(k)

  if (k > n_time) stop("k (", k, ") cannot exceed n_time (", n_time, ").")

  tt <- seq.int(0L, n_time - 1L)
  kk <- seq.int(0L, k - 1L)

  mat <- outer(
    tt, kk,
    function(t, j) cos(pi * (t + 0.5) * j / n_time)
  )

  if (norm == "ortho") {
    mat[, 1L] <- mat[, 1L] / sqrt(n_time)
    if (k > 1L) {
      mat[, 2L:k] <- mat[, 2L:k, drop = FALSE] * sqrt(2 / n_time)
    }
  }

  Matrix::Matrix(mat, sparse = FALSE)
}

#' Create a BasisHandle for a DCT temporal basis
#'
#' @param n_time Number of time points.
#' @param k      Number of components.
#' @param norm   Normalization passed to [build_dct_basis()].
#' @param id     Optional registry key; if NULL a deterministic id is derived
#'   from parameters.
#' @param label  Optional human-readable label.
#'
#' @return A \code{BasisHandle} object.
#' @export
dct_basis_handle <- function(n_time, k, norm = c("ortho", "none"),
                             id = NULL, label = NULL) {
  norm   <- match.arg(norm)
  n_time <- as.integer(n_time)
  k      <- as.integer(k)

  if (is.null(id)) {
    id <- sprintf("dct-%d-%d-%s", n_time, k, norm)
  }
  if (is.null(label)) {
    label <- sprintf("DCT(%d,%d,%s)", n_time, k, norm)
  }

  new("BasisHandle",
      id    = id,
      dim   = as.integer(c(n_time, k)),
      kind  = "dct",
      spec  = list(n_time = n_time, k = k, norm = norm),
      label = label)
}
