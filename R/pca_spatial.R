# Spatial PCA lift on clustered reductions

#' @include reduction.R
NULL

.pca_truncated_svd <- function(X, k, backend = c("auto", "svds", "svd")) {
  backend <- match.arg(backend)
  n_time <- nrow(X)
  n_vox <- ncol(X)
  k_max <- min(n_time, n_vox)
  k <- as.integer(min(k, k_max))

  if (k < 1L) {
    return(list(v = matrix(0, nrow = n_vox, ncol = 0L), d = numeric(0)))
  }

  if (backend == "auto") {
    backend <- if (requireNamespace("RSpectra", quietly = TRUE) && k < k_max) "svds" else "svd"
  }

  if (backend == "svds") {
    if (!requireNamespace("RSpectra", quietly = TRUE) || k >= k_max) {
      backend <- "svd"
    }
  }

  if (backend == "svds") {
    out <- RSpectra::svds(X, k = k, nu = 0, nv = k)
    list(v = out$v, d = out$d)
  } else {
    out <- svd(X, nu = 0, nv = k)
    list(v = out$v, d = out$d)
  }
}

.pca_flip_sign <- function(v) {
  if (ncol(v) < 1L) return(v)
  for (j in seq_len(ncol(v))) {
    idx <- which.max(abs(v[, j]))
    if (length(idx) && v[idx, j] < 0) {
      v[, j] <- -v[, j]
    }
  }
  v
}

#' Lift parcel/cluster-local PCA bases for ClusterReduction
#'
#' Computes PCA eigenvectors within each cluster and assembles a global
#' block-sparse loadings matrix (voxels x components). This is typically used
#' with `encode(..., spec_space_pca(...), reduction = ...)`.
#'
#' @param reduction ClusterReduction describing voxel-to-cluster map.
#' @param basis_spec PCA basis specification (from `basis_pca()`).
#' @param data Numeric matrix (time x voxels in mask order). Required.
#' @param center Logical; center voxels before PCA (default TRUE).
#' @param scale Logical; scale voxels before PCA (default FALSE).
#' @param offset Optional numeric vector of voxel means (length n_vox). If provided
#'   and `center = TRUE`, this is used instead of recomputing `colMeans(data)`.
#' @param backend SVD backend: "auto" (default), "svds" (RSpectra), or "svd" (base).
#' @param ... Unused.
#' @return A sparse Matrix (voxels x components) with attribute
#'   `fmrilatent.singular_values` giving per-component singular values.
#' @export
setMethod("lift", signature(reduction = "ClusterReduction", basis_spec = "spec_pca"),
  function(reduction, basis_spec, data = NULL,
           center = TRUE, scale = FALSE, offset = NULL,
           backend = c("auto", "svds", "svd"), ...) {
    backend <- match.arg(backend)
    if (is.null(data)) {
      stop("lift(ClusterReduction, spec_pca) requires 'data' (time x voxels).", call. = FALSE)
    }

    X <- as.matrix(data)
    map <- reduction@map
    ids <- reduction@cluster_ids
    n_vox <- length(map)
    if (ncol(X) != n_vox) {
      stop("data has ", ncol(X), " voxels, but reduction map has length ", n_vox, ".", call. = FALSE)
    }

    k_per_cluster <- as.integer(basis_spec$k %||% 3L)

    mu <- NULL
    if (isTRUE(center)) {
      if (!is.null(offset)) {
        if (!is.numeric(offset) || length(offset) != n_vox) {
          stop("offset must be a numeric vector of length n_vox when provided.", call. = FALSE)
        }
        mu <- as.numeric(offset)
      } else {
        mu <- colMeans(X)
      }
    }

    sig <- NULL
    if (isTRUE(scale)) {
      sig <- apply(X, 2, stats::sd)
      sig[!is.finite(sig) | sig == 0] <- 1
    }

    i_list <- list()
    j_list <- list()
    x_list <- list()
    d_list <- list()
    col_offset <- 0L

    for (cid in ids) {
      vox_idx <- which(map == cid)
      n_loc <- length(vox_idx)
      if (n_loc == 0L) next

      k_use <- min(k_per_cluster, n_loc, nrow(X))
      if (k_use < 1L) next

      X_loc <- X[, vox_idx, drop = FALSE]
      if (!is.null(mu)) {
        X_loc <- sweep(X_loc, 2, mu[vox_idx], "-")
      }
      if (!is.null(sig)) {
        X_loc <- sweep(X_loc, 2, sig[vox_idx], "/")
      }

      if (n_loc == 1L) {
        v <- matrix(1, nrow = 1L, ncol = 1L)
        d <- sqrt(sum(X_loc^2))
      } else {
        sv <- .pca_truncated_svd(X_loc, k = k_use, backend = backend)
        v <- .pca_flip_sign(sv$v)
        d <- sv$d
      }

      i_list[[length(i_list) + 1L]] <- rep.int(vox_idx, times = k_use)
      j_list[[length(j_list) + 1L]] <- rep.int(col_offset + seq_len(k_use), times = rep.int(n_loc, k_use))
      x_list[[length(x_list) + 1L]] <- as.vector(v[, seq_len(k_use), drop = FALSE])
      d_list[[length(d_list) + 1L]] <- as.numeric(d[seq_len(k_use)])

      col_offset <- col_offset + k_use
    }

    if (col_offset == 0L) {
      out <- Matrix::Matrix(0, nrow = n_vox, ncol = 0, sparse = TRUE)
      attr(out, "fmrilatent.singular_values") <- numeric(0)
      return(out)
    }

    out <- Matrix::sparseMatrix(
      i = unlist(i_list, use.names = FALSE),
      j = unlist(j_list, use.names = FALSE),
      x = unlist(x_list, use.names = FALSE),
      dims = c(n_vox, col_offset)
    )
    attr(out, "fmrilatent.singular_values") <- unlist(d_list, use.names = FALSE)
    out
  }
)
