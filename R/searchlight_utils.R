# Generic latent-space searchlight helpers

#' Compute local Gram matrices for neighborhoods
#'
#' @description
#' Given a loadings matrix L (voxels x k) and a list of voxel index sets
#' (neighborhoods), returns the per-neighborhood Gram matrices
#' `M_i = t(L_Vi) %*% L_Vi` without reconstructing voxels.
#'
#' @param loadings Matrix or LoadingsHandle (voxels x k) from a LatentNeuroVec.
#' @param neighborhoods List of integer vectors of voxel indices (mask order).
#' @param simplify Logical; if TRUE and all neighborhoods same size, returns an
#'   array k x k x n_neighborhood; otherwise a list of matrices.
#' @return List (or array) of Gram matrices.
#' @export
compute_local_gram <- function(loadings, neighborhoods, simplify = FALSE) {
  L <- loadings_mat(loadings)
  k <- ncol(L)
  out <- vector("list", length(neighborhoods))
  for (i in seq_along(neighborhoods)) {
    idx <- neighborhoods[[i]]
    Li <- L[idx, , drop = FALSE]
    out[[i]] <- Matrix::crossprod(Li)
  }
  if (simplify && length(out) > 0 && all(vapply(out, function(m) all(dim(m) == c(k, k)), logical(1)))) {
    arr <- array(NA_real_, dim = c(k, k, length(out)))
    for (i in seq_along(out)) arr[, , i] <- as.matrix(out[[i]])
    return(arr)
  }
  out
}

#' Apply a user-defined function in latent space over neighborhoods
#'
#' @description
#' Runs a user-supplied function `fun` for each neighborhood using only latent
#' quantities. `fun` is called with arguments `(B, L_V, M_V, idx, ...)`, where:
#'   - `B` is the basis matrix (time x k)
#'   - `L_V` is the loadings restricted to the neighborhood (|V| x k)
#'   - `M_V = t(L_V) %*% L_V` (k x k)
#'   - `idx` is the voxel indices of the neighborhood
#' `fun` should return any R object; results are collected in a list.
#'
#' @param basis Matrix or BasisHandle (time x k) from a LatentNeuroVec.
#' @param loadings Matrix or LoadingsHandle (voxels x k) from a LatentNeuroVec.
#' @param neighborhoods List of integer vectors of voxel indices (mask order).
#' @param fun Function(B, L_V, M_V, idx, ...) returning a result per neighborhood.
#' @param ... Passed through to `fun`.
#' @return List of results, one per neighborhood.
#' @export
latent_searchlight <- function(basis, loadings, neighborhoods, fun, ...) {
  B <- basis_mat(basis)
  L <- loadings_mat(loadings)
  out <- vector("list", length(neighborhoods))
  for (i in seq_along(neighborhoods)) {
    idx <- neighborhoods[[i]]
    L_V <- L[idx, , drop = FALSE]
    M_V <- Matrix::crossprod(L_V)
    out[[i]] <- fun(B, L_V, M_V, idx, ...)
  }
  out
}
