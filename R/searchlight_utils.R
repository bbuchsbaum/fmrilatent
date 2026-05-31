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
#' @keywords internal
compute_local_gram <- function(loadings, neighborhoods, simplify = FALSE) {
  L <- loadings_mat(loadings)
  neighborhoods <- .validate_searchlight_neighborhoods(
    neighborhoods,
    n_voxels = nrow(L),
    context = "compute_local_gram"
  )
  k <- ncol(L)
  out <- vector("list", length(neighborhoods))
  for (i in seq_along(neighborhoods)) {
    idx <- neighborhoods[[i]]
    Li <- L[idx, , drop = FALSE]
    out[[i]] <- Matrix::crossprod(Li)
  }
  if (simplify && length(out) > 0) {
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
  if (ncol(B) != ncol(L)) {
    .encoder_cli_abort(
      paste0("latent_searchlight requires ncol(basis) == ncol(loadings); got ",
             ncol(B), " and ", ncol(L), "."),
      class = "fmrilatent_error_dimension_mismatch"
    )
  }
  neighborhoods <- .validate_searchlight_neighborhoods(
    neighborhoods,
    n_voxels = nrow(L),
    context = "latent_searchlight"
  )
  if (!is.function(fun)) {
    .encoder_cli_abort("fun must be a function.",
                       class = "fmrilatent_error_invalid_argument")
  }
  out <- vector("list", length(neighborhoods))
  for (i in seq_along(neighborhoods)) {
    idx <- neighborhoods[[i]]
    L_V <- L[idx, , drop = FALSE]
    M_V <- Matrix::crossprod(L_V)
    out[[i]] <- fun(B, L_V, M_V, idx, ...)
  }
  out
}

.validate_searchlight_neighborhoods <- function(neighborhoods, n_voxels, context) {
  if (!is.list(neighborhoods)) {
    .encoder_cli_abort(
      paste0(context, " requires neighborhoods to be a list of voxel index vectors."),
      class = "fmrilatent_error_invalid_argument"
    )
  }
  lapply(seq_along(neighborhoods), function(i) {
    idx <- neighborhoods[[i]]
    if (!is.numeric(idx)) {
      .encoder_cli_abort(
        sprintf("%s neighborhood %d must be numeric indices.", context, i),
        class = "fmrilatent_error_invalid_index"
      )
    }
    if (length(idx) == 0L) {
      return(integer())
    }
    if (anyNA(idx) || any(!is.finite(idx)) || any(idx != floor(idx))) {
      .encoder_cli_abort(
        sprintf("%s neighborhood %d contains missing or non-integer indices.", context, i),
        class = "fmrilatent_error_invalid_index"
      )
    }
    idx <- as.integer(idx)
    if (any(idx < 1L | idx > n_voxels)) {
      .encoder_cli_abort(
        sprintf("%s neighborhood %d contains indices outside [1, %d].",
                context, i, n_voxels),
        class = "fmrilatent_error_invalid_index"
      )
    }
    idx
  })
}
