#' @include reduction.R
NULL

# Heat / diffusion graph wavelets (rgsp-backed)

#' Heat wavelet basis specification
#'
#' @param scales Numeric vector of heat diffusion scales.
#' @param order Polynomial approximation order.
#' @param threshold Threshold for small coefficients.
#' @return A `spec_heat_wavelet` object for `lift()`.
#' @export
basis_heat_wavelet <- function(scales = c(1, 2, 4, 8), order = 30, threshold = 1e-6) {
  structure(list(scales = scales, order = order, threshold = threshold),
            class = "spec_heat_wavelet")
}

#' Lift heat wavelets for clustered reduction
#' @param reduction ClusterReduction describing voxel-to-cluster map.
#' @param basis_spec Heat wavelet spec (basis_heat_wavelet()).
#' @param data Optional (unused) for API symmetry.
#' @param k_neighbors k for local graph building.
#' @param ... Additional arguments (unused).
#' @return Sparse Matrix (voxels x components).
#' @export
setMethod("lift", signature(reduction = "ClusterReduction", basis_spec = "spec_heat_wavelet"),
  function(reduction, basis_spec, data = NULL, k_neighbors = 6L, ...) {
    if (!requireNamespace("rgsp", quietly = TRUE)) {
      .encoder_cli_abort(
        "rgsp not installed; install bbuchsbaum/rgsp to use heat wavelets",
        class = "fmrilatent_error_missing_package"
      )
    }

    mask <- reduction@mask
    .build_sparse_dictionary_from_clusters(reduction, function(vox_idx, cid) {
      n_loc <- length(vox_idx)
      if (n_loc < 2L) {
        # Degenerate cluster: single voxel => wavelet is identity per scale
        W_list <- lapply(
          basis_spec$scales,
          function(...) Matrix::sparseMatrix(i = 1L, j = 1L, x = 1, dims = c(1, 1))
        )
      } else {
        g <- .build_voxel_knn_graph(mask, vox_idx, k_neighbors = k_neighbors)

        I_local <- diag(n_loc)
        W_list <- rgsp::wavelet_heat_transform(g, I_local,
                                               scales = basis_spec$scales,
                                               K = basis_spec$order)
      }

      blocks <- vector("list", length(W_list))
      for (s in seq_along(W_list)) {
        mat <- W_list[[s]]
        mat[abs(mat) < basis_spec$threshold] <- 0
        blocks[[s]] <- Matrix::Matrix(mat, sparse = TRUE)
      }
      do.call(cbind, blocks)
    })
  }
)

# make_cluster_reduction is now defined in R/reduction.R

#' Heat wavelet latent constructor (explicit basis)
#'
#' @param X Matrix time x voxels (mask order)
#' @param mask LogicalNeuroVol or 3D logical array
#' @param reduction ClusterReduction; if NULL, defaults to one cluster per voxel
#' @param spec heat wavelet spec (basis_heat_wavelet())
#' @param k_neighbors k for local graph building
#' @param label Optional label
#' @export
heat_wavelet_latent <- function(X, mask, reduction = NULL, spec = basis_heat_wavelet(),
                                k_neighbors = 6L, label = "") {
  mask_arr <- .extract_mask_array(mask, "heat_wavelet_latent")
  if (is.null(reduction)) {
    map <- seq_len(sum(mask_arr))
    reduction <- make_cluster_reduction(mask, map)
  }

  loadings <- lift(reduction, spec, k_neighbors = k_neighbors)
  meta <- list(family = "heat_wavelet", spec = spec, k_neighbors = k_neighbors)
  .make_latent_neurovector(X, mask, loadings, label = label, meta = meta,
                           location = "heat_wavelet_latent")
}
