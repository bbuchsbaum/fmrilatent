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
      stop("rgsp not installed; install bbuchsbaum/rgsp to use heat wavelets", call. = FALSE)
    }

    mask <- reduction@mask
    map <- reduction@map
    ids <- reduction@cluster_ids
    n_vox <- length(map)

    i_list <- list(); j_list <- list(); x_list <- list(); col_offset <- 0L

    for (cid in ids) {
      vox_idx <- which(map == cid)
      if (length(vox_idx) == 0L) next

      n_loc <- length(vox_idx)
      if (n_loc < 2L) {
        # Degenerate cluster: single voxel => wavelet is identity per scale
        W_list <- lapply(
          basis_spec$scales,
          function(...) Matrix::sparseMatrix(i = 1L, j = 1L, x = 1, dims = c(1, 1))
        )
      } else {
        g <- voxel_subset_to_gsp(mask, vox_idx, k_neighbors = k_neighbors)

        I_local <- diag(n_loc)
        W_list <- rgsp::wavelet_heat_transform(g, I_local,
                                               scales = basis_spec$scales,
                                               K = basis_spec$order)
      }

      for (s in seq_along(W_list)) {
        mat <- W_list[[s]]
        mat[abs(mat) < basis_spec$threshold] <- 0
        sp <- as(Matrix::Matrix(mat, sparse = TRUE), "TsparseMatrix")
        if (length(sp@x) > 0) {
          i_list[[length(i_list) + 1L]] <- vox_idx[sp@i + 1L]
          j_list[[length(j_list) + 1L]] <- col_offset + sp@j + 1L
          x_list[[length(x_list) + 1L]] <- sp@x
        }
        col_offset <- col_offset + ncol(mat)
      }
    }

    if (col_offset == 0L) {
      return(Matrix::Matrix(0, nrow = n_vox, ncol = 0, sparse = TRUE))
    }

    Matrix::sparseMatrix(
      i = unlist(i_list, use.names = FALSE),
      j = unlist(j_list, use.names = FALSE),
      x = unlist(x_list, use.names = FALSE),
      dims = c(n_vox, col_offset)
    )
  }
)

# Helper to build a ClusterReduction from a map/mask
make_cluster_reduction <- function(mask, map) {
  mask_vol <- if (inherits(mask, "LogicalNeuroVol")) mask else LogicalNeuroVol(mask, neuroim2::NeuroSpace(dim(mask)))
  new("ClusterReduction", mask = mask_vol, map = as.integer(map),
      cluster_ids = as.integer(sort(unique(map))), info = list())
}

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
  mask_arr <- .mask_to_array(mask, "heat_wavelet_latent")
  if (is.null(reduction)) {
    map <- seq_len(sum(mask_arr))
    reduction <- make_cluster_reduction(mask, map)
  }

  loadings <- lift(reduction, spec, k_neighbors = k_neighbors)
  basis <- as.matrix(X) %*% as.matrix(loadings)
  spc <- neuroim2::NeuroSpace(c(dim(mask_arr), nrow(X)))
  mask_vol <- if (inherits(mask, "LogicalNeuroVol")) mask else LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dim(mask_arr)))
  meta <- list(family = "heat_wavelet", spec = spec, k_neighbors = k_neighbors)
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask_vol, label = label, meta = meta)
}
