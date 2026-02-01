# Spatial Slepian lift on clustered reductions

#' Lift spatial Slepians for clustered reduction
#'
#' @param reduction ClusterReduction describing voxel-to-cluster map.
#' @param basis_spec Slepian basis specification (from `basis_slepian()`).
#' @param data Optional (unused) for API symmetry.
#' @param k_neighbors k for local graph building.
#' @param ... Additional arguments (unused).
#' @return Sparse Matrix (voxels x components), block-concatenated over clusters.
#' @export
setMethod("lift", signature(reduction = "ClusterReduction", basis_spec = "spec_slepian"),
  function(reduction, basis_spec, data = NULL, k_neighbors = 6L, ...) {
    mask <- reduction@mask
    map <- reduction@map
    ids <- reduction@cluster_ids
    n_vox <- length(map)

    k_per_cluster <- as.integer(basis_spec$k %||% 3L)
    type <- basis_spec$type %||% "laplacian"

    i_list <- list(); j_list <- list(); x_list <- list(); col_offset <- 0L

    for (cid in ids) {
      vox_idx <- which(map == cid)
      n_loc <- length(vox_idx)
      if (n_loc == 0L) next

      # Build local graph
      k_use <- min(k_per_cluster, n_loc)
      if (n_loc == 1L) {
        vecs <- Matrix::sparseMatrix(i = 1L, j = 1L, x = 1, dims = c(1, 1))
      } else if (n_loc < 3L) {
        # RSpectra requires dimension >= 3; use base eigen() for small clusters
        g <- voxel_subset_to_gsp(mask, vox_idx, k_neighbors = k_neighbors)
        L <- as.matrix(g$laplacian)
        eig <- eigen(L, symmetric = TRUE)
        ord <- order(eig$values)
        vecs <- Matrix::Matrix(eig$vectors[, ord[seq_len(k_use)], drop = FALSE], sparse = TRUE)
      } else {
        g <- voxel_subset_to_gsp(mask, vox_idx, k_neighbors = k_neighbors)
        L <- g$laplacian
        eig <- RSpectra::eigs(L, k = k_use, which = "SM")
        vecs <- Matrix::Matrix(eig$vectors, sparse = TRUE)
      }

      # robust sparse triplets
      trip <- Matrix::summary(as(vecs, "TsparseMatrix"))
      if (nrow(trip) > 0) {
        i_list[[length(i_list) + 1L]] <- vox_idx[trip$i]
        j_list[[length(j_list) + 1L]] <- col_offset + trip$j
        x_list[[length(x_list) + 1L]] <- trip$x
      }
      col_offset <- col_offset + ncol(vecs)
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

#' Slepian spatial latent constructor (explicit basis)
#'
#' @param X Matrix time x voxels (mask order).
#' @param mask LogicalNeuroVol or 3D logical array.
#' @param reduction ClusterReduction; if NULL, defaults to one cluster per voxel.
#' @param spec Slepian basis spec (basis_slepian()).
#' @param k_neighbors k for local graph building.
#' @param label Optional label.
#' @export
slepian_spatial_latent <- function(X, mask, reduction = NULL, spec = basis_slepian(),
                                   k_neighbors = 6L, label = "") {
  mask_arr <- .mask_to_array(mask, "slepian_spatial_latent")
  if (is.null(reduction)) {
    map <- seq_len(sum(mask_arr))
    reduction <- make_cluster_reduction(mask, map)
  }

  loadings <- lift(reduction, spec, k_neighbors = k_neighbors)
  basis <- as.matrix(X) %*% as.matrix(loadings)
  spc <- neuroim2::NeuroSpace(c(dim(mask_arr), nrow(X)))
  mask_vol <- if (inherits(mask, "LogicalNeuroVol")) mask else LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dim(mask_arr)))
  meta <- list(family = "slepian_spatial", spec = spec, k_neighbors = k_neighbors)
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask_vol, label = label, meta = meta)
}
