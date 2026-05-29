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
    k_per_cluster <- as.integer(basis_spec$k %||% 3L)

    .build_sparse_dictionary_from_clusters(reduction, function(vox_idx, cid) {
      n_loc <- length(vox_idx)

      # Build local graph
      k_use <- min(k_per_cluster, n_loc)
      if (n_loc == 1L) {
        vecs <- Matrix::sparseMatrix(i = 1L, j = 1L, x = 1, dims = c(1, 1))
      } else if (n_loc < 3L) {
        # RSpectra requires dimension >= 3; use base eigen() for small clusters
        g <- .build_voxel_knn_graph(mask, vox_idx, k_neighbors = k_neighbors)
        L <- as.matrix(g$laplacian)
        eig <- eigen(L, symmetric = TRUE)
        ord <- order(eig$values)
        vecs <- Matrix::Matrix(eig$vectors[, ord[seq_len(k_use)], drop = FALSE], sparse = TRUE)
      } else {
        g <- .build_voxel_knn_graph(mask, vox_idx, k_neighbors = k_neighbors)
        L <- g$laplacian
        eig <- RSpectra::eigs(L, k = k_use, which = "SM")
        vecs <- Matrix::Matrix(eig$vectors, sparse = TRUE)
      }
      vecs
    })
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
  mask_arr <- .extract_mask_array(mask, "slepian_spatial_latent")
  if (is.null(reduction)) {
    map <- seq_len(sum(mask_arr))
    reduction <- make_cluster_reduction(mask, map)
  }

  loadings <- lift(reduction, spec, k_neighbors = k_neighbors)
  meta <- list(family = "slepian_spatial", spec = spec, k_neighbors = k_neighbors)
  .make_latent_neurovector(X, mask, loadings, label = label, meta = meta,
                           location = "slepian_spatial_latent")
}
