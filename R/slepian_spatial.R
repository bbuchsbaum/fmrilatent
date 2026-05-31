# Spatial Slepian lift on clustered reductions

.slepian_spatial_laplacian_vectors <- function(L, k) {
  k <- .validate_positive_count(k, "k")
  n <- nrow(L)
  if (is.null(n) || n < 1L) {
    .encoder_cli_abort("L must have at least one row.",
                       class = "fmrilatent_error_dimension_mismatch")
  }
  k <- min(k, n)
  if (n < 3L || k >= n) {
    eig <- eigen(as.matrix(L), symmetric = TRUE)
    ord <- order(eig$values)
    keep <- ord[seq_len(k)]
    return(eig$vectors[, keep, drop = FALSE])
  }
  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    .encoder_cli_abort("RSpectra is required for spatial Slepian bases",
                       class = "fmrilatent_error_missing_dependency")
  }
  RSpectra::eigs(L, k = k, which = "SM")$vectors
}

#' Lift spatial Slepians for clustered reduction
#'
#' @param reduction ClusterReduction describing voxel-to-cluster map.
#' @param basis_spec Slepian basis specification (from `basis_slepian()`).
#' @param data Ignored for this graph-only spatial dictionary; accepted only
#'   for the shared `lift()` method signature.
#' @param k_neighbors k for local graph building.
#' @param ... Additional arguments (unused).
#' @return Sparse Matrix (voxels x components), block-concatenated over clusters.
#' @export
setMethod("lift", signature(reduction = "ClusterReduction", basis_spec = "spec_slepian"),
  function(reduction, basis_spec, data = NULL, k_neighbors = 6L, ...) {
    mask <- reduction@mask
    k_per_cluster <- .validate_positive_count(basis_spec$k %||% 3L, "basis_spec$k")
    k_neighbors <- .validate_positive_count(k_neighbors, "k_neighbors")

    .build_sparse_dictionary_from_clusters(reduction, function(vox_idx, cid) {
      n_loc <- length(vox_idx)

      # Build local graph
      k_use <- min(k_per_cluster, n_loc)
      if (n_loc == 1L) {
        vecs <- Matrix::sparseMatrix(i = 1L, j = 1L, x = 1, dims = c(1, 1))
      } else {
        g <- .build_voxel_knn_graph(mask, vox_idx, k_neighbors = k_neighbors)
        L <- g$laplacian
        vecs <- Matrix::Matrix(.slepian_spatial_laplacian_vectors(L, k_use), sparse = TRUE)
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
#' @return A `LatentNeuroVec` object.
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
