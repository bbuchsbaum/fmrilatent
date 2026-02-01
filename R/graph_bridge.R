# Graph bridge helpers (rgsp)

#' Convert voxel subset to an rgsp graph
#'
#' @param mask LogicalNeuroVol or 3D logical array.
#' @param voxel_indices Integer vector of voxel indices (mask order).
#' @param k_neighbors Number of nearest neighbours for k-NN graph.
#' @return `gsp_graph` object (from rgsp).
#' @export
voxel_subset_to_gsp <- function(mask, voxel_indices, k_neighbors = 6L) {
  if (!requireNamespace("rgsp", quietly = TRUE)) {
    stop("rgsp not installed; install bbuchsbaum/rgsp to use heat wavelets", call. = FALSE)
  }
  mask_arr <- .mask_to_array(mask, "voxel_subset_to_gsp")
  # derive voxel coordinates (mask order) directly to avoid dependency on indices()
  # Note: which() with arr.ind=TRUE on a 3D array returns a 3-col matrix
  idx_all <- which(mask_arr != 0, arr.ind = TRUE)
  coords <- idx_all[voxel_indices, , drop = FALSE]


  g <- rgsp::graph_knn(coords, k = k_neighbors, weight = "distance", sym = "union")
  g
}
