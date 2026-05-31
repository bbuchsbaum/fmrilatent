# Graph bridge helpers (rgsp)

#' Convert voxel subset to an rgsp graph
#'
#' @param mask LogicalNeuroVol or 3D logical array.
#' @param voxel_indices Integer vector of voxel indices (mask order).
#' @param k_neighbors Number of nearest neighbours for k-NN graph.
#' @return `gsp_graph` object (from rgsp).
#' @export
voxel_subset_to_gsp <- function(mask, voxel_indices, k_neighbors = 6L) {
  if (!.encoder_has_namespace("rgsp")) {
    .encoder_cli_abort(
      "rgsp not installed; install bbuchsbaum/rgsp to use heat wavelets",
      class = c("fmrilatent_error_missing_dependency", "fmrilatent_error"),
      call = rlang::caller_env()
    )
  }
  mask_arr <- .mask_to_array(mask, "voxel_subset_to_gsp")
  # derive voxel coordinates (mask order) directly to avoid dependency on indices()
  # Note: which() with arr.ind=TRUE on a 3D array returns a 3-col matrix
  idx_all <- which(mask_arr != 0, arr.ind = TRUE)
  coords <- idx_all[voxel_indices, , drop = FALSE]
  n_coords <- nrow(coords)
  if (n_coords < 2L) {
    .encoder_cli_abort(
      "voxel_subset_to_gsp requires at least 2 voxel coordinates.",
      class = "fmrilatent_error_invalid_graph",
      call = rlang::caller_env()
    )
  }
  if (length(k_neighbors) != 1L || is.na(k_neighbors) || k_neighbors < 1L) {
    .encoder_cli_abort(
      "k_neighbors must be a positive integer.",
      class = "fmrilatent_error_invalid_count",
      call = rlang::caller_env()
    )
  }
  k_use <- min(as.integer(k_neighbors), n_coords - 1L)

  g <- rgsp::graph_knn(coords, k = k_use, weight = "distance", sym = "union")
  g
}
