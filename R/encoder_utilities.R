# Shared utilities for encoder families

#' @include latent_handles.R
NULL

.encoder_cli_abort <- function(message, class = "fmrilatent_error_encoder", call = NULL) {
  cli::cli_abort(c(x = message), class = class, call = call)
}

.extract_mask_array <- function(mask, location = "unknown function") {
  mask_arr <- tryCatch(.mask_to_array(mask, location), error = function(e) NULL)
  if (is.null(mask_arr) || is.list(mask_arr)) {
    if (is.list(mask) && !is.null(mask$arr)) {
      mask_arr <- mask$arr
    } else {
      .encoder_cli_abort(
        sprintf("Cannot extract array from mask (%s)", location),
        class = "fmrilatent_error_mask_extract"
      )
    }
  }
  if (is.null(dim(mask_arr))) {
    .encoder_cli_abort(
      sprintf("mask must have dimensions (%s)", location),
      class = "fmrilatent_error_mask_extract"
    )
  }
  array(as.logical(mask_arr), dim = dim(mask_arr))
}

.mask_volume_from_array <- function(mask, mask_arr = NULL, location = "unknown function") {
  if (inherits(mask, "LogicalNeuroVol")) {
    return(mask)
  }
  mask_arr <- mask_arr %||% .extract_mask_array(mask, location)
  neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dim(mask_arr)))
}

.build_sparse_dictionary_from_clusters <- function(reduction, block_fun) {
  map <- reduction@map
  n_vox <- length(map)
  i_list <- list()
  j_list <- list()
  x_list <- list()
  col_offset <- 0L

  for (cid in reduction@cluster_ids) {
    vox_idx <- which(map == cid)
    if (length(vox_idx) == 0L) next

    block <- block_fun(vox_idx, cid)
    block <- Matrix::Matrix(block, sparse = TRUE)
    trip <- Matrix::summary(as(block, "TsparseMatrix"))
    if (nrow(trip) > 0L) {
      i_list[[length(i_list) + 1L]] <- vox_idx[trip$i]
      j_list[[length(j_list) + 1L]] <- col_offset + trip$j
      x_list[[length(x_list) + 1L]] <- trip$x
    }
    col_offset <- col_offset + ncol(block)
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

.build_voxel_knn_graph <- function(mask, voxel_indices, k_neighbors = 6L) {
  voxel_subset_to_gsp(mask, voxel_indices, k_neighbors = k_neighbors)
}

.build_cluster_centroid_graph <- function(reduction, k_neighbors = 6L) {
  if (!requireNamespace("rgsp", quietly = TRUE)) {
    .encoder_cli_abort(
      "rgsp not installed; install bbuchsbaum/rgsp to use diffusion wavelets",
      class = "fmrilatent_error_missing_package"
    )
  }

  mask_arr <- .extract_mask_array(reduction@mask, "build_cluster_centroid_graph")
  coords <- which(mask_arr, arr.ind = TRUE)
  ids <- reduction@cluster_ids
  grp <- factor(reduction@map, levels = ids)

  centers <- vapply(seq_len(ncol(coords)), function(d) {
    tapply(coords[, d], grp, mean)
  }, numeric(length(ids)))
  centers <- Matrix::Matrix(centers, sparse = FALSE)

  g <- rgsp::graph_knn(centers, k = k_neighbors, weight = "distance", sym = "union")
  row_normalize_sparse(g$adjacency)
}

.make_latent_neurovector <- function(X, mask, loadings, basis = NULL, offset = numeric(0),
                                     label = "", meta = list(), location = "encoder",
                                     expect_dense = FALSE) {
  mask_arr <- .extract_mask_array(mask, location)
  basis <- basis %||% (as.matrix(X) %*% as.matrix(loadings))
  spc <- .space_with_time_from_mask(mask, nrow(X), location)
  mask_vol <- .mask_volume_from_array(mask, mask_arr, location)
  LatentNeuroVec(
    basis = basis,
    loadings = loadings,
    space = spc,
    mask = mask_vol,
    offset = offset,
    label = label,
    meta = meta,
    expect_dense = expect_dense
  )
}
