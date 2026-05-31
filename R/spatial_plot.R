# Simple spatial atom plot helper

.plot_spatial_atom_display <- function(plot) {
  print(plot)
}

.plot_spatial_atom_idx <- function(idx, n_components) {
  if (length(idx) != 1L || is.na(idx) || !is.finite(idx) ||
      idx != floor(idx) || idx < 1L || idx > n_components) {
    .encoder_cli_abort(
      "idx out of range.",
      class = "fmrilatent_error_invalid_argument"
    )
  }
  as.integer(idx)
}

#' Plot a spatial atom (loading vector) on a mask
#'
#' @param loadings Matrix (voxels x k) or LoadingsHandle.
#' @param mask LogicalNeuroVol or logical array defining voxel order.
#' @param idx Component index (1-based).
#' @param main Optional title.
#' @return Invisibly, the 3D array plotted.
#' @export
plot_spatial_atom <- function(loadings, mask, idx = 1L, main = NULL) {
  if (inherits(loadings, "LoadingsHandle")) {
    idx <- .plot_spatial_atom_idx(idx, .latent_loadings_dim(loadings)[2L])
    loadings <- loadings_mat(loadings, j = idx)
    idx <- 1L
  }
  if (is.null(dim(loadings))) {
    .encoder_cli_abort(
      "loadings must be a matrix or LoadingsHandle.",
      class = "fmrilatent_error_invalid_type"
    )
  }
  if (ncol(loadings) < 1L) {
    .encoder_cli_abort(
      "loadings must have at least one component column.",
      class = "fmrilatent_error_dim"
    )
  }
  idx <- .plot_spatial_atom_idx(idx, ncol(loadings))
  mask_arr <- .mask_to_array(mask, "plot_spatial_atom")
  vec <- as.numeric(loadings[, idx])
  if (length(vec) != sum(mask_arr)) {
    .encoder_cli_abort(
      "loadings rows must match the number of TRUE voxels in mask.",
      class = "fmrilatent_error_dim"
    )
  }
  arr <- array(0, dim = dim(mask_arr))
  arr[mask_arr] <- vec
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    graphics::image(arr[, , 1, drop = TRUE], main = main %||% paste0("Atom ", idx))
    return(invisible(arr))
  }
  df <- expand.grid(x = seq_len(dim(arr)[1]), y = seq_len(dim(arr)[2]), z = seq_len(dim(arr)[3]))
  df$value <- as.vector(arr)
  # collapse across z with max intensity
  agg <- aggregate(value ~ x + y, data = df, FUN = function(v) v[which.max(abs(v))])
  p <- ggplot2::ggplot(agg, ggplot2::aes(x, y, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_distiller(palette = "RdBu", direction = -1) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = main %||% paste0("Atom ", idx), x = "X", y = "Y", fill = "Loading") +
    ggplot2::theme_minimal()
  .plot_spatial_atom_display(p)
  invisible(arr)
}
