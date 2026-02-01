# Simple spatial atom plot helper

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
    loadings <- loadings_mat(loadings)
  }
  if (is.null(dim(loadings))) stop("loadings must be a matrix or LoadingsHandle.", call. = FALSE)
  if (idx < 1 || idx > ncol(loadings)) stop("idx out of range.", call. = FALSE)
  mask_arr <- .mask_to_array(mask, "plot_spatial_atom")
  vec <- as.numeric(loadings[, idx])
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
  ggplot2::ggplot(agg, ggplot2::aes(x, y, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_distiller(palette = "RdBu", direction = -1) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = main %||% paste0("Atom ", idx), x = "X", y = "Y", fill = "Loading") +
    ggplot2::theme_minimal()
  invisible(arr)
}
