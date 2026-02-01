# DCT temporal basis + heat-wavelet spatial dictionary constructor

#' Create a LatentNeuroVec with heat-wavelet spatial dictionary
#'
#' Creates a template LatentNeuroVec with heat-wavelet spatial loadings.
#' The basis matrix is initialized to zeros and should be populated with
#' actual coefficients (e.g., via encoding fMRI data).
#'
#' @param n_time        Number of time points.
#' @param k_time        Ignored (kept for backwards compatibility). The number
#'   of components is determined by the heat-wavelet loadings.
#' @param mask          LogicalNeuroVol or logical array mask (3D).
#' @param cluster_map   Optional integer vector mapping voxels (mask order) to clusters.
#' @param reduction     Graph reduction object; if NULL, built via
#'   `make_cluster_reduction(mask, cluster_map)` with default one-cluster-per-voxel map.
#' @param hw_basis_spec Heat-wavelet basis spec; defaults to `basis_heat_wavelet()`.
#' @param offset        Optional voxel-wise offset (length n_vox).
#' @param label         Optional label.
#'
#' @return A \code{LatentNeuroVec} with placeholder basis matrix.
#' @export
latent_dct_heatwavelet <- function(n_time,
                                   k_time,
                                   mask,
                                   cluster_map   = NULL,
                                   reduction     = NULL,
                                   hw_basis_spec = NULL,
                                   offset        = numeric(0),
                                   label         = "DCT + heat-wavelet") {
  # Build mask volume and array first

  mask_arr <- .mask_to_array(mask, "latent_dct_heatwavelet")
  mask_vol <- if (inherits(mask, "LogicalNeuroVol")) {
    mask
  } else {
    neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dim(mask_arr)))
  }

  # Spatial reduction + heat-wavelet spec
  if (is.null(hw_basis_spec)) {
    hw_basis_spec <- basis_heat_wavelet()
  }
  if (is.null(reduction)) {
    if (is.null(cluster_map)) {
      cluster_map <- seq_len(sum(mask_arr))
    }
    reduction <- make_cluster_reduction(mask, cluster_map)
  }

  # Create loadings handle first - this determines the number of components

  L_handle <- heat_wavelet_loadings_handle(
    reduction  = reduction,
    basis_spec = hw_basis_spec,
    data       = NULL,
    label      = "heat-wavelet"
  )

  # Get dimensions from loadings
  n_vox <- .latent_loadings_dim(L_handle)[1L]
  k_components <- .latent_loadings_dim(L_handle)[2L]

  # Validate offset
  if (!length(offset)) {
    offset <- numeric(n_vox)
  } else if (length(offset) != n_vox) {
    stop("length(offset) (", length(offset),
         ") must equal n_vox (", n_vox, ").")
  }

  # Create explicit basis matrix with matching number of columns

  # Initialized to zeros - users should populate with actual coefficients
  basis0 <- Matrix::Matrix(
    0,
    nrow = as.integer(n_time),
    ncol = k_components,
    sparse = FALSE
  )

  # Build 4D space
  spc <- neuroim2::NeuroSpace(c(dim(mask_arr), as.integer(n_time)))

  LatentNeuroVec(
    basis    = basis0,
    loadings = L_handle,
    space    = spc,
    mask     = mask_vol,
    offset   = offset,
    label    = label,
    meta     = list(
      time_basis   = "template",
      time_k       = k_components,
      spatial_dict = "heat_wavelet"
    )
  )
}
