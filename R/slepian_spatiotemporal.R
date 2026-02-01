# Spatiotemporal Slepians via separable bases (implicit latent)

#' Spatiotemporal Slepian latent (implicit, separable)
#'
#' @param X Numeric matrix time x voxels (mask order).
#' @param mask LogicalNeuroVol or 3D logical array.
#' @param tr Repetition time (seconds).
#' @param bandwidth Half-bandwidth in Hz for temporal Slepians (default 0.1).
#' @param k_time Number of temporal Slepians; if NULL uses floor(2*NW)-1.
#' @param reduction ClusterReduction for spatial graph; if NULL, one cluster per voxel.
#' @param k_space Number of spatial Slepians per cluster (default 3).
#' @param k_neighbors k-NN for spatial graph (default 6).
#' @param label Optional label.
#' @return An `ImplicitLatent` with decoder using separable Slepians.
#' @export
slepian_spatiotemporal_latent <- function(X, mask, tr,
                                          bandwidth = 0.1,
                                          k_time = NULL,
                                          reduction = NULL,
                                          k_space = 3L,
                                          k_neighbors = 6L,
                                          label = "") {
  mask_arr <- .mask_to_array(mask, "slepian_spatiotemporal_latent")
  n_time <- nrow(X)
  if (is.null(n_time) || n_time < 1) stop("X must have time in rows")
  n_vox <- sum(mask_arr)

  # Temporal Slepians
  W <- bandwidth * tr
  NW <- n_time * W
  if (is.null(k_time)) k_time <- max(1L, floor(2 * NW) - 1L)
  B_t <- dpss_time_basis(n_time, tr = tr, bandwidth = bandwidth, k = k_time, backend = "tridiag")

  # Spatial reduction + Slepians
  if (is.null(reduction)) {
    map <- seq_len(n_vox)
    reduction <- make_cluster_reduction(mask, map)
  }
  spec_sp <- basis_slepian(k = k_space)
  L_s <- lift(reduction, spec_sp, k_neighbors = k_neighbors) # vox x k_space_total

  # Core tensor (kt x ks): project data into separable bases
  B_t_mat <- as.matrix(B_t)
  L_s_mat <- as.matrix(L_s)
  core <- crossprod(B_t_mat, X) %*% L_s_mat  # kt x ks

  decoder <- function(time_idx = NULL, roi_mask = NULL, ...) {
    t_sel <- if (is.null(time_idx)) seq_len(n_time) else as.integer(time_idx)
    B_sel <- B_t_mat[t_sel, , drop = FALSE]
    rec_mat <- B_sel %*% core %*% t(L_s_mat)  # time x vox

    if (!is.null(roi_mask)) {
      roi_idx <- which(as.logical(roi_mask))
      rec_mat <- rec_mat[, roi_idx, drop = FALSE]
    }
    rec_mat
  }

  meta <- list(
    family = "slepian_spatiotemporal",
    tr = tr,
    bandwidth = bandwidth,
    k_time = k_time,
    k_space = k_space,
    label = label
  )

  implicit_latent(coeff = list(core = core, B_t = B_t_mat, L_s = L_s_mat),
                  decoder = decoder,
                  meta = meta,
                  mask = mask_arr)
}
