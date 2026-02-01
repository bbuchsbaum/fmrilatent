#' Active-pencil wavelet latent (CDF 5/3)
#'
#' @param X Numeric matrix (time x voxels in mask order) or 4D array.
#' @param mask LogicalNeuroVol or 3D logical array.
#' @param levels_space Integer spatial lifting levels (default 2).
#' @param levels_time Integer temporal lifting levels (default 0 = none).
#' @param threshold Optional hard threshold after space transform.
#' @return An `ImplicitLatent` with decoder supporting `time_idx` and `roi_mask`.
#' @export
wavelet_active_latent <- function(X, mask, levels_space = 2L, levels_time = 0L, threshold = 0) {
  if (inherits(mask, "LogicalNeuroVol")) {
    mask_arr <- as.array(mask)
    coords <- which(mask_arr, arr.ind = TRUE)
  } else {
    mask_arr <- as.array(mask)
    coords <- which(mask_arr, arr.ind = TRUE)
    mask <- LogicalNeuroVol(mask_arr, NeuroSpace(dim(mask_arr)))
  }

  if (is.array(X) && length(dim(X)) == 4) {
    idx <- which(mask_arr)
    n_time <- dim(X)[4]
    X_mat <- matrix(0, nrow = n_time, ncol = length(idx))
    for (t in seq_len(n_time)) {
      X_mat[t, ] <- X[, , , t][idx]
    }
  } else {
    X_mat <- X
  }

  if (levels_time > 0) {
    X_mat <- cdf53_time_lift(X_mat, levels_time, TRUE)
  }

  coeffs_mat <- t(apply(X_mat, 1, function(vol_vec) {
    active_pencil_wavelet(vol_vec, coords, dim(mask_arr), levels_space, TRUE)
  }))

  if (threshold > 0) coeffs_mat[abs(coeffs_mat) < threshold] <- 0

  meta <- list(
    family = "wavelet_active",
    dims = dim(mask_arr),
    coords = coords,
    levels_space = levels_space,
    levels_time = levels_time,
    threshold = threshold,
    n_time = nrow(coeffs_mat)
  )

  decoder <- function(time_idx = NULL, roi_mask = NULL, ...) {
    # If time lifting was applied, we need full time inverse; subset afterward.
    need_full_time <- meta$levels_time > 0
    t_sel_full <- seq_len(meta$n_time)
    sub <- coeffs_mat[t_sel_full, , drop = FALSE]
    rec_mat <- t(apply(sub, 1, function(c_vec) {
      active_pencil_wavelet(c_vec, meta$coords, meta$dims, meta$levels_space, FALSE)
    }))
    if (meta$levels_time > 0) {
      rec_mat <- cdf53_time_lift(rec_mat, meta$levels_time, FALSE)
    }
    if (!is.null(time_idx)) {
      rec_mat <- rec_mat[as.integer(time_idx), , drop = FALSE]
    }
    if (!is.null(roi_mask)) {
      idx <- which(as.logical(roi_mask))
      rec_mat <- rec_mat[, idx, drop = FALSE]
    }
    rec_mat
  }

  implicit_latent(coeff = coeffs_mat, decoder = decoder, meta = meta, mask = mask_arr)
}
