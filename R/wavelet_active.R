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
  mask_arr <- .extract_mask_array(mask, "wavelet_active_latent")
  if (!any(mask_arr)) {
    .encoder_cli_abort(
      "wavelet_active_latent() requires a mask with at least one active voxel.",
      class = "fmrilatent_error_empty_mask"
    )
  }
  levels_space <- .validate_positive_count(levels_space, "levels_space")
  levels_time <- .validate_nonnegative_count(levels_time, "levels_time")
  threshold <- .validate_nonnegative_scalar(threshold, "threshold")
  coords <- which(mask_arr, arr.ind = TRUE)

  if (is.array(X) && length(dim(X)) == 4) {
    if (!identical(dim(X)[seq_len(3L)], dim(mask_arr))) {
      .encoder_cli_abort(
        "4D X spatial dimensions must match mask dimensions.",
        class = "fmrilatent_error_dim"
      )
    }
    idx <- which(mask_arr)
    n_time <- dim(X)[4]
    X_mat <- matrix(0, nrow = n_time, ncol = length(idx))
    for (t in seq_len(n_time)) {
      X_mat[t, ] <- X[, , , t][idx]
    }
  } else {
    X_mat <- X
  }
  X_mat <- as.matrix(X_mat)
  n_time <- nrow(X_mat)
  if (n_time < 1L) {
    .encoder_cli_abort(
      "wavelet_active_latent() requires at least one time point.",
      class = "fmrilatent_error_invalid_count"
    )
  }
  if (ncol(X_mat) != nrow(coords)) {
    .encoder_cli_abort(
      "X must have one column per active mask voxel.",
      class = "fmrilatent_error_dim"
    )
  }
  if (!.matrix_like_all_finite(X_mat)) {
    .encoder_cli_abort(
      "X must contain only finite values.",
      class = "fmrilatent_error_value"
    )
  }

  if (levels_time > 0) {
    X_mat <- cdf53_time_lift(X_mat, levels_time, TRUE)
  }

  coeffs_mat <- do.call(rbind, lapply(seq_len(nrow(X_mat)), function(row_idx) {
    active_pencil_wavelet(X_mat[row_idx, ], coords, dim(mask_arr), levels_space, TRUE)
  }))

  if (threshold > 0) coeffs_mat[abs(coeffs_mat) < threshold] <- 0

  meta <- list(
    family = "wavelet_active",
    dims = dim(mask_arr),
    coords = coords,
    levels_space = levels_space,
    levels_time = levels_time,
    threshold = threshold,
    n_time = n_time
  )

  decoder <- function(time_idx = NULL, roi_mask = NULL, ...) {
    t_sel_full <- seq_len(meta$n_time)
    sub <- coeffs_mat[t_sel_full, , drop = FALSE]
    rec_mat <- do.call(rbind, lapply(seq_len(nrow(sub)), function(row_idx) {
      active_pencil_wavelet(sub[row_idx, ], meta$coords, meta$dims, meta$levels_space, FALSE)
    }))
    if (meta$levels_time > 0) {
      rec_mat <- cdf53_time_lift(rec_mat, meta$levels_time, FALSE)
    }
    if (!is.null(time_idx)) {
      rec_mat <- rec_mat[as.integer(time_idx), , drop = FALSE]
    }
    if (!is.null(roi_mask)) {
      rec_mat <- roi_subset_columns(rec_mat, mask_arr, roi_mask)
    }
    rec_mat
  }

  implicit_latent(
    coeff = coeffs_mat,
    decoder = decoder,
    meta = meta,
    mask = if (inherits(mask, "LogicalNeuroVol")) mask else mask_arr
  )
}
