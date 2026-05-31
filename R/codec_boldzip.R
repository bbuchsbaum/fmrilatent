# Experimental BOLDZip-SR matrix codec
#
# This file holds the public orchestration core for the BOLDZip-SR codec.
# The codec was decomposed into thematic files (issue bd-01KSQQPCEQSHQ6WTVMCXZR1VME);
# helpers now live alongside the public functions they support:
#
#   R/codec_boldzip_validate.R    - scalar/matrix/basis validators + linalg
#                                   primitives (project, synthesize,
#                                   orthonormalize_columns,
#                                   canonicalize_eigenvectors, corr).
#   R/codec_boldzip_decompose.R   - temporal-spec materialization, carrier
#                                   bank (lag_signal, lagged_carrier_bank,
#                                   predict_texture, fit_sparse_texture),
#                                   events (pair_time, encode_events,
#                                   events_matrix), pairing (pair_indices,
#                                   row_cor).
#   R/codec_boldzip_quantize.R    - .boldzip_quantize_values(),
#                                   boldzip_quantization(), and the unified
#                                   noise-scale policy (.boldzip_noise_scale).
#   R/codec_boldzip_spatial.R     - boldzip_spatial_basis,
#                                   as_boldzip_spatial_basis (+ methods),
#                                   boldzip_graph_spatial_basis,
#                                   loadings_from_template.
#   R/codec_boldzip_diagnostics.R - evaluate_boldzip_sr, boldzip_sr_simulate,
#                                   boldzip_parcel_reconstruct,
#                                   boldzip_svd_reconstruct, compare_boldzip_sr,
#                                   boldzip_sr_payload_summary, and the
#                                   reliability/events settings constructors.
#
# Kept here: boldzip_sr_encode, boldzip_sr_decode, predict.BoldZipSR,
# as.matrix.BoldZipSR, as_implicit_latent.BoldZipSR, .boldzip_roi_from_mask.


#' Encode a matrix with experimental BOLDZip-SR compression
#'
#' @description
#' `boldzip_sr_encode()` is an experimental matrix-level implementation of
#' Split-Reliable Graph Carrier Compression. It expects a matrix with rows as
#' voxels or grayordinates and columns as time points. The implementation stores
#' temporally compressed carriers, sparse high-resolution texture loadings, and
#' sparse split-reliable residual events. It is intended as a research prototype,
#' not a production binary codec.
#'
#' @param X Numeric matrix with dimensions `voxels x time`.
#' @param spatial_basis Optional `BoldZipSRSpatialBasis` object or list with
#'   `phi_c` and `phi_d` matrices. If omitted, the detail basis is identity and
#'   carriers are learned from `X` directly.
#' @param k_carriers Number of carrier time courses to learn.
#' @param temporal_k Number of DCT coefficients per carrier. Defaults to
#'   `ceiling(n_time / 4)` when `temporal_spec` is `NULL`.
#' @param temporal_spec Optional temporal basis descriptor. May be a
#'   `SharedTemporalSpec`, `spec_time_dct`, `spec_time_bspline`, or numeric
#'   basis matrix with rows equal to time points. When supplied, it determines
#'   the temporal basis and `temporal_k`.
#' @param q_texture Maximum number of carrier loadings per detail atom.
#' @param texture_lags Integer vector of allowed carrier lags for texture
#'   loading fits. A positive lag uses `Z_k(t - lag)`.
#' @param reliability Reliability settings from [boldzip_reliability()].
#' @param quantization Quantization settings from [boldzip_quantization()].
#' @param events Event settings from [boldzip_events()].
#' @param center Whether to store and remove a voxel-wise mean before fitting.
#' @param label Optional label stored in metadata.
#' @return A `BoldZipSR` object.
#' @export
boldzip_sr_encode <- function(X,
                              spatial_basis = NULL,
                              k_carriers = 96L,
                              temporal_k = NULL,
                              temporal_spec = NULL,
                              q_texture = 2L,
                              texture_lags = 0L,
                              reliability = boldzip_reliability(),
                              quantization = boldzip_quantization(),
                              events = boldzip_events(),
                              center = TRUE,
                              label = NULL) {
  X <- .boldzip_validate_matrix(X, "X")
  n_voxels <- nrow(X)
  n_time <- ncol(X)
  if (n_time < 2L) {
    .boldzip_cli_abort("X must contain at least two time points for BOLDZip-SR encoding.")
  }
  k_carriers <- .boldzip_check_scalar_integer(k_carriers, "k_carriers")
  q_texture <- .boldzip_check_scalar_integer(q_texture, "q_texture")
  texture_lags <- sort(unique(as.integer(texture_lags)))
  if (length(texture_lags) == 0L || anyNA(texture_lags)) {
    .boldzip_cli_abort("texture_lags must be a non-empty integer vector.")
  }
  texture_lags <- texture_lags[abs(texture_lags) < n_time]
  if (length(texture_lags) == 0L) {
    .boldzip_cli_abort("texture_lags must contain at least one lag with abs(lag) < n_time.")
  }
  if (is.null(temporal_k)) {
    temporal_k <- max(1L, ceiling(n_time / 4))
  }
  if (is.null(temporal_spec)) {
    temporal_k <- min(.boldzip_check_scalar_integer(temporal_k, "temporal_k"), n_time)
    temporal_info <- list(
      basis = as.matrix(build_dct_basis(n_time, temporal_k, norm = "ortho")),
      temporal_k = temporal_k,
      spec = shared_temporal_spec(
        "dct",
        n_time = n_time,
        rank = temporal_k,
        params = list(norm = "ortho")
      ),
      label = "default_dct"
    )
  } else {
    temporal_info <- .boldzip_materialize_temporal_spec(temporal_spec, n_time)
    temporal_k <- temporal_info$temporal_k
  }

  if (is.null(spatial_basis)) {
    spatial_basis <- boldzip_spatial_basis(label = "identity_detail")
  } else if (!inherits(spatial_basis, "BoldZipSRSpatialBasis")) {
    spatial_basis <- as_boldzip_spatial_basis(spatial_basis)
  }
  phi_c <- .boldzip_validate_basis_matrix(spatial_basis$phi_c, "phi_c", n_voxels)
  phi_d <- .boldzip_validate_basis_matrix(spatial_basis$phi_d, "phi_d", n_voxels)
  .boldzip_validate_orthonormal_columns(phi_c, "phi_c")
  .boldzip_validate_orthonormal_columns(phi_d, "phi_d")
  .boldzip_validate_cross_orthogonal(phi_c, phi_d)

  mu <- if (isTRUE(center)) rowMeans(X) else rep(0, n_voxels)
  X_centered <- sweep(X, 1L, mu, "-")
  pairs <- .boldzip_pair_indices(n_time, reliability$split)

  y_carrier <- if (is.null(phi_c)) X_centered else .boldzip_project(phi_c, X_centered)
  k_use <- min(k_carriers, nrow(y_carrier), ncol(y_carrier))
  sv <- svd(y_carrier, nu = k_use, nv = k_use)
  singular_values <- sv$d[seq_len(k_use)]
  singular_tol <- max(dim(y_carrier)) * .Machine$double.eps *
    max(1, max(abs(y_carrier)))
  singular_values[singular_values <= singular_tol] <- 0
  carrier_u <- sv$u[, seq_len(k_use), drop = FALSE]
  z_raw <- diag(singular_values, nrow = k_use) %*%
    t(sv$v[, seq_len(k_use), drop = FALSE])

  carrier_reliability <- .boldzip_row_cor(
    z_raw[, pairs$first, drop = FALSE],
    z_raw[, pairs$second, drop = FALSE]
  )
  z_raw[carrier_reliability < reliability$min_temporal_reliability, ] <- 0

  temporal_basis <- temporal_info$basis
  theta <- z_raw %*% temporal_basis
  theta_reliability <- matrix(
    rep(carrier_reliability, times = ncol(theta)),
    nrow = nrow(theta),
    ncol = ncol(theta)
  )
  theta[] <- .boldzip_quantize_values(
    as.numeric(theta),
    reliability = as.numeric(theta_reliability),
    quantization = quantization,
    noise_scale = .boldzip_noise_scale(theta)
  )
  z_hat <- theta %*% t(temporal_basis)

  coarse_coef <- if (is.null(phi_c)) {
    matrix(0, nrow = 0L, ncol = n_time)
  } else {
    carrier_u %*% z_hat
  }
  coarse_recon <- if (is.null(phi_c)) {
    matrix(0, nrow = n_voxels, ncol = n_time)
  } else {
    .boldzip_synthesize(phi_c, coarse_coef)
  }

  detail_target <- X_centered - coarse_recon
  y_detail <- .boldzip_project(phi_d, detail_target)
  texture <- .boldzip_fit_sparse_texture(
    y_detail = y_detail,
    z = z_hat,
    pairs = pairs,
    q = q_texture,
    min_reliability = reliability$min_texture_reliability,
    quantization = quantization,
    lags = texture_lags
  )

  y_detail_pred <- .boldzip_predict_texture(texture$loadings, z_hat, nrow(y_detail))
  residual <- y_detail - y_detail_pred
  event_payload <- .boldzip_encode_events(
    residual = residual,
    pairs = pairs,
    split_method = reliability$split,
    events = events,
    quantization = quantization
  )

  out <- list(
    dimensions = c(voxels = n_voxels, time = n_time),
    mu = mu,
    spatial_basis = list(
      phi_c = phi_c,
      phi_d = phi_d,
      label = spatial_basis$label,
      basis_asset = spatial_basis$basis_asset %||% NULL,
      detail_identity = is.null(phi_d),
      coarse_identity = is.null(phi_c)
    ),
    temporal_basis = temporal_basis,
    temporal_spec = temporal_info$spec,
    temporal_label = temporal_info$label,
    carriers = list(
      u = carrier_u,
      theta = theta,
      reliability = carrier_reliability
    ),
    texture = list(
      loadings = texture$loadings,
      matrix = texture$matrix,
      matrix_index = texture$matrix_index
    ),
    events = event_payload,
    settings = list(
      reliability = reliability,
      quantization = quantization,
      events = events,
      k_carriers_requested = k_carriers,
      k_carriers = k_use,
      temporal_k = temporal_k,
      q_texture = q_texture,
      texture_lags = texture_lags,
      center = isTRUE(center)
    ),
    meta = list(
      label = label %||% "boldzip_sr",
      family = "boldzip_sr",
      version = 1L
    )
  )
  class(out) <- "BoldZipSR"
  out
}

#' Decode an experimental BOLDZip-SR object
#'
#' @param object A `BoldZipSR` object.
#' @param time_idx Optional integer time indices to return.
#' @param roi Optional integer or logical row subset to return.
#' @return Reconstructed matrix with rows as voxels/grayordinates and columns
#'   as time points.
#' @export
boldzip_sr_decode <- function(object, time_idx = NULL, roi = NULL) {
  if (!inherits(object, "BoldZipSR")) {
    .boldzip_cli_abort("object must be a BoldZipSR object.")
  }
  n_voxels <- object$dimensions[["voxels"]]
  n_time <- object$dimensions[["time"]]

  z_hat <- object$carriers$theta %*% t(object$temporal_basis)
  phi_c <- object$spatial_basis$phi_c
  phi_d <- object$spatial_basis$phi_d
  coarse <- if (is.null(phi_c)) {
    matrix(0, nrow = n_voxels, ncol = n_time)
  } else {
    phi_c %*% (object$carriers$u %*% z_hat)
  }

  n_detail <- if (is.null(phi_d)) n_voxels else ncol(phi_d)
  event_mat <- .boldzip_events_matrix(object$events, n_detail, n_time)
  y_detail <- .boldzip_predict_texture(object$texture$loadings, z_hat, n_detail) +
    as.matrix(event_mat)
  detail <- .boldzip_synthesize(phi_d, y_detail)
  recon <- coarse + detail
  recon <- sweep(recon, 1L, object$mu, "+")

  if (!is.null(time_idx)) {
    if (!is.numeric(time_idx) || anyNA(time_idx) ||
        any(time_idx != as.integer(time_idx)) ||
        any(time_idx < 1L | time_idx > n_time)) {
      .boldzip_cli_abort("time_idx must contain valid positive integer time indices.")
    }
    recon <- recon[, time_idx, drop = FALSE]
  }
  if (!is.null(roi)) {
    if (!is.numeric(roi) && !is.logical(roi)) {
      .boldzip_cli_abort("roi must be integer or logical.")
    }
    if (is.logical(roi)) {
      if (length(roi) != n_voxels || anyNA(roi)) {
        .boldzip_cli_abort("logical roi must have one non-missing value per row.")
      }
    } else if (anyNA(roi) ||
               any(roi != as.integer(roi)) ||
               any(roi < 1L | roi > n_voxels)) {
      .boldzip_cli_abort("roi must contain valid positive integer row indices.")
    }
    recon <- recon[roi, , drop = FALSE]
  }
  recon
}

#' @method as.matrix BoldZipSR
#' @export
as.matrix.BoldZipSR <- function(x, ...) {
  boldzip_sr_decode(x, ...)
}

.boldzip_roi_from_mask <- function(roi_mask, n_voxels, support = NULL,
                                   mask = NULL, context = "BoldZipSR") {
  if (is.null(roi_mask)) {
    return(NULL)
  }
  if (is.numeric(roi_mask)) {
    return(roi_mask)
  }
  if (!is.logical(roi_mask) || anyNA(roi_mask)) {
    .boldzip_cli_abort(context, " roi_mask must be numeric indices or a non-missing logical mask.")
  }
  if (length(roi_mask) == n_voxels && is.null(dim(roi_mask))) {
    return(roi_mask)
  }
  if (!is.null(mask)) {
    mask_arr <- as.logical(as.array(mask))
    if (!identical(dim(roi_mask), dim(mask_arr))) {
      .boldzip_cli_abort(context, " roi_mask dimensions must match mask dimensions.")
    }
    global_idx <- which(mask_arr)
    roi_global <- which(as.logical(roi_mask))
    return(which(global_idx %in% roi_global))
  }
  if (!is.null(support) && length(roi_mask) == length(support)) {
    return(roi_mask)
  }
  .boldzip_cli_abort(context, " roi_mask must have one logical value per decoded row.")
}

#' Predict from a BOLDZip-SR codec payload
#'
#' @param object A `BoldZipSR` object.
#' @param time_idx Optional integer time indices to return.
#' @param roi Optional integer or logical row subset to return.
#' @param ... Additional arguments passed to [boldzip_sr_decode()].
#' @return Reconstructed matrix with rows as voxels/grayordinates and columns
#'   as time points.
#' @export
predict.BoldZipSR <- function(object, time_idx = NULL, roi = NULL, ...) {
  boldzip_sr_decode(object, time_idx = time_idx, roi = roi, ...)
}

#' Coerce a BOLDZip-SR payload to an ImplicitLatent
#'
#' @param x A `BoldZipSR` object.
#' @param mask Optional logical 3D array or `LogicalNeuroVol` for volumetric
#'   wrapping.
#' @param domain Optional decoded output domain.
#' @param support Optional decoded support. Defaults to row indices.
#' @param ... Additional arguments ignored.
#' @return An `ImplicitLatent` whose decoder returns matrices as
#'   `time x voxels`, matching the rest of the implicit latent API.
#' @export
as_implicit_latent.BoldZipSR <- function(x, mask = NULL, domain = NULL,
                                         support = NULL, ...) {
  n_voxels <- x$dimensions[["voxels"]]
  if (is.null(support)) {
    support <- if (is.null(mask)) seq_len(n_voxels) else mask
  }
  decoder <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL, ...) {
    roi <- .boldzip_roi_from_mask(
      roi_mask,
      n_voxels = n_voxels,
      support = support,
      mask = mask,
      context = "as_implicit_latent.BoldZipSR"
    )
    t(boldzip_sr_decode(x, time_idx = time_idx, roi = roi))
  }
  out <- implicit_latent(
    coeff = x,
    decoder = decoder,
    meta = c(x$meta, list(
      source_class = "BoldZipSR",
      orientation = "time_x_voxels",
      codec_orientation = "voxels_x_time"
    )),
    mask = mask,
    domain = domain,
    support = support
  )
  out$basis_asset <- x$spatial_basis$basis_asset %||% NULL
  out
}
