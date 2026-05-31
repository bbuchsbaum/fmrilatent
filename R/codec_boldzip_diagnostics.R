# BOLDZip-SR codec: evaluation, simulation, baselines, and settings constructors
#
# Split out of R/codec_boldzip.R (see that file's header for the full map).
# Functions in this file:
#   boldzip_reliability
#   boldzip_events
#   boldzip_sr_payload_summary
#   evaluate_boldzip_sr
#   boldzip_sr_simulate
#   boldzip_parcel_reconstruct
#   boldzip_svd_reconstruct
#   compare_boldzip_sr


#' Construct split-half reliability settings for BOLDZip-SR
#'
#' @param split Split strategy. `"odd_even"` pairs adjacent odd/even time
#'   points; `"halves"` pairs the first and second half of the run.
#' @param min_texture_reliability Minimum held-out reliability required to keep
#'   a fine-detail texture loading.
#' @param min_temporal_reliability Minimum carrier reliability required before
#'   temporal coefficients are kept. Carriers below this threshold are zeroed.
#' @return A `BoldZipSRReliability` settings object.
#' @export
boldzip_reliability <- function(split = c("odd_even", "halves"),
                                min_texture_reliability = 0,
                                min_temporal_reliability = 0) {
  split <- match.arg(split)
  min_texture_reliability <- .boldzip_check_scalar_number(
    min_texture_reliability, "min_texture_reliability", min = 0, max = 1)
  min_temporal_reliability <- .boldzip_check_scalar_number(
    min_temporal_reliability, "min_temporal_reliability", min = 0, max = 1)
  structure(
    list(
      split = split,
      min_texture_reliability = min_texture_reliability,
      min_temporal_reliability = min_temporal_reliability
    ),
    class = "BoldZipSRReliability"
  )
}

#' Construct sparse innovation event settings for BOLDZip-SR
#'
#' @param max_events Maximum number of residual events to store.
#' @param threshold_sd Residual event threshold in robust standard deviations.
#' @param paired_fraction Minimum paired split amplitude agreement.
#' @return A `BoldZipSREvents` settings object.
#' @export
boldzip_events <- function(max_events = 256L, threshold_sd = 3,
                           paired_fraction = 0.25) {
  max_events <- .boldzip_check_scalar_integer(max_events, "max_events", min = 0L)
  threshold_sd <- .boldzip_check_scalar_number(threshold_sd, "threshold_sd", min = 0)
  paired_fraction <- .boldzip_check_scalar_number(
    paired_fraction, "paired_fraction", min = 0, max = 1)
  structure(
    list(
      max_events = max_events,
      threshold_sd = threshold_sd,
      paired_fraction = paired_fraction
    ),
    class = "BoldZipSREvents"
  )
}

#' Summarize an experimental BOLDZip-SR payload
#'
#' @param object A `BoldZipSR` object.
#' @return Data frame with component counts and approximate object bytes.
#' @export
boldzip_sr_payload_summary <- function(object) {
  if (!inherits(object, "BoldZipSR")) {
    .boldzip_cli_abort("object must be a BoldZipSR object.")
  }
  n_load <- nrow(object$texture$loadings)
  n_events <- nrow(object$events)
  n_texture_index <- if (is.null(object$texture$matrix_index)) {
    0L
  } else {
    nrow(object$texture$matrix_index) * 2L
  }
  rows <- data.frame(
    component = c("carriers_theta", "texture_loadings", "events", "baseline_mu",
                  "basis_metadata"),
    scalar_count = c(
      length(object$carriers$theta),
      n_load * 5L + n_texture_index,
      n_events * 5L,
      length(object$mu),
      length(object$temporal_basis) +
        length(object$spatial_basis$phi_c %||% numeric(0)) +
        length(object$spatial_basis$phi_d %||% numeric(0))
    ),
    bytes = c(
      as.numeric(utils::object.size(object$carriers$theta)),
      as.numeric(utils::object.size(object$texture$loadings)),
      as.numeric(utils::object.size(object$events)),
      as.numeric(utils::object.size(object$mu)),
      as.numeric(utils::object.size(object$temporal_basis)) +
        as.numeric(utils::object.size(object$spatial_basis))
    )
  )
  rows <- rbind(
    rows,
    data.frame(
      component = "total_object",
      scalar_count = sum(rows$scalar_count),
      bytes = as.numeric(utils::object.size(object))
    )
  )
  rownames(rows) <- NULL
  rows
}

#' Evaluate BOLDZip-SR reconstruction quality
#'
#' @param X Original matrix with rows as voxels/grayordinates and columns as
#'   time points.
#' @param object A `BoldZipSR` object or a reconstructed matrix.
#' @param reliability_weights Optional matrix or vector of reliability weights.
#' @return Named numeric vector with reconstruction metrics.
#' @export
evaluate_boldzip_sr <- function(X, object, reliability_weights = NULL) {
  X <- .boldzip_validate_matrix(X, "X")
  is_boldzip <- inherits(object, "BoldZipSR")
  X_hat <- if (inherits(object, "BoldZipSR")) {
    boldzip_sr_decode(object)
  } else {
    .boldzip_validate_matrix(object, "object")
  }
  if (!identical(dim(X), dim(X_hat))) {
    .boldzip_cli_abort("X and reconstruction must have identical dimensions.")
  }
  err <- X - X_hat
  mse <- mean(err * err)
  weighted_mse <- NA_real_
  if (!is.null(reliability_weights)) {
    if (is.null(dim(reliability_weights))) {
      w <- as.numeric(reliability_weights)
      if (length(w) != nrow(X)) {
        .boldzip_cli_abort("reliability_weights vector length must match nrow(X).")
      }
      w <- matrix(rep(w, times = ncol(X)), nrow = nrow(X), ncol = ncol(X))
    } else {
      w <- as.matrix(reliability_weights)
      if (!identical(dim(w), dim(X))) {
        .boldzip_cli_abort("reliability_weights matrix dimensions must match X.")
      }
    }
    if (any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) {
      .boldzip_cli_abort("reliability_weights must be finite non-negative weights with positive sum.")
    }
    weighted_mse <- sum(w * as.numeric(err)^2) / sum(w)
  }
  c(
    mse = mse,
    rmse = sqrt(mse),
    correlation = .boldzip_corr(X, X_hat),
    reliability_weighted_mse = weighted_mse,
    payload_scalars = if (is_boldzip) {
      psum <- boldzip_sr_payload_summary(object)
      psum$scalar_count[psum$component == "total_object"]
    } else {
      NA_real_
    }
  )
}

#' Simulate data with BOLDZip-SR carrier, texture, and event structure
#'
#' @param n_voxels Number of spatial rows.
#' @param n_time Number of time points.
#' @param k_carriers Number of latent carrier time courses.
#' @param q_texture Number of non-zero carrier loadings per voxel.
#' @param n_events Number of paired impulse events to add.
#' @param noise_sd Independent Gaussian noise standard deviation.
#' @param seed Optional random seed.
#' @return List containing `X`, `mu`, `carriers`, `texture`, and `events`.
#' @export
boldzip_sr_simulate <- function(n_voxels = 40L, n_time = 80L,
                                k_carriers = 3L, q_texture = 1L,
                                n_events = 8L, noise_sd = 0.05,
                                seed = NULL) {
  n_voxels <- .boldzip_check_scalar_integer(n_voxels, "n_voxels")
  n_time <- .boldzip_check_scalar_integer(n_time, "n_time", min = 4L)
  k_carriers <- .boldzip_check_scalar_integer(k_carriers, "k_carriers")
  q_texture <- .boldzip_check_scalar_integer(q_texture, "q_texture")
  n_events <- .boldzip_check_scalar_integer(n_events, "n_events", min = 0L)
  noise_sd <- .boldzip_check_scalar_number(noise_sd, "noise_sd", min = 0)
  if (!is.null(seed)) {
    seed <- .boldzip_check_scalar_integer(seed, "seed", min = 0L)
    set.seed(seed)
  }

  k_use <- min(k_carriers, n_time)
  q_use <- min(q_texture, k_use)
  tt <- seq_len(n_time)
  carrier_bank <- vapply(seq_len(k_use), function(k) {
    sin(2 * pi * k * tt / n_time) + 0.5 * cos(2 * pi * (k + 1L) * tt / n_time)
  }, numeric(n_time))
  carriers <- t(qr.Q(qr(carrier_bank))[, seq_len(k_use), drop = FALSE])

  texture <- matrix(0, nrow = n_voxels, ncol = k_use)
  for (idx in seq_len(n_voxels)) {
    active <- sample.int(k_use, q_use)
    texture[idx, active] <- stats::rnorm(q_use, sd = 1)
  }

  mu <- stats::rnorm(n_voxels, sd = 0.5)
  X <- texture %*% carriers
  events <- data.frame(
    atom = integer(), time = integer(), duration = integer(),
    amplitude = numeric(), shape_id = character(), reliability = numeric()
  )
  if (n_events > 0L) {
    event_atoms <- sample.int(n_voxels, n_events, replace = TRUE)
    odd_times <- seq.int(1L, n_time - 1L, by = 2L)
    event_times <- sample(odd_times, n_events, replace = TRUE)
    amps <- stats::rnorm(n_events, mean = 0, sd = 3)
    for (idx in seq_len(n_events)) {
      atom <- event_atoms[[idx]]
      time <- event_times[[idx]]
      amp <- amps[[idx]]
      X[atom, c(time, time + 1L)] <- X[atom, c(time, time + 1L)] + amp
      events <- rbind(
        events,
        data.frame(
          atom = rep.int(atom, 2L),
          time = c(time, time + 1L),
          duration = c(1L, 1L),
          amplitude = c(amp, amp),
          shape_id = c("impulse", "impulse"),
          reliability = c(1, 1)
        )
      )
    }
  }

  X <- sweep(X, 1L, mu, "+")
  if (noise_sd > 0) {
    X <- X + matrix(stats::rnorm(n_voxels * n_time, sd = noise_sd), nrow = n_voxels)
  }

  list(
    X = X,
    mu = mu,
    carriers = carriers,
    texture = texture,
    events = events
  )
}

#' Reconstruct a matrix from parcel-average time series
#'
#' @param X Numeric matrix with dimensions `voxels x time`.
#' @param parcels Integer, factor, or character vector with one parcel label per
#'   row of `X`.
#' @return Matrix with the same dimensions as `X`, expanded from parcel means.
#' @export
boldzip_parcel_reconstruct <- function(X, parcels) {
  X <- .boldzip_validate_matrix(X, "X")
  if (length(parcels) != nrow(X)) {
    .boldzip_cli_abort("parcels must have one label per row of X.")
  }
  parcels <- as.factor(parcels)
  out <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  for (level in levels(parcels)) {
    idx <- which(parcels == level)
    out[idx, ] <- rep(colMeans(X[idx, , drop = FALSE]), each = length(idx))
  }
  out
}

#' Reconstruct a matrix from a low-rank SVD baseline
#'
#' @param X Numeric matrix with dimensions `voxels x time`.
#' @param rank SVD rank.
#' @param center Whether to remove and restore row means.
#' @return Matrix with the same dimensions as `X`.
#' @export
boldzip_svd_reconstruct <- function(X, rank, center = TRUE) {
  X <- .boldzip_validate_matrix(X, "X")
  rank <- min(.boldzip_check_scalar_integer(rank, "rank"), min(dim(X)))
  mu <- if (isTRUE(center)) rowMeans(X) else rep(0, nrow(X))
  X_centered <- sweep(X, 1L, mu, "-")
  sv <- svd(X_centered, nu = rank, nv = rank)
  recon <- sv$u[, seq_len(rank), drop = FALSE] %*%
    diag(sv$d[seq_len(rank)], nrow = rank) %*%
    t(sv$v[, seq_len(rank), drop = FALSE])
  sweep(recon, 1L, mu, "+")
}

#' Compare BOLDZip-SR against simple reconstruction baselines
#'
#' @param X Original matrix with rows as voxels/grayordinates and columns as
#'   time points.
#' @param fit A `BoldZipSR` object.
#' @param parcels Optional parcel labels for a parcel-average baseline.
#' @param svd_ranks Optional integer vector of SVD ranks.
#' @return Data frame with method, MSE, RMSE, correlation, and scalar payload
#'   estimates.
#' @export
compare_boldzip_sr <- function(X, fit, parcels = NULL, svd_ranks = NULL) {
  X <- .boldzip_validate_matrix(X, "X")
  if (!inherits(fit, "BoldZipSR")) {
    .boldzip_cli_abort("fit must be a BoldZipSR object.")
  }

  rows <- list()
  metric <- evaluate_boldzip_sr(X, fit)
  rows[[length(rows) + 1L]] <- data.frame(
    method = "boldzip_sr",
    parameter = NA_character_,
    mse = metric[["mse"]],
    rmse = metric[["rmse"]],
    correlation = metric[["correlation"]],
    payload_scalars = metric[["payload_scalars"]]
  )

  if (!is.null(parcels)) {
    parcel_hat <- boldzip_parcel_reconstruct(X, parcels)
    metric <- evaluate_boldzip_sr(X, parcel_hat)
    rows[[length(rows) + 1L]] <- data.frame(
      method = "parcel",
      parameter = as.character(length(unique(parcels))),
      mse = metric[["mse"]],
      rmse = metric[["rmse"]],
      correlation = metric[["correlation"]],
      payload_scalars = length(unique(parcels)) * ncol(X)
    )
  }

  if (!is.null(svd_ranks)) {
    for (rank in as.integer(svd_ranks)) {
      svd_hat <- boldzip_svd_reconstruct(X, rank = rank)
      metric <- evaluate_boldzip_sr(X, svd_hat)
      rows[[length(rows) + 1L]] <- data.frame(
        method = "svd",
        parameter = as.character(rank),
        mse = metric[["mse"]],
        rmse = metric[["rmse"]],
        correlation = metric[["correlation"]],
        payload_scalars = rank * (nrow(X) + ncol(X) + 1L) + nrow(X)
      )
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
