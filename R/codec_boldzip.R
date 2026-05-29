# Experimental BOLDZip-SR matrix codec

.boldzip_check_scalar_integer <- function(x, name, min = 1L) {
  if (length(x) != 1L || is.na(x) || x != as.integer(x) || x < min) {
    stop(name, " must be a scalar integer >= ", min, ".", call. = FALSE)
  }
  as.integer(x)
}

.boldzip_check_scalar_number <- function(x, name, min = -Inf, max = Inf) {
  if (length(x) != 1L || !is.finite(x) || x < min || x > max) {
    stop(name, " must be a finite scalar in [", min, ", ", max, "].",
         call. = FALSE)
  }
  as.numeric(x)
}

.boldzip_pair_indices <- function(n_time, method = c("odd_even", "halves")) {
  method <- match.arg(method)
  n_time <- .boldzip_check_scalar_integer(n_time, "n_time", min = 2L)
  if (method == "odd_even") {
    n_pair <- floor(n_time / 2L)
    return(list(
      first = seq.int(1L, by = 2L, length.out = n_pair),
      second = seq.int(2L, by = 2L, length.out = n_pair)
    ))
  }

  n_pair <- floor(n_time / 2L)
  list(
    first = seq_len(n_pair),
    second = seq.int(n_pair + 1L, length.out = n_pair)
  )
}

.boldzip_row_cor <- function(a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  if (!identical(dim(a), dim(b))) {
    stop("paired matrices must have identical dimensions.", call. = FALSE)
  }
  a <- sweep(a, 1L, rowMeans(a), "-")
  b <- sweep(b, 1L, rowMeans(b), "-")
  numer <- rowSums(a * b)
  denom <- sqrt(rowSums(a * a) * rowSums(b * b))
  out <- numer / denom
  out[!is.finite(out)] <- 0
  pmax(0, out)
}

.boldzip_validate_matrix <- function(x, name) {
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    stop(name, " must be a numeric matrix.", call. = FALSE)
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  if (nrow(x) < 1L || ncol(x) < 2L) {
    stop(name, " must have at least one row and two columns.", call. = FALSE)
  }
  x
}

.boldzip_validate_basis_matrix <- function(x, name, n_voxels) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    stop(name, " must be a numeric matrix.", call. = FALSE)
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop(name, " must have at least one row and one column.", call. = FALSE)
  }
  if (nrow(x) != n_voxels) {
    stop(name, " must have ", n_voxels, " rows.", call. = FALSE)
  }
  x
}

.boldzip_validate_orthonormal_columns <- function(x, name, tol = 1e-8) {
  if (is.null(x)) {
    return(invisible(TRUE))
  }
  gram <- crossprod(x)
  if (!isTRUE(all.equal(gram, diag(ncol(x)), tolerance = tol))) {
    stop(name, " columns must be orthonormal for BOLDZip-SR spatial coding.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.boldzip_validate_cross_orthogonal <- function(phi_c, phi_d, tol = 1e-8) {
  if (is.null(phi_c) || is.null(phi_d)) {
    return(invisible(TRUE))
  }
  cross <- crossprod(phi_c, phi_d)
  if (!isTRUE(all.equal(cross, matrix(0, nrow = ncol(phi_c), ncol = ncol(phi_d)),
                        tolerance = tol))) {
    stop("phi_c and phi_d must be mutually orthogonal for BOLDZip-SR spatial coding.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.boldzip_validate_spatial_basis_input <- function(x, name) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    stop(name, " must be a numeric matrix.", call. = FALSE)
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop(name, " must have at least one row and one column.", call. = FALSE)
  }
  x
}

.boldzip_project <- function(phi, x) {
  if (is.null(phi)) {
    return(x)
  }
  crossprod(phi, x)
}

.boldzip_synthesize <- function(phi, coef) {
  if (is.null(phi)) {
    return(coef)
  }
  phi %*% coef
}

.boldzip_orthonormalize_columns <- function(x, tol = 1e-10) {
  x <- as.matrix(x)
  if (ncol(x) == 0L) {
    return(x)
  }
  q <- qr(x, tol = tol)
  rank <- q$rank
  if (rank < 1L) {
    stop("basis matrix must have at least one non-zero independent column.",
         call. = FALSE)
  }
  qr.Q(q)[, seq_len(rank), drop = FALSE]
}

.boldzip_materialize_temporal_spec <- function(spec, n_time) {
  validate_basis <- function(basis, label) {
    basis <- .boldzip_validate_basis_matrix(basis, label, n_time)
    gram <- crossprod(basis)
    if (!isTRUE(all.equal(gram, diag(ncol(basis)), tolerance = 1e-8))) {
      stop(label, " columns must be orthonormal for BOLDZip-SR temporal coding.",
           call. = FALSE)
    }
    basis
  }
  if (inherits(spec, "SharedTemporalSpec")) {
    if (spec$n_time != n_time) {
      stop("temporal_spec n_time must match the number of columns in X.",
           call. = FALSE)
    }
    basis <- validate_basis(materialize_shared_temporal_spec(spec), "temporal_spec")
    return(list(
      basis = basis,
      temporal_k = ncol(basis),
      spec = spec,
      label = paste0("shared_", spec$kind)
    ))
  }
  if (inherits(spec, "spec_time_dct")) {
    basis <- validate_basis(
      as.matrix(build_dct_basis(n_time, spec$k, norm = spec$norm)),
      "temporal_spec"
    )
    return(list(
      basis = basis,
      temporal_k = ncol(basis),
      spec = spec,
      label = "spec_time_dct"
    ))
  }
  if (inherits(spec, "spec_time_bspline")) {
    if (!isTRUE(spec$orthonormalize)) {
      stop("spec_time_bspline temporal_spec must use orthonormalize = TRUE for BOLDZip-SR.",
           call. = FALSE)
    }
    basis <- as.matrix(build_bspline_basis(
      n_time = n_time,
      k = spec$k,
      degree = spec$degree,
      include_intercept = spec$include_intercept,
      orthonormalize = spec$orthonormalize
    ))
    basis <- validate_basis(basis, "temporal_spec")
    return(list(
      basis = basis,
      temporal_k = ncol(basis),
      spec = spec,
      label = "spec_time_bspline"
    ))
  }
  if (is.matrix(spec) || inherits(spec, "Matrix")) {
    basis <- validate_basis(spec, "temporal_spec")
    return(list(
      basis = basis,
      temporal_k = ncol(basis),
      spec = shared_temporal_spec("custom", basis = basis),
      label = "matrix"
    ))
  }
  stop("temporal_spec must be NULL, a SharedTemporalSpec, a spec_time_dct, ",
       "a spec_time_bspline, or a numeric matrix.", call. = FALSE)
}

.boldzip_quantize_values <- function(values, reliability, quantization,
                                     noise_scale = NULL) {
  values <- as.numeric(values)
  if (length(values) == 0L) {
    return(values)
  }
  base_step <- quantization$base_step %||% 0
  if (!is.finite(base_step) || base_step <= 0) {
    return(values)
  }
  reliability <- rep(as.numeric(reliability), length.out = length(values))
  noise_scale <- rep(as.numeric(noise_scale %||% stats::sd(values)), length.out = length(values))
  noise_scale[!is.finite(noise_scale) | noise_scale <= 0] <- 1
  step <- base_step * noise_scale / sqrt(pmax(reliability, 0) + quantization$epsilon)
  step <- pmax(step, .Machine$double.eps)
  round(values / step) * step
}

.boldzip_lag_signal <- function(signal, lag) {
  signal <- as.numeric(signal)
  lag <- as.integer(lag)
  n <- length(signal)
  out <- numeric(n)
  if (lag == 0L) {
    return(signal)
  }
  if (abs(lag) >= n) {
    return(out)
  }
  if (lag > 0L) {
    out[(lag + 1L):n] <- signal[seq_len(n - lag)]
  } else {
    lead <- abs(lag)
    out[seq_len(n - lead)] <- signal[(lead + 1L):n]
  }
  out
}

.boldzip_lagged_carrier_bank <- function(z, lags) {
  lags <- sort(unique(as.integer(lags)))
  rows <- vector("list", nrow(z) * length(lags))
  carrier <- integer(length(rows))
  lag <- integer(length(rows))
  pos <- 0L
  for (carrier_idx in seq_len(nrow(z))) {
    for (lag_idx in seq_along(lags)) {
      pos <- pos + 1L
      rows[[pos]] <- .boldzip_lag_signal(z[carrier_idx, ], lags[[lag_idx]])
      carrier[[pos]] <- carrier_idx
      lag[[pos]] <- lags[[lag_idx]]
    }
  }
  list(
    bank = do.call(rbind, rows),
    carrier = carrier,
    lag = lag
  )
}

.boldzip_predict_texture <- function(loadings, z, n_detail) {
  pred <- matrix(0, nrow = n_detail, ncol = ncol(z))
  if (nrow(loadings) == 0L) {
    return(pred)
  }
  for (idx in seq_len(nrow(loadings))) {
    pred[loadings$atom[[idx]], ] <- pred[loadings$atom[[idx]], ] +
      loadings$amplitude[[idx]] *
        .boldzip_lag_signal(z[loadings$carrier[[idx]], ], loadings$lag[[idx]])
  }
  pred
}

.boldzip_fit_sparse_texture <- function(y_detail, z, pairs, q,
                                        min_reliability, quantization,
                                        lags = 0L) {
  n_detail <- nrow(y_detail)
  n_carriers <- nrow(z)
  q <- min(as.integer(q), n_carriers)
  if (q < 1L || n_detail == 0L || n_carriers == 0L) {
    return(list(
      loadings = data.frame(
        atom = integer(), carrier = integer(), amplitude = numeric(),
        lag = integer(), reliability = numeric()
      ),
      matrix = Matrix::Matrix(0, nrow = n_detail, ncol = n_carriers, sparse = TRUE),
      matrix_index = data.frame(
        carrier = seq_len(n_carriers),
        lag = rep.int(0L, n_carriers)
      )
    ))
  }

  lagged <- .boldzip_lagged_carrier_bank(z, lags)
  z_bank <- lagged$bank
  train <- pairs$first
  validate <- pairs$second
  z_train <- z_bank[, train, drop = FALSE]
  z_validate <- z_bank[, validate, drop = FALSE]
  z_norm <- sqrt(rowSums(z_train * z_train))
  z_norm[z_norm <= .Machine$double.eps] <- NA_real_

  atom_out <- integer()
  carrier_out <- integer()
  amplitude_out <- numeric()
  lag_out <- integer()
  reliability_out <- numeric()

  for (idx in seq_len(n_detail)) {
    y_train <- y_detail[idx, train]
    y_validate <- y_detail[idx, validate]
    y_norm <- sqrt(sum(y_train * y_train))
    if (!is.finite(y_norm) || y_norm <= .Machine$double.eps ||
        stats::sd(y_validate) <= .Machine$double.eps) {
      next
    }

    score <- as.numeric(z_train %*% y_train) / (z_norm * y_norm)
    score[!is.finite(score)] <- 0
    active <- integer()
    for (candidate in order(abs(score), decreasing = TRUE)) {
      if (abs(score[[candidate]]) <= 0) {
        break
      }
      if (lagged$carrier[[candidate]] %in% lagged$carrier[active]) {
        next
      }
      active <- c(active, candidate)
      if (length(active) >= q) {
        break
      }
    }
    active <- active[abs(score[active]) > 0]
    if (length(active) == 0L) {
      next
    }

    design <- t(z_train[active, , drop = FALSE])
    coef <- tryCatch(
      qr.solve(design, y_train),
      error = function(e) rep(0, length(active))
    )
    coef[!is.finite(coef)] <- 0
    pred_validate <- as.numeric(coef %*% z_validate[active, , drop = FALSE])
    rel <- .boldzip_row_cor(
      matrix(pred_validate, nrow = 1L),
      matrix(y_validate, nrow = 1L)
    )
    if (rel < min_reliability) {
      next
    }

    coef <- .boldzip_quantize_values(
      coef,
      reliability = rel,
      quantization = quantization,
      noise_scale = stats::sd(y_train - as.numeric(coef %*% z_train[active, , drop = FALSE]))
    )
    keep <- abs(coef) > .Machine$double.eps
    if (!any(keep)) {
      next
    }

    atom_out <- c(atom_out, rep.int(idx, sum(keep)))
    carrier_out <- c(carrier_out, lagged$carrier[active[keep]])
    amplitude_out <- c(amplitude_out, coef[keep])
    lag_out <- c(lag_out, lagged$lag[active[keep]])
    reliability_out <- c(reliability_out, rep.int(rel, sum(keep)))
  }

  matrix_index <- data.frame(
    carrier = lagged$carrier,
    lag = lagged$lag
  )
  if (length(atom_out) == 0L) {
    loading_mat <- Matrix::Matrix(
      0,
      nrow = n_detail,
      ncol = nrow(matrix_index),
      sparse = TRUE
    )
  } else {
    matrix_col <- match(
      paste(carrier_out, lag_out, sep = ":"),
      paste(matrix_index$carrier, matrix_index$lag, sep = ":")
    )
    loading_mat <- Matrix::sparseMatrix(
      i = atom_out,
      j = matrix_col,
      x = amplitude_out,
      dims = c(n_detail, nrow(matrix_index))
    )
  }

  list(
    loadings = data.frame(
      atom = atom_out,
      carrier = carrier_out,
      amplitude = amplitude_out,
      lag = lag_out,
      reliability = reliability_out
    ),
    matrix = loading_mat,
    matrix_index = matrix_index
  )
}

.boldzip_pair_time <- function(t, n_time, method) {
  if (method == "odd_even") {
    if (t %% 2L == 1L && t < n_time) {
      return(t + 1L)
    }
    if (t %% 2L == 0L) {
      return(t - 1L)
    }
    return(NA_integer_)
  }
  half <- floor(n_time / 2L)
  if (t <= half && t + half <= n_time) {
    return(t + half)
  }
  if (t > half && t - half >= 1L) {
    return(t - half)
  }
  NA_integer_
}

.boldzip_encode_events <- function(residual, pairs, split_method, events, quantization) {
  max_events <- events$max_events %||% 0L
  max_events <- as.integer(max_events)
  if (max_events <= 0L || length(residual) == 0L) {
    return(data.frame(
      atom = integer(), time = integer(), duration = integer(),
      amplitude = numeric(), shape_id = character(), reliability = numeric()
    ))
  }

  n_time <- ncol(residual)
  atom <- integer()
  time <- integer()
  amplitude <- numeric()
  reliability <- numeric()

  threshold_sd <- events$threshold_sd %||% 3
  paired_fraction <- events$paired_fraction %||% 0.25

  for (idx in seq_len(nrow(residual))) {
    row <- residual[idx, ]
    scale <- stats::mad(row, center = stats::median(row), constant = 1.4826)
    if (!is.finite(scale) || scale <= .Machine$double.eps) {
      scale <- stats::sd(row)
    }
    if (!is.finite(scale) || scale <= .Machine$double.eps) {
      next
    }
    threshold <- threshold_sd * scale
    candidates <- which(abs(row) >= threshold)
    if (length(candidates) == 0L) {
      next
    }
    for (tt in candidates) {
      pair_t <- .boldzip_pair_time(tt, n_time, split_method)
      if (is.na(pair_t)) {
        next
      }
      if (sign(row[tt]) != sign(row[pair_t])) {
        next
      }
      rel <- min(abs(row[tt]), abs(row[pair_t])) / max(abs(row[tt]), abs(row[pair_t]))
      if (!is.finite(rel) || rel < paired_fraction) {
        next
      }
      atom <- c(atom, idx)
      time <- c(time, tt)
      amplitude <- c(amplitude, row[tt])
      reliability <- c(reliability, rel)
    }
  }

  if (length(atom) == 0L) {
    return(data.frame(
      atom = integer(), time = integer(), duration = integer(),
      amplitude = numeric(), shape_id = character(), reliability = numeric()
    ))
  }

  ord <- order(abs(amplitude) * reliability, decreasing = TRUE)
  ord <- ord[seq_len(min(max_events, length(ord)))]
  amplitude <- .boldzip_quantize_values(
    amplitude[ord],
    reliability = reliability[ord],
    quantization = quantization,
    noise_scale = stats::sd(amplitude)
  )

  data.frame(
    atom = atom[ord],
    time = time[ord],
    duration = rep.int(1L, length(ord)),
    amplitude = amplitude,
    shape_id = rep.int("impulse", length(ord)),
    reliability = reliability[ord]
  )
}

.boldzip_events_matrix <- function(events, n_detail, n_time) {
  if (nrow(events) == 0L) {
    return(Matrix::Matrix(0, nrow = n_detail, ncol = n_time, sparse = TRUE))
  }
  Matrix::sparseMatrix(
    i = events$atom,
    j = events$time,
    x = events$amplitude,
    dims = c(n_detail, n_time)
  )
}

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

#' Construct reliability-aware quantization settings for BOLDZip-SR
#'
#' @param base_step Base quantization step. Set to 0 to disable quantization.
#' @param epsilon Small positive value used in reliability-shaped step sizes.
#' @return A `BoldZipSRQuantization` settings object.
#' @export
boldzip_quantization <- function(base_step = 0, epsilon = 1e-6) {
  base_step <- .boldzip_check_scalar_number(base_step, "base_step", min = 0)
  epsilon <- .boldzip_check_scalar_number(epsilon, "epsilon", min = .Machine$double.eps)
  structure(
    list(base_step = base_step, epsilon = epsilon),
    class = "BoldZipSRQuantization"
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

#' Build a matrix spatial basis descriptor for BOLDZip-SR
#'
#' @param phi_c Optional coarse basis matrix with rows equal to voxels.
#' @param phi_d Optional detail basis matrix with rows equal to voxels. If
#'   `NULL`, the detail basis is the identity basis.
#' @param label Optional label stored in metadata.
#' @param basis_asset Optional source template or shared-basis object used to
#'   build this descriptor.
#' @return A `BoldZipSRSpatialBasis` object.
#' @export
boldzip_spatial_basis <- function(phi_c = NULL, phi_d = NULL, label = NULL,
                                  basis_asset = NULL) {
  if (!is.null(phi_c)) {
    phi_c <- .boldzip_validate_spatial_basis_input(phi_c, "phi_c")
  }
  if (!is.null(phi_d)) {
    phi_d <- .boldzip_validate_spatial_basis_input(phi_d, "phi_d")
  }
  if (!is.null(phi_c) && !is.null(phi_d) && nrow(phi_c) != nrow(phi_d)) {
    stop("phi_c and phi_d must have the same number of rows.", call. = FALSE)
  }
  structure(
    list(
      phi_c = phi_c,
      phi_d = phi_d,
      label = label %||% "matrix",
      basis_asset = basis_asset
    ),
    class = "BoldZipSRSpatialBasis"
  )
}

.boldzip_loadings_from_template <- function(x) {
  tryCatch(template_loadings(x), error = function(e) NULL)
}

#' Coerce a spatial object to a BOLDZip-SR spatial basis
#'
#' @description
#' This helper lets the standalone BOLDZip-SR codec consume matrix-like shared
#' basis assets without registering BOLDZip as an `encode()` family. Matrix and
#' template inputs are used as the detail basis by default and are
#' orthonormalized because BOLDZip projection currently uses the transpose as
#' the analysis operator.
#'
#' @param x A `BoldZipSRSpatialBasis`, matrix-like object, shared reference, or
#'   template object supporting [template_loadings()].
#' @param ... Additional arguments reserved for methods. The default method
#'   accepts `label`, `role`, and `orthonormalize`.
#' @return A `BoldZipSRSpatialBasis` object.
#' @export
as_boldzip_spatial_basis <- function(x, ...) {
  UseMethod("as_boldzip_spatial_basis")
}

#' @export
as_boldzip_spatial_basis.BoldZipSRSpatialBasis <- function(x, ...) {
  x
}

#' @export
as_boldzip_spatial_basis.SharedReference <- function(x, ...) {
  out <- as_boldzip_spatial_basis(resolve_shared_reference(x), ...)
  out$basis_asset <- x
  out
}

#' @export
as_boldzip_spatial_basis.default <- function(x, label = NULL,
                                             role = c("detail", "coarse"),
                                             orthonormalize = TRUE, ...) {
  role <- match.arg(role)
  if (is.list(x) && (any(c("phi_c", "phi_d") %in% names(x)))) {
    args <- x[intersect(names(x), c("phi_c", "phi_d", "label", "basis_asset"))]
    return(do.call(boldzip_spatial_basis, args))
  }

  loadings <- if (is.matrix(x) || inherits(x, "Matrix")) {
    x
  } else {
    .boldzip_loadings_from_template(x)
  }
  if (is.null(loadings)) {
    stop("x must be a BoldZipSRSpatialBasis, matrix-like object, SharedReference, ",
         "or template object with template_loadings().", call. = FALSE)
  }
  phi <- .boldzip_validate_spatial_basis_input(loadings, "template_loadings(x)")
  if (isTRUE(orthonormalize)) {
    phi <- .boldzip_orthonormalize_columns(phi)
  }
  label <- label %||% tryCatch(template_meta(x)$family, error = function(e) NULL) %||%
    class(x)[[1L]]
  if (identical(role, "coarse")) {
    return(boldzip_spatial_basis(phi_c = phi, label = label, basis_asset = x))
  }
  boldzip_spatial_basis(phi_d = phi, label = label, basis_asset = x)
}

.boldzip_canonicalize_eigenvectors <- function(vectors) {
  for (idx in seq_len(ncol(vectors))) {
    pivot <- which.max(abs(vectors[, idx]))
    if (length(pivot) && vectors[pivot, idx] < 0) {
      vectors[, idx] <- -vectors[, idx]
    }
  }
  vectors
}

#' Build a spectral graph spatial basis for BOLDZip-SR
#'
#' @description
#' Constructs an experimental graph-adapted spatial basis from a symmetric
#' adjacency matrix. Low-frequency graph Laplacian eigenvectors become coarse
#' carrier support functions and the following higher-frequency eigenvectors
#' become detail atoms. This is a materialized small-to-medium graph MVP, not a
#' production graph-wavelet operator.
#'
#' @param adjacency Symmetric non-negative adjacency matrix.
#' @param n_coarse Number of low-frequency eigenvectors to use as `phi_c`.
#' @param n_detail Number of following eigenvectors to use as `phi_d`. Defaults
#'   to all remaining eigenvectors.
#' @param normalized Whether to use the normalized graph Laplacian.
#' @param label Optional label stored in metadata.
#' @return A `BoldZipSRSpatialBasis` object.
#' @export
boldzip_graph_spatial_basis <- function(adjacency,
                                        n_coarse = 8L,
                                        n_detail = NULL,
                                        normalized = TRUE,
                                        label = NULL) {
  adjacency <- .boldzip_validate_matrix(adjacency, "adjacency")
  if (nrow(adjacency) != ncol(adjacency)) {
    stop("adjacency must be square.", call. = FALSE)
  }
  if (any(adjacency < 0)) {
    stop("adjacency must be non-negative.", call. = FALSE)
  }
  adjacency <- 0.5 * (adjacency + t(adjacency))
  diag(adjacency) <- 0
  n <- nrow(adjacency)
  n_coarse <- min(.boldzip_check_scalar_integer(n_coarse, "n_coarse", min = 0L), n)
  if (is.null(n_detail)) {
    n_detail <- n - n_coarse
  }
  n_detail <- min(.boldzip_check_scalar_integer(n_detail, "n_detail", min = 0L),
                  n - n_coarse)
  if (n_coarse + n_detail < 1L) {
    stop("n_coarse + n_detail must be at least 1.", call. = FALSE)
  }

  degree <- rowSums(adjacency)
  if (isTRUE(normalized)) {
    inv_sqrt <- ifelse(degree > 0, 1 / sqrt(degree), 0)
    laplacian <- diag(n) - (inv_sqrt * adjacency) * rep(inv_sqrt, each = n)
  } else {
    laplacian <- diag(degree, nrow = n) - adjacency
  }
  laplacian <- 0.5 * (laplacian + t(laplacian))
  eig <- eigen(laplacian, symmetric = TRUE)
  ord <- order(eig$values, decreasing = FALSE)
  vectors <- .boldzip_canonicalize_eigenvectors(eig$vectors[, ord, drop = FALSE])

  phi_c <- if (n_coarse > 0L) {
    vectors[, seq_len(n_coarse), drop = FALSE]
  } else {
    NULL
  }
  detail_start <- n_coarse + 1L
  phi_d <- if (n_detail > 0L) {
    vectors[, seq.int(detail_start, length.out = n_detail), drop = FALSE]
  } else {
    NULL
  }

  out <- boldzip_spatial_basis(
    phi_c = phi_c,
    phi_d = phi_d,
    label = label %||% "spectral_graph"
  )
  out$graph <- list(
    normalized = isTRUE(normalized),
    eigenvalues = eig$values[ord][seq_len(n_coarse + n_detail)],
    n_coarse = n_coarse,
    n_detail = n_detail
  )
  out
}

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
  k_carriers <- .boldzip_check_scalar_integer(k_carriers, "k_carriers")
  q_texture <- .boldzip_check_scalar_integer(q_texture, "q_texture")
  texture_lags <- sort(unique(as.integer(texture_lags)))
  if (length(texture_lags) == 0L || anyNA(texture_lags)) {
    stop("texture_lags must be a non-empty integer vector.", call. = FALSE)
  }
  texture_lags <- texture_lags[abs(texture_lags) < n_time]
  if (length(texture_lags) == 0L) {
    stop("texture_lags must contain at least one lag with abs(lag) < n_time.",
         call. = FALSE)
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
    noise_scale = stats::sd(z_raw)
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
    stop("object must be a BoldZipSR object.", call. = FALSE)
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
      stop("time_idx must contain valid positive integer time indices.",
           call. = FALSE)
    }
    recon <- recon[, time_idx, drop = FALSE]
  }
  if (!is.null(roi)) {
    if (!is.numeric(roi) && !is.logical(roi)) {
      stop("roi must be integer or logical.", call. = FALSE)
    }
    if (is.logical(roi)) {
      if (length(roi) != n_voxels || anyNA(roi)) {
        stop("logical roi must have one non-missing value per row.",
             call. = FALSE)
      }
    } else if (anyNA(roi) ||
               any(roi != as.integer(roi)) ||
               any(roi < 1L | roi > n_voxels)) {
      stop("roi must contain valid positive integer row indices.",
           call. = FALSE)
    }
    recon <- recon[roi, , drop = FALSE]
  }
  recon
}

#' Summarize an experimental BOLDZip-SR payload
#'
#' @param object A `BoldZipSR` object.
#' @return Data frame with component counts and approximate object bytes.
#' @export
boldzip_sr_payload_summary <- function(object) {
  if (!inherits(object, "BoldZipSR")) {
    stop("object must be a BoldZipSR object.", call. = FALSE)
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

.boldzip_corr <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (stats::sd(x) <= .Machine$double.eps || stats::sd(y) <= .Machine$double.eps) {
    return(NA_real_)
  }
  stats::cor(x, y)
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
    stop("X and reconstruction must have identical dimensions.", call. = FALSE)
  }
  err <- X - X_hat
  mse <- mean(err * err)
  weighted_mse <- NA_real_
  if (!is.null(reliability_weights)) {
    w <- as.numeric(reliability_weights)
    w <- rep(w, length.out = length(err))
    if (any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) {
      stop("reliability_weights must be finite non-negative weights with positive sum.",
           call. = FALSE)
    }
    weighted_mse <- sum(w * as.numeric(err)^2) / sum(w)
  }
  c(
    mse = mse,
    rmse = sqrt(mse),
    correlation = .boldzip_corr(X, X_hat),
    reliability_weighted_mse = weighted_mse,
    payload_scalars = if (is_boldzip) {
      summary <- boldzip_sr_payload_summary(object)
      summary$scalar_count[summary$component == "total_object"]
    } else {
      NA_real_
    }
  )
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
    stop(context, " roi_mask must be numeric indices or a non-missing logical mask.",
         call. = FALSE)
  }
  if (length(roi_mask) == n_voxels && is.null(dim(roi_mask))) {
    return(roi_mask)
  }
  if (!is.null(mask)) {
    mask_arr <- as.logical(as.array(mask))
    if (!identical(dim(roi_mask), dim(mask_arr))) {
      stop(context, " roi_mask dimensions must match mask dimensions.",
           call. = FALSE)
    }
    global_idx <- which(mask_arr)
    roi_global <- which(as.logical(roi_mask))
    return(which(global_idx %in% roi_global))
  }
  if (!is.null(support) && length(roi_mask) == length(support)) {
    return(roi_mask)
  }
  stop(context, " roi_mask must have one logical value per decoded row.",
       call. = FALSE)
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
    stop("parcels must have one label per row of X.", call. = FALSE)
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
    stop("fit must be a BoldZipSR object.", call. = FALSE)
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
