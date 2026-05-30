# BOLDZip-SR codec: temporal-spec materialization, carrier bank, texture, events, pairing
#
# Split out of R/codec_boldzip.R (see that file's header for the full map).
# Functions in this file:
#   .boldzip_pair_indices
#   .boldzip_row_cor
#   .boldzip_materialize_temporal_spec
#   .boldzip_lag_signal
#   .boldzip_lagged_carrier_bank
#   .boldzip_predict_texture
#   .boldzip_fit_sparse_texture
#   .boldzip_pair_time
#   .boldzip_encode_events
#   .boldzip_events_matrix


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
      noise_scale = .boldzip_noise_scale(y_train - as.numeric(coef %*% z_train[active, , drop = FALSE]))
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
    noise_scale = .boldzip_noise_scale(amplitude)
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
