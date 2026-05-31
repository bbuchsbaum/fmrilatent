manual_boldzip_lag <- function(signal, lag) {
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

manual_boldzip_texture <- function(loadings, z, n_detail) {
  out <- matrix(0, nrow = n_detail, ncol = ncol(z))
  if (nrow(loadings) == 0L) {
    return(out)
  }
  for (idx in seq_len(nrow(loadings))) {
    out[loadings$atom[[idx]], ] <- out[loadings$atom[[idx]], ] +
      loadings$amplitude[[idx]] *
        manual_boldzip_lag(z[loadings$carrier[[idx]], ], loadings$lag[[idx]])
  }
  out
}

manual_boldzip_events <- function(events, n_detail, n_time) {
  out <- matrix(0, nrow = n_detail, ncol = n_time)
  if (nrow(events) == 0L) {
    return(out)
  }
  for (idx in seq_len(nrow(events))) {
    out[events$atom[[idx]], events$time[[idx]]] <-
      out[events$atom[[idx]], events$time[[idx]]] + events$amplitude[[idx]]
  }
  out
}

manual_boldzip_decode <- function(fit) {
  n_voxels <- fit$dimensions[["voxels"]]
  n_time <- fit$dimensions[["time"]]
  z_hat <- fit$carriers$theta %*% t(fit$temporal_basis)
  phi_c <- fit$spatial_basis$phi_c
  phi_d <- fit$spatial_basis$phi_d
  coarse <- if (is.null(phi_c)) {
    matrix(0, nrow = n_voxels, ncol = n_time)
  } else {
    phi_c %*% (fit$carriers$u %*% z_hat)
  }
  n_detail <- if (is.null(phi_d)) n_voxels else ncol(phi_d)
  detail_coef <- manual_boldzip_texture(fit$texture$loadings, z_hat, n_detail) +
    manual_boldzip_events(fit$events, n_detail, n_time)
  detail <- if (is.null(phi_d)) detail_coef else phi_d %*% detail_coef
  sweep(coarse + detail, 1L, fit$mu, "+")
}

expect_boldzip_spatial_basis_contract <- function(fit, tol = 1e-8) {
  phi_c <- fit$spatial_basis$phi_c
  phi_d <- fit$spatial_basis$phi_d
  if (!is.null(phi_c)) {
    expect_equal(crossprod(phi_c), diag(ncol(phi_c)), tolerance = tol)
  }
  if (!is.null(phi_d)) {
    expect_equal(crossprod(phi_d), diag(ncol(phi_d)), tolerance = tol)
  }
  if (!is.null(phi_c) && !is.null(phi_d)) {
    expect_equal(
      crossprod(phi_c, phi_d),
      matrix(0, nrow = ncol(phi_c), ncol = ncol(phi_d)),
      tolerance = tol
    )
  }
}

expect_boldzip_payload_contract <- function(fit) {
  expect_s3_class(fit, "BoldZipSR")
  n_voxels <- fit$dimensions[["voxels"]]
  n_time <- fit$dimensions[["time"]]
  k <- fit$settings$k_carriers
  n_lags <- length(fit$settings$texture_lags)

  expect_equal(length(fit$mu), n_voxels)
  expect_equal(dim(fit$temporal_basis), c(n_time, fit$settings$temporal_k))
  expect_equal(crossprod(fit$temporal_basis), diag(fit$settings$temporal_k),
               tolerance = 1e-8)
  expect_equal(dim(fit$carriers$theta), c(k, fit$settings$temporal_k))
  expect_equal(length(fit$carriers$reliability), k)
  expect_true(all(is.finite(fit$carriers$theta)))
  expect_true(all(is.finite(fit$carriers$reliability)))
  expect_true(all(fit$carriers$reliability >= 0 & fit$carriers$reliability <= 1))

  if (!is.null(fit$spatial_basis$phi_c)) {
    expect_equal(nrow(fit$spatial_basis$phi_c), n_voxels)
  }
  n_detail <- if (is.null(fit$spatial_basis$phi_d)) {
    n_voxels
  } else {
    expect_equal(nrow(fit$spatial_basis$phi_d), n_voxels)
    ncol(fit$spatial_basis$phi_d)
  }
  expect_boldzip_spatial_basis_contract(fit)

  expect_equal(nrow(fit$texture$matrix), n_detail)
  expect_equal(ncol(fit$texture$matrix), k * n_lags)
  expect_equal(nrow(fit$texture$matrix_index), ncol(fit$texture$matrix))
  expect_equal(sort(unique(fit$texture$matrix_index$lag)),
               sort(fit$settings$texture_lags))
  expect_equal(as.integer(table(fit$texture$matrix_index$carrier)),
               rep(n_lags, k))

  if (nrow(fit$texture$loadings) > 0L) {
    loadings <- fit$texture$loadings
    expect_true(all(loadings$atom >= 1L & loadings$atom <= n_detail))
    expect_true(all(loadings$carrier >= 1L & loadings$carrier <= k))
    expect_true(all(loadings$lag %in% fit$settings$texture_lags))
    expect_true(all(is.finite(loadings$amplitude)))
    expect_true(all(loadings$reliability >=
                      fit$settings$reliability$min_texture_reliability - 1e-12))

    matrix_col <- match(
      paste(loadings$carrier, loadings$lag, sep = ":"),
      paste(fit$texture$matrix_index$carrier, fit$texture$matrix_index$lag,
            sep = ":")
    )
    ref <- Matrix::sparseMatrix(
      i = loadings$atom,
      j = matrix_col,
      x = loadings$amplitude,
      dims = dim(fit$texture$matrix)
    )
    expect_equal(as.matrix(fit$texture$matrix), as.matrix(ref),
                 tolerance = 1e-12)
  }

  if (nrow(fit$events) > 0L) {
    expect_true(all(fit$events$atom >= 1L & fit$events$atom <= n_detail))
    expect_true(all(fit$events$time >= 1L & fit$events$time <= n_time))
    expect_true(all(fit$events$duration >= 1L))
    expect_true(all(is.finite(fit$events$amplitude)))
    expect_true(all(fit$events$reliability >=
                      fit$settings$events$paired_fraction - 1e-12))
    expect_true(all(fit$events$reliability <= 1))
  }

  summary <- boldzip_sr_payload_summary(fit)
  total <- summary$component == "total_object"
  expect_equal(summary$scalar_count[total],
               sum(summary$scalar_count[!total]))
  expect_true(all(is.finite(boldzip_sr_decode(fit))))
}

test_that("BOLDZip-SR payload algebra exactly defines decoder output", {
  set.seed(301)
  sim <- boldzip_sr_simulate(
    n_voxels = 9L,
    n_time = 28L,
    k_carriers = 2L,
    q_texture = 1L,
    n_events = 4L,
    noise_sd = 0.01,
    seed = 301
  )
  fit <- boldzip_sr_encode(
    sim$X,
    k_carriers = 2L,
    temporal_spec = shared_temporal_spec("dct", n_time = 28L, rank = 16L),
    q_texture = 2L,
    texture_lags = -1:1,
    reliability = boldzip_reliability(min_texture_reliability = 0),
    events = boldzip_events(max_events = 8L, threshold_sd = 2.5)
  )

  expect_boldzip_payload_contract(fit)
  manual <- manual_boldzip_decode(fit)
  expect_equal(boldzip_sr_decode(fit), manual, tolerance = 1e-12)
  expect_equal(
    boldzip_sr_decode(fit, time_idx = c(2L, 7L, 11L), roi = c(1L, 5L, 9L)),
    manual[c(1L, 5L, 9L), c(2L, 7L, 11L), drop = FALSE],
    tolerance = 1e-12
  )
})

test_that("BOLDZip-SR metrics match independent reference calculations", {
  set.seed(302)
  X <- matrix(rnorm(7 * 18), nrow = 7L)
  fit <- boldzip_sr_encode(
    X,
    k_carriers = 3L,
    temporal_k = 8L,
    q_texture = 2L,
    events = boldzip_events(max_events = 3L)
  )
  X_hat <- boldzip_sr_decode(fit)
  weights <- matrix(seq(0.2, 1.5, length.out = length(X)), nrow = nrow(X))
  err <- X - X_hat
  ref <- c(
    mse = mean(err * err),
    rmse = sqrt(mean(err * err)),
    correlation = stats::cor(as.numeric(X), as.numeric(X_hat)),
    reliability_weighted_mse =
      sum(as.numeric(weights) * as.numeric(err)^2) / sum(weights)
  )
  metrics <- evaluate_boldzip_sr(X, fit, reliability_weights = weights)

  expect_equal(metrics[names(ref)], ref, tolerance = 1e-12)
  expect_equal(
    evaluate_boldzip_sr(X, X_hat)[c("mse", "rmse", "correlation")],
    metrics[c("mse", "rmse", "correlation")],
    tolerance = 1e-12
  )
})

test_that("BOLDZip-SR reliability weight vectors are row aligned", {
  X <- matrix(seq_len(12), nrow = 3L)
  X_hat <- X
  X_hat[1L, ] <- X_hat[1L, ] + 1
  X_hat[3L, ] <- X_hat[3L, ] - 2
  weights <- c(10, 1, 2)
  err <- X - X_hat
  weight_matrix <- matrix(rep(weights, times = ncol(X)), nrow = nrow(X))

  metrics <- evaluate_boldzip_sr(X, X_hat, reliability_weights = weights)

  expect_equal(
    metrics[["reliability_weighted_mse"]],
    sum(weight_matrix * err^2) / sum(weight_matrix),
    tolerance = 1e-12
  )
})

test_that("BOLDZip-SR advertised low-budget codec preserves clean modeled signal compactly", {
  sim <- boldzip_sr_simulate(
    n_voxels = 60L,
    n_time = 48L,
    k_carriers = 2L,
    q_texture = 1L,
    n_events = 0L,
    noise_sd = 0,
    seed = 303
  )
  fit <- boldzip_sr_encode(
    sim$X,
    k_carriers = 2L,
    temporal_k = 16L,
    q_texture = 1L,
    reliability = boldzip_reliability(min_texture_reliability = 0),
    events = boldzip_events(max_events = 0L)
  )
  metrics <- evaluate_boldzip_sr(sim$X, fit)

  expect_boldzip_payload_contract(fit)
  expect_lt(metrics[["payload_scalars"]], length(sim$X))
  expect_gt(metrics[["correlation"]], 0.98)
  expect_lt(metrics[["rmse"]], 0.15)
})

test_that("BOLDZip-SR clips impossible carrier requests but records the request", {
  set.seed(305)
  X <- matrix(rnorm(3L * 10L), nrow = 3L)
  fit <- boldzip_sr_encode(
    X,
    k_carriers = 99L,
    temporal_k = 6L,
    q_texture = 2L,
    events = boldzip_events(max_events = 0L)
  )

  expect_boldzip_payload_contract(fit)
  expect_equal(fit$settings$k_carriers_requested, 99L)
  expect_equal(fit$settings$k_carriers, 3L)
  expect_equal(dim(fit$carriers$theta), c(3L, 6L))
  expect_equal(length(fit$carriers$reliability), 3L)
})

test_that("BOLDZip-SR is equivariant to baseline shifts and positive scaling", {
  n_vox <- 7L
  n_time <- 30L
  tt <- seq_len(n_time)
  X <- outer(seq(0.4, 1.6, length.out = n_vox),
             sin(2 * pi * tt / n_time) + 0.2 * cos(4 * pi * tt / n_time))
  shift <- seq(-5, 5, length.out = n_vox)
  common <- list(
    k_carriers = 1L,
    temporal_k = n_time,
    q_texture = 1L,
    reliability = boldzip_reliability(min_texture_reliability = 0),
    events = boldzip_events(max_events = 0L)
  )

  fit <- do.call(boldzip_sr_encode, c(list(X = X), common))
  shifted <- do.call(boldzip_sr_encode, c(list(X = sweep(X, 1L, shift, "+")),
                                          common))
  scaled <- do.call(boldzip_sr_encode, c(list(X = 3.25 * X), common))

  expect_equal(
    boldzip_sr_decode(shifted),
    sweep(boldzip_sr_decode(fit), 1L, shift, "+"),
    tolerance = 1e-10
  )
  expect_equal(
    boldzip_sr_decode(scaled),
    3.25 * boldzip_sr_decode(fit),
    tolerance = 1e-10
  )
})

test_that("BOLDZip-SR temporal reliability threshold zeroes unstable carriers", {
  n_time <- 20L
  z <- rep(c(1, -1), length.out = n_time)
  phi_c <- diag(3L)[, 1L, drop = FALSE]
  X <- phi_c %*% matrix(z, nrow = 1L)
  basis <- boldzip_spatial_basis(phi_c = phi_c)
  common <- list(
    spatial_basis = basis,
    k_carriers = 1L,
    temporal_k = n_time,
    q_texture = 1L,
    center = FALSE,
    events = boldzip_events(max_events = 0L)
  )

  permissive <- do.call(
    boldzip_sr_encode,
    c(list(X = X),
      common,
      list(reliability = boldzip_reliability(min_temporal_reliability = 0)))
  )
  strict <- do.call(
    boldzip_sr_encode,
    c(list(X = X),
      common,
      list(reliability = boldzip_reliability(min_temporal_reliability = 0.01)))
  )

  expect_boldzip_payload_contract(permissive)
  expect_boldzip_payload_contract(strict)
  expect_lt(max(strict$carriers$reliability), 0.01)
  expect_gt(sum(abs(permissive$carriers$theta)), 0)
  expect_equal(strict$carriers$theta, matrix(0, nrow = 1L, ncol = n_time),
               tolerance = 1e-12)
  expect_lt(
    evaluate_boldzip_sr(X, permissive)[["mse"]],
    evaluate_boldzip_sr(X, strict)[["mse"]]
  )
})

test_that("BOLDZip-SR texture sparsity budget recovers mixed detail atoms", {
  n_time <- 24L
  tt <- seq_len(n_time)
  z1 <- sin(2 * pi * tt / n_time)
  z2 <- cos(2 * pi * tt / n_time)
  X <- rbind(z1, z2, z1 + 0.5 * z2)
  basis <- boldzip_spatial_basis(
    phi_c = diag(3L)[, 1:2, drop = FALSE],
    phi_d = diag(3L)[, 3L, drop = FALSE]
  )
  common <- list(
    spatial_basis = basis,
    k_carriers = 2L,
    temporal_k = n_time,
    center = FALSE,
    reliability = boldzip_reliability(min_texture_reliability = 0),
    events = boldzip_events(max_events = 0L)
  )

  one_link <- do.call(boldzip_sr_encode, c(list(X = X), common,
                                           list(q_texture = 1L)))
  two_links <- do.call(boldzip_sr_encode, c(list(X = X), common,
                                            list(q_texture = 2L)))

  expect_boldzip_payload_contract(one_link)
  expect_boldzip_payload_contract(two_links)
  expect_equal(nrow(one_link$texture$loadings), 1L)
  expect_equal(nrow(two_links$texture$loadings), 2L)
  expect_gt(evaluate_boldzip_sr(X, one_link)[["mse"]], 1e-3)
  expect_lt(evaluate_boldzip_sr(X, two_links)[["mse"]], 1e-20)
})

test_that("BOLDZip-SR fuzzed small payloads satisfy shape and index invariants", {
  for (seed in 401:408) {
    set.seed(seed)
    n_vox <- sample(4:9, 1L)
    n_time <- sample(seq(10L, 18L, by = 2L), 1L)
    X <- matrix(rnorm(n_vox * n_time), nrow = n_vox)
    k <- sample(1:min(3L, n_vox), 1L)
    temporal_k <- sample(2:min(n_time, 8L), 1L)
    lags <- sort(sample((-2):2, sample(1:3, 1L)))
    fit <- boldzip_sr_encode(
      X,
      k_carriers = k,
      temporal_spec = shared_temporal_spec("dct", n_time = n_time,
                                           rank = temporal_k),
      q_texture = sample(1:2, 1L),
      texture_lags = lags,
      quantization = boldzip_quantization(base_step = sample(c(0, 0.05), 1L)),
      events = boldzip_events(max_events = sample(0:4, 1L)),
      center = sample(c(TRUE, FALSE), 1L)
    )

    expect_boldzip_payload_contract(fit)
    full <- boldzip_sr_decode(fit)
    roi <- sort(sample(seq_len(n_vox), min(3L, n_vox)))
    t_idx <- sort(sample(seq_len(n_time), min(4L, n_time)))
    expect_equal(
      boldzip_sr_decode(fit, time_idx = t_idx, roi = roi),
      full[roi, t_idx, drop = FALSE],
      tolerance = 1e-12
    )
  }
})

test_that("BOLDZip-SR handles degenerate constant data as baseline-only signal", {
  baseline <- seq(-2, 2, length.out = 6L)
  X <- matrix(rep(baseline, times = 12L), nrow = 6L)
  fit <- boldzip_sr_encode(
    X,
    k_carriers = 3L,
    temporal_k = 4L,
    q_texture = 2L,
    events = boldzip_events(max_events = 4L)
  )

  expect_boldzip_payload_contract(fit)
  expect_equal(boldzip_sr_decode(fit), X, tolerance = 1e-12)
  expect_true(all(fit$carriers$reliability == 0))
  expect_equal(nrow(fit$texture$loadings), 0L)
  expect_equal(nrow(fit$events), 0L)
})

test_that("BOLDZip-SR rejects non-orthonormal temporal specs explicitly", {
  X <- matrix(rnorm(5 * 12), nrow = 5L)
  expect_error(
    boldzip_sr_encode(
      X,
      temporal_spec = matrix(1, nrow = 12L, ncol = 2L)
    ),
    "orthonormal"
  )
  expect_error(
    boldzip_sr_encode(
      X,
      temporal_spec = spec_time_bspline(k = 5L, orthonormalize = FALSE)
    ),
    "orthonormalize = TRUE"
  )
})

test_that("BOLDZip-SR rejects supplied spatial bases that violate projection contracts", {
  X4 <- matrix(rnorm(4L * 12L), nrow = 4L)
  expect_error(
    boldzip_sr_encode(
      X4,
      spatial_basis = boldzip_spatial_basis(
        phi_d = matrix(c(1, 1, 0, 0), nrow = 4L)
      )
    ),
    "phi_d columns must be orthonormal"
  )

  X3 <- matrix(rnorm(3L * 12L), nrow = 3L)
  expect_error(
    boldzip_sr_encode(
      X3,
      spatial_basis = boldzip_spatial_basis(
        phi_c = diag(3L)[, 1L, drop = FALSE],
        phi_d = diag(3L)[, 1L, drop = FALSE]
      )
    ),
    "mutually orthogonal"
  )
})

test_that("BOLDZip-SR spatial coercion orthonormalizes near-collinear template inputs", {
  raw <- cbind(
    c(1, 1, 0, 0),
    c(1, 1, 0, 0) + 1e-13,
    c(0, 0, 1, 1)
  )
  expect_warning(
    basis <- as_boldzip_spatial_basis(raw),
    "dropped 1 linearly dependent column"
  )

  expect_s3_class(basis, "BoldZipSRSpatialBasis")
  expect_equal(nrow(basis$phi_d), 4L)
  expect_equal(ncol(basis$phi_d), 2L)
  expect_equal(crossprod(basis$phi_d), diag(2L), tolerance = 1e-10)

  z <- sin(2 * pi * seq_len(16L) / 16L)
  X <- basis$phi_d[, 1L, drop = FALSE] %*% matrix(z, nrow = 1L)
  fit <- boldzip_sr_encode(
    X,
    spatial_basis = basis,
    k_carriers = 1L,
    temporal_k = 16L,
    q_texture = 1L,
    center = FALSE,
    events = boldzip_events(max_events = 0L)
  )
  expect_lt(evaluate_boldzip_sr(X, fit)[["mse"]], 1e-10)
})

test_that("BOLDZip-SR event encoding requires paired split support", {
  residual <- matrix(0, nrow = 2L, ncol = 12L)
  residual[1L, c(2L, 8L)] <- c(4, 3)
  residual[2L, c(3L, 9L)] <- c(5, -5)
  pairs <- fmrilatent:::.boldzip_pair_indices(12L, method = "halves")

  events <- fmrilatent:::.boldzip_encode_events(
    residual = residual,
    pairs = pairs,
    split_method = "halves",
    events = boldzip_events(
      max_events = 10L,
      threshold_sd = 1,
      paired_fraction = 0.5
    ),
    quantization = boldzip_quantization(base_step = 0)
  )

  expect_equal(events$atom, c(1L, 1L))
  expect_setequal(events$time, c(2L, 8L))
  expect_equal(events$reliability, c(0.75, 0.75), tolerance = 1e-12)
  expect_true(all(sign(events$amplitude) > 0))
})

test_that("BOLDZip-SR half split pairing drops odd trailing time points", {
  expect_equal(fmrilatent:::.boldzip_pair_time(1L, 5L, "halves"), 3L)
  expect_equal(fmrilatent:::.boldzip_pair_time(3L, 5L, "halves"), 1L)
  expect_equal(fmrilatent:::.boldzip_pair_time(5L, 5L, "halves"), NA_integer_)
  expect_equal(fmrilatent:::.boldzip_pair_time(4L, 4L, "halves"), 2L)
})

test_that("evaluate_boldzip_sr rejects invalid reliability weights", {
  X <- matrix(rnorm(4L * 10L), nrow = 4L)

  expect_error(
    evaluate_boldzip_sr(X, X, reliability_weights = rep(1, ncol(X))),
    "length must match nrow\\(X\\)"
  )
  expect_error(
    evaluate_boldzip_sr(X, X, reliability_weights = matrix(1, nrow = nrow(X), ncol = ncol(X) - 1L)),
    "matrix dimensions must match X"
  )
  expect_error(
    evaluate_boldzip_sr(X, X, reliability_weights = c(1, -1, 1, 1)),
    "finite non-negative"
  )
  expect_error(
    evaluate_boldzip_sr(X, X, reliability_weights = rep(0, nrow(X))),
    "positive sum"
  )
  expect_error(
    evaluate_boldzip_sr(X, X, reliability_weights = c(1, Inf, 1, 1)),
    "finite non-negative"
  )
})

test_that("BOLDZip-SR representative encode/decode stays within a smoke runtime envelope", {
  skip_on_cran()
  set.seed(304)
  sim <- boldzip_sr_simulate(
    n_voxels = 35L,
    n_time = 40L,
    k_carriers = 3L,
    q_texture = 1L,
    n_events = 4L,
    noise_sd = 0.02,
    seed = 304
  )
  elapsed <- system.time({
    fit <- boldzip_sr_encode(
      sim$X,
      k_carriers = 3L,
      temporal_k = 12L,
      q_texture = 2L,
      events = boldzip_events(max_events = 8L)
    )
    X_hat <- boldzip_sr_decode(fit)
  })[["elapsed"]]

  expect_equal(dim(X_hat), dim(sim$X))
  expect_true(all(is.finite(X_hat)))
  expect_lt(elapsed, 5)
})
