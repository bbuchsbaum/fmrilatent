test_that("BOLDZip-SR reconstructs a deterministic rank-one texture model", {
  n_vox <- 8L
  n_time <- 24L
  tt <- seq_len(n_time)
  carrier <- sin(2 * pi * tt / n_time)
  texture <- seq(0.25, 2, length.out = n_vox)
  mu <- seq(-1, 1, length.out = n_vox)
  X <- texture %*% t(carrier)
  X <- sweep(X, 1L, mu, "+")

  fit <- boldzip_sr_encode(
    X,
    k_carriers = 1L,
    temporal_k = n_time,
    q_texture = 1L,
    reliability = boldzip_reliability(min_texture_reliability = 0.1),
    events = boldzip_events(max_events = 0L)
  )
  X_hat <- boldzip_sr_decode(fit)

  expect_s3_class(fit, "BoldZipSR")
  expect_equal(dim(X_hat), dim(X))
  expect_lt(mean((X - X_hat)^2), 1e-10)
  expect_equal(nrow(fit$texture$loadings), n_vox)
  expect_true(all(is.finite(fit$carriers$reliability)))
})

test_that("BOLDZip-SR payload summary accounts for all core components", {
  set.seed(10)
  X <- matrix(rnorm(6 * 18), nrow = 6L)
  fit <- boldzip_sr_encode(
    X,
    k_carriers = 2L,
    temporal_k = 6L,
    q_texture = 1L,
    events = boldzip_events(max_events = 3L)
  )
  summary <- boldzip_sr_payload_summary(fit)

  expect_true(all(c("carriers_theta", "texture_loadings", "events",
                    "baseline_mu", "basis_metadata", "total_object") %in%
                    summary$component))
  expect_true(all(summary$scalar_count >= 0))
  expect_true(all(summary$bytes > 0))
  expect_equal(
    summary$scalar_count[summary$component == "carriers_theta"],
    length(fit$carriers$theta)
  )
})

test_that("BOLDZip-SR residual events improve paired impulse reconstruction", {
  n_vox <- 5L
  n_time <- 30L
  tt <- seq_len(n_time)
  carrier <- sin(2 * pi * tt / n_time)
  X <- seq(0.5, 1.5, length.out = n_vox) %*% t(carrier)
  X[3L, c(13L, 14L)] <- X[3L, c(13L, 14L)] + 8

  common <- list(
    X = X,
    k_carriers = 1L,
    temporal_k = 4L,
    q_texture = 1L,
    reliability = boldzip_reliability(min_texture_reliability = 0),
    center = FALSE
  )
  no_events <- do.call(
    boldzip_sr_encode,
    c(common, list(events = boldzip_events(max_events = 0L)))
  )
  with_events <- do.call(
    boldzip_sr_encode,
    c(common, list(events = boldzip_events(max_events = 8L, threshold_sd = 2)))
  )

  expect_gt(nrow(with_events$events), 0)
  expect_lt(
    evaluate_boldzip_sr(X, with_events)[["mse"]],
    evaluate_boldzip_sr(X, no_events)[["mse"]]
  )
})

test_that("BOLDZip-SR can recover lagged texture loadings", {
  n_time <- 32L
  carrier <- sin(2 * pi * seq_len(n_time) / n_time)
  lagged <- c(0, 0, carrier[seq_len(n_time - 2L)])
  X <- rbind(carrier, 0.75 * lagged)
  basis <- boldzip_spatial_basis(
    phi_c = matrix(c(1, 0), nrow = 2L),
    phi_d = matrix(c(0, 1), nrow = 2L)
  )

  fit <- boldzip_sr_encode(
    X,
    spatial_basis = basis,
    k_carriers = 1L,
    temporal_k = n_time,
    q_texture = 1L,
    texture_lags = 0:3,
    center = FALSE,
    events = boldzip_events(max_events = 0L)
  )

  expect_true(any(fit$texture$loadings$lag == 2L))
  expect_lt(evaluate_boldzip_sr(X, fit)[["mse"]], 1e-10)
})

test_that("BOLDZip-SR supports supplied coarse and detail bases", {
  n_vox <- 6L
  n_time <- 20L
  phi_c <- matrix(1 / sqrt(n_vox), nrow = n_vox, ncol = 1L)
  phi_d <- qr.Q(qr(cbind(phi_c, diag(n_vox))))[, -1L, drop = FALSE]
  z <- sin(2 * pi * seq_len(n_time) / n_time)
  X <- phi_c %*% matrix(3 * z, nrow = 1L) +
    phi_d[, 1L, drop = FALSE] %*% matrix(0.8 * z, nrow = 1L)

  fit <- boldzip_sr_encode(
    X,
    spatial_basis = boldzip_spatial_basis(phi_c = phi_c, phi_d = phi_d,
                                          label = "orthogonal_toy"),
    k_carriers = 1L,
    temporal_k = n_time,
    q_texture = 1L,
    events = boldzip_events(max_events = 0L),
    center = FALSE
  )
  X_hat <- boldzip_sr_decode(fit)

  expect_equal(fit$spatial_basis$label, "orthogonal_toy")
  expect_equal(dim(X_hat), dim(X))
  expect_lt(mean((X - X_hat)^2), 1e-10)
})

test_that("BOLDZip-SR builds and uses a spectral graph spatial basis", {
  adjacency <- matrix(0, nrow = 6L, ncol = 6L)
  for (idx in seq_len(5L)) {
    adjacency[idx, idx + 1L] <- 1
    adjacency[idx + 1L, idx] <- 1
  }
  basis <- boldzip_graph_spatial_basis(
    adjacency,
    n_coarse = 2L,
    n_detail = 3L,
    normalized = TRUE
  )
  gram <- crossprod(cbind(basis$phi_c, basis$phi_d))

  expect_s3_class(basis, "BoldZipSRSpatialBasis")
  expect_equal(dim(basis$phi_c), c(6L, 2L))
  expect_equal(dim(basis$phi_d), c(6L, 3L))
  expect_equal(gram, diag(5L), tolerance = 1e-10)
  expect_true(all(diff(basis$graph$eigenvalues) >= -1e-10))

  tt <- seq_len(24L)
  X <- basis$phi_c[, 1L, drop = FALSE] %*% matrix(sin(2 * pi * tt / 24), nrow = 1L) +
    basis$phi_d[, 1L, drop = FALSE] %*% matrix(0.5 * sin(2 * pi * tt / 24), nrow = 1L)
  fit <- boldzip_sr_encode(
    X,
    spatial_basis = basis,
    k_carriers = 1L,
    temporal_k = 24L,
    q_texture = 1L,
    center = FALSE,
    events = boldzip_events(max_events = 0L)
  )
  expect_lt(evaluate_boldzip_sr(X, fit)[["mse"]], 1e-10)
  expect_error(boldzip_graph_spatial_basis(adjacency[1:5, ], n_coarse = 1L), "square")
  expect_error(boldzip_graph_spatial_basis(-adjacency, n_coarse = 1L), "non-negative")
})

test_that("BOLDZip-SR validation rejects malformed inputs and settings", {
  expect_error(boldzip_sr_encode(c(1, 2, 3)), "numeric matrix")
  expect_error(boldzip_sr_encode(matrix(1, nrow = 2, ncol = 1)), "two columns")
  expect_error(boldzip_sr_encode(matrix(c(1, Inf), nrow = 1)), "finite")
  expect_error(boldzip_reliability(min_texture_reliability = 2), "finite scalar")
  expect_error(boldzip_quantization(base_step = -0.1), "finite scalar")
  expect_error(boldzip_events(max_events = -1L), "scalar integer")
  expect_error(
    boldzip_sr_encode(
      matrix(rnorm(4 * 8), nrow = 4L),
      spatial_basis = boldzip_spatial_basis(phi_d = matrix(1, nrow = 3L))
    ),
    "phi_d must have 4 rows"
  )
})

test_that("BOLDZip-SR quantization keeps finite decodable payloads", {
  set.seed(20)
  X <- matrix(rnorm(7 * 22), nrow = 7L)
  fit <- boldzip_sr_encode(
    X,
    k_carriers = 3L,
    temporal_k = 8L,
    q_texture = 2L,
    quantization = boldzip_quantization(base_step = 0.05),
    events = boldzip_events(max_events = 4L)
  )
  X_hat <- boldzip_sr_decode(fit)
  metrics <- evaluate_boldzip_sr(X, fit)

  expect_equal(dim(X_hat), dim(X))
  expect_true(all(is.finite(X_hat)))
  expect_true(all(is.finite(metrics[c("mse", "rmse", "payload_scalars")])))
})

test_that("BOLDZip-SR reconstruction is equivariant to voxel permutations", {
  n_vox <- 7L
  n_time <- 26L
  carrier <- sin(2 * pi * seq_len(n_time) / n_time)
  texture <- seq(0.4, 1.8, length.out = n_vox)
  X <- texture %*% t(carrier)
  perm <- c(4L, 1L, 7L, 2L, 6L, 3L, 5L)

  fit <- boldzip_sr_encode(
    X,
    k_carriers = 1L,
    temporal_k = n_time,
    q_texture = 1L,
    center = FALSE,
    events = boldzip_events(max_events = 0L)
  )
  fit_perm <- boldzip_sr_encode(
    X[perm, , drop = FALSE],
    k_carriers = 1L,
    temporal_k = n_time,
    q_texture = 1L,
    center = FALSE,
    events = boldzip_events(max_events = 0L)
  )

  expect_equal(
    boldzip_sr_decode(fit_perm),
    boldzip_sr_decode(fit)[perm, , drop = FALSE],
    tolerance = 1e-10
  )
})

test_that("BOLDZip-SR decoding subsets are consistent with full reconstruction", {
  set.seed(21)
  X <- matrix(rnorm(8 * 24), nrow = 8L)
  fit <- boldzip_sr_encode(
    X,
    k_carriers = 3L,
    temporal_k = 10L,
    q_texture = 2L,
    events = boldzip_events(max_events = 5L)
  )
  full <- boldzip_sr_decode(fit)

  expect_equal(
    boldzip_sr_decode(fit, time_idx = c(2L, 5L, 7L), roi = c(1L, 4L, 8L)),
    full[c(1L, 4L, 8L), c(2L, 5L, 7L), drop = FALSE]
  )
  expect_error(boldzip_sr_decode(fit, time_idx = 0L), "time_idx")
  expect_error(boldzip_sr_decode(fit, roi = 0L), "roi")
})

test_that("BOLDZip-SR temporal budget has monotone reconstruction quality on smooth signals", {
  n_vox <- 6L
  n_time <- 36L
  tt <- seq_len(n_time)
  carrier <- sin(2 * pi * tt / n_time) + 0.25 * cos(4 * pi * tt / n_time)
  X <- seq(0.5, 1.5, length.out = n_vox) %*% t(carrier)

  low_budget <- boldzip_sr_encode(
    X,
    k_carriers = 1L,
    temporal_k = 2L,
    q_texture = 1L,
    center = FALSE,
    events = boldzip_events(max_events = 0L)
  )
  high_budget <- boldzip_sr_encode(
    X,
    k_carriers = 1L,
    temporal_k = n_time,
    q_texture = 1L,
    center = FALSE,
    events = boldzip_events(max_events = 0L)
  )

  expect_lte(
    evaluate_boldzip_sr(X, high_budget)[["mse"]],
    evaluate_boldzip_sr(X, low_budget)[["mse"]] + 1e-12
  )
})

test_that("BOLDZip-SR handles extreme but finite magnitudes without non-finite output", {
  n_time <- 18L
  base <- sin(2 * pi * seq_len(n_time) / n_time)
  X <- rbind(
    1e-8 * base,
    1e8 * base,
    rep(1e4, n_time),
    -1e6 * base
  )

  fit <- boldzip_sr_encode(
    X,
    k_carriers = 2L,
    temporal_k = 8L,
    q_texture = 1L,
    quantization = boldzip_quantization(base_step = 0.01),
    events = boldzip_events(max_events = 2L)
  )
  X_hat <- boldzip_sr_decode(fit)
  metrics <- evaluate_boldzip_sr(X, fit)

  expect_true(all(is.finite(X_hat)))
  expect_true(all(is.finite(metrics[c("mse", "rmse")])))
  expect_equal(dim(X_hat), dim(X))
})

test_that("evaluate_boldzip_sr accepts raw reconstruction matrices", {
  X <- matrix(rnorm(4 * 10), nrow = 4L)
  metrics <- evaluate_boldzip_sr(X, X)

  expect_equal(metrics[["mse"]], 0)
  expect_true(is.na(metrics[["payload_scalars"]]))
})

test_that("BOLDZip-SR simulation harness is deterministic and shaped correctly", {
  sim1 <- boldzip_sr_simulate(
    n_voxels = 9L,
    n_time = 28L,
    k_carriers = 2L,
    q_texture = 1L,
    n_events = 4L,
    noise_sd = 0,
    seed = 99
  )
  sim2 <- boldzip_sr_simulate(
    n_voxels = 9L,
    n_time = 28L,
    k_carriers = 2L,
    q_texture = 1L,
    n_events = 4L,
    noise_sd = 0,
    seed = 99
  )

  expect_equal(sim1$X, sim2$X)
  expect_equal(dim(sim1$X), c(9L, 28L))
  expect_equal(dim(sim1$carriers), c(2L, 28L))
  expect_equal(dim(sim1$texture), c(9L, 2L))
  expect_equal(nrow(sim1$events), 8L)
  expect_true(all(rowSums(sim1$texture != 0) == 1L))
})

test_that("parcel and SVD baselines reconstruct expected reference cases", {
  X <- matrix(seq_len(6 * 5), nrow = 6L)
  parcels <- c(1, 1, 2, 2, 3, 3)
  parcel_hat <- boldzip_parcel_reconstruct(X, parcels)

  expect_equal(dim(parcel_hat), dim(X))
  expect_equal(parcel_hat[1L, ], colMeans(X[1:2, ]))
  expect_equal(parcel_hat[2L, ], colMeans(X[1:2, ]))
  expect_error(boldzip_parcel_reconstruct(X, parcels[-1L]), "one label per row")

  rank_one <- matrix(seq_len(6), nrow = 6L) %*% matrix(seq_len(5), nrow = 1L)
  svd_hat <- boldzip_svd_reconstruct(rank_one, rank = 1L, center = FALSE)
  expect_lt(mean((rank_one - svd_hat)^2), 1e-10)
  expect_error(boldzip_svd_reconstruct(X, rank = 0L), "rank must be")
})

test_that("SVD baseline error is monotone as rank budget increases", {
  set.seed(22)
  X <- matrix(rnorm(8 * 10), nrow = 8L)
  rank1 <- boldzip_svd_reconstruct(X, rank = 1L, center = TRUE)
  rank3 <- boldzip_svd_reconstruct(X, rank = 3L, center = TRUE)

  expect_lte(
    mean((X - rank3)^2),
    mean((X - rank1)^2) + 1e-12
  )
})

test_that("compare_boldzip_sr returns baseline metrics and payload estimates", {
  sim <- boldzip_sr_simulate(
    n_voxels = 10L,
    n_time = 30L,
    k_carriers = 2L,
    n_events = 2L,
    noise_sd = 0.01,
    seed = 7
  )
  fit <- boldzip_sr_encode(
    sim$X,
    k_carriers = 2L,
    temporal_k = 12L,
    q_texture = 1L,
    events = boldzip_events(max_events = 4L)
  )
  cmp <- compare_boldzip_sr(
    sim$X,
    fit,
    parcels = rep(1:5, each = 2L),
    svd_ranks = c(1L, 2L)
  )

  expect_equal(cmp$method, c("boldzip_sr", "parcel", "svd", "svd"))
  expect_true(all(is.finite(cmp$mse)))
  expect_true(all(is.finite(cmp$payload_scalars)))
  expect_equal(cmp$payload_scalars[cmp$method == "parcel"], 5L * ncol(sim$X))
  expect_error(compare_boldzip_sr(sim$X, sim$X), "fit must be")
})
