# Tests for the unified noise-scale quantization policy and decode/reliability
# behavior (issue bd-01KSQQPCEQSHQ6WTVMCXZR1VME).

test_that(".boldzip_noise_scale is the canonical scale (sd of values)", {
  set.seed(1)
  v <- rnorm(50, sd = 3)
  expect_equal(fmrilatent:::.boldzip_noise_scale(v), stats::sd(v))
  # Coerces to numeric and matches sd on the coerced vector.
  iv <- as.integer(c(1L, 5L, 9L, 2L))
  expect_equal(fmrilatent:::.boldzip_noise_scale(iv), stats::sd(as.numeric(iv)))
  expect_equal(fmrilatent:::.boldzip_noise_scale(7), 0)
  expect_equal(fmrilatent:::.boldzip_noise_scale(numeric()), 0)
  expect_true(is.na(fmrilatent:::.boldzip_corr(1, 1)))
})

test_that("quantization step follows the reliability-shaped formula", {
  quant <- boldzip_quantization(base_step = 0.5, epsilon = 1e-6)
  values <- c(1.234, -2.345, 0.5, 3.21)
  reliability <- c(0.9, 0.1, 0.5, 0.0)
  noise_scale <- 2.0

  out <- fmrilatent:::.boldzip_quantize_values(
    values, reliability = reliability, quantization = quant,
    noise_scale = noise_scale
  )

  step <- quant$base_step * noise_scale /
    sqrt(pmax(reliability, 0) + quant$epsilon)
  step <- pmax(step, .Machine$double.eps)
  expected <- round(values / step) * step
  expect_equal(out, expected)

  # More reliable coefficients get a finer (smaller) step than less reliable
  # ones, given the same noise_scale and base_step.
  rel_steps <- quant$base_step * noise_scale /
    sqrt(pmax(c(0.9, 0.1), 0) + quant$epsilon)
  expect_lt(rel_steps[1], rel_steps[2])
})

test_that("default noise_scale routes through .boldzip_noise_scale", {
  quant <- boldzip_quantization(base_step = 0.3)
  values <- c(0.7, -1.4, 2.2, -0.9, 1.1)
  reliability <- rep(0.5, length(values))

  # Omitting noise_scale must equal supplying .boldzip_noise_scale(values).
  default_out <- fmrilatent:::.boldzip_quantize_values(
    values, reliability = reliability, quantization = quant
  )
  explicit_out <- fmrilatent:::.boldzip_quantize_values(
    values, reliability = reliability, quantization = quant,
    noise_scale = fmrilatent:::.boldzip_noise_scale(values)
  )
  expect_identical(default_out, explicit_out)
})

test_that("event quantization noise scale uses the trimmed event amplitudes", {
  skip_if_not(exists("local_mocked_bindings", envir = asNamespace("testthat")),
              "testthat local_mocked_bindings unavailable")
  recorded_noise_scale <- NULL
  local_mocked_bindings(
    .boldzip_quantize_values = function(values, reliability, quantization, noise_scale = NULL) {
      recorded_noise_scale <<- noise_scale
      values
    },
    .package = "fmrilatent"
  )

  residual <- matrix(c(100, 100, 1, 1), nrow = 1L)
  events <- fmrilatent:::.boldzip_encode_events(
    residual = residual,
    pairs = NULL,
    split_method = "odd_even",
    events = boldzip_events(max_events = 1L, threshold_sd = 0, paired_fraction = 0),
    quantization = boldzip_quantization(base_step = 0.25)
  )

  expect_equal(nrow(events), 1L)
  expect_equal(recorded_noise_scale, 0)
})

test_that("disabled quantization (base_step <= 0) returns values unchanged", {
  quant <- boldzip_quantization(base_step = 0)
  values <- c(1.5, -2.7, 0.33)
  out <- fmrilatent:::.boldzip_quantize_values(
    values, reliability = rep(1, 3), quantization = quant
  )
  expect_identical(out, values)
})

test_that("non-finite / non-positive noise_scale is treated as 1", {
  quant <- boldzip_quantization(base_step = 0.25, epsilon = 1e-6)
  values <- c(1.0, -1.0, 0.5)
  reliability <- rep(1, length(values))

  out_bad <- fmrilatent:::.boldzip_quantize_values(
    values, reliability = reliability, quantization = quant, noise_scale = NA_real_
  )
  out_one <- fmrilatent:::.boldzip_quantize_values(
    values, reliability = reliability, quantization = quant, noise_scale = 1
  )
  expect_equal(out_bad, out_one)
})

test_that("boldzip_sr_simulate validates seed before set.seed", {
  expect_error(
    boldzip_sr_simulate(seed = c(1L, 2L)),
    class = "fmrilatent_error_boldzip"
  )
  expect_error(
    boldzip_sr_simulate(seed = -1L),
    class = "fmrilatent_error_boldzip"
  )
})

test_that("decode is consistent with reliability-driven carrier suppression", {
  sim <- boldzip_sr_simulate(n_voxels = 24L, n_time = 48L, seed = 42L)
  X <- sim$X
  fit <- boldzip_sr_encode(X)

  # Reconstruction must stay finite and correctly shaped.
  recon <- boldzip_sr_decode(fit)
  expect_true(all(is.finite(recon)))
  expect_equal(dim(recon), unname(fit$dimensions[c("voxels", "time")]))
  expect_equal(dim(recon), dim(X))

  # Carrier reliability is bounded in [0, 1].
  rel <- fit$carriers$reliability
  expect_true(all(rel >= 0 & rel <= 1))

  # A high temporal reliability floor zeroes carriers whose split-half
  # reliability falls below the floor, so the stored carrier coefficients
  # (theta) carry no more energy than with no floor.
  fit_floor <- boldzip_sr_encode(
    X,
    reliability = boldzip_reliability(min_temporal_reliability = 0.99)
  )
  recon_floor <- boldzip_sr_decode(fit_floor)
  expect_true(all(is.finite(recon_floor)))
  expect_equal(dim(recon_floor), dim(X))
  expect_lte(sum(fit_floor$carriers$theta^2), sum(fit$carriers$theta^2) + 1e-8)

  # Any carrier below the floor must have been zeroed in theta.
  below <- fit_floor$carriers$reliability < 0.99
  if (any(below)) {
    expect_true(all(abs(fit_floor$carriers$theta[below, , drop = FALSE]) == 0))
  }
})
