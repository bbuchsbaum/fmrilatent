library(testthat)

test_that("BOLDZip lagged carrier bank matches hand-shifted carriers", {
  z <- rbind(
    c(10, 20, 30, 40),
    c(1, 3, 5, 7)
  )

  bank <- fmrilatent:::.boldzip_lagged_carrier_bank(z, lags = c(1L, -1L, 0L))

  expect_equal(bank$carrier, c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_equal(bank$lag, c(-1L, 0L, 1L, -1L, 0L, 1L))
  expect_equal(
    bank$bank,
    rbind(
      c(20, 30, 40, 0),
      c(10, 20, 30, 40),
      c(0, 10, 20, 30),
      c(3, 5, 7, 0),
      c(1, 3, 5, 7),
      c(0, 1, 3, 5)
    )
  )
})

test_that("BOLDZip sparse texture fit recovers exact lagged support and coefficient", {
  z <- matrix(c(-3, 1, -3, -2, 2, -1, 3, -2, -1, -1), nrow = 1L)
  target <- 3 * fmrilatent:::.boldzip_lag_signal(z[1, ], lag = 2L)
  y_detail <- matrix(target, nrow = 1L)

  fit <- fmrilatent:::.boldzip_fit_sparse_texture(
    y_detail = y_detail,
    z = z,
    pairs = fmrilatent:::.boldzip_pair_indices(ncol(z), "halves"),
    q = 1L,
    min_reliability = 0.99,
    quantization = boldzip_quantization(base_step = 0),
    lags = -3:3
  )

  expect_equal(nrow(fit$loadings), 1L)
  expect_equal(fit$loadings$atom, 1L)
  expect_equal(fit$loadings$carrier, 1L)
  expect_equal(fit$loadings$lag, 2L)
  expect_equal(fit$loadings$amplitude, 3, tolerance = 1e-12)
  expect_equal(fit$loadings$reliability, 1, tolerance = 1e-12)
  expect_equal(Matrix::nnzero(fit$matrix), 1L)

  predicted <- fmrilatent:::.boldzip_predict_texture(
    fit$loadings,
    z = z,
    n_detail = nrow(y_detail)
  )
  expect_equal(predicted, y_detail, tolerance = 1e-12)
})
