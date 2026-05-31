library(testthat)

make_boldzip_error_fit <- function() {
  set.seed(1001)
  X <- matrix(rnorm(6L * 8L), nrow = 6L)
  boldzip_sr_encode(
    X,
    k_carriers = 2L,
    temporal_k = 2L,
    q_texture = 1L,
    reliability = boldzip_reliability(
      min_temporal_reliability = 0,
      min_texture_reliability = 0
    ),
    events = boldzip_events(max_events = 0L)
  )
}

test_that("BOLDZip scalar and matrix validators raise codec-classed errors", {
  expect_error(
    boldzip_sr_encode(c(1, 2, 3)),
    class = "fmrilatent_error_boldzip"
  )
  expect_equal(
    dim(fmrilatent:::.boldzip_validate_matrix(matrix(1, nrow = 2L, ncol = 1L), "X")),
    c(2L, 1L)
  )
  expect_error(
    boldzip_sr_encode(matrix(1, nrow = 2L, ncol = 1L)),
    "at least two time points",
    class = "fmrilatent_error_boldzip"
  )
  expect_error(
    boldzip_reliability(min_texture_reliability = 2),
    class = "fmrilatent_error_boldzip"
  )
  expect_error(
    boldzip_events(max_events = -1L),
    class = "fmrilatent_error_boldzip"
  )
})

test_that("BOLDZip spatial coercion and graph basis errors are classed", {
  expect_error(
    boldzip_spatial_basis(
      phi_c = matrix(1, nrow = 2L, ncol = 1L),
      phi_d = matrix(1, nrow = 3L, ncol = 1L)
    ),
    class = "fmrilatent_error_boldzip"
  )
  expect_error(
    as_boldzip_spatial_basis(list(not_a_basis = TRUE)),
    class = "fmrilatent_error_boldzip"
  )
  expect_error(
    boldzip_graph_spatial_basis(matrix(1, nrow = 2L, ncol = 3L), n_coarse = 1L),
    class = "fmrilatent_error_boldzip"
  )
})

test_that("BOLDZip temporal-spec and linalg errors are classed", {
  expect_error(
    fmrilatent:::.boldzip_materialize_temporal_spec(
      matrix(1, nrow = 3L, ncol = 4L),
      n_time = 3L
    ),
    regexp = "more columns than rows",
    class = "fmrilatent_error_boldzip"
  )
  expect_error(
    fmrilatent:::.boldzip_materialize_temporal_spec(
      spec_time_bspline(k = 5L, orthonormalize = FALSE),
      n_time = 8L
    ),
    class = "fmrilatent_error_boldzip"
  )
  expect_error(
    fmrilatent:::.boldzip_orthonormalize_columns(matrix(0, nrow = 2L, ncol = 1L)),
    class = "fmrilatent_error_boldzip"
  )
})

test_that("BOLDZip orthonormalization warns when rank drops", {
  x <- cbind(c(1, 0, 0), c(1, 0, 0), c(0, 1, 0))
  expect_warning(
    out <- fmrilatent:::.boldzip_orthonormalize_columns(x),
    "dropped 1 linearly dependent column"
  )
  expect_equal(dim(out), c(3L, 2L))
})

test_that("BOLDZip correlation treats near-zero variance as undefined", {
  nearly_constant <- 1 + c(-1e-9, 0, 1e-9)
  expect_true(is.na(fmrilatent:::.boldzip_corr(nearly_constant, c(1, 2, 3))))
})

test_that("BOLDZip decode and ROI errors are classed", {
  fit <- make_boldzip_error_fit()
  expect_error(
    boldzip_sr_decode(list()),
    class = "fmrilatent_error_boldzip"
  )
  expect_error(
    boldzip_sr_decode(fit, time_idx = 0L),
    class = "fmrilatent_error_boldzip"
  )
  expect_error(
    boldzip_sr_decode(fit, roi = "bad"),
    class = "fmrilatent_error_boldzip"
  )
  expect_error(
    fmrilatent:::.boldzip_roi_from_mask(c(TRUE, NA), n_voxels = 2L),
    class = "fmrilatent_error_boldzip"
  )
})

test_that("BOLDZip diagnostics and baselines raise codec-classed errors", {
  X <- matrix(0, nrow = 2L, ncol = 2L)
  expect_error(
    boldzip_sr_payload_summary(list()),
    class = "fmrilatent_error_boldzip"
  )
  expect_error(
    evaluate_boldzip_sr(X, matrix(0, nrow = 3L, ncol = 2L)),
    class = "fmrilatent_error_boldzip"
  )
  expect_error(
    evaluate_boldzip_sr(X, X, reliability_weights = -1),
    class = "fmrilatent_error_boldzip"
  )
  expect_error(
    boldzip_parcel_reconstruct(X, parcels = 1L),
    class = "fmrilatent_error_boldzip"
  )
  expect_error(
    compare_boldzip_sr(X, fit = X),
    class = "fmrilatent_error_boldzip"
  )
})
