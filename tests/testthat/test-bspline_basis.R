library(testthat)

expected_bspline_matrix <- function(n_time,
                                    k,
                                    degree = 3L,
                                    include_intercept = FALSE,
                                    boundary_knots = c(0, 1)) {
  n_int <- max(0L, k - degree - as.integer(include_intercept))
  knots <- NULL
  if (n_int > 0L) {
    knots <- seq(boundary_knots[1L], boundary_knots[2L],
                 length.out = n_int + 2L)[-c(1L, n_int + 2L)]
  }
  mat <- as.matrix(splines::bs(
    seq_len(n_time) / n_time,
    df = k,
    degree = degree,
    knots = knots,
    Boundary.knots = boundary_knots,
    intercept = include_intercept
  ))
  matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat))
}

test_that("build_bspline_basis generates exactly k non-intercept columns", {
  B <- as.matrix(build_bspline_basis(
    n_time = 20L,
    k = 6L,
    degree = 3L,
    include_intercept = FALSE,
    orthonormalize = FALSE
  ))
  expected <- expected_bspline_matrix(
    n_time = 20L,
    k = 6L,
    degree = 3L,
    include_intercept = FALSE
  )

  expect_equal(dim(B), c(20L, 6L))
  expect_equal(unname(B), unname(expected), tolerance = 1e-12)
  expect_true(all(colSums(abs(B)) > 0))
})

test_that("build_bspline_basis generates exactly k intercept columns", {
  B <- as.matrix(build_bspline_basis(
    n_time = 20L,
    k = 6L,
    degree = 3L,
    include_intercept = TRUE,
    orthonormalize = FALSE
  ))
  expected <- expected_bspline_matrix(
    n_time = 20L,
    k = 6L,
    degree = 3L,
    include_intercept = TRUE
  )

  expect_equal(dim(B), c(20L, 6L))
  expect_equal(unname(B), unname(expected), tolerance = 1e-12)
  expect_true(all(colSums(abs(B)) > 0))
})

test_that("build_bspline_basis does not fabricate padded QR directions", {
  B <- as.matrix(build_bspline_basis(
    n_time = 20L,
    k = 6L,
    degree = 3L,
    include_intercept = FALSE,
    orthonormalize = TRUE
  ))

  expect_equal(dim(B), c(20L, 6L))
  expect_equal(crossprod(B), diag(6L), tolerance = 1e-10)
})

test_that("build_bspline_basis rejects explicit knots with wrong column count", {
  expect_error(
    build_bspline_basis(
      n_time = 20L,
      k = 6L,
      degree = 3L,
      knots = c(0.25, 0.5),
      include_intercept = FALSE,
      orthonormalize = FALSE
    ),
    class = "fmrilatent_error_dimension_mismatch"
  )
})

test_that("build_bspline_basis rejects invalid counts before spline construction", {
  expect_error(
    build_bspline_basis(n_time = 0L, k = 4L),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    build_bspline_basis(n_time = 8L, k = 0L),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    build_bspline_basis(n_time = 8L, k = 4.5),
    class = "fmrilatent_error_invalid_count"
  )
})

test_that("bspline_basis_handle rejects invalid dimensions immediately", {
  expect_error(
    bspline_basis_handle(n_time = 0L, k = 4L),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    bspline_basis_handle(n_time = 8L, k = 0L),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    bspline_basis_handle(n_time = 4L, k = 5L),
    class = "fmrilatent_error_invalid_count"
  )
})
