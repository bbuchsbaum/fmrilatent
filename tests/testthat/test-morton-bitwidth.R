# Tests that the Morton (Z-order) bit-width guards agree between the R and C++
# implementations. The shared policy is 10 bits per axis (max dimension 1024,
# max 0-based coordinate 1023). A volume with a dimension > 1024 along any axis
# must be REJECTED identically by the R reference encoder and the C++ encoder,
# rather than being silently accepted by one and rejected by the other.

test_that("shared Morton bit-width constants are as documented", {
  expect_identical(fmrilatent:::.morton_max_bits_per_axis, 10L)
  expect_identical(fmrilatent:::.morton_max_coord, 1023L)
  expect_match(fmrilatent:::.morton_overflow_msg, "max 10 bits / 1024 per axis")
})

test_that("C++ Morton ordering accepts dim 1024 and rejects dim 1025", {
  skip_if_not(fmrilatent:::use_haar_rcpp())
  # Sparse masks keep the allocation tiny; only the bit-width guard is exercised.
  mk <- function(d) {
    a <- array(FALSE, dim = c(d, 1L, 1L))
    a[1L, 1L, 1L] <- TRUE
    a
  }
  # 1024 -> exactly 10 bits -> accepted.
  expect_silent(fmrilatent:::get_morton_ordered_indices_rcpp(mk(1024L)))
  # 1025 -> 11 bits -> rejected with the shared message.
  expect_error(
    fmrilatent:::get_morton_ordered_indices_rcpp(mk(1025L)),
    "too large for Morton codes \\(max 10 bits / 1024 per axis\\)"
  )
})

test_that("R reference Morton ordering rejects dim 1025 with shared message", {
  # Force the pure-R path so we exercise the R guard, not the C++ one.
  mk <- function(d) {
    a <- array(FALSE, dim = c(d, 1L, 1L))
    a[1L, 1L, 1L] <- TRUE
    a[d, 1L, 1L] <- TRUE
    a
  }
  testthat::local_mocked_bindings(use_haar_rcpp = function() FALSE)
  # 1024 -> accepted.
  expect_silent(fmrilatent:::get_morton_ordered_indices(mk(1024L)))
  # 1025 -> rejected; same human-readable text and classed condition as C++.
  err <- expect_error(
    fmrilatent:::get_morton_ordered_indices(mk(1025L)),
    "too large for Morton codes \\(max 10 bits / 1024 per axis\\)"
  )
  expect_s3_class(err, "fmrilatent_error_morton_overflow")
})

test_that("R and C++ agree at the dim-1024/1025 boundary", {
  skip_if_not(fmrilatent:::use_haar_rcpp())
  mk <- function(d) {
    a <- array(FALSE, dim = c(d, 1L, 1L))
    a[1L, 1L, 1L] <- TRUE
    a[d, 1L, 1L] <- TRUE
    a
  }
  cpp_status <- function(d) {
    tryCatch({
      fmrilatent:::get_morton_ordered_indices_rcpp(mk(d))
      "accept"
    }, error = function(e) "reject")
  }
  r_status <- function(d) {
    testthat::with_mocked_bindings(
      tryCatch({
        fmrilatent:::get_morton_ordered_indices(mk(d))
        "accept"
      }, error = function(e) "reject"),
      use_haar_rcpp = function() FALSE
    )
  }
  for (d in c(1024L, 1025L)) {
    expect_identical(cpp_status(d), r_status(d),
                     info = sprintf("dim = %d", d))
  }
})
