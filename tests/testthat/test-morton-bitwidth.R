# Tests that the Morton (Z-order) bit-width guards agree between the R and C++
# implementations. The shared policy is 21 bits per axis (0-based coordinate
# in [0, 2097151]); both 64-bit Morton encoders mask inputs to 21 bits, so a
# coordinate >= 2^21 must be rejected (not silently truncated).

test_that("Morton cap constant matches across R", {
  expect_identical(fmrilatent:::.morton_max_coord, 2097151L)
})

test_that("haar_octree_order accepts at the cap and rejects above it", {
  cap <- fmrilatent:::.morton_max_coord  # 0-based max; 1-based = cap + 1
  # 1-based coordinate exactly at the per-axis ceiling (0-based == cap).
  at_cap <- matrix(c(cap + 1L, 1L, 1L), nrow = 1L)
  expect_silent(fmrilatent:::haar_octree_order(at_cap))

  # One past the ceiling must be rejected with the shared message.
  over <- matrix(c(cap + 2L, 1L, 1L), nrow = 1L)
  err <- expect_error(
    fmrilatent:::haar_octree_order(over),
    "21-bit Morton limit \\(max 2097151 per axis\\)"
  )
  # Classed-condition contract from abort_haar().
  expect_s3_class(err, "fmrilatent_haar_error")
  expect_s3_class(err, "fmrilatent_error")
})

test_that("haar_octree_order_cpp rejects 0-based coord above the cap", {
  cap <- fmrilatent:::.morton_max_coord
  # At the cap (0-based) is fine.
  expect_silent(
    fmrilatent:::haar_octree_order_cpp(cap, 0L, 0L, cap + 1L, 1L, 1L)
  )
  # Just above the cap is rejected with the identical wording.
  expect_error(
    fmrilatent:::haar_octree_order_cpp(cap + 1L, 0L, 0L, cap + 2L, 1L, 1L),
    "21-bit Morton limit \\(max 2097151 per axis\\)"
  )
})

test_that("active_morton_order_cpp rejects coord above the cap", {
  cap <- fmrilatent:::.morton_max_coord
  ok <- matrix(c(cap, 0L, 0L), nrow = 1L)
  expect_silent(fmrilatent:::active_morton_order_cpp(ok))

  over <- matrix(c(cap + 1L, 0L, 0L), nrow = 1L)
  expect_error(
    fmrilatent:::active_morton_order_cpp(over),
    "21-bit Morton limit \\(max 2097151 per axis\\)"
  )
})
