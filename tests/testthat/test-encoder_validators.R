test_that("count validators reject malformed counts with stable classes", {
  invalid_counts <- list(NA_integer_, Inf, c(1L, 2L), 1.5)

  for (value in invalid_counts) {
    expect_error(
      fmrilatent:::.validate_positive_count(value, "k"),
      class = "fmrilatent_error_invalid_count"
    )
    expect_error(
      fmrilatent:::.validate_nonnegative_count(value, "levels"),
      class = "fmrilatent_error_invalid_count"
    )
  }

  expect_error(
    fmrilatent:::.validate_positive_count(0L, "k"),
    "positive integer",
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    fmrilatent:::.validate_nonnegative_count(-1L, "levels"),
    "non-negative integer",
    class = "fmrilatent_error_invalid_count"
  )
  expect_equal(fmrilatent:::.validate_positive_count(3, "k"), 3L)
  expect_equal(fmrilatent:::.validate_nonnegative_count(0, "levels"), 0L)
})

test_that("scalar validators reject malformed finite scalars with stable classes", {
  invalid_scalars <- list(NA_real_, Inf, c(1, 2))

  for (value in invalid_scalars) {
    expect_error(
      fmrilatent:::.validate_positive_scalar(value, "sigma0"),
      class = "fmrilatent_error_invalid_scalar"
    )
    expect_error(
      fmrilatent:::.validate_nonnegative_scalar(value, "threshold"),
      class = "fmrilatent_error_invalid_scalar"
    )
  }

  expect_error(
    fmrilatent:::.validate_positive_scalar(0, "sigma0"),
    "positive finite number",
    class = "fmrilatent_error_invalid_scalar"
  )
  expect_error(
    fmrilatent:::.validate_nonnegative_scalar(-0.1, "threshold"),
    "non-negative finite number",
    class = "fmrilatent_error_invalid_scalar"
  )
  expect_equal(fmrilatent:::.validate_positive_scalar(2L, "sigma0"), 2)
  expect_equal(fmrilatent:::.validate_nonnegative_scalar(0, "threshold"), 0)
})

test_that("flag validator accepts only scalar logical flags", {
  expect_true(fmrilatent:::.validate_flag_scalar(TRUE, "center"))
  expect_false(fmrilatent:::.validate_flag_scalar(FALSE, "center"))

  for (value in list(NA, c(TRUE, FALSE), 1, "TRUE", NULL)) {
    expect_error(
      fmrilatent:::.validate_flag_scalar(value, "center"),
      class = "fmrilatent_error_invalid_flag"
    )
  }
})

test_that("HRBF parameter validator coerces good values and rejects bad ones", {
  params <- fmrilatent:::.validate_hrbf_params(list(
    sigma0 = 2L,
    levels = 1,
    radius_factor = 2.5,
    num_extra_fine_levels = 0,
    seed = 11,
    kernel_type = "gaussian",
    kernel_type_fine_levels = "wendland_c4"
  ))

  expect_equal(params$sigma0, 2)
  expect_equal(params$levels, 1L)
  expect_equal(params$radius_factor, 2.5)
  expect_equal(params$num_extra_fine_levels, 0L)
  expect_equal(params$seed, 11L)
  expect_equal(params$kernel_type_fine_levels, "wendland_c4")

  expect_error(
    fmrilatent:::.validate_hrbf_params("not-list"),
    "params must be a list",
    class = "fmrilatent_error_invalid_params"
  )
  expect_error(
    fmrilatent:::.validate_hrbf_params(list(sigma0 = 0)),
    class = "fmrilatent_error_invalid_scalar"
  )
  expect_error(
    fmrilatent:::.validate_hrbf_params(list(levels = -1)),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    fmrilatent:::.validate_hrbf_params(list(seed = 1.5)),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    fmrilatent:::.validate_hrbf_params(list(kernel_type_fine_levels = "bad")),
    "kernel_type_fine_levels",
    class = "fmrilatent_error_invalid_params"
  )
})
