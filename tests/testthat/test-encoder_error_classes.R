test_that("encoder validators use structured cli errors", {
  expect_error(
    fmrilatent:::.validate_positive_count(0, "k"),
    "positive integer",
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    fmrilatent:::.validate_positive_count(.Machine$integer.max + 1, "k"),
    "positive integer",
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    fmrilatent:::.validate_positive_count(4.5, "k"),
    "positive integer",
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    fmrilatent:::.validate_nonnegative_count(.Machine$integer.max + 1, "levels"),
    "non-negative integer",
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    fmrilatent:::.validate_nonnegative_count(1.5, "levels"),
    "non-negative integer",
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    fmrilatent:::.validate_hrbf_params(list(kernel_type = "bad")),
    "kernel_type",
    class = "fmrilatent_error_invalid_params"
  )
})

test_that("encoder-facing constructors expose stable error classes", {
  expect_error(
    fmrilatent:::build_dct_basis(n_time = 2L, k = 3L),
    "cannot exceed n_time",
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    spec_space_parcel(list()),
    "ParcelBasisTemplate",
    class = "fmrilatent_error_invalid_parcel_template"
  )
  expect_error(
    hrbf_reconstruct_partial(matrix(1, nrow = 2, ncol = 1), array(TRUE, c(2, 1, 1)), list()),
    "voxel_idx",
    class = "fmrilatent_error_invalid_index"
  )
})

test_that("encode core validation helpers expose stable error classes", {
  expect_error(
    fmrilatent:::.run_lengths_from_info(5L, c(2L, 2L)),
    "run lengths must sum",
    class = "fmrilatent_error_dim"
  )
  expect_error(
    fmrilatent:::.normalize_penalty_matrix(c(1, 2), 3L, context = "spatial_penalty"),
    "vector must have length 3",
    class = "fmrilatent_error_dim"
  )
  expect_error(
    fmrilatent:::.normalize_penalty_matrix(matrix(0, nrow = 2L, ncol = 2L), 3L,
                                           context = "spatial_penalty"),
    "must have dimensions 3x3",
    class = "fmrilatent_error_dim"
  )
  expect_error(
    latent_factory("awpt", x = matrix(0, nrow = 1L, ncol = 1L),
                   mask = array(TRUE, dim = c(1L, 1L, 1L))),
    "does not support AWPT",
    class = "fmrilatent_error_unsupported_operation"
  )
})

test_that("encode dispatch uses structured cli errors", {
  expect_error(
    encode(list(a = 1), spec_time_dct(k = 3), mask = array(TRUE, dim = c(1, 1, 1))),
    class = "fmrilatent_error_unsupported_encode_input"
  )
  expect_error(
    encode_spec(matrix(1, 1, 1), structure(list(), class = "spec_unknown")),
    class = "fmrilatent_error_unsupported_spec"
  )
})
