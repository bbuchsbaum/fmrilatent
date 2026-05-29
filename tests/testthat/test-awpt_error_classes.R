# Structured-condition coverage for awpt.R (bd-01KSQQNYTQCW36WNEYPYMX8188).
# The AWPT error paths were converted from bare stop() to classed
# .encoder_cli_abort() conditions; assert the classes so downstream callers can
# branch on them programmatically. These triggers need neither rgsp nor
# neurosurf so they run unconditionally.

library(testthat)

test_that("basis_awpt_wavelet flags missing custom weights with a classed error", {
  expect_error(
    basis_awpt_wavelet(penalty_rule = "custom"),
    class = "fmrilatent_error_invalid_awpt_spec"
  )
})

test_that("awpt_mean_conductance rejects an empty conductance list", {
  expect_error(
    awpt_mean_conductance(list()),
    class = "fmrilatent_error_invalid_conductance"
  )
})

test_that("awpt_mean_conductance rejects mismatched conductance dimensions", {
  expect_error(
    awpt_mean_conductance(list(diag(2), diag(3))),
    class = "fmrilatent_error_dimension_mismatch"
  )
})

test_that("awpt_mean_conductance rejects out-of-range shrinkage", {
  expect_error(
    awpt_mean_conductance(list(diag(2)), shrinkage = 2),
    class = "fmrilatent_error_invalid_argument"
  )
})

test_that("awpt_basis_template rejects a non-reduction parcellation", {
  expect_error(
    awpt_basis_template(42),
    class = "fmrilatent_error_invalid_reduction"
  )
})

test_that("classed AWPT errors still carry the original message text", {
  # The message text is preserved through the cli_abort conversion, so callers
  # matching on substrings continue to work alongside class-based handling.
  expect_error(
    awpt_mean_conductance(list(diag(2), diag(3))),
    "identical dimensions"
  )
})
