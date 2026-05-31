# Structured-condition coverage for transport_latent.R
# (bd-01KSQQNYTQCW36WNEYPYMX8188). Bare stop()/warning() calls were converted to
# classed .encoder_cli_abort()/.encoder_cli_warn() conditions; the
# validate_portable_linear_map() re-raise preserves the classed inner condition.
# Triggers below need neither rgsp nor neurosurf.

library(testthat)

test_that("validate_portable_linear_map flags missing contract fields (via re-raise)", {
  # .normalize_linear_map() aborts with a class; validate_portable_linear_map()
  # re-raises it with stop(res), so the class must survive to the caller.
  expect_error(
    validate_portable_linear_map(list(forward = function(x) x)),
    class = "fmrilatent_error_invalid_argument"
  )
  expect_error(
    validate_portable_linear_map(list(forward = function(x) x)),
    "missing required fields"
  )
  expect_error(
    validate_portable_linear_map(list()),
    "missing required fields: forward, adjoint_apply, n_source, n_target"
  )
  # error = FALSE still returns FALSE rather than signalling.
  expect_false(
    validate_portable_linear_map(list(forward = function(x) x), error = FALSE)
  )
})

test_that("validate_portable_linear_map rejects an unsupported contract object", {
  expect_error(
    validate_portable_linear_map(42),
    class = "fmrilatent_error_unsupported_operation"
  )
})

test_that("validate_portable_linear_map rejects a non-positive dimension", {
  expect_error(
    validate_portable_linear_map(
      list(forward = function(x) x, adjoint_apply = function(x) x,
           n_source = 0L, n_target = 2L)
    ),
    class = "fmrilatent_error_invalid_argument"
  )
})

test_that("strict map composition warns with a classed condition on domain mismatch", {
  m1 <- list(forward = function(x) x, adjoint_apply = function(x) x,
             n_source = 2L, n_target = 2L, target_domain_id = "A")
  m2 <- list(forward = function(x) x, adjoint_apply = function(x) x,
             n_source = 2L, n_target = 2L, source_domain_id = "B")
  expect_warning(
    fmrilatent:::.compose_linear_maps(m1, m2, strict = TRUE),
    class = "fmrilatent_warning_domain_mismatch"
  )
})
