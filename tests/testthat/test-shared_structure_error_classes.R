# Structured-condition coverage for shared_structure.R
# (bd-01KSQQNYTQCW36WNEYPYMX8188). Bare stop() calls were converted to classed
# .encoder_cli_abort() conditions; validators re-raise the classed inner
# condition rather than flattening it. These triggers are dependency-free.

library(testthat)

test_that("shared_reference rejects a non-function materialize with a classed error", {
  expect_error(
    shared_reference(materialize = 42),
    class = "fmrilatent_error_invalid_argument"
  )
})

test_that("shared_component_contract flags a support/loadings mismatch", {
  loadings <- matrix(rnorm(6), nrow = 3, ncol = 2)
  expect_error(
    shared_component_contract(loadings, support = c("a", "b")),
    class = "fmrilatent_error_dimension_mismatch"
  )
  # Message text is preserved through the cli_abort conversion.
  expect_error(
    shared_component_contract(loadings, support = c("a", "b")),
    "support cardinality 2 does not match loadings rows 3"
  )
})

test_that("group_delta_loadings flags mismatched group/delta dimensions", {
  group <- matrix(0, 2, 2)
  delta <- matrix(0, 1, 2)
  expect_error(
    group_delta_loadings(group, delta),
    class = "fmrilatent_error_dimension_mismatch"
  )
})

test_that("shared_temporal_spec rejects rank greater than n_time", {
  expect_error(
    shared_temporal_spec("dct", n_time = 3L, rank = 5L),
    class = "fmrilatent_error_invalid_argument"
  )
})

test_that("validators re-raise the classed inner condition (not a flattened error)", {
  bad <- structure(
    list(id = "x", contract_version = 1L, family = "f",
         n_features = 3L, n_components = 2L, support = c("a", "b", "c", "d")),
    class = "SharedComponentContract"
  )
  expect_error(
    validate_shared_component_contract(bad, error = TRUE),
    class = "fmrilatent_error_dimension_mismatch"
  )
  # error = FALSE still returns FALSE rather than signalling.
  expect_false(validate_shared_component_contract(bad, error = FALSE))
})

test_that("render_shared_events stays robust to glue braces in a user shape_id", {
  dict <- shared_event_dictionary(shapes = list(impulse = 1), n_time = 4L)
  events <- data.frame(
    atom = 1L, time = 1L, amplitude = 1,
    shape_id = "a{b}", stringsAsFactors = FALSE
  )
  # A "{" in user data must not crash cli's interpolation; it must surface as
  # the intended classed error with the offending id echoed verbatim.
  expect_error(
    render_shared_events(dict, events, n_atoms = 1L),
    class = "fmrilatent_error_invalid_events"
  )
  expect_error(
    render_shared_events(dict, events, n_atoms = 1L),
    "a\\{b\\}"
  )
})
