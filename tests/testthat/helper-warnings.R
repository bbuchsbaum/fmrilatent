without_warnings_matching <- function(expr, pattern) {
  withCallingHandlers(
    force(expr),
    warning = function(w) {
      if (grepl(pattern, conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

without_dense_basis_warning <- function(expr) {
  without_warnings_matching(
    expr,
    "^Input 'basis' is dense [(]100% non-zero[)]; storing as dense dgeMatrix[.]$"
  )
}

without_dct_fallback_warnings <- function(expr) {
  without_warnings_matching(
    expr,
    "^encode_spec[.]spec_time_dct (was singular|ridge solve failed)"
  )
}
