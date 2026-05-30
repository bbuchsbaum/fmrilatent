#' Validate that a set of basis matrices share identical structure and values
#'
#' @param bases A list of materialized basis matrices. The first element is the
#'   reference; every other element must match it.
#' @param labels Optional character vector of component labels used in error
#'   messages (defaults to \code{block1}, \code{block2}, ...).
#' @param tolerance Numeric tolerance passed to \code{\link[base]{all.equal}}.
#' @param dim_msg A \code{sprintf} template (single \code{\%s} for the label)
#'   used when dimensions differ.
#' @param value_msg A \code{sprintf} template (single \code{\%s}) used when
#'   values differ.
#' @return \code{TRUE} if all bases match, otherwise a character vector of
#'   error messages.
#' @keywords internal
#' @noRd
.validate_shared_basis <- function(bases,
                                   labels = NULL,
                                   tolerance = 1e-8,
                                   dim_msg = "Component '%s' dimensions do not match the reference.",
                                   value_msg = "Component '%s' must match the shared reference basis.") {
  n <- length(bases)
  if (n <= 1L) {
    return(TRUE)
  }

  ref <- as.matrix(bases[[1L]])
  ref_dim <- dim(ref)
  if (is.null(labels)) {
    labels <- paste0("block", seq_len(n))
  }

  errors <- character()
  for (i in seq_len(n)) {
    cur <- as.matrix(bases[[i]])
    if (!identical(dim(cur), ref_dim)) {
      errors <- c(errors, sprintf(dim_msg, labels[[i]]))
      next
    }
    if (!isTRUE(all.equal(ref, cur, tolerance = tolerance,
                          check.attributes = FALSE))) {
      errors <- c(errors, sprintf(value_msg, labels[[i]]))
    }
  }
  if (length(errors) == 0L) TRUE else errors
}
