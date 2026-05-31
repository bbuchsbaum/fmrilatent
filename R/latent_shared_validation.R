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
  for (i in seq.int(2L, n)) {
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

.explicit_latent_basis_dim <- function(x) {
  if (methods::is(x, "LatentNeuroVec") || methods::is(x, "LatentNeuroSurfaceVector")) {
    return(.latent_basis_dim(x@basis))
  }
  if (methods::is(x, "BilatLatentNeuroSurfaceVector")) {
    return(.explicit_latent_basis_dim(x@left))
  }
  if (methods::is(x, "BlockLatentNeuroVector")) {
    if (length(x@blocks) < 1L) {
      .encoder_cli_abort(
        "BlockLatentNeuroVector must contain at least one block.",
        class = "fmrilatent_error_value"
      )
    }
    return(.explicit_latent_basis_dim(x@blocks[[1L]]))
  }
  .latent_basis_dim(basis(x))
}

.validate_shared_basis_dims <- function(dims,
                                        labels = NULL,
                                        dim_msg = "Component '%s' dimensions do not match the reference.") {
  n <- length(dims)
  if (n <= 1L) {
    return(TRUE)
  }
  ref_dim <- as.integer(dims[[1L]])
  if (is.null(labels)) {
    labels <- paste0("block", seq_len(n))
  }

  errors <- character()
  for (i in seq_len(n)) {
    cur_dim <- as.integer(dims[[i]])
    if (!identical(cur_dim, ref_dim)) {
      errors <- c(errors, sprintf(dim_msg, labels[[i]]))
    }
  }
  if (length(errors) == 0L) TRUE else errors
}
