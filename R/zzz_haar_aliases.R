# Aliases and meta helpers for Haar implicit latent

#' Test if object is a Haar latent representation
#' @param x Object to test
#' @return Logical indicating if x is a Haar latent
#' @export
is_haar_latent <- function(x) {
  inherits(x, "HaarLatent") || (is_implicit_latent(x) && identical(x$meta$family, "haar"))
}

#' Get metadata from Haar latent object
#' @param x A Haar latent object
#' @return Metadata list or NULL
#' @export
haar_meta <- function(x) {
  if (inherits(x, "HaarLatent") || is_implicit_latent(x)) {
    m <- x$meta %||% NULL
    if (!is.null(m) && identical(m$family, "haar")) return(m)
  }
  NULL
}

#' Convert to HaarLatent class
#' @param obj Object to convert
#' @return Object with HaarLatent class added
#' @export
as_haar_latent <- function(obj) {
  if (inherits(obj, "HaarLatent")) return(obj)
  if (is_implicit_latent(obj) && identical(obj$meta$family, "haar")) {
    class(obj) <- unique(c("HaarLatent", class(obj)))
    return(obj)
  }
  stop("Object is not Haar implicit latent", call. = FALSE)
}

