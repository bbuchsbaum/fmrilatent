# BOLDZip-SR codec: spatial basis descriptors and graph spatial basis
#
# Split out of R/codec_boldzip.R (see that file's header for the full map).
# Functions in this file:
#   boldzip_spatial_basis
#   .boldzip_loadings_from_template
#   as_boldzip_spatial_basis
#   as_boldzip_spatial_basis.BoldZipSRSpatialBasis
#   as_boldzip_spatial_basis.SharedReference
#   as_boldzip_spatial_basis.default
#   boldzip_graph_spatial_basis


#' Build a matrix spatial basis descriptor for BOLDZip-SR
#'
#' @param phi_c Optional coarse basis matrix with rows equal to voxels.
#' @param phi_d Optional detail basis matrix with rows equal to voxels. If
#'   `NULL`, the detail basis is the identity basis.
#' @param label Optional label stored in metadata.
#' @param basis_asset Optional source template or shared-basis object used to
#'   build this descriptor.
#' @return A `BoldZipSRSpatialBasis` object.
#' @export
boldzip_spatial_basis <- function(phi_c = NULL, phi_d = NULL, label = NULL,
                                  basis_asset = NULL) {
  if (!is.null(phi_c)) {
    phi_c <- .boldzip_validate_spatial_basis_input(phi_c, "phi_c")
  }
  if (!is.null(phi_d)) {
    phi_d <- .boldzip_validate_spatial_basis_input(phi_d, "phi_d")
  }
  if (!is.null(phi_c) && !is.null(phi_d) && nrow(phi_c) != nrow(phi_d)) {
    .boldzip_cli_abort("phi_c and phi_d must have the same number of rows.")
  }
  structure(
    list(
      phi_c = phi_c,
      phi_d = phi_d,
      label = label %||% "matrix",
      basis_asset = basis_asset
    ),
    class = "BoldZipSRSpatialBasis"
  )
}

.boldzip_loadings_from_template <- function(x) {
  tryCatch(template_loadings(x), error = function(e) NULL)
}

#' Coerce a spatial object to a BOLDZip-SR spatial basis
#'
#' @description
#' This helper lets the standalone BOLDZip-SR codec consume matrix-like shared
#' basis assets without registering BOLDZip as an `encode()` family. Matrix and
#' template inputs are used as the detail basis by default and are
#' orthonormalized because BOLDZip projection currently uses the transpose as
#' the analysis operator.
#'
#' @param x A `BoldZipSRSpatialBasis`, matrix-like object, shared reference, or
#'   template object supporting [template_loadings()].
#' @param ... Additional arguments reserved for methods. The default method
#'   accepts `label`, `role`, and `orthonormalize`.
#' @return A `BoldZipSRSpatialBasis` object.
#' @export
as_boldzip_spatial_basis <- function(x, ...) {
  UseMethod("as_boldzip_spatial_basis")
}

#' @export
as_boldzip_spatial_basis.BoldZipSRSpatialBasis <- function(x, ...) {
  x
}

#' @export
as_boldzip_spatial_basis.SharedReference <- function(x, ...) {
  out <- as_boldzip_spatial_basis(resolve_shared_reference(x), ...)
  out$basis_asset <- x
  out
}

#' @export
as_boldzip_spatial_basis.default <- function(x, label = NULL,
                                             role = c("detail", "coarse"),
                                             orthonormalize = TRUE, ...) {
  role <- match.arg(role)
  if (is.list(x) && (any(c("phi_c", "phi_d") %in% names(x)))) {
    args <- x[intersect(names(x), c("phi_c", "phi_d", "label", "basis_asset"))]
    return(do.call(boldzip_spatial_basis, args))
  }

  loadings <- if (is.matrix(x) || inherits(x, "Matrix")) {
    x
  } else {
    .boldzip_loadings_from_template(x)
  }
  if (is.null(loadings)) {
    .boldzip_cli_abort(
      "x must be a BoldZipSRSpatialBasis, matrix-like object, SharedReference, ",
      "or template object with template_loadings()."
    )
  }
  phi <- .boldzip_validate_spatial_basis_input(loadings, "template_loadings(x)")
  if (isTRUE(orthonormalize)) {
    phi <- .boldzip_orthonormalize_columns(phi)
  }
  label <- label %||% tryCatch(template_meta(x)$family, error = function(e) NULL) %||%
    class(x)[[1L]]
  if (identical(role, "coarse")) {
    return(boldzip_spatial_basis(phi_c = phi, label = label, basis_asset = x))
  }
  boldzip_spatial_basis(phi_d = phi, label = label, basis_asset = x)
}

#' Build a spectral graph spatial basis for BOLDZip-SR
#'
#' @description
#' Constructs an experimental graph-adapted spatial basis from a symmetric
#' adjacency matrix. Low-frequency graph Laplacian eigenvectors become coarse
#' carrier support functions and the following higher-frequency eigenvectors
#' become detail atoms. This is a materialized small-to-medium graph MVP, not a
#' production graph-wavelet operator.
#'
#' @param adjacency Symmetric non-negative adjacency matrix.
#' @param n_coarse Number of low-frequency eigenvectors to use as `phi_c`.
#' @param n_detail Number of following eigenvectors to use as `phi_d`. Defaults
#'   to all remaining eigenvectors.
#' @param normalized Whether to use the normalized graph Laplacian.
#' @param label Optional label stored in metadata.
#' @return A `BoldZipSRSpatialBasis` object.
#' @export
boldzip_graph_spatial_basis <- function(adjacency,
                                        n_coarse = 8L,
                                        n_detail = NULL,
                                        normalized = TRUE,
                                        label = NULL) {
  adjacency <- .boldzip_validate_matrix(adjacency, "adjacency")
  if (nrow(adjacency) != ncol(adjacency)) {
    .boldzip_cli_abort("adjacency must be square.")
  }
  if (any(adjacency < 0)) {
    .boldzip_cli_abort("adjacency must be non-negative.")
  }
  adjacency <- 0.5 * (adjacency + t(adjacency))
  diag(adjacency) <- 0
  n <- nrow(adjacency)
  n_coarse <- min(.boldzip_check_scalar_integer(n_coarse, "n_coarse", min = 0L), n)
  if (is.null(n_detail)) {
    n_detail <- n - n_coarse
  }
  n_detail <- min(.boldzip_check_scalar_integer(n_detail, "n_detail", min = 0L),
                  n - n_coarse)
  if (n_coarse + n_detail < 1L) {
    .boldzip_cli_abort("n_coarse + n_detail must be at least 1.")
  }

  degree <- rowSums(adjacency)
  if (isTRUE(normalized)) {
    inv_sqrt <- ifelse(degree > 0, 1 / sqrt(degree), 0)
    laplacian <- diag(n) - (inv_sqrt * adjacency) * rep(inv_sqrt, each = n)
  } else {
    laplacian <- diag(degree, nrow = n) - adjacency
  }
  laplacian <- 0.5 * (laplacian + t(laplacian))
  eig <- eigen(laplacian, symmetric = TRUE)
  ord <- order(eig$values, decreasing = FALSE)
  vectors <- .boldzip_canonicalize_eigenvectors(eig$vectors[, ord, drop = FALSE])

  phi_c <- if (n_coarse > 0L) {
    vectors[, seq_len(n_coarse), drop = FALSE]
  } else {
    NULL
  }
  detail_start <- n_coarse + 1L
  phi_d <- if (n_detail > 0L) {
    vectors[, seq.int(detail_start, length.out = n_detail), drop = FALSE]
  } else {
    NULL
  }

  out <- boldzip_spatial_basis(
    phi_c = phi_c,
    phi_d = phi_d,
    label = label %||% "spectral_graph"
  )
  out$graph <- list(
    normalized = isTRUE(normalized),
    eigenvalues = eig$values[ord][seq_len(n_coarse + n_detail)],
    n_coarse = n_coarse,
    n_detail = n_detail
  )
  out
}
