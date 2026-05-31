# Handle materialization functions for LatentNeuroVec
# Extracted from latent_neurovector.R for modularity

#' @include latent_handles.R all_class.R all_generic.R bspline_basis.R dct_basis.R slepian_handles.R reduction.R
#' @importFrom methods setGeneric setMethod
#' @importFrom Matrix Matrix
NULL

# --- Internal generics for operator-aware access ---

#' @keywords internal
setGeneric(
  "basis_mat",
  function(x, i = NULL, j = NULL, ...) standardGeneric("basis_mat")
)

#' @keywords internal
setGeneric(
  "loadings_mat",
  function(x, i = NULL, j = NULL, ...) standardGeneric("loadings_mat")
)

# --- Methods for Matrix class ---

.subset_mat <- function(x, i = NULL, j = NULL, ...) {
  i <- i %||% seq_len(nrow(x))
  j <- j %||% seq_len(ncol(x))
  x[i, j, drop = FALSE]
}

setMethod(
  "basis_mat", "Matrix",
  .subset_mat
)

setMethod(
  "loadings_mat", "Matrix",
  .subset_mat
)

# --- Methods for base matrix class ---

setMethod(
  "basis_mat", "matrix",
  .subset_mat
)

setMethod(
  "loadings_mat", "matrix",
  .subset_mat
)

# --- Materialization from spec ---

#' Materialize a BasisHandle into a concrete matrix
#' @keywords internal
#' @noRd
materialize_basis_from_spec <- function(handle, i = NULL, j = NULL) {
  if (identical(handle@kind, "explicit") && !is.null(handle@spec$matrix)) {
    mat <- handle@spec$matrix
    i <- i %||% seq_len(nrow(mat))
    j <- j %||% seq_len(ncol(mat))
    return(Matrix::Matrix(mat[i, j, drop = FALSE], sparse = FALSE))
  }
  obj <- .handle_kind_materializer("basis", handle@kind)(handle)
  basis_mat(obj, i = i, j = j)
}

#' Materialize a LoadingsHandle into a concrete matrix
#' @keywords internal
#' @noRd
materialize_loadings_from_spec <- function(handle, i = NULL, j = NULL) {
  if (identical(handle@kind, "explicit") && !is.null(handle@spec$matrix)) {
    mat <- handle@spec$matrix
    i <- i %||% seq_len(nrow(mat))
    j <- j %||% seq_len(ncol(mat))
    return(Matrix::Matrix(mat[i, j, drop = FALSE], sparse = FALSE))
  }
  obj <- .handle_kind_materializer("loadings", handle@kind)(handle)
  loadings_mat(obj, i = i, j = j)
}

register_handle_kind("dct", function(handle) {
  spec <- handle@spec
  build_dct_basis(
    n_time = spec$n_time,
    k      = spec$k,
    norm   = spec$norm %||% "ortho"
  )
}, type = "basis")

register_handle_kind("slepian_temporal", function(handle) {
  spec <- handle@spec
  Matrix::Matrix(dpss_time_basis(
    n_time   = spec$n_time,
    tr       = spec$tr,
    bandwidth = spec$bandwidth,
    k        = spec$k,
    backend  = spec$backend %||% "tridiag"
  ), sparse = FALSE)
}, type = "basis")

register_handle_kind("bspline", function(handle) {
  spec <- handle@spec
  build_bspline_basis(
    n_time = spec$n_time,
    k      = spec$k,
    degree = spec$degree %||% 3L,
    knots  = spec$knots %||% NULL,
    boundary_knots = spec$boundary_knots %||% NULL,
    include_intercept = spec$include_intercept %||% FALSE,
    orthonormalize = spec$orthonormalize %||% TRUE
  )
}, type = "basis")

register_handle_kind("lifted", function(handle) {
  spec <- handle@spec
  lift(spec$reduction, spec$basis_spec, data = spec$data)
}, type = "basis")

register_handle_kind("explicit", function(handle) {
  spec <- handle@spec
  if (!is.null(spec$matrix)) {
    Matrix::Matrix(spec$matrix, sparse = FALSE)
  } else {
    .encoder_cli_abort("BasisHandle(kind = 'explicit') requires spec$matrix.",
                       class = "fmrilatent_error_missing_argument", call = rlang::caller_env())
  }
}, type = "basis")

register_handle_kind("lifted", function(handle) {
  spec <- handle@spec
  lift(spec$reduction, spec$basis_spec, data = spec$data,
       k_neighbors = spec$k_neighbors %||% 6L)
}, type = "loadings")

register_handle_kind("slepian_spatial", function(handle) {
  spec <- handle@spec
  lift(spec$reduction, spec$basis_spec, data = spec$data,
       k_neighbors = spec$k_neighbors %||% 6L)
}, type = "loadings")

register_handle_kind("explicit", function(handle) {
  spec <- handle@spec
  if (!is.null(spec$matrix)) {
    Matrix::Matrix(spec$matrix, sparse = FALSE)
  } else {
    .encoder_cli_abort("LoadingsHandle(kind = 'explicit') requires spec$matrix.",
                       class = "fmrilatent_error_missing_argument", call = rlang::caller_env())
  }
}, type = "loadings")

# --- Methods for Handle classes ---

setMethod(
  "basis_mat", "BasisHandle",
  function(x, i = NULL, j = NULL, ...) {
    fingerprint <- .latent_handle_fingerprint(x)
    obj <- .latent_get_matrix(x@id, type = "basis", fingerprint = fingerprint)
    if (is.null(obj)) {
      if (!is.null(i) || !is.null(j)) {
        return(materialize_basis_from_spec(x, i = i, j = j))
      }
      obj <- materialize_basis_from_spec(x)
      .latent_register_matrix(
        x@id, obj, type = "basis", overwrite = FALSE,
        fingerprint = fingerprint
      )
    }
    basis_mat(obj, i = i, j = j, ...)
  }
)

setMethod(
  "loadings_mat", "LoadingsHandle",
  function(x, i = NULL, j = NULL, ...) {
    fingerprint <- .latent_handle_fingerprint(x)
    obj <- .latent_get_matrix(x@id, type = "loadings", fingerprint = fingerprint)
    if (is.null(obj)) {
      if (!is.null(i) || !is.null(j)) {
        return(materialize_loadings_from_spec(x, i = i, j = j))
      }
      obj <- materialize_loadings_from_spec(x)
      .latent_register_matrix(
        x@id, obj, type = "loadings", overwrite = FALSE,
        fingerprint = fingerprint
      )
    }
    loadings_mat(obj, i = i, j = j, ...)
  }
)

#' Coerce latent handles to dense matrices
#'
#' @param x A \code{BasisHandle} or \code{LoadingsHandle}.
#' @param ... Optional row/column subsetting arguments passed to
#'   \code{basis_mat()} or \code{loadings_mat()}.
#' @return A base dense matrix.
#' @keywords internal
#' @name as.matrix-handle-methods
NULL

#' @export
#' @rdname as.matrix-handle-methods
setMethod(
  "as.matrix", "BasisHandle",
  function(x, ...) {
    as.matrix(basis_mat(x, ...))
  }
)

#' @export
#' @rdname as.matrix-handle-methods
setMethod(
  "as.matrix", "LoadingsHandle",
  function(x, ...) {
    as.matrix(loadings_mat(x, ...))
  }
)
