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

setMethod(
  "basis_mat", "Matrix",
  function(x, i = NULL, j = NULL, ...) {
    i <- i %||% seq_len(nrow(x))
    j <- j %||% seq_len(ncol(x))
    x[i, j, drop = FALSE]
  }
)

setMethod(
  "loadings_mat", "Matrix",
  function(x, i = NULL, j = NULL, ...) {
    i <- i %||% seq_len(nrow(x))
    j <- j %||% seq_len(ncol(x))
    x[i, j, drop = FALSE]
  }
)

# --- Methods for base matrix class ---

setMethod(
  "basis_mat", "matrix",
  function(x, i = NULL, j = NULL, ...) {
    i <- i %||% seq_len(nrow(x))
    j <- j %||% seq_len(ncol(x))
    x[i, j, drop = FALSE]
  }
)

setMethod(
  "loadings_mat", "matrix",
  function(x, i = NULL, j = NULL, ...) {
    i <- i %||% seq_len(nrow(x))
    j <- j %||% seq_len(ncol(x))
    x[i, j, drop = FALSE]
  }
)

# --- Materialization from spec ---

#' Materialize a BasisHandle into a concrete matrix
#' @keywords internal
#' @noRd
materialize_basis_from_spec <- function(handle) {
  kind <- handle@kind
  spec <- handle@spec

  if (identical(kind, "dct")) {
    build_dct_basis(
      n_time = spec$n_time,
      k      = spec$k,
      norm   = spec$norm %||% "ortho"
    )
  } else if (identical(kind, "slepian_temporal")) {
    Matrix::Matrix(dpss_time_basis(
      n_time   = spec$n_time,
      tr       = spec$tr,
      bandwidth = spec$bandwidth,
      k        = spec$k,
      backend  = spec$backend %||% "tridiag"
    ), sparse = FALSE)
  } else if (identical(kind, "bspline")) {
    build_bspline_basis(
      n_time = spec$n_time,
      k      = spec$k,
      degree = spec$degree %||% 3L,
      knots  = spec$knots %||% NULL,
      boundary_knots = spec$boundary_knots %||% NULL,
      include_intercept = spec$include_intercept %||% FALSE,
      orthonormalize = spec$orthonormalize %||% TRUE
    )
  } else if (identical(kind, "lifted")) {
    lift(spec$reduction, spec$basis_spec, data = spec$data)
  } else if (identical(kind, "explicit")) {
    if (!is.null(spec$matrix)) {
      Matrix::Matrix(spec$matrix)
    } else {
      stop("BasisHandle(kind = 'explicit') requires spec$matrix.")
    }
  } else {
    stop("Unknown BasisHandle kind: ", kind)
  }
}

#' Materialize a LoadingsHandle into a concrete matrix
#' @keywords internal
#' @noRd
materialize_loadings_from_spec <- function(handle) {
  kind <- handle@kind
  spec <- handle@spec

  if (identical(kind, "lifted")) {
    lift(spec$reduction, spec$basis_spec, data = spec$data,
         k_neighbors = spec$k_neighbors %||% 6L)
  } else if (identical(kind, "slepian_spatial")) {
    lift(spec$reduction, spec$basis_spec, data = spec$data,
         k_neighbors = spec$k_neighbors %||% 6L)
  } else if (identical(kind, "explicit")) {
    if (!is.null(spec$matrix)) {
      Matrix::Matrix(spec$matrix)
    } else {
      stop("LoadingsHandle(kind = 'explicit') requires spec$matrix.")
    }
  } else {
    stop("Unknown LoadingsHandle kind: ", kind)
  }
}

# --- Methods for Handle classes ---

setMethod(
  "basis_mat", "BasisHandle",
  function(x, i = NULL, j = NULL, ...) {
    obj <- .latent_get_matrix(x@id, type = "basis")
    if (is.null(obj)) {
      obj <- materialize_basis_from_spec(x)
      .latent_register_matrix(x@id, obj, type = "basis", overwrite = FALSE)
    }
    basis_mat(obj, i = i, j = j, ...)
  }
)

setMethod(
  "loadings_mat", "LoadingsHandle",
  function(x, i = NULL, j = NULL, ...) {
    obj <- .latent_get_matrix(x@id, type = "loadings")
    if (is.null(obj)) {
      obj <- materialize_loadings_from_spec(x)
      .latent_register_matrix(x@id, obj, type = "loadings", overwrite = FALSE)
    }
    loadings_mat(obj, i = i, j = j, ...)
  }
)
