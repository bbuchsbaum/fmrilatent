#' @include all_generic.R latent_neurovec_materialize.R
#' @importFrom Matrix Matrix
#' @importFrom methods setMethod setValidity
NULL

.require_neurosurf_surface_latent <- function(context = "LatentNeuroSurfaceVector") {
  if (!requireNamespace("neurosurf", quietly = TRUE)) {
    stop(context, " requires the 'neurosurf' package.", call. = FALSE)
  }
}

#' Construct a LatentNeuroSurfaceVector
#'
#' @param basis Temporal basis matrix (\code{n_time x k}) or \code{BasisHandle}.
#' @param loadings Surface loadings matrix (\code{n_support x k}) or \code{LoadingsHandle}.
#' @param geometry A \code{neurosurf::SurfaceGeometry} or \code{neurosurf::SurfaceSet}.
#' @param support Surface support as vertex indices or a logical vector over all vertices.
#' @param offset Optional numeric vector of length \code{n_support}.
#' @param label Optional label.
#' @param meta Optional metadata list.
#' @return A \code{LatentNeuroSurfaceVector}.
#' @export
LatentNeuroSurfaceVector <- function(basis, loadings, geometry, support = NULL,
                                     offset = NULL, label = "", meta = list()) {
  .require_neurosurf_surface_latent()
  if (!(methods::is(geometry, "SurfaceGeometry") || methods::is(geometry, "SurfaceSet"))) {
    stop("'geometry' must be a neurosurf::SurfaceGeometry or neurosurf::SurfaceSet object.",
         call. = FALSE)
  }
  if (!is.matrix(basis) && !inherits(basis, c("Matrix", "BasisHandle"))) {
    stop("'basis' must be a matrix, Matrix, or BasisHandle object")
  }
  if (!is.matrix(loadings) && !inherits(loadings, c("Matrix", "LoadingsHandle"))) {
    stop("'loadings' must be a matrix, Matrix, or LoadingsHandle object")
  }

  n_nodes <- length(neurosurf::nodes(geometry))
  if (is.null(support)) {
    support <- seq_len(n_nodes)
  }
  support <- .normalize_surface_support(
    support,
    domain = geometry,
    context = "LatentNeuroSurfaceVector support"
  )

  if (is.null(offset) || length(offset) == 0L) {
    offset <- numeric(0)
  } else {
    offset <- as.numeric(offset)
  }

  b_dim <- .latent_basis_dim(basis)
  l_dim <- .latent_loadings_dim(loadings)

  if (b_dim[2L] != l_dim[2L]) {
    stop("'basis' and 'loadings' must have the same number of columns (components)")
  }
  if (l_dim[1L] != length(support)) {
    stop("'loadings' must have ", length(support), " rows to match the support cardinality.")
  }
  if (length(offset) > 0L && length(offset) != l_dim[1L]) {
    stop("'offset' length must match number of rows in 'loadings'")
  }
  if (is.matrix(basis) && !is(basis, "Matrix")) {
    density_basis <- sum(basis != 0) / length(basis)
    basis <- Matrix::Matrix(basis, sparse = !(density_basis > 0.5))
  }
  if (is.matrix(loadings) && !is(loadings, "Matrix")) {
    density_loadings <- sum(loadings != 0) / length(loadings)
    loadings <- Matrix::Matrix(loadings, sparse = !(density_loadings > 0.5))
  }
  if (!is.list(meta)) stop("'meta' must be a list")

  new("LatentNeuroSurfaceVector",
      basis = basis,
      loadings = loadings,
      offset = offset,
      geometry = geometry,
      support = as.integer(support),
      label = label,
      meta = meta)
}

setMethod("basis_mat", "LatentNeuroSurfaceVector",
          function(x, i = NULL, j = NULL, ...) basis_mat(x@basis, i = i, j = j, ...))

setMethod("loadings_mat", "LatentNeuroSurfaceVector",
          function(x, i = NULL, j = NULL, ...) loadings_mat(x@loadings, i = i, j = j, ...))

#' @export
#' @rdname basis-methods
setMethod("basis", "LatentNeuroSurfaceVector", function(x) basis_mat(x))

#' @export
#' @rdname loadings-methods
setMethod("loadings", "LatentNeuroSurfaceVector", function(x) loadings_mat(x))

#' @export
#' @rdname offset-methods
setMethod("offset", "LatentNeuroSurfaceVector", function(object) object@offset)

#' @export
#' @rdname latent_meta
setMethod("latent_meta", "LatentNeuroSurfaceVector", function(x, ...) x@meta %||% list())

#' @export
#' @rdname latent_domain
setMethod("latent_domain", "LatentNeuroSurfaceVector", function(x, ...) x@geometry)

#' @export
#' @rdname latent_support
setMethod("latent_support", "LatentNeuroSurfaceVector", function(x, ...) x@support)

#' @export
#' @rdname is_explicit_latent
setMethod("is_explicit_latent", "LatentNeuroSurfaceVector", function(x, ...) TRUE)

#' @export
#' @rdname reconstruct_matrix
setMethod("reconstruct_matrix", "LatentNeuroSurfaceVector",
          function(x, time_idx = NULL, roi_mask = NULL, ...) {
            mat <- as.matrix(x)
            if (!is.null(time_idx)) {
              mat <- mat[as.integer(time_idx), , drop = FALSE]
            }
            .subset_reconstruction_by_support(
              mat,
              support = x@support,
              roi_mask = roi_mask,
              domain = x@geometry,
              context = "reconstruct_matrix.LatentNeuroSurfaceVector"
            )
          })

#' @export
#' @rdname reconstruct_array
setMethod("reconstruct_array", "LatentNeuroSurfaceVector",
          function(x, time_idx = NULL, roi_mask = NULL, ...) {
            stop("reconstruct_array() is not defined for surface latent objects. ",
                 "Use reconstruct_matrix() plus wrap_decoded().", call. = FALSE)
          })

#' @export
#' @rdname wrap_decoded
setMethod("wrap_decoded", "LatentNeuroSurfaceVector",
          function(x, values, time_idx = NULL, space = c("native", "template"), ...) {
            space <- match.arg(space)
            if (space != "native") {
              stop("wrap_decoded() for LatentNeuroSurfaceVector currently supports only native-space wrapping.",
                   call. = FALSE)
            }
            .wrap_decoded_surface(values, x@geometry, x@support,
                                  context = "wrap_decoded.LatentNeuroSurfaceVector")
          })

#' @export
#' @rdname as.matrix-LatentNeuroVec-method
setMethod("as.matrix", "LatentNeuroSurfaceVector",
          function(x, ...) {
            B <- basis_mat(x)
            L <- loadings_mat(x)
            mat <- B %*% t(L)
            if (length(x@offset) > 0) {
              mat <- sweep(mat, 2, x@offset, "+")
            }
            as.matrix(mat)
          })

#' @export
#' @rdname show-methods
setMethod("show", "LatentNeuroSurfaceVector",
          function(object) {
            cat("LatentNeuroSurfaceVector\n")
            cat("  Time points:", .latent_basis_dim(object@basis)[1L], "\n")
            cat("  Components:", .latent_basis_dim(object@basis)[2L], "\n")
            cat("  Support vertices:", length(object@support), "\n")
            cat("  Geometry:", class(object@geometry)[1], "\n")
            invisible(object)
          })

#' @export
#' @rdname coef_time
setMethod("coef_time", "LatentNeuroSurfaceVector",
          function(x, coordinates = c("analysis", "raw"), ...) {
            .explicit_coef_time(x, as.matrix(basis(x)), coordinates = coordinates)
          })

#' @export
#' @rdname coef_metric
setMethod("coef_metric", "LatentNeuroSurfaceVector",
          function(x, coordinates = c("raw", "analysis"), ...) {
            .explicit_coef_metric(x, ncol(basis(x)), coordinates = coordinates)
          })

#' @export
#' @rdname analysis_transform
setMethod("analysis_transform", "LatentNeuroSurfaceVector",
          function(x, ...) .explicit_latent_analysis_transform(x, ncol(basis(x))))

#' @export
#' @rdname basis_asset
setMethod("basis_asset", "LatentNeuroSurfaceVector", function(x, ...) .explicit_latent_basis_asset(x))

.validate_LatentNeuroSurfaceVector <- function(object) {
  errors <- character()
  if (!inherits(object@basis, c("Matrix", "matrix", "BasisHandle"))) {
    errors <- c(errors, "Slot @basis must be Matrix, matrix, or BasisHandle.")
  }
  if (!inherits(object@loadings, c("Matrix", "matrix", "LoadingsHandle"))) {
    errors <- c(errors, "Slot @loadings must be Matrix, matrix, or LoadingsHandle.")
  }
  if (!is.numeric(object@offset)) {
    errors <- c(errors, "Slot @offset must be numeric.")
  }
  if (!is.integer(object@support)) {
    errors <- c(errors, "Slot @support must be integer.")
  }
  if (!is.character(object@label) || length(object@label) != 1L) {
    errors <- c(errors, "Slot @label must be a single character string.")
  }
  if (!is.list(object@meta)) {
    errors <- c(errors, "Slot @meta must be a list.")
  }
  if (length(errors) == 0L) {
    if (!(methods::is(object@geometry, "SurfaceGeometry") || methods::is(object@geometry, "SurfaceSet"))) {
      errors <- c(errors, "Slot @geometry must be a neurosurf SurfaceGeometry or SurfaceSet.")
    } else {
      n_nodes <- length(neurosurf::nodes(object@geometry))
      if (any(object@support < 1L | object@support > n_nodes)) {
        errors <- c(errors, "Slot @support contains invalid vertex indices.")
      }
    }
    b_dim <- .latent_basis_dim(object@basis)
    l_dim <- .latent_loadings_dim(object@loadings)
    if (b_dim[2L] != l_dim[2L]) {
      errors <- c(errors, "Basis and loadings component counts do not match.")
    }
    if (l_dim[1L] != length(object@support)) {
      errors <- c(errors, "Loadings row count does not match support cardinality.")
    }
    if (length(object@offset) > 0L && length(object@offset) != l_dim[1L]) {
      errors <- c(errors, "Offset length does not match loadings row count.")
    }
  }
  if (length(errors) == 0L) TRUE else errors
}

setValidity("LatentNeuroSurfaceVector", .validate_LatentNeuroSurfaceVector)
