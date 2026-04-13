#' @include all_generic.R latent_surface_vector.R
#' @importFrom methods setMethod setValidity new
NULL

#' Construct a bilateral surface latent object
#'
#' @param left Left \code{LatentNeuroSurfaceVector}.
#' @param right Right \code{LatentNeuroSurfaceVector}.
#' @param label Optional label.
#' @param meta Optional metadata list.
#' @return A \code{BilatLatentNeuroSurfaceVector}.
#' @export
BilatLatentNeuroSurfaceVector <- function(left, right, label = "", meta = list()) {
  if (!methods::is(left, "LatentNeuroSurfaceVector")) {
    stop("'left' must be a LatentNeuroSurfaceVector.", call. = FALSE)
  }
  if (!methods::is(right, "LatentNeuroSurfaceVector")) {
    stop("'right' must be a LatentNeuroSurfaceVector.", call. = FALSE)
  }
  if (!identical(dim(basis(left)), dim(basis(right)))) {
    stop("left and right basis matrices must have identical dimensions.", call. = FALSE)
  }
  if (!is.list(meta)) {
    stop("'meta' must be a list.", call. = FALSE)
  }
  new("BilatLatentNeuroSurfaceVector",
      left = left,
      right = right,
      label = label,
      meta = meta)
}

#' @export
#' @rdname latent_meta
setMethod("latent_meta", "BilatLatentNeuroSurfaceVector",
          function(x, ...) {
            utils::modifyList(
              list(
                family = "bilat_surface_explicit",
                left = latent_meta(x@left),
                right = latent_meta(x@right)
              ),
              x@meta %||% list()
            )
          })

#' @export
#' @rdname latent_domain
setMethod("latent_domain", "BilatLatentNeuroSurfaceVector",
          function(x, ...) {
            structure(
              list(left = latent_domain(x@left), right = latent_domain(x@right)),
              class = "BilateralSurfaceDomain"
            )
          })

#' @export
#' @rdname latent_support
setMethod("latent_support", "BilatLatentNeuroSurfaceVector",
          function(x, ...) {
            list(left = latent_support(x@left), right = latent_support(x@right))
          })

#' @export
#' @rdname is_explicit_latent
setMethod("is_explicit_latent", "BilatLatentNeuroSurfaceVector", function(x, ...) TRUE)

#' @export
#' @rdname basis-methods
setMethod("basis", "BilatLatentNeuroSurfaceVector", function(x) basis(x@left))

#' @export
#' @rdname loadings-methods
setMethod("loadings", "BilatLatentNeuroSurfaceVector",
          function(x) {
            Matrix::Matrix(
              rbind(as.matrix(loadings(x@left)), as.matrix(loadings(x@right))),
              sparse = FALSE
            )
          })

#' @export
#' @rdname offset-methods
setMethod("offset", "BilatLatentNeuroSurfaceVector",
          function(object) c(offset(object@left), offset(object@right)))

#' @export
#' @rdname reconstruct_matrix
setMethod("reconstruct_matrix", "BilatLatentNeuroSurfaceVector",
          function(x, time_idx = NULL, roi_mask = NULL, ...) {
            left_roi <- if (is.list(roi_mask)) roi_mask$left %||% NULL else NULL
            right_roi <- if (is.list(roi_mask)) roi_mask$right %||% NULL else NULL
            cbind(
              reconstruct_matrix(x@left, time_idx = time_idx, roi_mask = left_roi, ...),
              reconstruct_matrix(x@right, time_idx = time_idx, roi_mask = right_roi, ...)
            )
          })

#' @export
#' @rdname reconstruct_array
setMethod("reconstruct_array", "BilatLatentNeuroSurfaceVector",
          function(x, time_idx = NULL, roi_mask = NULL, ...) {
            stop("reconstruct_array() is not defined for bilateral surface latent objects. ",
                 "Use reconstruct_matrix() plus wrap_decoded().", call. = FALSE)
          })

#' @export
#' @rdname wrap_decoded
setMethod("wrap_decoded", "BilatLatentNeuroSurfaceVector",
          function(x, values, time_idx = NULL, space = c("native", "template"), ...) {
            space <- match.arg(space)
            if (space != "native") {
              stop("wrap_decoded() for BilatLatentNeuroSurfaceVector currently supports only native-space wrapping.",
                   call. = FALSE)
            }
            if (!requireNamespace("neurosurf", quietly = TRUE)) {
              stop("wrap_decoded() requires the 'neurosurf' package.", call. = FALSE)
            }
            if (is.atomic(values) && is.null(dim(values))) {
              n_left <- length(latent_support(x@left))
              n_right <- length(latent_support(x@right))
              if (length(values) != n_left + n_right) {
                stop("wrap_decoded() vector length does not match bilateral support cardinality.",
                     call. = FALSE)
              }
              left_obj <- wrap_decoded(x@left, values[seq_len(n_left)], space = "native")
              right_obj <- wrap_decoded(x@right, values[seq.int(n_left + 1L, n_left + n_right)],
                                        space = "native")
              return(new("BilatNeuroSurfaceVector",
                         left = neurosurf::NeuroSurfaceVector(
                           neurosurf::geometry(left_obj), neuroim2::indices(left_obj), Matrix::Matrix(left_obj@data, ncol = 1)
                         ),
                         right = neurosurf::NeuroSurfaceVector(
                           neurosurf::geometry(right_obj), neuroim2::indices(right_obj), Matrix::Matrix(right_obj@data, ncol = 1)
                         )))
            }
            values <- as.matrix(values)
            n_left <- length(latent_support(x@left))
            n_right <- length(latent_support(x@right))
            if (ncol(values) != n_left + n_right) {
              stop("wrap_decoded() matrix column count does not match bilateral support cardinality.",
                   call. = FALSE)
            }
            left_vals <- values[, seq_len(n_left), drop = FALSE]
            right_vals <- values[, seq.int(n_left + 1L, n_left + n_right), drop = FALSE]
            new("BilatNeuroSurfaceVector",
                left = wrap_decoded(x@left, left_vals, space = "native"),
                right = wrap_decoded(x@right, right_vals, space = "native"))
          })

#' @export
#' @rdname as.matrix-LatentNeuroVec-method
setMethod("as.matrix", "BilatLatentNeuroSurfaceVector",
          function(x, ...) reconstruct_matrix(x, ...))

#' @export
#' @rdname show-methods
setMethod("show", "BilatLatentNeuroSurfaceVector",
          function(object) {
            cat("BilatLatentNeuroSurfaceVector\n")
            cat("  Time points:", nrow(basis(object)), "\n")
            cat("  Components:", ncol(basis(object)), "\n")
            cat("  Left support:", length(latent_support(object)$left), "\n")
            cat("  Right support:", length(latent_support(object)$right), "\n")
            invisible(object)
          })

#' @export
#' @rdname coef_time
setMethod("coef_time", "BilatLatentNeuroSurfaceVector",
          function(x, coordinates = c("analysis", "raw"), ...) as.matrix(basis(x)))

#' @export
#' @rdname coef_metric
setMethod("coef_metric", "BilatLatentNeuroSurfaceVector",
          function(x, coordinates = c("raw", "analysis"), ...) diag(ncol(basis(x))))

#' @export
#' @rdname analysis_transform
setMethod("analysis_transform", "BilatLatentNeuroSurfaceVector",
          function(x, ...) .transport_identity_transform(ncol(basis(x))))

#' @export
#' @rdname basis_asset
setMethod("basis_asset", "BilatLatentNeuroSurfaceVector", function(x, ...) NULL)

#' @export
#' @rdname decoder
setMethod("decoder", "BilatLatentNeuroSurfaceVector",
          function(x, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), ...) {
            space <- match.arg(space)
            coordinates <- match.arg(coordinates)
            if (space == "template") {
              warning("BilatLatentNeuroSurfaceVector has no separate template domain; returning the stored surface decoder.",
                      call. = FALSE)
            }
            .latent_loadings_map(x)
          })

#' @export
#' @rdname decode_coefficients
setMethod("decode_coefficients", "BilatLatentNeuroSurfaceVector",
          function(x, gamma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"),
                   wrap = c("none", "auto"), ...) {
            .decode_coefficients_via_decoder(x, gamma, space = space,
                                             coordinates = coordinates,
                                             wrap = wrap, ...)
          })

#' @export
#' @rdname decode_covariance
setMethod("decode_covariance", "BilatLatentNeuroSurfaceVector",
          function(x, Sigma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), diag_only = TRUE, ...) {
            map <- decoder(x, space = space, coordinates = coordinates, ...)
            Sigma <- .as_square_matrix(Sigma, map$n_source, context = "Sigma")
            if (isTRUE(diag_only)) {
              .project_covariance_diag(map, Sigma)
            } else {
              D <- .materialize_linear_map(map)
              D %*% Sigma %*% t(D)
            }
          })

#' @export
#' @rdname project_effect
setMethod("project_effect", "BilatLatentNeuroSurfaceVector",
          function(x, gamma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), ...) {
            decode_coefficients(x, gamma, space = space, coordinates = coordinates, ...)
          })

#' @export
#' @rdname project_vcov
setMethod("project_vcov", "BilatLatentNeuroSurfaceVector",
          function(x, Sigma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), diag_only = TRUE, ...) {
            decode_covariance(x, Sigma, space = space, coordinates = coordinates,
                              diag_only = diag_only, ...)
          })

.validate_BilatLatentNeuroSurfaceVector <- function(object) {
  errors <- character()
  if (!methods::is(object@left, "LatentNeuroSurfaceVector")) {
    errors <- c(errors, "Slot @left must be a LatentNeuroSurfaceVector.")
  }
  if (!methods::is(object@right, "LatentNeuroSurfaceVector")) {
    errors <- c(errors, "Slot @right must be a LatentNeuroSurfaceVector.")
  }
  if (length(errors) == 0L) {
    if (!identical(dim(basis(object@left)), dim(basis(object@right)))) {
      errors <- c(errors, "Left and right basis matrices must have identical dimensions.")
    }
  }
  if (length(errors) == 0L) TRUE else errors
}

setValidity("BilatLatentNeuroSurfaceVector", .validate_BilatLatentNeuroSurfaceVector)
