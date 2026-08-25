#' @include all_generic.R latent_surface_vector.R latent_shared_validation.R
#' @importFrom methods setMethod setValidity new
NULL

.bilat_neurosurfacevector_as_matrix <- function(x, ...) {
  left_obj <- neurosurf::left(x)
  right_obj <- neurosurf::right(x)
  mat_left <- Matrix::Matrix(as.matrix(left_obj), sparse = FALSE)
  mat_right <- Matrix::Matrix(as.matrix(right_obj), sparse = FALSE)
  out <- rbind(mat_left, mat_right)
  attr(out, "coords") <- rbind(
    neuroim2::coords(left_obj),
    neuroim2::coords(right_obj)
  )
  attr(out, "indices") <- c(
    neuroim2::indices(left_obj),
    neuroim2::indices(right_obj)
  )
  attr(out, "hemi") <- c(
    rep(1, nrow(mat_left)),
    rep(2, nrow(mat_right))
  )
  out
}

.bilat_neurosurfacevector_registry <- new.env(parent = emptyenv())
.bilat_neurosurfacevector_registry$as_matrix <- FALSE

.require_bilat_neurosurfacevector_class <- function(context) {
  caller <- rlang::caller_env()
  tryCatch(
    methods::getClass("BilatNeuroSurfaceVector", where = asNamespace("neurosurf")),
    error = function(e) {
      .encoder_cli_abort(
        paste0(context, " requires neurosurf::BilatNeuroSurfaceVector."),
        class = "fmrilatent_error_missing_dependency",
        call = caller
      )
    }
  )
}

.register_bilat_neurosurfacevector_methods <- function(context) {
  .require_bilat_neurosurfacevector_class(context)
  if (!isTRUE(.bilat_neurosurfacevector_registry$as_matrix)) {
    methods::setMethod("as.matrix", "BilatNeuroSurfaceVector",
                       .bilat_neurosurfacevector_as_matrix)
    .bilat_neurosurfacevector_registry$as_matrix <- TRUE
  }
  invisible(TRUE)
}

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
    .encoder_cli_abort("'left' must be a LatentNeuroSurfaceVector.",
                       class = "fmrilatent_error_type", call = rlang::caller_env())
  }
  if (!methods::is(right, "LatentNeuroSurfaceVector")) {
    .encoder_cli_abort("'right' must be a LatentNeuroSurfaceVector.",
                       class = "fmrilatent_error_type", call = rlang::caller_env())
  }
  basis_check <- .validate_shared_basis(
    list(as.matrix(basis(left)), as.matrix(basis(right))),
    labels = c("left", "right"),
    tolerance = 1e-8,
    dim_msg = "left and right basis matrices must have identical dimensions.",
    value_msg = "left and right basis matrices must be equal (the %s basis differs from the left basis)."
  )
  if (!isTRUE(basis_check)) {
    .encoder_cli_abort(paste(basis_check, collapse = "\n"),
                       class = "fmrilatent_error_dim", call = rlang::caller_env())
  }
  if (!is.list(meta)) {
    .encoder_cli_abort("'meta' must be a list.",
                       class = "fmrilatent_error_type", call = rlang::caller_env())
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
#' @rdname basis-methods
setMethod("basis", "BilatLatentNeuroSurfaceVector", function(x) basis(x@left))

#' @export
#' @rdname loadings-methods
setMethod("loadings", "BilatLatentNeuroSurfaceVector",
          function(x) {
            rbind(loadings(x@left), loadings(x@right))
          })

#' @export
#' @rdname offset-methods
setMethod("offset", "BilatLatentNeuroSurfaceVector",
          function(object, ...) c(offset(object@left), offset(object@right)))

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
            .encoder_cli_abort(
              paste0("reconstruct_array() is not defined for bilateral surface latent objects. ",
                     "Use reconstruct_matrix() plus wrap_decoded()."),
              class = "fmrilatent_error_unsupported_operation", call = rlang::caller_env())
          })

#' @export
#' @rdname wrap_decoded
setMethod("wrap_decoded", "BilatLatentNeuroSurfaceVector",
          function(x, values, time_idx = NULL, space = c("native", "template"), ...) {
            space <- match.arg(space)
            if (space != "native") {
              .encoder_cli_abort(
                "wrap_decoded() for BilatLatentNeuroSurfaceVector currently supports only native-space wrapping.",
                class = "fmrilatent_error_unsupported_operation", call = rlang::caller_env())
            }
            if (!requireNamespace("neurosurf", quietly = TRUE)) {
              .encoder_cli_abort("wrap_decoded() requires the 'neurosurf' package.",
                         class = "fmrilatent_error_missing_dependency", call = rlang::caller_env())
            }
            .register_bilat_neurosurfacevector_methods("wrap_decoded()")
            if (is.atomic(values) && is.null(dim(values))) {
              n_left <- length(latent_support(x@left))
              n_right <- length(latent_support(x@right))
              if (length(values) != n_left + n_right) {
                .encoder_cli_abort(
                  "wrap_decoded() vector length does not match bilateral support cardinality.",
                  class = "fmrilatent_error_dim", call = rlang::caller_env())
              }
              values <- matrix(values, nrow = 1L)
            }
            values <- as.matrix(values)
            n_left <- length(latent_support(x@left))
            n_right <- length(latent_support(x@right))
            if (ncol(values) != n_left + n_right) {
              .encoder_cli_abort(
                "wrap_decoded() matrix column count does not match bilateral support cardinality.",
                class = "fmrilatent_error_dim", call = rlang::caller_env())
            }
            left_vals <- values[, seq_len(n_left), drop = FALSE]
            right_vals <- values[, seq.int(n_left + 1L, n_left + n_right), drop = FALSE]
            methods::new(
              "BilatNeuroSurfaceVector",
              left = wrap_decoded(x@left, left_vals, space = "native"),
              right = wrap_decoded(x@right, right_vals, space = "native")
            )
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
          function(x, coordinates = c("analysis", "raw"), ...) diag(ncol(basis(x))))

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
              .encoder_cli_warn(
                "BilatLatentNeuroSurfaceVector has no separate template domain; returning the stored surface decoder.",
                class = "fmrilatent_warning_no_template_domain", call = rlang::caller_env())
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
    handle_backed <- methods::is(object@left@basis, "BasisHandle") ||
      methods::is(object@right@basis, "BasisHandle")
    basis_check <- if (handle_backed) {
      .validate_shared_basis_dims(
        list(
          .explicit_latent_basis_dim(object@left),
          .explicit_latent_basis_dim(object@right)
        ),
        labels = c("left", "right"),
        dim_msg = "Left and right basis matrices must have identical dimensions."
      )
    } else {
      .validate_shared_basis(
        list(as.matrix(object@left@basis), as.matrix(object@right@basis)),
        labels = c("left", "right"),
        tolerance = 1e-8,
        dim_msg = "Left and right basis matrices must have identical dimensions.",
        value_msg = "Left and right basis matrices must be equal (the %s basis differs from the left basis)."
      )
    }
    if (!isTRUE(basis_check)) {
      errors <- c(errors, basis_check)
    }
  }
  if (length(errors) == 0L) TRUE else errors
}

setValidity("BilatLatentNeuroSurfaceVector", .validate_BilatLatentNeuroSurfaceVector)
