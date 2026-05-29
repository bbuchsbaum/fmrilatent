# Generic implicit latent container for basis-free decoders

#' @include all_generic.R
NULL

.wrap_decoded_volume <- function(values, mask_like, context = "wrap_decoded") {
  has_space <- inherits(mask_like, "LogicalNeuroVol")
  mask_vol <- if (has_space) {
    mask_like
  } else {
    neuroim2::LogicalNeuroVol(as.logical(mask_like), neuroim2::NeuroSpace(dim(mask_like)))
  }
  mask_arr <- as.array(mask_vol)
  card <- sum(mask_arr)
  space_3d <- neuroim2::space(mask_vol)

  if (is.atomic(values) && is.null(dim(values))) {
    if (length(values) != card) {
      .encoder_cli_abort(
        paste0(context, " vector length ", length(values),
               " does not match support cardinality ", card, "."),
        class = "fmrilatent_error_dimension_mismatch"
      )
    }
    out <- numeric(prod(dim(mask_arr)))
    out[which(mask_arr)] <- values
    arr <- array(out, dim = dim(mask_arr))
    if (has_space) {
      return(neuroim2::DenseNeuroVol(arr, space_3d))
    }
    return(arr)
  }

  values <- as.matrix(values)
  if (ncol(values) != card) {
    .encoder_cli_abort(
      paste0(context, " matrix column count ", ncol(values),
             " does not match support cardinality ", card, "."),
      class = "fmrilatent_error_dimension_mismatch"
    )
  }

  n_time <- nrow(values)
  arr <- array(0, dim = c(dim(mask_arr), n_time))
  fill_idx <- which(mask_arr)
  for (t in seq_len(n_time)) {
    slice <- numeric(prod(dim(mask_arr)))
    slice[fill_idx] <- values[t, ]
    arr[, , , t] <- array(slice, dim = dim(mask_arr))
  }
  if (has_space) {
    space_4d <- neuroim2::NeuroSpace(
      c(dim(mask_arr), n_time),
      spacing = neuroim2::spacing(space_3d),
      origin = neuroim2::origin(space_3d)
    )
    return(neuroim2::DenseNeuroVec(arr, space_4d))
  }
  arr
}

.is_surface_domain <- function(domain) {
  if (is.null(domain) || !requireNamespace("neurosurf", quietly = TRUE)) {
    return(FALSE)
  }
  methods::is(domain, "SurfaceGeometry") || methods::is(domain, "SurfaceSet")
}

.normalize_surface_support <- function(support, domain = NULL, context = "surface support") {
  if (is.null(support)) {
    .encoder_cli_abort(paste0(context, " is required for surface-like domains."),
                       class = "fmrilatent_error_missing_argument")
  }
  if (is.logical(support)) {
    support <- which(as.logical(support))
  }
  support <- as.integer(support)
  if (length(support) == 0L || anyNA(support) || any(support < 1L)) {
    .encoder_cli_abort(paste0(context, " must contain positive vertex indices."),
                       class = "fmrilatent_error_invalid_support")
  }
  if (.is_surface_domain(domain)) {
    n_nodes <- length(neurosurf::nodes(domain))
    if (any(support > n_nodes)) {
      .encoder_cli_abort(
        paste0(context, " contains indices beyond the domain node count ", n_nodes, "."),
        class = "fmrilatent_error_invalid_support"
      )
    }
  }
  support
}

.wrap_decoded_surface <- function(values, domain, support, context = "wrap_decoded") {
  if (!.is_surface_domain(domain)) {
    .encoder_cli_abort(paste0(context, " requires a SurfaceGeometry or SurfaceSet domain."),
                       class = "fmrilatent_error_invalid_geometry")
  }
  support <- .normalize_surface_support(support, domain = domain, context = context)
  if (is.atomic(values) && is.null(dim(values))) {
    if (length(values) != length(support)) {
      .encoder_cli_abort(
        paste0(context, " vector length ", length(values),
               " does not match support cardinality ", length(support), "."),
        class = "fmrilatent_error_dimension_mismatch"
      )
    }
    return(neurosurf::NeuroSurface(domain, support, as.numeric(values)))
  }

  values <- as.matrix(values)
  if (ncol(values) != length(support)) {
    .encoder_cli_abort(
      paste0(context, " matrix column count ", ncol(values),
             " does not match support cardinality ", length(support), "."),
      class = "fmrilatent_error_dimension_mismatch"
    )
  }

  neurosurf::NeuroSurfaceVector(domain, support, t(values))
}

#' Construct an ImplicitLatent object
#'
#' @param coeff Arbitrary coefficient payload (list or matrix) needed by decoder.
#' @param decoder Function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL) returning matrix.
#' @param meta List metadata; must include `family` string.
#' @param mask Logical 3D array (or LogicalNeuroVol) describing volumetric support.
#' @param domain Optional decoded output domain. For non-volumetric latent objects,
#'   supply a domain such as a \code{neurosurf::SurfaceGeometry} together with
#'   \code{support}.
#' @param support Optional decoded output support. For volumetric latent objects this
#'   is usually derived from \code{mask}; for surface-like latent objects it should be
#'   a vector of vertex indices (or a logical mask over vertices).
#' @return An object of class `ImplicitLatent`.
#' @export
implicit_latent <- function(coeff, decoder, meta, mask = NULL, domain = NULL, support = NULL) {
  if (is.null(meta$family)) {
    .encoder_cli_abort("meta$family required for implicit_latent",
                       class = "fmrilatent_error_invalid_argument")
  }
  if (!is.function(decoder)) {
    .encoder_cli_abort(
      paste0("decoder must be a function, got ", paste(class(decoder), collapse = "/")),
      class = "fmrilatent_error_invalid_argument"
    )
  }
  fmls <- names(formals(decoder))
  has_dots <- "..." %in% fmls
  missing <- setdiff(c("time_idx", "roi_mask"), fmls)
  if (length(missing) > 0L && !has_dots) {
    .encoder_cli_abort(
      paste0("decoder must accept formals 'time_idx' and 'roi_mask' ",
             "(or use '...'). Missing: ", paste(missing, collapse = ", ")),
      class = "fmrilatent_error_invalid_argument"
    )
  }
  if (is.null(mask) && is.null(support)) {
    .encoder_cli_abort(
      "implicit_latent requires either a volumetric mask or an explicit support.",
      class = "fmrilatent_error_missing_argument"
    )
  }
  if (!is.null(mask) && is.null(support)) {
    support <- mask
  }
  if (is.null(domain) && !is.null(mask)) {
    domain <- if (inherits(mask, "LogicalNeuroVol")) {
      neuroim2::space(mask)
    } else {
      neuroim2::NeuroSpace(dim(mask))
    }
  }
  structure(list(
    coeff = coeff,
    decoder = decoder,
    meta = meta,
    mask = mask,
    domain = domain,
    support = support
  ),
            class = "ImplicitLatent")
}

#' Coerce an object to an ImplicitLatent decoder
#'
#' @param x Object to coerce.
#' @param ... Additional arguments passed to methods.
#' @return An object of class `ImplicitLatent`.
#' @export
as_implicit_latent <- function(x, ...) {
  UseMethod("as_implicit_latent")
}

#' @export
as_implicit_latent.ImplicitLatent <- function(x, ...) {
  x
}

#' @export
as_implicit_latent.default <- function(x, ...) {
  .encoder_cli_abort(
    paste0("No as_implicit_latent method for class: ", paste(class(x), collapse = ",")),
    class = "fmrilatent_error_unsupported_operation"
  )
}

#' Predict method for ImplicitLatent
#' @param object ImplicitLatent object
#' @param roi_mask Optional ROI mask
#' @param time_idx Optional time indices
#' @param levels_keep Optional levels to keep
#' @param ... Additional arguments (unused)
#' @return Matrix of predicted values
#' @export
predict.ImplicitLatent <- function(object, roi_mask = NULL, time_idx = NULL,
                                   levels_keep = NULL, ...) {
  object$decoder(time_idx = time_idx, roi_mask = roi_mask, levels_keep = levels_keep)
}

#' Test if object is an ImplicitLatent
#' @param x Object to test
#' @return Logical indicating if x is an ImplicitLatent
#' @export
is_implicit_latent <- function(x) inherits(x, "ImplicitLatent")

#' Get metadata from ImplicitLatent object
#' @param x An ImplicitLatent object
#' @return Metadata list or NULL
#' @export
implicit_meta <- function(x) {
  if (!is_implicit_latent(x)) return(NULL)
  x$meta %||% NULL
}

#' Reconstruct an ImplicitLatent as a matrix
#'
#' @method as.matrix ImplicitLatent
#' @param x An \code{ImplicitLatent} object.
#' @param time_idx Optional integer time indices to keep.
#' @param roi_mask Optional logical ROI mask for spatial subsetting.
#' @param ... Additional arguments passed to the decoder.
#' @return A numeric matrix with rows = time and columns = voxels within the
#'   requested mask support.
#' @export
as.matrix.ImplicitLatent <- function(x, time_idx = NULL, roi_mask = NULL, ...) {
  reconstruct_matrix(x, time_idx = time_idx, roi_mask = roi_mask, ...)
}

#' Reconstruct an ImplicitLatent as an array
#'
#' @method as.array ImplicitLatent
#' @param x An \code{ImplicitLatent} object.
#' @param time_idx Optional integer time indices to keep.
#' @param roi_mask Optional logical ROI mask; voxels outside the ROI are zero.
#' @param ... Additional arguments passed to the decoder.
#' @return A numeric array with dimensions \code{c(x, y, z, time)}.
#' @export
as.array.ImplicitLatent <- function(x, time_idx = NULL, roi_mask = NULL, ...) {
  reconstruct_array(x, time_idx = time_idx, roi_mask = roi_mask, ...)
}

#' @export
#' @rdname mask-methods
setMethod("mask", "ImplicitLatent", function(x) {
  if (is.null(x$mask)) {
    .encoder_cli_abort("mask() is not defined for non-volumetric ImplicitLatent objects.",
                       class = "fmrilatent_error_unsupported_operation")
  }
  if (inherits(x$mask, "LogicalNeuroVol")) {
    x$mask
  } else {
    neuroim2::LogicalNeuroVol(as.logical(x$mask), neuroim2::NeuroSpace(dim(x$mask)))
  }
})

#' @export
#' @rdname latent_meta
setMethod("latent_meta", "ImplicitLatent", function(x, ...) x$meta %||% list())

#' @export
#' @rdname latent_domain
setMethod("latent_domain", "ImplicitLatent",
          function(x, ...) {
            if (!is.null(x$domain)) {
              return(x$domain)
            }
            if (!is.null(x$mask)) {
              return(neuroim2::space(mask(x)))
            }
            NULL
          })

#' @export
#' @rdname latent_support
setMethod("latent_support", "ImplicitLatent",
          function(x, ...) x$support %||% x$mask %||% NULL)

#' @export
#' @rdname is_explicit_latent
setMethod("is_explicit_latent", "ImplicitLatent", function(x, ...) FALSE)

#' @export
#' @rdname reconstruct_matrix
setMethod("reconstruct_matrix", "ImplicitLatent",
          function(x, time_idx = NULL, roi_mask = NULL, ...) {
            predict.ImplicitLatent(x, time_idx = time_idx, roi_mask = roi_mask, ...)
          })

#' @export
#' @rdname reconstruct_array
setMethod("reconstruct_array", "ImplicitLatent",
          function(x, time_idx = NULL, roi_mask = NULL, ...) {
            if (is.null(x$mask)) {
              .encoder_cli_abort(
                paste0("reconstruct_array() is only defined for volumetric ImplicitLatent objects. ",
                       "Use reconstruct_matrix() plus wrap_decoded() for other domains."),
                class = "fmrilatent_error_unsupported_operation"
              )
            }
            mask_arr <- as.array(mask(x))
            roi_arr <- .normalize_roi_mask(mask_arr, roi_mask, "reconstruct_array.ImplicitLatent")
            rec <- reconstruct_matrix(x, time_idx = time_idx, roi_mask = roi_arr, ...)
            fill_mask <- if (is.null(roi_arr)) mask_arr else (mask_arr & roi_arr)
            .wrap_decoded_volume(rec, fill_mask, context = "reconstruct_array.ImplicitLatent")
          })

.wrap_decoded_by_support <- function(values, support, domain, context) {
  if (inherits(support, "LogicalNeuroVol") ||
      (is.logical(support) && length(dim(support)) == 3L)) {
    return(.wrap_decoded_volume(values, support, context = context))
  }
  if (.is_surface_domain(domain)) {
    return(.wrap_decoded_surface(values, domain, support, context = context))
  }
  .encoder_cli_abort(paste0(context, " has no wrapper for this domain/support combination."),
                     class = "fmrilatent_error_unsupported_operation")
}

#' @export
#' @rdname wrap_decoded
setMethod("wrap_decoded", "ImplicitLatent",
          function(x, values, time_idx = NULL, space = c("native", "template"), ...) {
            space <- match.arg(space)
            if (space == "native") {
              return(.wrap_decoded_by_support(
                values,
                support = latent_support(x),
                domain = latent_domain(x),
                context = "wrap_decoded.ImplicitLatent[native]"
              ))
            }
            # Template-space wrapping is only meaningful for transport-backed
            # latent objects, which carry a basis_asset whose template_support()
            # / template_domain() describe the template-side sample layout.
            if (!.is_transport_latent(x)) {
              .encoder_cli_abort(
                paste0("wrap_decoded() does not support space = \"template\" for ",
                       "non-transport latent objects (no basis asset to describe ",
                       "template-side support)."),
                class = "fmrilatent_error_unsupported_operation"
              )
            }
            asset <- basis_asset(x)
            if (is.null(asset)) {
              .encoder_cli_abort(
                paste0("wrap_decoded(space = \"template\") requires a basis_asset ",
                       "on the transport latent object."),
                class = "fmrilatent_error_missing_argument"
              )
            }
            t_support <- tryCatch(template_support(asset), error = function(e) NULL)
            t_domain <- tryCatch(template_domain(asset), error = function(e) NULL)
            if (is.null(t_support)) {
              .encoder_cli_abort(
                paste0("wrap_decoded(space = \"template\") requires the basis ",
                       "asset to implement template_support()."),
                class = "fmrilatent_error_missing_argument"
              )
            }
            .wrap_decoded_by_support(
              values,
              support = t_support,
              domain = t_domain,
              context = "wrap_decoded.ImplicitLatent[template]"
            )
          })

#' Subset reconstruction matrix columns by ROI mask
#'
#' Extracts columns from a reconstruction matrix that correspond to voxels
#' inside an ROI mask. This is a shared utility used by decoder functions
#' to handle ROI subsetting consistently.
#'
#' @param rec_mat Numeric matrix (time x voxels-in-mask) to subset.
#' @param mask_arr Logical array (the full brain mask).
#' @param roi_mask Logical array (the ROI mask), or \code{NULL} to return
#'   \code{rec_mat} unchanged.
#' @return The subsetted matrix, or \code{rec_mat} unchanged when
#'   \code{roi_mask} is \code{NULL}.
#' @export
#' @examples
#' mask <- array(c(TRUE, TRUE, FALSE, TRUE), dim = c(2, 2, 1))
#' rec  <- matrix(1:9, nrow = 3, ncol = 3)  # 3 time x 3 masked voxels
#' roi  <- array(c(TRUE, FALSE, FALSE, TRUE), dim = c(2, 2, 1))
#' roi_subset_columns(rec, mask, roi)  # keeps columns 1 and 3
roi_subset_columns <- function(rec_mat, mask_arr, roi_mask = NULL) {
  roi_arr <- .normalize_roi_mask(mask_arr, roi_mask, "roi_subset_columns")
  if (is.null(roi_arr)) return(rec_mat)
  global_idx <- which(as.logical(mask_arr))
  roi_global <- which(roi_arr)
  col_keep   <- which(global_idx %in% roi_global)
  rec_mat[, col_keep, drop = FALSE]
}
