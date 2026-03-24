# Generic implicit latent container for basis-free decoders

#' @include all_generic.R
NULL

#' Construct an ImplicitLatent object
#'
#' @param coeff Arbitrary coefficient payload (list or matrix) needed by decoder.
#' @param decoder Function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL) returning matrix.
#' @param meta List metadata; must include `family` string.
#' @param mask Logical 3D array (or LogicalNeuroVol) describing voxel support.
#' @return An object of class `ImplicitLatent`.
#' @export
implicit_latent <- function(coeff, decoder, meta, mask) {
  if (is.null(meta$family)) stop("meta$family required for implicit_latent", call. = FALSE)
  if (!is.function(decoder)) {
    stop("decoder must be a function, got ", paste(class(decoder), collapse = "/"),
         call. = FALSE)
  }
  fmls <- names(formals(decoder))
  has_dots <- "..." %in% fmls
  missing <- setdiff(c("time_idx", "roi_mask"), fmls)
  if (length(missing) > 0L && !has_dots) {
    stop("decoder must accept formals 'time_idx' and 'roi_mask' ",
         "(or use '...'). Missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  structure(list(coeff = coeff, decoder = decoder, meta = meta, mask = mask),
            class = "ImplicitLatent")
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

#' @export
as.matrix.ImplicitLatent <- function(x, time_idx = NULL, roi_mask = NULL, ...) {
  reconstruct_matrix(x, time_idx = time_idx, roi_mask = roi_mask, ...)
}

#' @export
as.array.ImplicitLatent <- function(x, time_idx = NULL, roi_mask = NULL, ...) {
  reconstruct_array(x, time_idx = time_idx, roi_mask = roi_mask, ...)
}

#' @export
#' @rdname mask-methods
setMethod("mask", "ImplicitLatent", function(x) {
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
            mask_arr <- as.array(mask(x))
            rec <- reconstruct_matrix(x, time_idx = time_idx, roi_mask = roi_mask, ...)
            n_time <- nrow(rec)
            arr <- array(0, dim = c(dim(mask_arr), n_time))
            fill_mask <- if (is.null(roi_mask)) mask_arr else (mask_arr & as.logical(roi_mask))
            fill_idx <- which(fill_mask)
            for (t in seq_len(n_time)) {
              slice <- numeric(prod(dim(mask_arr)))
              slice[fill_idx] <- rec[t, ]
              arr[, , , t] <- array(slice, dim = dim(mask_arr))
            }
            arr
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
  if (is.null(roi_mask)) return(rec_mat)
  global_idx <- which(as.logical(mask_arr))
  roi_global <- which(as.logical(roi_mask))
  col_keep   <- which(global_idx %in% roi_global)
  rec_mat[, col_keep, drop = FALSE]
}
