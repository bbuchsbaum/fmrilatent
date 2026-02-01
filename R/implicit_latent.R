# Generic implicit latent container for basis-free decoders

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

