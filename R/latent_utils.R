# Package-wide utility primitives loaded before handle and class definitions.

# Utility: safe NULL coalesce
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Convert mask to array
#'
#' Safely converts a \code{LogicalNeuroVol} or array-like mask to a plain
#' logical array, with informative error messages on failure.
#'
#' @param mask A \code{LogicalNeuroVol} or logical array.
#' @param location Character string used in error messages to identify the caller.
#' @return A logical array.
#' @export
mask_to_array <- function(mask, location = "unknown function") {
  result <- tryCatch(
    as.array(mask),
    error = function(e) {
      .encoder_cli_abort(
        paste0(sprintf("In %s: mask must be array-like or LogicalNeuroVol. ", location),
               "Underlying error: ", conditionMessage(e)),
        class = "fmrilatent_error_invalid_mask", call = rlang::caller_env()
      )
    }
  )
  if (is.null(result)) {
    .encoder_cli_abort(
      paste0(sprintf("In %s: mask must be array-like or LogicalNeuroVol ", location),
             "(conversion returned NULL)."),
      class = "fmrilatent_error_invalid_mask", call = rlang::caller_env()
    )
  }
  result
}

# Internal alias so existing callers don't break
.mask_to_array <- function(mask, location = "unknown function") {
  mask_to_array(mask, location)
}

.normalize_roi_mask <- function(mask_arr, roi_mask = NULL, location = "unknown function") {
  if (is.null(roi_mask)) {
    return(NULL)
  }

  roi_arr <- .mask_to_array(roi_mask, location)
  if (!identical(dim(roi_arr), dim(mask_arr))) {
    .encoder_cli_abort(
      paste0("roi_mask dimensions (", paste(dim(roi_arr), collapse = " x "),
             ") do not match mask dimensions (", paste(dim(mask_arr), collapse = " x "), ")."),
      class = "fmrilatent_error_dim", call = rlang::caller_env()
    )
  }

  roi_arr <- array(as.logical(roi_arr), dim = dim(mask_arr))
  if (any(roi_arr & !as.logical(mask_arr))) {
    .encoder_cli_abort("roi_mask contains voxels outside the object mask.",
                       class = "fmrilatent_error_value", call = rlang::caller_env())
  }

  roi_arr
}

.space_with_time_from_mask <- function(mask, n_time, location = "unknown function") {
  mask_arr <- .mask_to_array(mask, location)
  if (inherits(mask, "LogicalNeuroVol")) {
    mask_space <- neuroim2::space(mask)
    return(neuroim2::NeuroSpace(
      c(dim(mask_arr), as.integer(n_time)),
      spacing = neuroim2::spacing(mask_space),
      origin = neuroim2::origin(mask_space)
    ))
  }

  neuroim2::NeuroSpace(c(dim(mask_arr), as.integer(n_time)))
}

.assert_template_mask_match <- function(mask, template_mask, location = "unknown function") {
  supplied_mask_arr <- .mask_to_array(mask, location)
  template_mask_arr <- as.array(template_mask)

  if (!identical(supplied_mask_arr, template_mask_arr)) {
    .encoder_cli_abort(
      "mask does not match the template mask. Shared templates require identical voxel support and ordering.",
      class = "fmrilatent_error_dim", call = rlang::caller_env()
    )
  }

  if (inherits(mask, "LogicalNeuroVol")) {
    if (!isTRUE(all.equal(neuroim2::space(mask), neuroim2::space(template_mask)))) {
      .encoder_cli_abort("mask space does not match the template mask space.",
                         class = "fmrilatent_error_dim", call = rlang::caller_env())
    }
  }

  invisible(template_mask_arr)
}
