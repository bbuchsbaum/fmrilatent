# Internal argument validators used by spec constructors and encoders.
#
# These all follow the same pattern: validate the scalar shape/sign,
# return a coerced canonical value, or cli::cli_abort() with a clear message.
# Used directly by R/encode_spec.R and indirectly by encoder families
# in awpt.R, hrbf.R, parcel_basis.R, etc.

#' @include encode.R
NULL

.validate_positive_count <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 1 ||
      x != floor(x) || x > .Machine$integer.max) {
    .encoder_cli_abort(
      paste0(name, " must be a positive integer."),
      class = "fmrilatent_error_invalid_count"
    )
  }
  as.integer(x)
}

.validate_optional_positive_count <- function(x, name) {
  if (is.null(x)) {
    return(NULL)
  }
  .validate_positive_count(x, name)
}

.validate_auto_positive_count <- function(x, n, default = NULL, name = "k",
                                          n_name = "n") {
  n <- .validate_positive_count(n, n_name)
  if (is.null(x)) {
    x <- default %||% n
  }
  x <- .validate_positive_count(x, name)
  if (x > n) {
    .encoder_cli_abort(
      paste0(name, " (", x, ") cannot exceed ", n_name, " (", n, ")."),
      class = "fmrilatent_error_invalid_count"
    )
  }
  x
}

.validate_nonnegative_count <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0 ||
      x != floor(x) || x > .Machine$integer.max) {
    .encoder_cli_abort(
      paste0(name, " must be a non-negative integer."),
      class = "fmrilatent_error_invalid_count"
    )
  }
  as.integer(x)
}

.validate_positive_scalar <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    .encoder_cli_abort(
      paste0(name, " must be a positive finite number."),
      class = "fmrilatent_error_invalid_scalar"
    )
  }
  as.numeric(x)
}

.validate_nonnegative_scalar <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
    .encoder_cli_abort(
      paste0(name, " must be a non-negative finite number."),
      class = "fmrilatent_error_invalid_scalar"
    )
  }
  as.numeric(x)
}

.validate_flag_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .encoder_cli_abort(
      paste0(name, " must be TRUE or FALSE."),
      class = "fmrilatent_error_invalid_flag"
    )
  }
  isTRUE(x)
}

.validate_hrbf_params <- function(params) {
  if (!is.list(params)) {
    .encoder_cli_abort("params must be a list.", class = "fmrilatent_error_invalid_params")
  }
  params_clean <- params
  param_names <- names(params_clean)
  if ("sigma0" %in% param_names) {
    params_clean$sigma0 <- .validate_positive_scalar(params_clean[["sigma0"]], "params$sigma0")
  }
  if ("levels" %in% param_names) {
    params_clean$levels <- .validate_nonnegative_count(params_clean[["levels"]], "params$levels")
  }
  if ("radius_factor" %in% param_names) {
    params_clean$radius_factor <- .validate_positive_scalar(params_clean[["radius_factor"]], "params$radius_factor")
  }
  if ("num_extra_fine_levels" %in% param_names) {
    params_clean$num_extra_fine_levels <- .validate_nonnegative_count(
      params_clean[["num_extra_fine_levels"]],
      "params$num_extra_fine_levels"
    )
  }
  if ("seed" %in% param_names) {
    params_clean$seed <- .validate_nonnegative_count(params_clean[["seed"]], "params$seed")
  }
  if ("kernel_type" %in% param_names &&
      !params_clean[["kernel_type"]] %in% c("gaussian", "wendland_c4", "wendland_c6")) {
    .encoder_cli_abort(
      "params$kernel_type must be one of: gaussian, wendland_c4, wendland_c6.",
      class = "fmrilatent_error_invalid_params"
    )
  }
  if ("kernel_type_fine_levels" %in% param_names &&
      !params_clean[["kernel_type_fine_levels"]] %in% c("gaussian", "wendland_c4", "wendland_c6")) {
    .encoder_cli_abort(
      "params$kernel_type_fine_levels must be one of: gaussian, wendland_c4, wendland_c6.",
      class = "fmrilatent_error_invalid_params"
    )
  }
  params_clean
}
