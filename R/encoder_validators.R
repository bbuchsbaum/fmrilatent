# Internal argument validators used by spec constructors and encoders.
#
# These all follow the same pattern: validate the scalar shape/sign,
# return a coerced canonical value, or stop() with a clear message.
# Used directly by R/encode_spec.R and indirectly by encoder families
# in awpt.R, hrbf.R, parcel_basis.R, etc.

#' @include encode.R
NULL

.validate_positive_count <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 1 || !isTRUE(all.equal(x, round(x)))) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  as.integer(round(x))
}

.validate_nonnegative_count <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0 || !isTRUE(all.equal(x, round(x)))) {
    stop(name, " must be a non-negative integer.", call. = FALSE)
  }
  as.integer(round(x))
}

.validate_positive_scalar <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    stop(name, " must be a positive finite number.", call. = FALSE)
  }
  as.numeric(x)
}

.validate_nonnegative_scalar <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
    stop(name, " must be a non-negative finite number.", call. = FALSE)
  }
  as.numeric(x)
}

.validate_flag_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  isTRUE(x)
}

.validate_hrbf_params <- function(params) {
  if (!is.list(params)) {
    stop("params must be a list.", call. = FALSE)
  }
  params_clean <- params
  if (!is.null(params_clean$sigma0)) {
    params_clean$sigma0 <- .validate_positive_scalar(params_clean$sigma0, "params$sigma0")
  }
  if (!is.null(params_clean$levels)) {
    params_clean$levels <- .validate_nonnegative_count(params_clean$levels, "params$levels")
  }
  if (!is.null(params_clean$radius_factor)) {
    params_clean$radius_factor <- .validate_positive_scalar(params_clean$radius_factor, "params$radius_factor")
  }
  if (!is.null(params_clean$num_extra_fine_levels)) {
    params_clean$num_extra_fine_levels <- .validate_nonnegative_count(
      params_clean$num_extra_fine_levels,
      "params$num_extra_fine_levels"
    )
  }
  if (!is.null(params_clean$seed)) {
    params_clean$seed <- .validate_nonnegative_count(params_clean$seed, "params$seed")
  }
  if (!is.null(params_clean$kernel_type) &&
      !params_clean$kernel_type %in% c("gaussian", "wendland_c4", "wendland_c6")) {
    stop("params$kernel_type must be one of: gaussian, wendland_c4, wendland_c6.", call. = FALSE)
  }
  if (!is.null(params_clean$kernel_type_fine_levels) &&
      !params_clean$kernel_type_fine_levels %in% c("gaussian", "wendland_c4", "wendland_c6")) {
    stop(
      "params$kernel_type_fine_levels must be one of: gaussian, wendland_c4, wendland_c6.",
      call. = FALSE
    )
  }
  params_clean
}
