#' @include RcppExports.R encoder_utilities.R
NULL

# Classed condition wrappers for generated Rcpp entry points. RcppExports.R is
# generated, so keep the package-specific condition policy in this file.

.cpp_boundary_raw <- local({
  list(
    cdf53_time_lift = cdf53_time_lift,
    active_pencil_wavelet = active_pencil_wavelet,
    get_morton_ordered_indices_rcpp = get_morton_ordered_indices_rcpp,
    precompute_haar_scalings_rcpp = precompute_haar_scalings_rcpp,
    forward_lift_rcpp = forward_lift_rcpp,
    inverse_lift_rcpp = inverse_lift_rcpp,
    get_valid_finest_blocks_rcpp = get_valid_finest_blocks_rcpp,
    hrbf_atoms_rcpp = hrbf_atoms_rcpp,
    generate_dpss_basis_rcpp = generate_dpss_basis_rcpp,
    generate_dpss_tridiag_rcpp = generate_dpss_tridiag_rcpp
  )
})

.cpp_boundary_error_class <- function(message) {
  if (grepl("too large for Morton codes", message, fixed = TRUE)) {
    return(c("fmrilatent_error_morton_overflow", "fmrilatent_error_cpp_boundary"))
  }
  if (grepl("failed", message, fixed = TRUE)) {
    return(c("fmrilatent_error_cpp_boundary", "fmrilatent_error_computation"))
  }
  c("fmrilatent_error_cpp_boundary", "fmrilatent_error_invalid_argument")
}

.cpp_boundary_call <- function(expr, call = rlang::caller_env()) {
  tryCatch(
    expr,
    error = function(e) {
      if (inherits(e, "fmrilatent_error_cpp_boundary")) {
        stop(e)
      }
      msg <- conditionMessage(e)
      .encoder_cli_abort(
        msg,
        class = .cpp_boundary_error_class(msg),
        call = call
      )
    }
  )
}

cdf53_time_lift <- function(X, levels, forward) {
  .cpp_boundary_call(
    .cpp_boundary_raw$cdf53_time_lift(X, levels, forward),
    call = rlang::caller_env()
  )
}

active_pencil_wavelet <- function(data_voxels, coords, dims, levels, forward) {
  .cpp_boundary_call(
    .cpp_boundary_raw$active_pencil_wavelet(data_voxels, coords, dims, levels, forward),
    call = rlang::caller_env()
  )
}

get_morton_ordered_indices_rcpp <- function(mask) {
  .cpp_boundary_call(
    .cpp_boundary_raw$get_morton_ordered_indices_rcpp(mask),
    call = rlang::caller_env()
  )
}

precompute_haar_scalings_rcpp <- function(mask, levels) {
  .cpp_boundary_call(
    .cpp_boundary_raw$precompute_haar_scalings_rcpp(mask, levels),
    call = rlang::caller_env()
  )
}

forward_lift_rcpp <- function(data_morton, mask_flat_morton, mask_dims, levels, scalings) {
  .cpp_boundary_call(
    .cpp_boundary_raw$forward_lift_rcpp(data_morton, mask_flat_morton, mask_dims, levels, scalings),
    call = rlang::caller_env()
  )
}

inverse_lift_rcpp <- function(root_coeff, detail_vecs, mask_flat_morton, mask_dims, levels, scalings) {
  .cpp_boundary_call(
    .cpp_boundary_raw$inverse_lift_rcpp(root_coeff, detail_vecs, mask_flat_morton, mask_dims, levels, scalings),
    call = rlang::caller_env()
  )
}

get_valid_finest_blocks_rcpp <- function(mask) {
  .cpp_boundary_call(
    .cpp_boundary_raw$get_valid_finest_blocks_rcpp(mask),
    call = rlang::caller_env()
  )
}

hrbf_atoms_rcpp <- function(mask_xyz_world, centres_xyz_world, sigma_vec_mm, kernel_type, value_threshold = 1e-8) {
  .cpp_boundary_call(
    .cpp_boundary_raw$hrbf_atoms_rcpp(
      mask_xyz_world,
      centres_xyz_world,
      sigma_vec_mm,
      kernel_type,
      value_threshold
    ),
    call = rlang::caller_env()
  )
}

generate_dpss_basis_rcpp <- function(n, NW, k) {
  .cpp_boundary_call(
    .cpp_boundary_raw$generate_dpss_basis_rcpp(n, NW, k),
    call = rlang::caller_env()
  )
}

generate_dpss_tridiag_rcpp <- function(n, NW, k) {
  .cpp_boundary_call(
    .cpp_boundary_raw$generate_dpss_tridiag_rcpp(n, NW, k),
    call = rlang::caller_env()
  )
}
