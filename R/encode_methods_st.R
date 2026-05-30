# Spatiotemporal encode_spec methods

#' @include encode.R
NULL

#' @exportS3Method
encode_spec.spec_st <- function(x, spec, mask, reduction, materialize, label, ...) {
  materialize <- .resolve_materialize(
    materialize,
    supported = c("matrix", "handle"),
    default = "matrix",
    context = "encode_spec.spec_st"
  )
  mask_arr <- .extract_mask_array(mask, "encode_spec.spec_st")
  n_time <- nrow(x)

  if (inherits(spec$time, "spec_time_slepian")) {
    k_time <- spec$time$k %||% .slepian_default_k(n_time, spec$time$tr, spec$time$bandwidth)
    B_t <- basis_mat(slepian_temporal_handle(n_time = n_time, tr = spec$time$tr,
                                             bandwidth = spec$time$bandwidth, k = k_time,
                                             backend = spec$time$backend))
  } else if (inherits(spec$time, "spec_time_bspline")) {
    k_time <- spec$time$k
    B_t <- as.matrix(build_bspline_basis(
      n_time = n_time,
      k = k_time,
      degree = spec$time$degree,
      include_intercept = spec$time$include_intercept,
      orthonormalize = spec$time$orthonormalize
    ))
    if (!spec$time$orthonormalize) {
      q <- qr(B_t)
      B_t <- qr.Q(q)
    }
  } else if (inherits(spec$time, "spec_time_dct")) {
    k_time <- spec$time$k
    if (k_time > n_time) {
      .encoder_cli_abort(
        paste0("spec_time_dct$k (", k_time, ") cannot exceed n_time (", n_time, ")."),
        class = "fmrilatent_error_invalid_count"
      )
    }
    B_t <- as.matrix(build_dct_basis(n_time = n_time, k = k_time, norm = spec$time$norm))
  } else {
    .encoder_cli_abort(
      paste0("Unsupported spec_st$time class: ", class(spec$time)[1]),
      class = "fmrilatent_error_unsupported_spec"
    )
  }

  if (is.null(reduction)) reduction <- make_cluster_reduction(mask, seq_len(sum(mask_arr)))
  if (inherits(spec$space, "spec_space_slepian")) {
    L_s <- as.matrix(lift(reduction, basis_slepian(k = spec$space$k), k_neighbors = spec$space$k_neighbors))
    B_atoms <- t(L_s)
  } else if (inherits(spec$space, "spec_space_hrbf")) {
    params <- spec$space$params %||% list()
    B_atoms <- as.matrix(hrbf_generate_basis(params, mask))
    L_s <- t(B_atoms)
  } else if (inherits(spec$space, "spec_space_wavelet_active")) {
    .encoder_cli_abort(
      "spec_st with space = spec_space_wavelet_active not yet supported; use spec_space_wavelet_active directly via encode().",
      class = "fmrilatent_error_unsupported_spec"
    )
  } else {
    .encoder_cli_abort(
      paste0("Unsupported spec_st$space class: ", class(spec$space)[1]),
      class = "fmrilatent_error_unsupported_spec"
    )
  }

  C_atoms <- .atom_coefficients(x, B_atoms, context = "encode_spec.spec_st spatial projection")

  core <- .basis_coefficients(
    C_atoms,
    B_t,
    context = "encode_spec.spec_st temporal projection"
  )

  decoder <- function(time_idx = NULL, roi_mask = NULL, ...) {
    t_sel <- if (is.null(time_idx)) seq_len(n_time) else as.integer(time_idx)
    B_sel <- B_t[t_sel, , drop = FALSE]
    rec <- B_sel %*% core %*% t(L_s)
    if (!is.null(roi_mask)) {
      global_idx <- which(as.logical(mask_arr))
      roi_global <- which(as.logical(roi_mask))
      col_keep <- which(global_idx %in% roi_global)
      rec <- rec[, col_keep, drop = FALSE]
    }
    rec
  }

  meta <- list(
    family = "st_separable",
    time = class(spec$time)[1],
    space = class(spec$space)[1],
    k_time = ncol(B_t),
    k_space = ncol(L_s),
    label = label
  )

  implicit_latent(
    coeff = list(core = core, B_t = B_t, L_s = L_s),
    decoder = decoder,
    meta = meta,
    mask = mask_arr
  )
}
