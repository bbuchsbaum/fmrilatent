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
    k_time <- .validate_auto_positive_count(spec$time$k, n_time,
                                            default = min(5L, n_time),
                                            name = "spec_time_bspline$k",
                                            n_name = "n_time")
    B_t <- as.matrix(build_bspline_basis(
      n_time = n_time,
      k = k_time,
      degree = spec$time$degree,
      include_intercept = spec$time$include_intercept,
      orthonormalize = spec$time$orthonormalize
    ))
  } else if (inherits(spec$time, "spec_time_dct")) {
    k_time <- .validate_auto_positive_count(spec$time$k, n_time, default = n_time,
                                            name = "spec_time_dct$k",
                                            n_name = "n_time")
    B_t <- as.matrix(build_dct_basis(n_time = n_time, k = k_time, norm = spec$time$norm))
  } else {
    .encoder_cli_abort(
      paste0("Unsupported spec_st$time class: ", class(spec$time)[1]),
      class = "fmrilatent_error_unsupported_spec"
    )
  }

  if (is.null(reduction)) reduction <- make_cluster_reduction(mask, seq_len(sum(mask_arr)))
  if (inherits(spec$space, "spec_space_slepian")) {
    k_space <- spec$space$k %||% 3L
    k_space <- .validate_positive_count(k_space, "spec_space_slepian$k")
    L_s <- as.matrix(lift(reduction, basis_slepian(k = k_space), k_neighbors = spec$space$k_neighbors))
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
      rec <- roi_subset_columns(rec, mask_arr, roi_mask)
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

  core_mode <- spec$core_mode %||% "auto"
  if (identical(core_mode, "explicit")) {
    # An explicit LatentNeuroVec is feasible whenever both separable factors
    # are materializable dense matrices, which is the case for every supported
    # spec_st time/space combination that reaches this point (wavelet_active
    # space already aborts above). The separable decoder output
    #   rec = B_t %*% core %*% t(L_s)
    # factors exactly as (B_t %*% core) %*% t(L_s), so we store
    # basis = B_t %*% core and loadings = L_s with no offset.
    factors_explicit <- is.numeric(as.matrix(B_t)) &&
      is.numeric(as.matrix(core)) && is.numeric(as.matrix(L_s))
    if (factors_explicit) {
      basis <- as.matrix(B_t) %*% as.matrix(core)
      loadings <- .loadings_for_materialize(
        Matrix::Matrix(as.matrix(L_s), sparse = FALSE),
        materialize,
        id_prefix = "st-loadings",
        label = "st-loadings"
      )
      spc <- .space_with_time_from_mask(mask, n_time, "encode_spec.spec_st")
      meta$core_mode <- "explicit"
      return(LatentNeuroVec(
        basis = Matrix::Matrix(basis, sparse = FALSE),
        loadings = loadings,
        space = spc,
        mask = mask,
        offset = numeric(0),
        label = label,
        meta = meta
      ))
    }
    .encoder_cli_warn(
      paste0("spec_st core_mode = 'explicit' requested but the separable ",
             "factors are not both materializable; returning an ImplicitLatent."),
      class = "fmrilatent_warning_explicit_core_unavailable"
    )
  }

  implicit_latent(
    coeff = list(core = core, B_t = B_t, L_s = L_s),
    decoder = decoder,
    meta = meta,
    mask = if (inherits(mask, "LogicalNeuroVol")) mask else mask_arr
  )
}
