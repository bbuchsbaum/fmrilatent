# Transport/operator encode entry points

#' @include encode.R
NULL

#' Encode data using a shared basis asset and subject field operator
#'
#' @param x Numeric matrix (time x target samples) or a \code{NeuroVec}.
#' @param template Shared basis asset providing \code{basis_decoder()}.
#' @param field_operator Subject field operator mapping template field
#'   samples to observed native samples.
#' @param observation_operator Legacy alias for \code{field_operator}.
#' @details
#' The external field-operator contract is intentionally narrow. `fmrilatent`
#' consumes an operator-like object with:
#' \describe{
#'   \item{\code{n_source}, \code{n_target}}{Source and target dimensions.}
#'   \item{\code{source_domain_id}, \code{target_domain_id}}{Stable domain identifiers.}
#'   \item{\code{forward(z)}}{Applies the operator to template field samples.}
#'   \item{\code{adjoint_apply(y)}}{Applies the adjoint map.}
#'   \item{\code{provenance$target_mask}}{Optional target-domain mask when
#'   the caller does not pass \code{mask} explicitly.}
#' }
#'
#' On the main quadratic and sparse transport/AWPT paths, coefficient recovery
#' is matrix-free: `fmrilatent` uses the operator's forward and adjoint methods
#' rather than materializing the full subject decoder.
#' @param mask Optional volumetric target mask for the field-operator
#'   target domain.
#' @param domain Optional non-volumetric target domain, for example a surface
#'   geometry.
#' @param support Optional target support aligned to \code{domain}, for example
#'   vertex indices on a surface.
#'   At least one of \code{mask} or \code{support} must be available either
#'   explicitly or via \code{field_operator$provenance}.
#' @param lambda Ridge penalty strength.
#' @param spatial_lambda Strength of the spatial coefficient penalty.
#' @param spatial_penalty Optional coefficient-space roughness matrix or diagonal weights.
#' @param temporal_lambda Strength of temporal smoothing.
#' @param temporal_order Difference order used for temporal smoothing.
#' @param sparse_lambda Strength of optional sparse coefficient shrinkage.
#' @param sparse_mode Sparse penalty mode. Use \code{"group_l2"} for atom-wise group shrinkage.
#' @param max_iter Maximum iterations for sparse AWPT optimization.
#' @param tol Relative convergence tolerance for sparse AWPT optimization.
#' @param center Logical; if \code{TRUE}, center target samples before solving.
#' @param run_info Optional run metadata carried on the resulting latent object.
#' @param label Optional label stored in metadata.
#' @param units Optional declared [latent_units_record()] captured on the
#'   encoded object.
#' @param ... Reserved for future extensions.
#' @return A \code{TransportLatent} object.
#' @export
encode_operator <- function(x, template, field_operator = NULL, observation_operator = NULL, mask = NULL,
                            domain = NULL, support = NULL,
                            lambda = 0, center = TRUE, run_info = NULL,
                            spatial_lambda = lambda,
                            spatial_penalty = NULL,
                            temporal_lambda = 0,
                            temporal_order = 1L,
                            sparse_lambda = 0,
                            sparse_mode = c("none", "group_l2", "lasso"),
                            max_iter = 200L,
                            tol = 1e-6,
                            label = "", units = NULL, ...) {
  sparse_mode <- match.arg(sparse_mode)
  field_operator <- .resolve_field_operator(
    field_operator = field_operator,
    observation_operator = observation_operator,
    context = "encode_operator() field_operator"
  )
  obs_map <- .normalize_field_operator_map(field_operator)
  target_info <- .resolve_transport_target_support(
    mask = mask,
    domain = domain,
    support = support,
    field_operator = obs_map,
    location = "encode_operator"
  )
  mask_use <- target_info$mask
  support_use <- target_info$support
  domain_use <- target_info$domain
  if (inherits(x, "NeuroVec")) {
    if (is.null(mask_use)) {
      .encoder_cli_abort(
        "NeuroVec inputs currently require a volumetric target mask for encode_operator().",
        class = "fmrilatent_error_missing_target_mask"
      )
    }
    X <- t(neuroim2::series(x, mask_use != 0))
  } else {
    X <- as.matrix(x)
  }

  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    .encoder_cli_abort(
      "lambda must be a single non-negative finite number.",
      class = "fmrilatent_error_invalid_scalar"
    )
  }
  if (!is.numeric(spatial_lambda) || length(spatial_lambda) != 1L ||
      !is.finite(spatial_lambda) || spatial_lambda < 0) {
    .encoder_cli_abort(
      "spatial_lambda must be a single non-negative finite number.",
      class = "fmrilatent_error_invalid_scalar"
    )
  }
  if (!is.numeric(temporal_lambda) || length(temporal_lambda) != 1L ||
      !is.finite(temporal_lambda) || temporal_lambda < 0) {
    .encoder_cli_abort(
      "temporal_lambda must be a single non-negative finite number.",
      class = "fmrilatent_error_invalid_scalar"
    )
  }
  if (!is.numeric(sparse_lambda) || length(sparse_lambda) != 1L ||
      !is.finite(sparse_lambda) || sparse_lambda < 0) {
    .encoder_cli_abort(
      "sparse_lambda must be a single non-negative finite number.",
      class = "fmrilatent_error_invalid_scalar"
    )
  }

  basis_map <- .normalize_basis_decoder_map(template)
  D_map <- .compose_linear_maps(basis_map, obs_map, context = "subject decoder")
  analysis_transform_use <- .template_asset_analysis_transform(template) %||%
    .transport_identity_transform(D_map$n_source)

  if (ncol(X) != D_map$n_target) {
    .encoder_cli_abort(
      paste0("x has ", ncol(X), " target samples but the composed decoder has ",
             D_map$n_target, "."),
      class = "fmrilatent_error_dimension_mismatch"
    )
  }

  offset <- numeric(0)
  X_proj <- X
  if (isTRUE(center)) {
    offset <- colMeans(X)
    X_proj <- sweep(X, 2L, offset, "-")
  }

  penalty_use <- spatial_penalty
  if (is.null(penalty_use) && spatial_lambda > 0) {
    penalty_use <- template_roughness(template, coordinates = "analysis") %||%
      diag(D_map$n_source)
  }
  coeff_analysis <- if (sparse_lambda > 0 && sparse_mode != "none") {
    .solve_transport_coefficients_sparse_matrix_free(
      X = X_proj,
      decoder_map = D_map,
      spatial_lambda = spatial_lambda,
      spatial_penalty = penalty_use,
      temporal_lambda = temporal_lambda,
      temporal_order = temporal_order,
      run_info = run_info,
      sparse_lambda = sparse_lambda,
      sparse_mode = sparse_mode,
      max_iter = max_iter,
      tol = tol
    )
  } else {
    .solve_transport_coefficients_matrix_free(
      X = X_proj,
      decoder_map = D_map,
      spatial_lambda = spatial_lambda,
      spatial_penalty = penalty_use,
      temporal_lambda = temporal_lambda,
      temporal_order = temporal_order,
      run_info = run_info,
      max_iter = max_iter,
      tol = tol
    )
  }
  coeff_raw <- t(analysis_transform_use$to_raw(t(coeff_analysis)))

  value <- transport_latent(
    coeff_raw = coeff_raw,
    coeff_analysis = coeff_analysis,
    basis_asset = template,
    field_operator = field_operator,
    mask = mask_use,
    domain = domain_use,
    support = support_use,
    analysis_transform = analysis_transform_use,
    offset = offset,
    run_info = run_info,
    meta = list(
      family = "transport",
      label = label,
      lambda = lambda,
      spatial_lambda = spatial_lambda,
      temporal_lambda = temporal_lambda,
      temporal_order = temporal_order,
      sparse_lambda = sparse_lambda,
      sparse_mode = sparse_mode,
      method = if (is_awpt_template(template)) "AWPT" else NULL,
      centered = isTRUE(center),
      target_mask_source = target_info$source
    )
  )
  .with_latent_units(value, units)
}

#' Encode data using transport-backed latent semantics
#'
#' @param x Numeric matrix (time x target samples) or a \code{NeuroVec}.
#' @param basis_asset Shared basis asset.
#' @param field_operator Subject field operator. See
#'   \code{\link{encode_operator}()} for the required contract.
#' @param observation_operator Legacy alias for \code{field_operator}.
#' @param mask Optional volumetric target mask for the field-operator target domain.
#' @param domain Optional non-volumetric target domain.
#' @param support Optional target support aligned to \code{domain}.
#' @param lambda Ridge penalty strength.
#' @param spatial_lambda Strength of the spatial coefficient penalty.
#' @param spatial_penalty Optional coefficient-space roughness matrix or diagonal weights.
#' @param temporal_lambda Strength of temporal smoothing.
#' @param temporal_order Difference order used for temporal smoothing.
#' @param sparse_lambda Strength of optional sparse coefficient shrinkage.
#' @param sparse_mode Sparse penalty mode. Use \code{"group_l2"} for atom-wise group shrinkage.
#' @param max_iter Maximum iterations for sparse AWPT optimization.
#' @param tol Relative convergence tolerance for sparse AWPT optimization.
#' @param center Logical; if \code{TRUE}, center target samples before solving.
#' @param run_info Optional run metadata carried on the resulting latent object.
#' @param label Optional label stored in metadata.
#' @param units Optional declared [latent_units_record()] captured on the
#'   encoded object.
#' @param ... Reserved for future extensions.
#' @return A \code{TransportLatent} object.
#' @export
encode_transport <- function(x, basis_asset, field_operator = NULL, observation_operator = NULL, mask = NULL,
                             domain = NULL, support = NULL,
                             lambda = 0, center = TRUE, run_info = NULL,
                             spatial_lambda = lambda,
                             spatial_penalty = NULL,
                             temporal_lambda = 0,
                             temporal_order = 1L,
                             sparse_lambda = 0,
                             sparse_mode = c("none", "group_l2", "lasso"),
                             max_iter = 200L,
                             tol = 1e-6,
                             label = "", units = NULL, ...) {
  encode_operator(
    x = x,
    template = basis_asset,
    field_operator = field_operator,
    observation_operator = observation_operator,
    mask = mask,
    domain = domain,
    support = support,
    lambda = lambda,
    spatial_lambda = spatial_lambda,
    spatial_penalty = spatial_penalty,
    temporal_lambda = temporal_lambda,
    temporal_order = temporal_order,
    sparse_lambda = sparse_lambda,
    sparse_mode = sparse_mode,
    max_iter = max_iter,
    tol = tol,
    center = center,
    run_info = run_info,
    label = label,
    units = units,
    ...
  )
}

#' Encode data using the AWPT model
#'
#' @param x Numeric matrix (time x target samples) or a \code{NeuroVec}.
#' @param basis_asset An \code{AWPTBasisTemplate}.
#' @param field_operator Subject field operator. See
#'   \code{\link{encode_operator}()} for the required contract.
#' @param observation_operator Legacy alias for \code{field_operator}.
#' @param mask Optional volumetric target mask for the field-operator target domain.
#' @param domain Optional non-volumetric target domain.
#' @param support Optional target support aligned to \code{domain}.
#' @param spatial_lambda Strength of the anatomical roughness penalty.
#' @param temporal_lambda Strength of temporal smoothing.
#' @param temporal_order Temporal difference order for smoothing.
#' @param sparse_lambda Strength of optional sparse atom selection.
#' @param sparse_mode Sparse penalty mode. Use \code{"group_l2"} for atom-wise group shrinkage.
#' @param max_iter Maximum iterations for sparse AWPT optimization.
#' @param tol Relative convergence tolerance for sparse AWPT optimization.
#' @param center Logical; if \code{TRUE}, center target samples before solving.
#' @param run_info Optional run metadata; \code{run_lengths} control temporal blocks.
#' @param label Optional label stored in metadata.
#' @param units Optional declared [latent_units_record()] captured on the
#'   encoded object.
#' @param ... Reserved for future extensions.
#' @details
#' AWPT does not expose a separate ridge penalty. The returned metadata records
#' \code{lambda = 0}; the anatomical roughness penalty is recorded separately as
#' \code{spatial_lambda}.
#' @return A \code{TransportLatent} object with AWPT metadata.
#' @export
encode_awpt <- function(x, basis_asset, field_operator = NULL, observation_operator = NULL, mask = NULL,
                        domain = NULL, support = NULL,
                        spatial_lambda = 0, temporal_lambda = 0,
                        temporal_order = 1L,
                        sparse_lambda = 0,
                        sparse_mode = c("none", "group_l2", "lasso"),
                        max_iter = 200L,
                        tol = 1e-6,
                        center = TRUE,
                        run_info = NULL, label = "", units = NULL, ...) {
  sparse_mode <- match.arg(sparse_mode)
  if (!is_awpt_template(basis_asset)) {
    .encoder_cli_abort(
      "basis_asset must be an AWPT basis template.",
      class = "fmrilatent_error_invalid_basis_asset"
    )
  }
  encode_operator(
    x = x,
    template = basis_asset,
    field_operator = field_operator,
    observation_operator = observation_operator,
    mask = mask,
    domain = domain,
    support = support,
    lambda = 0,
    spatial_lambda = spatial_lambda,
    spatial_penalty = if (spatial_lambda > 0) {
      template_roughness(basis_asset, coordinates = "analysis")
    } else {
      NULL
    },
    temporal_lambda = temporal_lambda,
    temporal_order = temporal_order,
    sparse_lambda = sparse_lambda,
    sparse_mode = sparse_mode,
    max_iter = max_iter,
    tol = tol,
    center = center,
    run_info = run_info,
    label = label,
    units = units,
    ...
  )
}
