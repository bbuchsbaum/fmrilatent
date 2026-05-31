# High-level encoding API and specs

.encoder_cli_warn <- function(message, class = "fmrilatent_warning_encoder", call = NULL) {
  warning(
    structure(
      list(message = message, call = call),
      class = unique(c(class, "fmrilatent_warning_encoder", "warning", "condition"))
    )
  )
}

# Robust Gram matrix solve with ridge/SVD fallback diagnostics.
.robust_gram_solve <- function(gram, rhs, ridge = 1e-8, context = "Gram solve") {
  gram_mat <- as.matrix(gram)
  rhs <- as.matrix(rhs)
  result <- tryCatch(
    solve(gram_mat, rhs),
    error = function(e) NULL
  )
  if (is.null(result)) {
    gram_ridge <- gram_mat
    if (ridge > 0) {
      .encoder_cli_warn(
        paste0(context, " was singular; retrying with ridge = ", ridge, "."),
        class = "fmrilatent_warning_gram_ridge"
      )
      diag(gram_ridge) <- diag(gram_ridge) + ridge
    }
    result <- tryCatch(
      solve(gram_ridge, rhs),
      error = function(e) {
        .encoder_cli_warn(
          paste0(context, " ridge solve failed; using SVD pseudoinverse fallback."),
          class = "fmrilatent_warning_gram_svd_fallback"
        )
        svd_g <- svd(gram_ridge)
        max_d <- if (length(svd_g$d)) max(svd_g$d) else 0
        tol <- max(dim(gram_ridge)) * max_d * .Machine$double.eps
        pos <- svd_g$d > tol
        if (any(pos)) {
          if (sum(pos) < length(svd_g$d)) {
            .encoder_cli_warn(
              paste0(context, " is rank deficient under SVD fallback (rank ",
                     sum(pos), " of ", length(svd_g$d), ")."),
              class = "fmrilatent_warning_gram_rank_deficient"
            )
          }
          svd_g$v[, pos, drop = FALSE] %*% (t(svd_g$u[, pos, drop = FALSE]) / svd_g$d[pos]) %*% rhs
        } else {
          .encoder_cli_warn(
            paste0(context, " collapsed to rank zero; returning zero coefficients."),
            class = "fmrilatent_warning_gram_rank_collapse"
          )
          matrix(0, nrow = ncol(gram_ridge), ncol = ncol(rhs))
        }
      }
    )
  }
  result
}

.basis_coefficients <- function(x, basis, ridge = 1e-8, context = "basis projection") {
  x_mat <- as.matrix(x)
  basis_mat <- as.matrix(basis)
  gram <- crossprod(basis_mat)
  rhs <- crossprod(basis_mat, x_mat)
  .robust_gram_solve(gram, rhs, ridge = ridge, context = context)
}

.temporal_loadings_from_basis <- function(x, basis, ridge = 1e-8, context = "temporal projection") {
  coeff <- .basis_coefficients(x, basis, ridge = ridge, context = context)
  Matrix::Matrix(t(coeff), sparse = FALSE)
}

# Default Slepian/DPSS component count from the time-bandwidth product.
# Mirrors the classic 2*N*W*tr - 1 rule, floored at 1. Shared by the
# temporal and spatiotemporal slepian encoders.
.slepian_default_k <- function(n_time, tr, bandwidth) {
  max(1L, min(as.integer(n_time), floor(2 * n_time * bandwidth * tr) - 1L))
}

# Least-squares coefficients of x (time x vox) onto a spatial atom dictionary
# (atoms x vox), returned as a time x atoms matrix. The Gram is formed over the
# voxel axis (atoms x atoms), so this complements .basis_coefficients (which
# forms the Gram over rows). Used by the HRBF and spec_st spatial atom paths.
.atom_coefficients <- function(x, atoms, ridge = 1e-8, context = "atom projection") {
  gram <- as.matrix(atoms %*% Matrix::t(atoms))
  rhs <- t(as.matrix(x) %*% Matrix::t(atoms))
  t(.robust_gram_solve(gram, rhs, ridge = ridge, context = context))
}

# Shared column-mean centering for encoders that build a per-voxel offset.
# Given x (time x vox) and a `center` flag, returns the offset and the centered
# matrix. When `center` is FALSE the offset is numeric(0) (the package
# convention for "no offset") and x is returned unchanged. When `center` is
# TRUE the offset is the supplied `offset` (validated for length) or the column
# means of x, and x is centered by subtracting it. Single source of truth for
# the column-mean centering used by the PCA and template-projection paths.
.encode_center <- function(x, center = FALSE, offset = NULL,
                           context = "encode centering") {
  x_mat <- as.matrix(x)
  if (!isTRUE(center)) {
    return(list(offset = numeric(0), x_centered = x_mat))
  }
  off <- if (!is.null(offset)) {
    off <- as.numeric(offset)
    if (length(off) != ncol(x_mat)) {
      .encoder_cli_abort(
        sprintf("%s: offset must have length %d (one per column).",
                context, ncol(x_mat)),
        class = "fmrilatent_error_invalid_argument"
      )
    }
    off
  } else {
    colMeans(x_mat)
  }
  list(offset = off, x_centered = sweep(x_mat, 2L, off, "-"))
}

.resolve_materialize <- function(materialize,
                                 supported = c("handle", "matrix"),
                                 default = "handle",
                                 context = "encode") {
  materialize <- match.arg(materialize, c("handle", "auto", "matrix"))
  if (identical(materialize, "auto")) {
    materialize <- default
  }
  if (!(materialize %in% supported)) {
    .encoder_cli_abort(
      paste0(context, " does not support materialize = \"", materialize, "\"."),
      class = "fmrilatent_error_unsupported_materialize"
    )
  }
  materialize
}

.loadings_for_materialize <- function(loadings, materialize, id_prefix, label = "loadings") {
  L <- Matrix::Matrix(loadings)
  if (!identical(materialize, "handle")) {
    return(L)
  }
  id <- .latent_handle_id(
    id_prefix,
    list(kind = "explicit_loadings", dim = dim(L), matrix = L)
  )
  handle <- new(
    "LoadingsHandle",
    id = id,
    dim = as.integer(dim(L)),
    kind = "explicit",
    spec = list(matrix = L),
    label = label
  )
  .latent_register_matrix(
    id, L, type = "loadings", overwrite = FALSE,
    fingerprint = .latent_handle_fingerprint(handle)
  )
  handle
}

.difference_penalty_matrix <- function(n, order = 1L) {
  n <- as.integer(n)
  order <- as.integer(order)
  if (order < 1L || n <= 1L) {
    return(matrix(0, nrow = n, ncol = n))
  }
  D <- diag(n)
  for (idx in seq_len(order)) {
    D <- D[-1L, , drop = FALSE] - D[-nrow(D), , drop = FALSE]
  }
  crossprod(D)
}

# Template helpers (.template_*, .symmetric_matrix_factor,
# .analysis_transform_from_metric, .analysis_loadings_from_transform,
# .transform_quadratic_form) live in R/encode_template.R.

.run_lengths_from_info <- function(n_time, run_info = NULL) {
  if (is.null(run_info)) {
    return(as.integer(n_time))
  }
  if (is.numeric(run_info) && is.null(dim(run_info))) {
    vals <- as.integer(run_info)
  } else if (is.list(run_info) && !is.null(run_info$run_lengths)) {
    vals <- as.integer(run_info$run_lengths)
  } else {
    vals <- as.integer(n_time)
  }
  if (sum(vals) != n_time) {
    .encoder_cli_abort(
      "run lengths must sum to the number of time points.",
      class = "fmrilatent_error_dim",
      call = rlang::caller_env()
    )
  }
  vals
}

.block_temporal_penalty <- function(n_time, temporal_order = 1L, run_info = NULL) {
  run_lengths <- .run_lengths_from_info(n_time, run_info)
  blocks <- lapply(run_lengths, function(len) .difference_penalty_matrix(len, order = temporal_order))
  if (length(blocks) == 1L) {
    return(blocks[[1L]])
  }
  as.matrix(Matrix::bdiag(blocks))
}

.normalize_penalty_matrix <- function(penalty, k, context = "penalty") {
  if (is.null(penalty)) {
    return(matrix(0, nrow = k, ncol = k))
  }
  if (is.vector(penalty) && is.null(dim(penalty))) {
    if (length(penalty) != k) {
      .encoder_cli_abort(
        paste0(context, " vector must have length ", k, "."),
        class = "fmrilatent_error_dim",
        call = rlang::caller_env()
      )
    }
    return(diag(as.numeric(penalty), nrow = k))
  }
  penalty <- as.matrix(penalty)
  if (!identical(dim(penalty), c(k, k))) {
    .encoder_cli_abort(
      paste0(context, " must have dimensions ", k, "x", k, "."),
      class = "fmrilatent_error_dim",
      call = rlang::caller_env()
    )
  }
  penalty
}

.largest_symmetric_eigenvalue <- function(mat) {
  mat <- as.matrix(mat)
  if (nrow(mat) == 0L) return(0)
  vals <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
  max(Re(vals), 0)
}

# Transport solvers (.solve_transport_coefficients*, .cg_transport_quadratic,
# .estimate_transport_lipschitz, .awpt_objective_matrix_free, plus the
# .apply_transport_*/.transport_* helpers and .frobenius_inner) live in
# R/encode_transport_solve.R.

# Validators (.validate_positive_count, .validate_nonnegative_count,
# .validate_positive_scalar, .validate_nonnegative_scalar,
# .validate_flag_scalar, .validate_hrbf_params) live in
# R/encoder_validators.R.

# Spec constructors (spec_time_*, spec_space_*, spec_st) live in
# R/encode_spec.R and are loaded by Collate order.

# --- Encode generic -----------------------------------------------------------

#' Encode data into a latent representation using a spec
#'
#' @param x Matrix (time x voxels in mask order).
#' @param spec Standard encode spec object created by `spec_time_*`,
#'   `spec_space_*`, or `spec_st`. AWPT specs created by
#'   [basis_awpt_wavelet()] describe shared templates and are intentionally not
#'   accepted by `encode()`; use [encode_awpt()] or [encode_operator()] for
#'   transport-backed AWPT fits.
#' @param mask LogicalNeuroVol or logical array (required for spatial pieces).
#' @param reduction Optional GraphReduction (for spatial specs).
#' @param materialize "handle", "matrix", or "auto" (default "handle").
#' @param label Optional label.
#' @param ... Additional arguments passed to methods.
#' @return The return class depends on the spec family:
#'   \describe{
#'     \item{Explicit spatial families}{`spec_space_slepian`,
#'       `spec_space_heat`, `spec_space_hrbf`, `spec_space_pca`, and
#'       `spec_space_wavelet_active` return a [LatentNeuroVec], which is a
#'       concrete `ExplicitLatent` (the virtual marker defined at
#'       `R/all_class.R`). It stores an explicit `basis x loadings + offset`
#'       factorization.}
#'     \item{Explicit temporal families}{`spec_time_slepian`,
#'       `spec_time_dct`, and `spec_time_bspline` likewise return a
#'       [LatentNeuroVec] (`ExplicitLatent`).}
#'     \item{Spatiotemporal (`spec_st`)}{**Always** returns an
#'       `ImplicitLatent` (a decoder-only / separable representation), even
#'       when both the time and space factors are fully explicit dense
#'       bases. `ImplicitLatent` is *not* an `ExplicitLatent`.}
#'     \item{Transport (`spec_transport` / AWPT encoders)}{return a
#'       `TransportLatent`, which is also *not* an `ExplicitLatent`.}
#'   }
#'   In short: `ExplicitLatent` is the virtual S4 marker inherited by
#'   [LatentNeuroVec]; `ImplicitLatent` and `TransportLatent` are S3
#'   classes that deliberately do not inherit it.
#' @section Dispatch model:
#'   For standard in-mask matrix encoders, `encode()` routes to the S3 generic
#'   [encode_spec()], which dispatches on the spec class and builds the latent
#'   object directly. Transport-backed AWPT is the explicit exception: an
#'   AWPT basis spec is used to build a shared template, while the subject fit
#'   also requires a `basis_asset` and `field_operator`. Those assets are not
#'   part of the standard `encode_spec()` signature, so AWPT enters through
#'   [encode_awpt()] or [encode_operator()] instead of [encode()].
#' @section Offset and centering contract:
#'   A [LatentNeuroVec] reconstructs its data as
#'   `basis %*% t(loadings) + offset`, where `offset` is a per-voxel vector
#'   (length = number of in-mask voxels) added back after the low-rank term.
#'   The offset is owned by the encoder's *lift* step, which is the single
#'   place that decides whether and how to center:
#'   \describe{
#'     \item{Families that populate `offset`}{Only `spec_space_pca` produces a
#'       non-empty offset, and only when `center = TRUE` (the default): the
#'       per-voxel column means are removed before the PCA and stored in the
#'       `offset` slot so reconstruction restores them. The mean removal is
#'       performed exactly once, inside `lift.basis_pca` (see the
#'       `fmrilatent.mean_scores` attribute it returns); the encode caller does
#'       not re-center. With `center = FALSE`, PCA stores `offset = numeric(0)`.}
#'     \item{Families that never center}{All other explicit families
#'       (`spec_space_slepian`, `spec_space_heat`, `spec_space_hrbf`,
#'       `spec_time_slepian`, `spec_time_dct`, `spec_time_bspline`) store
#'       `offset = numeric(0)`, i.e. no offset.}
#'   }
#'   By convention `offset = numeric(0)` means "no offset" and reconstruction
#'   treats it as a zero vector; a populated `offset` always has one entry per
#'   in-mask voxel. The shared `.encode_center()` helper is the single
#'   implementation of column-mean centering used by the offset-producing
#'   paths.
#' @export
encode <- function(x, spec, mask, reduction = NULL,
                   materialize = c("auto", "handle", "matrix"),
                   label = "", ...) {
  UseMethod("encode")
}

#' @export
encode.default <- function(x, spec, mask, reduction = NULL,
                           materialize = c("auto", "handle", "matrix"),
                           label = "", ...) {
  .encoder_cli_abort(
    paste0("No encode method for class: ", paste(class(x), collapse = ",")),
    class = "fmrilatent_error_unsupported_encode_input"
  )
}

#' @export
encode.matrix <- function(x, spec, mask, reduction = NULL,
                          materialize = c("auto", "handle", "matrix"),
                          label = "", ...) {
  materialize <- match.arg(materialize)
  encode_spec(
    x, spec,
    mask = mask,
    reduction = reduction,
    materialize = materialize,
    label = label,
    ...
  )
}

#' @export
encode.NeuroVec <- function(x, spec, mask, reduction = NULL,
                            materialize = c("handle", "auto", "matrix"),
                            label = "", ...) {
  materialize <- match.arg(materialize)
  X <- t(neuroim2::series(x, mask != 0))  # series returns voxels x time, transpose to time x voxels
  encode(X, spec, mask = mask, reduction = reduction,
         materialize = materialize, label = label, ...)
}

# --- Factory helper -----------------------------------------------------------

#' Simple factory to build a spec and encode in one call
#'
#' @param family Character scalar naming one of the standard `encode()`
#'   families. See **Accepted family names** for the canonical names and
#'   supported aliases.
#' @param x Data matrix (time x voxels).
#' @param mask Mask (required for spatial families).
#' @param reduction Optional GraphReduction for spatial specs.
#' @param ... Passed to spec constructors and encode().
#' @param materialize "handle", "matrix", or "auto" (default "handle").
#' @param label Optional label for the resulting object.
#' @return The class follows the same per-family contract as [encode()]:
#'   explicit spatial families and explicit temporal families return a
#'   [LatentNeuroVec] (a concrete `ExplicitLatent`); the spatiotemporal
#'   families (`st_slepian`, `st_bspline_hrbf`) build a
#'   `spec_st` and therefore always return an `ImplicitLatent`. See the
#'   `@return` section of [encode()] for the full taxonomy.
#'
#' @section Accepted family names:
#' Canonical names are listed first; aliases in parentheses are accepted for
#' compatibility.
#' \describe{
#'   \item{Temporal}{`time_dct` (`dct_time`), `time_slepian`
#'     (`slepian_time`).}
#'   \item{Spatial}{`space_slepian` (`slepian_space`), `space_pca`
#'     (`pca_space`), `space_parcel` (`parcel_space`), `space_heat`
#'     (`heat_space`), `space_hrbf` (`hrbf_space`), `space_wavelet_active`
#'     (`wavelet_active`), and `hierarchical`.}
#'   \item{Spatiotemporal}{`st` (requires explicit `time` and `space` specs),
#'     `st_slepian` (`slepian_st`), and `st_bspline_hrbf`
#'     (`bspline_hrbf_st`).}
#' }
#'
#' AWPT is intentionally not a `latent_factory()` family because it requires a
#' shared `basis_asset` and a subject `field_operator`; use [encode_awpt()] or
#' [encode_operator()] for AWPT subject fitting.
#' @export
latent_factory <- function(family, x, mask, reduction = NULL, ..., materialize = "auto", label = "") {
  if (identical(family, "awpt")) {
    .encoder_cli_abort(
      paste0("latent_factory() does not support AWPT because AWPT requires a shared ",
             "basis_asset and subject field_operator. Use encode_awpt() or ",
             "encode_operator() instead."),
      class = "fmrilatent_error_unsupported_operation",
      call = rlang::caller_env()
    )
  }
  family_alias <- c(
    dct_time = "time_dct",
    slepian_time = "time_slepian",
    slepian_space = "space_slepian",
    pca_space = "space_pca",
    parcel_space = "space_parcel",
    heat_space = "space_heat",
    hrbf_space = "space_hrbf",
    wavelet_active = "space_wavelet_active",
    slepian_st = "st_slepian",
    bspline_hrbf_st = "st_bspline_hrbf"
  )
  choices <- c(
    "time_dct", "time_slepian", "space_slepian", "space_pca",
    "space_parcel", "space_heat", "space_hrbf", "space_wavelet_active",
    "st", "st_slepian", "st_bspline_hrbf", "hierarchical",
    names(family_alias)
  )
  family <- match.arg(family, choices)
  family_canonical <- unname(family_alias[family])
  if (!is.na(family_canonical)) {
    family <- family_canonical
  }
  spec <- switch(
    family,
    time_dct = spec_time_dct(...),
    time_slepian = spec_time_slepian(...),
    space_slepian = spec_space_slepian(...),
    space_pca = spec_space_pca(...),
    space_heat = spec_space_heat(...),
    space_hrbf = spec_space_hrbf(...),
    st = {
      args <- list(...)
      if (is.null(args$time) || is.null(args$space)) {
        .encoder_cli_abort(
          "latent_factory('st') requires explicit 'time' and 'space' specs.",
          class = "fmrilatent_error_missing_argument",
          call = rlang::caller_env()
        )
      }
      spec_st(time = args$time, space = args$space,
              core_mode = args$core_mode %||% "auto")
    },
    st_bspline_hrbf = {
      args <- list(...)
      time_spec <- args$time %||% spec_time_bspline(k = args$k_time %||% args$k %||% 5L,
                                                    degree = args$degree %||% 3L,
                                                    include_intercept = args$include_intercept %||% FALSE,
                                                    orthonormalize = TRUE)
      space_spec <- args$space %||% spec_space_hrbf(params = args$params %||% list())
      spec_st(time = time_spec, space = space_spec)
    },
    space_parcel = spec_space_parcel(...),
    space_wavelet_active = spec_space_wavelet_active(...),
    hierarchical = spec_hierarchical_template(...),
    st_slepian = {
      args <- list(...)
      time_spec <- args$time %||% spec_time_slepian(
        tr = args$tr, bandwidth = args$bandwidth %||% 0.1,
        k = args$k_time %||% args$k)
      space_spec <- args$space %||% spec_space_slepian(
        k = args$k_space %||% args$k %||% 3L,
        k_neighbors = args$k_neighbors %||% 6L)
      spec_st(time = time_spec, space = space_spec)
    }
  )
  encode(x, spec, mask = mask, reduction = reduction, materialize = materialize, label = label)
}
#' Dispatch standard encoding based on spec type
#'
#' @param x Data matrix.
#' @param spec Spec object.
#' @param ... Additional arguments passed to methods.
#' @return Encoded representation.
#' @details
#' `encode_spec()` is the S3 dispatch seam for standard `encode()` specs:
#' temporal, spatial, hierarchical, parcel, and spatiotemporal specs. AWPT's
#' `basis_awpt_wavelet()` object is a template-construction spec, not a complete
#' subject-encoding spec, because AWPT fitting additionally needs a shared
#' template/basis asset and subject field or observation operators. AWPT
#' therefore uses the parallel transport API [encode_awpt()] /
#' [encode_operator()].
#' @export
encode_spec <- function(x, spec, ...) UseMethod("encode_spec", spec)

#' @exportS3Method
encode_spec.default <- function(x, spec, ...) {
  .encoder_cli_abort(
    paste0("Unknown spec class: ", paste(class(spec), collapse = ",")),
    class = "fmrilatent_error_unsupported_spec"
  )
}

#' @exportS3Method
encode_spec.spec_awpt_wavelet <- function(x, spec, mask, reduction, materialize, label, ...) {
  .encoder_cli_abort(
    "basis_awpt_wavelet() creates an AWPT template spec, not a complete encode() spec. AWPT subject fitting uses encode_awpt() or encode_operator() because it requires a shared basis_asset and subject field_operator. See ?encode_awpt for details.",
    class = "fmrilatent_error_unsupported_spec"
  )
}
