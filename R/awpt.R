# AWPT basis assets and helpers

#' @include reduction.R heat_wavelet.R encode.R
#' @importFrom Matrix crossprod Cholesky Diagonal
#' @importFrom methods setMethod
NULL

#' AWPT wave-packet basis specification
#'
#' @param scales Numeric vector of anatomical wave-packet scales.
#' @param order Polynomial approximation order for the underlying heat-wavelet construction.
#' @param threshold Threshold for small coefficients.
#' @param k_neighbors Graph neighborhood parameter.
#' @param penalty_rule Rule used to convert scales into roughness weights.
#' @param custom_weights Optional explicit weights matching \code{scales}.
#' @return A \code{spec_awpt_wavelet} object.
#' @export
basis_awpt_wavelet <- function(scales = c(1, 2, 4, 8), order = 30L, threshold = 1e-6,
                               k_neighbors = 6L,
                               penalty_rule = c("inverse_scale", "inverse_scale_sq",
                                                "scale", "none", "custom"),
                               custom_weights = NULL) {
  penalty_rule <- match.arg(penalty_rule)
  if (penalty_rule == "custom" && (is.null(custom_weights) || length(custom_weights) != length(scales))) {
    stop("custom_weights must be supplied with one value per scale when penalty_rule = 'custom'.",
         call. = FALSE)
  }
  structure(
    list(
      scales = as.numeric(scales),
      order = as.integer(order),
      threshold = threshold,
      k_neighbors = as.integer(k_neighbors),
      penalty_rule = penalty_rule,
      custom_weights = custom_weights
    ),
    class = "spec_awpt_wavelet"
  )
}

.coerce_awpt_reduction <- function(parcellation) {
  if (inherits(parcellation, "ClusterReduction")) {
    parcellation
  } else if (inherits(parcellation, "ClusteredNeuroVol")) {
    as_cluster_reduction(parcellation)
  } else {
    stop("parcellation must be a ClusterReduction or ClusteredNeuroVol.", call. = FALSE)
  }
}

.awpt_scale_weights <- function(scales, penalty_rule, custom_weights = NULL) {
  switch(
    penalty_rule,
    inverse_scale = 1 / scales,
    inverse_scale_sq = 1 / (scales ^ 2),
    scale = scales,
    none = rep(0, length(scales)),
    custom = as.numeric(custom_weights)
  )
}

.awpt_atom_metadata <- function(reduction, scales, scale_weights) {
  counts <- tabulate(match(reduction@map, reduction@cluster_ids), nbins = length(reduction@cluster_ids))
  rows <- vector("list", length(reduction@cluster_ids) * length(scales))
  idx <- 1L
  col_offset <- 0L
  for (cid_idx in seq_along(reduction@cluster_ids)) {
    n_loc <- counts[[cid_idx]]
    if (n_loc == 0L) next
    cid <- reduction@cluster_ids[[cid_idx]]
    for (scale_idx in seq_along(scales)) {
      col_ids <- seq.int(col_offset + 1L, col_offset + n_loc)
      rows[[idx]] <- data.frame(
        col_id = col_ids,
        cluster_id = cid,
        scale = scales[[scale_idx]],
        scale_index = scale_idx,
        roughness_weight = scale_weights[[scale_idx]],
        stringsAsFactors = FALSE
      )
      col_offset <- col_offset + n_loc
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows[seq_len(idx - 1L)])
}

.awpt_as_square_matrix <- function(x, n, context) {
  x <- as.matrix(x)
  if (!identical(dim(x), c(n, n))) {
    stop(context, " must have dimensions ", n, "x", n, ".", call. = FALSE)
  }
  x
}

.awpt_graph_laplacian <- function(conductance, n_field) {
  W <- .awpt_as_square_matrix(conductance, n_field, context = "conductance")
  if (!isTRUE(all.equal(W, t(W), tolerance = 1e-8))) {
    stop("conductance must be symmetric.", call. = FALSE)
  }
  diag(rowSums(W)) - W
}

.awpt_enforce_symmetric <- function(mat, tol = 1e-8, context = "matrix") {
  mat <- as.matrix(mat)
  if (!identical(dim(mat)[1], dim(mat)[2])) {
    stop(context, " must be square.", call. = FALSE)
  }
  sym <- 0.5 * (mat + t(mat))
  if (!isTRUE(all.equal(sym, t(sym), tolerance = tol))) {
    stop(context, " could not be symmetrized stably.", call. = FALSE)
  }
  sym
}

.awpt_psd_project <- function(mat, tol = 1e-8) {
  eig <- eigen(mat, symmetric = TRUE)
  vals <- pmax(Re(eig$values), tol)
  eig$vectors %*% (vals * t(eig$vectors))
}

.awpt_log_euclidean_mean <- function(mats, tol = 1e-8) {
  logm_list <- lapply(mats, function(mat) {
    eig <- eigen(mat, symmetric = TRUE)
    vals <- pmax(Re(eig$values), tol)
    eig$vectors %*% (log(vals) * t(eig$vectors))
  })
  mean_log <- Reduce(`+`, logm_list) / length(logm_list)
  eig <- eigen(mean_log, symmetric = TRUE)
  eig$vectors %*% (exp(Re(eig$values)) * t(eig$vectors))
}

.awpt_roughness_from_scale <- function(reduction, basis_spec, n_coeff) {
  scale_weights <- .awpt_scale_weights(
    basis_spec$scales,
    basis_spec$penalty_rule,
    custom_weights = basis_spec$custom_weights
  )
  atoms <- .awpt_atom_metadata(reduction, basis_spec$scales, scale_weights)
  if (nrow(atoms) != n_coeff) {
    stop("AWPT atom metadata does not align with the basis rank.", call. = FALSE)
  }
  list(
    roughness = Matrix::Diagonal(x = atoms$roughness_weight),
    atoms = atoms,
    source = "scale_surrogate",
    field_operator = NULL
  )
}

.awpt_roughness_from_field_operator <- function(loadings, reduction, basis_spec,
                                                field_operator, source) {
  B <- as.matrix(loadings)
  n_field <- nrow(B)
  n_coeff <- ncol(B)
  L_field <- .awpt_as_square_matrix(field_operator, n_field, context = "anatomical_operator")
  Q <- crossprod(B, L_field %*% B)
  scale_weights <- diag(Q)
  atoms <- .awpt_atom_metadata(reduction, basis_spec$scales, rep(NA_real_, length(basis_spec$scales)))
  if (nrow(atoms) != n_coeff) {
    atoms <- data.frame(
      col_id = seq_len(n_coeff),
      cluster_id = NA_integer_,
      scale = NA_real_,
      scale_index = NA_integer_,
      roughness_weight = scale_weights,
      stringsAsFactors = FALSE
    )
  } else {
    atoms$roughness_weight <- scale_weights
  }
  list(
    roughness = Matrix::Matrix(Q, sparse = FALSE),
    atoms = atoms,
    source = source,
    field_operator = Matrix::Matrix(L_field, sparse = FALSE)
  )
}

.awpt_roughness_from_coefficient <- function(loadings, reduction, basis_spec, coefficient_roughness) {
  n_coeff <- ncol(loadings)
  Q <- .awpt_as_square_matrix(coefficient_roughness, n_coeff, context = "coefficient_roughness")
  scale_weights <- diag(Q)
  atoms <- .awpt_atom_metadata(reduction, basis_spec$scales, rep(NA_real_, length(basis_spec$scales)))
  if (nrow(atoms) != n_coeff) {
    atoms <- data.frame(
      col_id = seq_len(n_coeff),
      cluster_id = NA_integer_,
      scale = NA_real_,
      scale_index = NA_integer_,
      roughness_weight = scale_weights,
      stringsAsFactors = FALSE
    )
  } else {
    atoms$roughness_weight <- scale_weights
  }
  list(
    roughness = Matrix::Matrix(Q, sparse = FALSE),
    atoms = atoms,
    source = "coefficient_operator",
    field_operator = NULL
  )
}

.awpt_require_neurosurf <- function(context = "surface AWPT") {
  if (!requireNamespace("neurosurf", quietly = TRUE)) {
    stop(context, " requires the 'neurosurf' package.", call. = FALSE)
  }
}

.awpt_normalize_surface_support <- function(support, geometry, context = "surface support") {
  if (is.null(support)) {
    stop(context, " is required.", call. = FALSE)
  }
  if (is.logical(support)) {
    support <- which(as.logical(support))
  }
  support <- as.integer(support)
  if (length(support) == 0L || anyNA(support) || any(support < 1L)) {
    stop(context, " must contain positive vertex indices.", call. = FALSE)
  }
  n_nodes <- length(neurosurf::nodes(geometry))
  if (any(support > n_nodes)) {
    stop(context, " contains indices beyond the geometry node count ", n_nodes, ".",
         call. = FALSE)
  }
  support
}

.surface_awpt_atom_metadata <- function(centers, support, scales) {
  rows <- vector("list", length(scales))
  col_offset <- 0L
  for (scale_idx in seq_along(scales)) {
    n_loc <- length(centers)
    rows[[scale_idx]] <- data.frame(
      col_id = seq.int(col_offset + 1L, col_offset + n_loc),
      center_id = support[centers],
      center_index = centers,
      scale = scales[[scale_idx]],
      scale_index = scale_idx,
      stringsAsFactors = FALSE
    )
    col_offset <- col_offset + n_loc
  }
  do.call(rbind, rows)
}

.surface_awpt_filters <- function(eigenvalues, scales) {
  filters <- vector("list", length(scales))
  prev <- rep(1, length(eigenvalues))
  for (idx in seq_along(scales)) {
    cur <- exp(-scales[[idx]] * eigenvalues)
    filters[[idx]] <- if (idx < length(scales)) {
      prev - cur
    } else {
      cur
    }
    prev <- cur
  }
  filters
}

.surface_awpt_loadings <- function(geometry, support, basis_spec,
                                   centers = NULL, threshold = 1e-6) {
  .awpt_require_neurosurf("awpt_surface_basis_template")
  support <- .awpt_normalize_surface_support(
    support,
    geometry = geometry,
    context = "awpt_surface_basis_template support"
  )
  L_full <- as.matrix(neurosurf::laplacian(geometry))
  L_support <- L_full[support, support, drop = FALSE]
  eig <- eigen(0.5 * (L_support + t(L_support)), symmetric = TRUE)
  vals <- pmax(Re(eig$values), 0)
  vecs <- Re(eig$vectors)

  if (is.null(centers)) {
    centers <- seq_along(support)
  } else {
    if (is.logical(centers)) {
      centers <- which(centers)
    }
    centers <- as.integer(centers)
    if (anyNA(centers) || any(centers < 1L) || any(centers > length(support))) {
      stop("centers must index the support vertices.", call. = FALSE)
    }
  }

  filters <- .surface_awpt_filters(vals, basis_spec$scales)
  atoms <- lapply(seq_along(filters), function(idx) {
    F <- vecs %*% (filters[[idx]] * t(vecs))
    atom_block <- F[, centers, drop = FALSE]
    atom_block[abs(atom_block) < threshold] <- 0
    atom_block
  })
  list(
    loadings = Matrix::Matrix(do.call(cbind, atoms), sparse = FALSE),
    field_operator = Matrix::Matrix(L_support, sparse = FALSE),
    atoms = .surface_awpt_atom_metadata(centers, support, basis_spec$scales)
  )
}

#' Average subject conductance matrices on a shared template graph
#'
#' @param conductances List of symmetric conductance matrices on the same template graph.
#' @param method Averaging rule.
#' @param shrinkage Optional shrinkage toward the isotropic identity.
#' @param enforce_psd Logical; if \code{TRUE}, project the result to the PSD cone.
#' @param tol Numerical tolerance for SPD operations.
#' @return A symmetric averaged conductance matrix.
#' @export
awpt_mean_conductance <- function(conductances,
                                  method = c("log_euclidean", "arithmetic"),
                                  shrinkage = 0,
                                  enforce_psd = TRUE,
                                  tol = 1e-8) {
  method <- match.arg(method)
  if (!is.list(conductances) || length(conductances) == 0L) {
    stop("conductances must be a non-empty list of matrices.", call. = FALSE)
  }
  mats <- lapply(seq_along(conductances), function(idx) {
    .awpt_enforce_symmetric(conductances[[idx]], tol = tol,
                            context = paste0("conductance[[", idx, "]]"))
  })
  dims <- lapply(mats, dim)
  if (!all(vapply(dims, identical, logical(1), dims[[1L]]))) {
    stop("All conductance matrices must have identical dimensions.", call. = FALSE)
  }
  n <- nrow(mats[[1L]])
  if (!is.numeric(shrinkage) || length(shrinkage) != 1L || !is.finite(shrinkage) ||
      shrinkage < 0 || shrinkage > 1) {
    stop("shrinkage must be a single number in [0, 1].", call. = FALSE)
  }
  mean_mat <- switch(
    method,
    arithmetic = Reduce(`+`, mats) / length(mats),
    log_euclidean = .awpt_log_euclidean_mean(mats, tol = tol)
  )
  if (shrinkage > 0) {
    baseline <- diag(mean(diag(mean_mat)), nrow = n)
    mean_mat <- (1 - shrinkage) * mean_mat + shrinkage * baseline
  }
  mean_mat <- .awpt_enforce_symmetric(mean_mat, tol = tol, context = "mean conductance")
  if (isTRUE(enforce_psd)) {
    mean_mat <- .awpt_psd_project(mean_mat, tol = tol)
  }
  Matrix::Matrix(mean_mat, sparse = FALSE)
}

#' Build an AWPT field operator from a conductance matrix
#'
#' @param conductance Symmetric conductance matrix on the template graph.
#' @param normalize Laplacian normalization convention.
#' @param tol Numerical tolerance.
#' @return A field-space roughness operator.
#' @export
awpt_operator_from_conductance <- function(conductance,
                                           normalize = c("none", "sym", "rw"),
                                           tol = 1e-8) {
  normalize <- match.arg(normalize)
  W <- .awpt_enforce_symmetric(conductance, tol = tol, context = "conductance")
  d <- rowSums(W)
  L <- diag(d) - W
  if (normalize == "none") {
    return(Matrix::Matrix(L, sparse = FALSE))
  }
  d_safe <- pmax(d, tol)
  if (normalize == "rw") {
    return(Matrix::Matrix(diag(1 / d_safe) %*% L, sparse = FALSE))
  }
  inv_sqrt <- diag(1 / sqrt(d_safe))
  Matrix::Matrix(inv_sqrt %*% L %*% inv_sqrt, sparse = FALSE)
}

#' Build an AWPT field operator from subject conductance summaries
#'
#' @param conductances List of subject conductance matrices on a shared template graph.
#' @param mean_method Averaging rule for the conductance mean.
#' @param normalize Laplacian normalization convention.
#' @param shrinkage Optional shrinkage toward an isotropic baseline before Laplacian construction.
#' @param enforce_psd Logical; if \code{TRUE}, project the mean conductance to the PSD cone.
#' @param tol Numerical tolerance.
#' @return A list with \code{conductance_mean} and \code{operator}.
#' @export
awpt_operator_from_subject_conductances <- function(conductances,
                                                    mean_method = c("log_euclidean", "arithmetic"),
                                                    normalize = c("none", "sym", "rw"),
                                                    shrinkage = 0,
                                                    enforce_psd = TRUE,
                                                    tol = 1e-8) {
  mean_method <- match.arg(mean_method)
  normalize <- match.arg(normalize)
  conductance_mean <- awpt_mean_conductance(
    conductances = conductances,
    method = mean_method,
    shrinkage = shrinkage,
    enforce_psd = enforce_psd,
    tol = tol
  )
  list(
    conductance_mean = conductance_mean,
    operator = awpt_operator_from_conductance(conductance_mean, normalize = normalize, tol = tol)
  )
}

#' Build an AWPT basis template
#'
#' @param parcellation A \code{ClusterReduction} or \code{ClusteredNeuroVol}.
#' @param basis_spec An AWPT basis specification created by \code{basis_awpt_wavelet()}.
#' @param loadings Optional explicit template loadings matrix. When supplied,
#'   fmrilatent skips wave-packet lifting and uses these loadings directly as
#'   the decoder basis \eqn{B}.
#' @param anatomical_operator Optional field-space roughness operator on the
#'   template domain. When supplied, the coefficient roughness is computed as
#'   \eqn{Q = B^T L B}. In the current v1 implementation this affects the
#'   roughness penalty, not the basis construction itself.
#' @param conductance Optional symmetric field-space conductance matrix. When
#'   supplied, fmrilatent builds the corresponding graph Laplacian and then
#'   forms \eqn{Q = B^T L B}. As with \code{anatomical_operator}, this shapes
#'   the v1 roughness model rather than directly adapting the lifted basis.
#' @param coefficient_roughness Optional coefficient-space roughness matrix.
#'   This bypasses field-space construction and is stored directly as \eqn{Q}.
#' @param center Logical; if \code{TRUE}, center data before projection.
#' @param ridge Small positive ridge added to the Gram diagonal if needed.
#' @param label Optional label stored in metadata.
#' @param ... Additional arguments passed to the underlying lift path.
#' @return An \code{AWPTBasisTemplate}.
#' @export
awpt_basis_template <- function(parcellation,
                                basis_spec = basis_awpt_wavelet(),
                                loadings = NULL,
                                anatomical_operator = NULL,
                                conductance = NULL,
                                coefficient_roughness = NULL,
                                center = FALSE,
                                ridge = 1e-8,
                                label = "awpt_wavelet",
                                ...) {
  supplied_sources <- sum(!vapply(
    list(anatomical_operator, conductance, coefficient_roughness),
    is.null,
    logical(1)
  ))
  if (supplied_sources > 1L) {
    stop("Supply at most one of anatomical_operator, conductance, or coefficient_roughness.",
         call. = FALSE)
  }
  reduction <- .coerce_awpt_reduction(parcellation)
  if (is.null(loadings)) {
    lift_spec <- basis_heat_wavelet(
      scales = basis_spec$scales,
      order = basis_spec$order,
      threshold = basis_spec$threshold
    )
    loadings <- lift(
      reduction,
      lift_spec,
      k_neighbors = basis_spec$k_neighbors,
      ...
    )
  } else {
    loadings <- Matrix::Matrix(as.matrix(loadings), sparse = FALSE)
    n_vox <- sum(as.array(reduction@mask))
    if (nrow(loadings) != n_vox) {
      stop("Explicit loadings must have ", n_vox, " rows to match the reduction mask.",
           call. = FALSE)
    }
  }

  G <- Matrix::crossprod(loadings)
  gram_factor <- tryCatch(
    Matrix::Cholesky(G, perm = TRUE),
    error = function(e) {
      G_ridge <- G + ridge * Matrix::Diagonal(n = ncol(G))
      Matrix::Cholesky(G_ridge, perm = TRUE)
    }
  )

  roughness_info <- if (!is.null(coefficient_roughness)) {
    .awpt_roughness_from_coefficient(loadings, reduction, basis_spec, coefficient_roughness)
  } else if (!is.null(anatomical_operator)) {
    .awpt_roughness_from_field_operator(
      loadings = loadings,
      reduction = reduction,
      basis_spec = basis_spec,
      field_operator = anatomical_operator,
      source = "anatomical_operator"
    )
  } else if (!is.null(conductance)) {
    .awpt_roughness_from_field_operator(
      loadings = loadings,
      reduction = reduction,
      basis_spec = basis_spec,
      field_operator = .awpt_graph_laplacian(conductance, nrow(loadings)),
      source = "conductance_laplacian"
    )
  } else {
    .awpt_roughness_from_scale(reduction, basis_spec, ncol(loadings))
  }

  structure(
    list(
      loadings = loadings,
      gram_factor = gram_factor,
      roughness = roughness_info$roughness,
      field_operator = roughness_info$field_operator,
      reduction = reduction,
      basis_spec = basis_spec,
      atoms = roughness_info$atoms,
      center = isTRUE(center),
      meta = list(
        family = "awpt_wavelet",
        method = "AWPT",
        label = label,
        k = ncol(loadings),
        ridge = ridge,
        scales = basis_spec$scales,
        order = basis_spec$order,
        threshold = basis_spec$threshold,
        k_neighbors = basis_spec$k_neighbors,
        penalty_rule = basis_spec$penalty_rule,
        roughness_source = roughness_info$source
      )
    ),
    class = "AWPTBasisTemplate"
  )
}

#' Build a surface AWPT basis template
#'
#' @param geometry A \code{neurosurf::SurfaceGeometry} or \code{neurosurf::SurfaceSet}.
#' @param basis_spec An AWPT basis specification created by \code{basis_awpt_wavelet()}.
#' @param support Surface support as vertex indices or a logical vector over all vertices.
#' @param loadings Optional explicit surface decoder loadings.
#' @param centers Optional support-local center indices used for automatic wave-packet construction.
#' @param anatomical_operator Optional field-space roughness operator on the supported surface domain.
#' @param conductance Optional conductance matrix on the supported surface graph.
#' @param coefficient_roughness Optional coefficient-space roughness matrix.
#' @param measure Optional support-aligned weighting or mass information.
#' @param center Logical; if \code{TRUE}, center data before projection.
#' @param ridge Small positive ridge added to the Gram diagonal if needed.
#' @param label Optional label stored in metadata.
#' @return A \code{SurfaceAWPTBasisTemplate}.
#' @export
awpt_surface_basis_template <- function(geometry,
                                        basis_spec = basis_awpt_wavelet(),
                                        support = NULL,
                                        loadings = NULL,
                                        centers = NULL,
                                        anatomical_operator = NULL,
                                        conductance = NULL,
                                        coefficient_roughness = NULL,
                                        measure = NULL,
                                        center = FALSE,
                                        ridge = 1e-8,
                                        label = "surface_awpt_wavelet") {
  .awpt_require_neurosurf("awpt_surface_basis_template")
  if (!(methods::is(geometry, "SurfaceGeometry") || methods::is(geometry, "SurfaceSet"))) {
    stop("geometry must be a neurosurf::SurfaceGeometry or neurosurf::SurfaceSet.",
         call. = FALSE)
  }
  supplied_sources <- sum(!vapply(
    list(anatomical_operator, conductance, coefficient_roughness),
    is.null,
    logical(1)
  ))
  if (supplied_sources > 1L) {
    stop("Supply at most one of anatomical_operator, conductance, or coefficient_roughness.",
         call. = FALSE)
  }

  n_nodes <- length(neurosurf::nodes(geometry))
  if (is.null(support)) {
    support <- seq_len(n_nodes)
  }
  support <- .awpt_normalize_surface_support(
    support,
    geometry = geometry,
    context = "awpt_surface_basis_template support"
  )

  if (is.null(loadings)) {
    auto <- .surface_awpt_loadings(
      geometry = geometry,
      support = support,
      basis_spec = basis_spec,
      centers = centers,
      threshold = basis_spec$threshold
    )
    loadings <- auto$loadings
    field_operator_default <- auto$field_operator
    atoms <- auto$atoms
  } else {
    loadings <- Matrix::Matrix(as.matrix(loadings), sparse = FALSE)
    if (nrow(loadings) != length(support)) {
      stop("Explicit loadings must have ", length(support),
           " rows to match the surface support cardinality.", call. = FALSE)
    }
    field_operator_default <- NULL
    atoms <- data.frame(
      col_id = seq_len(ncol(loadings)),
      center_id = NA_integer_,
      center_index = NA_integer_,
      scale = NA_real_,
      scale_index = NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  field_operator_use <- if (!is.null(coefficient_roughness)) {
    NULL
  } else if (!is.null(anatomical_operator)) {
    .awpt_as_square_matrix(anatomical_operator, nrow(loadings), "anatomical_operator")
  } else if (!is.null(conductance)) {
    .awpt_graph_laplacian(conductance, nrow(loadings))
  } else {
    field_operator_default
  }

  roughness <- if (!is.null(coefficient_roughness)) {
    .awpt_as_square_matrix(coefficient_roughness, ncol(loadings), "coefficient_roughness")
  } else if (!is.null(field_operator_use)) {
    crossprod(as.matrix(loadings), as.matrix(field_operator_use) %*% as.matrix(loadings))
  } else {
    {
      sw <- .awpt_scale_weights(
        basis_spec$scales,
        basis_spec$penalty_rule,
        custom_weights = basis_spec$custom_weights
      )
      expected_len <- length(sw) * length(if (!is.null(centers)) centers else seq_along(support))
      if (ncol(loadings) != expected_len && ncol(loadings) %% length(sw) != 0L) {
        stop("Scale-weight vector length (", length(sw),
             ") does not evenly divide the number of loadings columns (",
             ncol(loadings), ").", call. = FALSE)
      }
      Diagonal(x = rep_len(sw, ncol(loadings)))
    }
  }

  G <- Matrix::crossprod(loadings)
  gram_factor <- tryCatch(
    withCallingHandlers(
      Matrix::Cholesky(G, perm = TRUE),
      warning = function(w) {
        if (grepl("rank deficient|not positive definite", conditionMessage(w), ignore.case = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    ),
    error = function(e) {
      Matrix::Cholesky(G + ridge * Matrix::Diagonal(n = ncol(G)), perm = TRUE)
    }
  )

  structure(
    list(
      geometry = geometry,
      support = support,
      loadings = loadings,
      roughness = Matrix::Matrix(as.matrix(roughness), sparse = FALSE),
      field_operator = if (is.null(field_operator_use)) NULL else Matrix::Matrix(as.matrix(field_operator_use), sparse = FALSE),
      measure = measure,
      gram_factor = gram_factor,
      basis_spec = basis_spec,
      atoms = atoms,
      center = isTRUE(center),
      meta = list(
        family = "surface_awpt_wavelet",
        method = "AWPT",
        label = label,
        k = ncol(loadings),
        ridge = ridge,
        scales = basis_spec$scales,
        order = basis_spec$order,
        threshold = basis_spec$threshold,
        penalty_rule = basis_spec$penalty_rule,
        roughness_source = if (!is.null(coefficient_roughness)) {
          "coefficient_operator"
        } else if (!is.null(anatomical_operator)) {
          "anatomical_operator"
        } else if (!is.null(conductance)) {
          "conductance_laplacian"
        } else if (!is.null(field_operator_default)) {
          "surface_laplacian"
        } else {
          "scale_surrogate"
        }
      )
    ),
    class = c("SurfaceAWPTBasisTemplate", "SurfaceBasisTemplate")
  )
}

is_awpt_template <- function(x) {
  inherits(x, "AWPTBasisTemplate") || inherits(x, "SurfaceAWPTBasisTemplate")
}

#' @export
#' @rdname template_loadings
setMethod("template_loadings", "SurfaceAWPTBasisTemplate", function(x, ...) x$loadings)

#' @export
#' @rdname template_mask
setMethod("template_mask", "SurfaceAWPTBasisTemplate",
          function(x, ...) {
            stop("SurfaceAWPTBasisTemplate has no volumetric mask; use template_support().",
                 call. = FALSE)
          })

#' @export
#' @rdname template_meta
setMethod("template_meta", "SurfaceAWPTBasisTemplate", function(x, ...) x$meta %||% list())

#' @export
#' @rdname basis_decoder
setMethod("basis_decoder", "SurfaceAWPTBasisTemplate",
          function(template, ...) {
            decoder_map <- .template_custom_field(template, "decoder_map") %||%
              template_meta(template)$decoder_map %||% NULL
            if (!is.null(decoder_map)) {
              return(.normalize_linear_map(decoder_map, context = "basis decoder"))
            }
            payload <- .template_coordinate_payload(
              raw_loadings = template_loadings(template),
              measure = template_measure(template),
              analysis_transform = .template_custom_field(template, "analysis_transform") %||%
                template_meta(template)$analysis_transform %||% NULL,
              default_measure = "null"
            )
            B <- payload$analysis_loadings
            .linear_map_from_matrix(
              B,
              source_domain_id = paste0("latent:", template_meta(template)$family %||% "surface_awpt"),
              target_domain_id = digest::digest(list(
                support = template_support(template),
                meta = template_meta(template)
              )),
              provenance = list(
                basis_asset_class = "SurfaceAWPTBasisTemplate",
                basis_family = template_meta(template)$family %||% "surface_awpt_wavelet",
                basis_id = digest::digest(B),
                target_support = template_support(template)
              )
            )
          })

#' @export
#' @rdname template_rank
setMethod("template_rank", "SurfaceAWPTBasisTemplate",
          function(template, ...) ncol(template_loadings(template)))

#' @export
#' @rdname template_domain
setMethod("template_domain", "SurfaceAWPTBasisTemplate",
          function(template, ...) template$geometry)

#' @export
#' @rdname template_support
setMethod("template_support", "SurfaceAWPTBasisTemplate",
          function(template, ...) template$support)

#' @export
#' @rdname template_measure
setMethod("template_measure", "SurfaceAWPTBasisTemplate",
          function(template, ...) template$measure %||% NULL)

#' @export
#' @rdname template_roughness
setMethod("template_roughness", "SurfaceAWPTBasisTemplate",
          function(template, coordinates = c("analysis", "raw"), ...) {
            coordinates <- match.arg(coordinates)
            .transform_quadratic_form(
              template$roughness %||% NULL,
              .template_coordinate_payload(
                raw_loadings = template_loadings(template),
                measure = template_measure(template),
                analysis_transform = .template_custom_field(template, "analysis_transform") %||%
                  template_meta(template)$analysis_transform %||% NULL,
                default_measure = "null"
              )$analysis_transform,
              coordinates = coordinates
            )
          })

#' @export
#' @rdname template_project
setMethod("template_project", signature(x = "SurfaceAWPTBasisTemplate", data = "ANY"),
          function(x, data, ...) {
            X <- as.matrix(data)
            if (ncol(X) != length(template_support(x))) {
              stop("data must have ", length(template_support(x)),
                   " columns to match the surface support cardinality.", call. = FALSE)
            }
            .template_projection_payload(
              data = X,
              raw_loadings = template_loadings(x),
              measure = template_measure(x),
              center = isTRUE(x$center),
              analysis_transform = .template_custom_field(x, "analysis_transform") %||%
                template_meta(x)$analysis_transform %||% NULL,
              default_measure = "null"
            )
          })

#' @export
#' @rdname save_template
setMethod("save_template", signature(template = "SurfaceAWPTBasisTemplate"),
          function(template, file, compress = TRUE, ...) {
            saveRDS(template, file = file, compress = compress)
            invisible(file)
          })

#' @export
print.SurfaceAWPTBasisTemplate <- function(x, ...) {
  cat("SurfaceAWPTBasisTemplate\n")
  cat("  Method: AWPT (surface)\n")
  cat("  Atoms:", ncol(x$loadings), "\n")
  cat("  Support vertices:", length(x$support), "\n")
  cat("  Scales:", paste(x$meta$scales, collapse = ", "), "\n")
  cat("  Roughness:", x$meta$roughness_source %||% "unknown", "\n")
  cat("  Center at encode:", x$center, "\n")
  invisible(x)
}

#' @export
print.AWPTBasisTemplate <- function(x, ...) {
  cat("AWPTBasisTemplate\n")
  cat("  Method: AWPT\n")
  cat("  Atoms:", ncol(x$loadings), "\n")
  cat("  Voxels:", nrow(x$loadings), "\n")
  cat("  Scales:", paste(x$meta$scales, collapse = ", "), "\n")
  cat("  Roughness:", x$meta$roughness_source %||% "unknown", "\n")
  cat("  Center at encode:", x$center, "\n")
  invisible(x)
}

#' @export
#' @rdname template_loadings
setMethod("template_loadings", "AWPTBasisTemplate", function(x, ...) x$loadings)

#' @export
#' @rdname template_mask
setMethod("template_mask", "AWPTBasisTemplate", function(x, ...) x$reduction@mask)

#' @export
#' @rdname template_meta
setMethod("template_meta", "AWPTBasisTemplate", function(x, ...) x$meta %||% list())

#' @export
#' @rdname basis_decoder
setMethod("basis_decoder", "AWPTBasisTemplate",
          function(template, ...) {
            decoder_map <- .template_custom_field(template, "decoder_map") %||%
              template_meta(template)$decoder_map %||% NULL
            if (!is.null(decoder_map)) {
              return(.normalize_linear_map(decoder_map, context = "basis decoder"))
            }
            payload <- .template_coordinate_payload(
              raw_loadings = template_loadings(template),
              measure = template_measure(template),
              analysis_transform = .template_custom_field(template, "analysis_transform") %||%
                template_meta(template)$analysis_transform %||% NULL,
              default_measure = "unit"
            )
            B <- payload$analysis_loadings
            .linear_map_from_matrix(
              B,
              source_domain_id = "latent:awpt",
              target_domain_id = digest::digest(list(
                mask = as.array(template_mask(template)),
                meta = template_meta(template)
              )),
              provenance = list(
                basis_asset_class = "AWPTBasisTemplate",
                basis_family = "awpt_wavelet",
                basis_id = digest::digest(B)
              )
            )
          })

#' @export
#' @rdname template_rank
setMethod("template_rank", "AWPTBasisTemplate",
          function(template, ...) ncol(template_loadings(template)))

#' @export
#' @rdname template_domain
setMethod("template_domain", "AWPTBasisTemplate",
          function(template, ...) neuroim2::space(template_mask(template)))

#' @export
#' @rdname template_support
setMethod("template_support", "AWPTBasisTemplate",
          function(template, ...) template_mask(template))

#' @export
#' @rdname template_measure
setMethod("template_measure", "AWPTBasisTemplate",
          function(template, ...) {
            template_meta(template)$measure %||% rep(1, nrow(template_loadings(template)))
          })

#' @export
#' @rdname template_roughness
setMethod("template_roughness", "AWPTBasisTemplate",
          function(template, coordinates = c("analysis", "raw"), ...) {
            coordinates <- match.arg(coordinates)
            .transform_quadratic_form(
              template$roughness,
              .template_coordinate_payload(
                raw_loadings = template_loadings(template),
                measure = template_measure(template),
                analysis_transform = .template_custom_field(template, "analysis_transform") %||%
                  template_meta(template)$analysis_transform %||% NULL,
                default_measure = "unit"
              )$analysis_transform,
              coordinates = coordinates
            )
          })

#' @export
#' @rdname template_project
setMethod("template_project", signature(x = "AWPTBasisTemplate", data = "ANY"),
          function(x, data, ...) {
            X <- as.matrix(data)
            .template_projection_payload(
              data = X,
              raw_loadings = template_loadings(x),
              measure = template_measure(x),
              center = isTRUE(x$center),
              analysis_transform = .template_custom_field(x, "analysis_transform") %||%
                template_meta(x)$analysis_transform %||% NULL,
              default_measure = "unit"
            )
          })

#' @export
#' @rdname save_template
setMethod("save_template", signature(template = "AWPTBasisTemplate"),
          function(template, file, compress = "xz", ...) {
            saveRDS(template, file = file, compress = compress)
            invisible(normalizePath(file, winslash = "/", mustWork = FALSE))
          })
