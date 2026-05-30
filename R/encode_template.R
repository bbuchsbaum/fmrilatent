# Template helpers shared by hierarchical, parcel, AWPT, and surface
# encoders. These functions resolve template measure semantics, build
# analysis transforms from raw metrics, project data into analysis
# coordinates, and rotate quadratic forms between raw / analysis bases.
#
# External callers: R/parcel_basis.R, R/awpt.R, plus encode_spec.* methods
# in R/encode.R itself. Two `.transport_*` helpers used here
# (.transport_identity_transform, .transport_raw_metric) live in
# R/transport_latent.R and resolve at call time.

#' @include encode.R
NULL

.template_custom_field <- function(template, field) {
  if (isS4(template) || !is.list(template)) {
    return(NULL)
  }
  template[[field]] %||% NULL
}

.template_measure_resolve <- function(measure, n_support,
                                      default = c("unit", "null"),
                                      context = "template_measure") {
  default <- match.arg(default)
  if (is.null(measure)) {
    if (default == "unit") {
      return(rep(1, n_support))
    }
    return(NULL)
  }

  if (is.atomic(measure) && is.null(dim(measure))) {
    measure <- as.numeric(measure)
    if (length(measure) == 1L) {
      return(rep(measure, n_support))
    }
    if (length(measure) != n_support) {
      stop(context, " vector must have length ", n_support, ".", call. = FALSE)
    }
    return(measure)
  }

  measure <- as.matrix(measure)
  if (!identical(dim(measure), c(n_support, n_support))) {
    stop(context, " matrix must have dimensions ",
         n_support, "x", n_support, ".", call. = FALSE)
  }
  measure
}

.template_weighted_crossprod <- function(loadings, measure = NULL) {
  loadings <- as.matrix(loadings)
  if (is.null(measure)) {
    return(crossprod(loadings))
  }
  if (is.atomic(measure) && is.null(dim(measure))) {
    return(crossprod(loadings, loadings * as.numeric(measure)))
  }
  crossprod(loadings, as.matrix(measure) %*% loadings)
}

.template_weighted_right_projection <- function(data, loadings, measure = NULL) {
  data <- as.matrix(data)
  loadings <- as.matrix(loadings)
  if (is.null(measure)) {
    return(data %*% loadings)
  }
  if (is.atomic(measure) && is.null(dim(measure))) {
    return((data * rep(as.numeric(measure), each = nrow(data))) %*% loadings)
  }
  data %*% as.matrix(measure) %*% loadings
}

.symmetric_matrix_factor <- function(mat, tol = 1e-10, context = "matrix factor") {
  mat <- 0.5 * (as.matrix(mat) + t(as.matrix(mat)))
  if (nrow(mat) != ncol(mat)) {
    stop(context, " requires a square matrix.", call. = FALSE)
  }

  if (isTRUE(all.equal(mat, diag(nrow(mat)), tolerance = tol))) {
    return(diag(nrow(mat)))
  }

  chol_try <- try(chol(mat), silent = TRUE)
  if (!inherits(chol_try, "try-error")) {
    return(chol_try)
  }

  eig <- eigen(mat, symmetric = TRUE)
  vals <- pmax(Re(eig$values), tol)
  eig$vectors %*% diag(sqrt(vals), nrow = length(vals))
}

.analysis_transform_from_metric <- function(raw_metric, tol = 1e-10) {
  raw_metric <- 0.5 * (as.matrix(raw_metric) + t(as.matrix(raw_metric)))
  k <- nrow(raw_metric)
  if (!identical(dim(raw_metric), c(k, k))) {
    stop("raw_metric must be square.", call. = FALSE)
  }
  if (isTRUE(all.equal(raw_metric, diag(k), tolerance = tol))) {
    return(.transport_identity_transform(k))
  }

  factor <- .symmetric_matrix_factor(raw_metric, tol = tol,
                                     context = "analysis transform")
  list(
    type = "metric_factor",
    dim = as.integer(k),
    raw_metric = raw_metric,
    matrix = factor,
    to_analysis = function(data) factor %*% as.matrix(data),
    to_raw = function(data) solve(factor, as.matrix(data))
  )
}

.analysis_loadings_from_transform <- function(raw_loadings, transform) {
  raw_loadings <- as.matrix(raw_loadings)
  mat <- transform$matrix %||% diag(ncol(raw_loadings))
  raw_loadings %*% solve(as.matrix(mat))
}

.template_coordinate_payload <- function(raw_loadings, measure = NULL,
                                         analysis_transform = NULL,
                                         default_measure = c("unit", "null"),
                                         tol = 1e-10) {
  raw_loadings <- as.matrix(raw_loadings)
  measure_use <- .template_measure_resolve(
    measure,
    n_support = nrow(raw_loadings),
    default = match.arg(default_measure),
    context = "template measure"
  )

  transform <- analysis_transform
  if (is.null(transform)) {
    raw_metric <- .template_weighted_crossprod(raw_loadings, measure_use)
    transform <- .analysis_transform_from_metric(raw_metric, tol = tol)
  } else {
    raw_metric <- .transport_raw_metric(transform, ncol(raw_loadings))
  }

  list(
    raw_loadings = raw_loadings,
    analysis_loadings = .analysis_loadings_from_transform(raw_loadings, transform),
    measure = measure_use,
    analysis_transform = transform,
    raw_metric = raw_metric
  )
}

.template_projection_payload <- function(data, raw_loadings, measure = NULL,
                                         center = FALSE, offset = NULL,
                                         analysis_transform = NULL,
                                         default_measure = c("unit", "null"),
                                         tol = 1e-10) {
  centered <- .encode_center(data, center = center, offset = offset,
                             context = ".template_projection_payload")
  offset_out <- centered$offset
  X_proj <- centered$x_centered

  payload <- .template_coordinate_payload(
    raw_loadings = raw_loadings,
    measure = measure,
    analysis_transform = analysis_transform,
    default_measure = default_measure,
    tol = tol
  )

  coeff_analysis <- .template_weighted_right_projection(
    X_proj,
    payload$analysis_loadings,
    payload$measure
  )
  coeff_raw <- t(payload$analysis_transform$to_raw(t(coeff_analysis)))

  list(
    coefficients = Matrix::Matrix(coeff_analysis, sparse = FALSE),
    raw_coefficients = Matrix::Matrix(coeff_raw, sparse = FALSE),
    offset = offset_out,
    analysis_transform = payload$analysis_transform,
    analysis_loadings = Matrix::Matrix(payload$analysis_loadings, sparse = FALSE),
    raw_metric = payload$raw_metric,
    measure = payload$measure
  )
}

.template_default_measure_mode <- function(template) {
  if (inherits(template, "SurfaceBasisTemplate") ||
      inherits(template, "SurfaceAWPTBasisTemplate")) {
    "null"
  } else {
    "unit"
  }
}

.template_asset_analysis_transform <- function(template) {
  transform <- .template_custom_field(template, "analysis_transform") %||%
    tryCatch(template_meta(template)$analysis_transform %||% NULL,
             error = function(e) NULL)
  if (!is.null(transform)) {
    return(transform)
  }

  raw_loadings <- tryCatch(template_loadings(template), error = function(e) NULL)
  if (is.null(raw_loadings)) {
    return(NULL)
  }

  measure <- tryCatch(template_measure(template), error = function(e) NULL)
  .template_coordinate_payload(
    raw_loadings = raw_loadings,
    measure = measure,
    default_measure = .template_default_measure_mode(template)
  )$analysis_transform
}

.transform_quadratic_form <- function(mat, transform, coordinates = c("analysis", "raw")) {
  coordinates <- match.arg(coordinates)
  if (is.null(mat)) {
    return(NULL)
  }
  mat <- as.matrix(mat)
  if (coordinates == "raw") {
    return(mat)
  }

  factor <- transform$matrix %||% diag(nrow(mat))
  solve(t(as.matrix(factor)), mat) %*% solve(as.matrix(factor))
}
