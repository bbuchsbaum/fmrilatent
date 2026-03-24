# Shared parcel basis templates and projection-only encoding
#
# Builds a reusable spatial dictionary (loadings) from any lift()-compatible
# spec (default: Laplacian eigenvectors via basis_slepian), caches the Gram
# factorization, and provides a project-only encode path for applying the
# template to new subjects.

#' @include reduction.R encode.R
#' @importFrom Matrix crossprod Cholesky solve Matrix
#' @importFrom methods is
NULL

# --- ParcelBasisTemplate -----------------------------------------------------

#' Build a shared parcel basis template
#'
#' Computes a spatial dictionary within each parcel using \code{\link{lift}()}
#' and caches the Gram factorization for efficient projection. The resulting
#' template can be reused across subjects via \code{\link{spec_space_parcel}()}.
#'
#' @param parcellation A \code{\link{ClusterReduction}} or a
#'   \code{ClusteredNeuroVol} (coerced via \code{\link{as_cluster_reduction}()}).
#' @param basis_spec A basis specification for \code{\link{lift}()}.
#'   Default is \code{basis_slepian(k = 5)} which computes the k smallest
#'   Laplacian eigenvectors of the voxel adjacency graph within each parcel.
#'   Use \code{basis_pca(k)} with \code{data} for data-driven bases.
#' @param data Optional numeric matrix (time x voxels, mask order). Required
#'   for data-driven specs like \code{basis_pca()}.
#' @param center Logical; if \code{TRUE} (default), center voxel time series
#'   both when building PCA templates and when projecting new subjects. The
#'   encoded object stores subject-specific voxel means in
#'   \code{LatentNeuroVec@offset} so reconstruction preserves voxel means.
#'   For geometric bases (Laplacian, Slepian) this controls whether subjects
#'   are centered prior to projection.
#' @param ridge Small positive scalar added to the Gram diagonal if Cholesky
#'   fails (default 1e-8).
#' @param ... Additional arguments passed to \code{lift()} (e.g.,
#'   \code{k_neighbors} for graph-based specs). For PCA templates, projection
#'   must remain replayable at encode time, so preprocessing arguments such as
#'   \code{offset} and \code{scale} are not supported here.
#'
#' @return A \code{ParcelBasisTemplate} object (S3) with components:
#'   \describe{
#'     \item{loadings}{Sparse \code{Matrix} (voxels x atoms), block-diagonal by parcel.}
#'     \item{gram_factor}{Cached Cholesky factorization of \eqn{L^T L}.}
#'     \item{reduction}{\code{ClusterReduction} used.}
#'     \item{basis_spec}{The spec that produced the loadings.}
#'     \item{center}{Whether centering is applied when building/projecting.}
#'     \item{meta}{List with \code{family}, \code{k}, \code{ridge},
#'       \code{label_map}, and \code{cluster_map}.}
#'   }
#'
#' @details
#' The default \code{basis_slepian(k)} computes the k smallest eigenvectors
#' of the graph Laplacian built from voxel coordinates within each parcel.
#' These smooth spatial functions form a data-independent dictionary suitable
#' for projecting any subject's data.
#'
#' For data-driven bases, pass \code{basis_pca(k)} with \code{data =} a
#' training matrix (e.g., group-average data). The resulting PCA loadings
#' are then fixed and reused for each subject.
#'
#' @examples
#' \dontrun{
#' # Geometric (Laplacian) basis -- no training data needed
#' atlas <- load_atlas("schaefer_400")
#' tmpl <- parcel_basis_template(atlas, basis_slepian(k = 8))
#' lvec <- encode(bold, spec_space_parcel(tmpl))
#'
#' # Data-driven shared PCA basis
#' tmpl_pca <- parcel_basis_template(atlas, basis_pca(k = 5), data = group_bold)
#' lvec <- encode(subj_bold, spec_space_parcel(tmpl_pca))
#' }
#'
#' @seealso \code{\link{spec_space_parcel}}, \code{\link{as_cluster_reduction}},
#'   \code{\link{lift}}, \code{\link{basis_slepian}}, \code{\link{basis_pca}}
#' @export
parcel_basis_template <- function(parcellation,
                                  basis_spec = basis_slepian(k = 5),
                                  data = NULL,
                                  center = TRUE,
                                  ridge = 1e-8,
                                  ...) {
  extra_args <- list(...)

  # Coerce parcellation
  reduction <- if (inherits(parcellation, "ClusterReduction")) {
    parcellation
  } else if (inherits(parcellation, "ClusteredNeuroVol")) {
    as_cluster_reduction(parcellation)
  } else {
    stop("parcellation must be a ClusterReduction or ClusteredNeuroVol.", call. = FALSE)
  }

  if (inherits(basis_spec, "spec_pca")) {
    if ("center" %in% names(extra_args)) {
      stop("Use parcel_basis_template(center = ...) rather than passing center in ... for PCA templates.",
           call. = FALSE)
    }
    if ("offset" %in% names(extra_args)) {
      stop("PCA parcel templates do not support a custom offset; encode-time projection must remain replayable.",
           call. = FALSE)
    }
    if ("scale" %in% names(extra_args) && isTRUE(extra_args$scale)) {
      stop("PCA parcel templates do not support scale = TRUE; projection only replays centering.",
           call. = FALSE)
    }
    extra_args$center <- isTRUE(center)
  }

  # Build per-parcel spatial dictionary via lift()
  loadings <- do.call(
    lift,
    c(
      list(reduction = reduction, basis_spec = basis_spec, data = data),
      extra_args
    )
  )

  if (ncol(loadings) == 0L) {
    stop("lift() returned zero columns; check basis_spec and parcellation.", call. = FALSE)
  }

  # Cache Gram factorization: G = L^T L
  G <- Matrix::crossprod(loadings)
  gram_factor <- tryCatch(
    Matrix::Cholesky(G, perm = TRUE),
    error = function(e) {
      # Add ridge and retry
      G_ridge <- G
      Matrix::diag(G_ridge) <- Matrix::diag(G_ridge) + ridge
      tryCatch(
        Matrix::Cholesky(G_ridge, perm = TRUE),
        error = function(e2) {
          stop("Gram matrix factorization failed even with ridge = ", ridge,
               ". The loadings may be rank-deficient.", call. = FALSE)
        }
      )
    }
  )

  structure(
    list(
      loadings = loadings,
      gram_factor = gram_factor,
      reduction = reduction,
      basis_spec = basis_spec,
      center = isTRUE(center),
      meta = list(
        family = class(basis_spec)[[1L]],
        k = ncol(loadings),
        ridge = ridge,
        label_map = reduction@info$label_map,
        cluster_map = reduction@info$cluster_map
      )
    ),
    class = "ParcelBasisTemplate"
  )
}

#' Print method for ParcelBasisTemplate
#' @param x A ParcelBasisTemplate object.
#' @param ... Ignored.
#' @export
print.ParcelBasisTemplate <- function(x, ...) {
  cat("ParcelBasisTemplate\n")
  cat("  Basis family:", x$meta$family, "\n")
  cat("  Atoms:", x$meta$k, "\n")
  cat("  Voxels:", nrow(x$loadings), "\n")
  cat("  Parcels:", length(x$reduction@cluster_ids), "\n")
  cat("  Center at encode:", x$center, "\n")
  if (!is.null(x$meta$label_map)) {
    cat("  Label map: available\n")
  }
  invisible(x)
}

#' @export
#' @rdname template_loadings
setMethod("template_loadings", "ParcelBasisTemplate", function(x, ...) x$loadings)

#' @export
#' @rdname template_mask
setMethod("template_mask", "ParcelBasisTemplate", function(x, ...) x$reduction@mask)

#' @export
#' @rdname template_meta
setMethod("template_meta", "ParcelBasisTemplate", function(x, ...) x$meta %||% list())

#' @export
#' @rdname template_project
setMethod("template_project", signature(x = "ParcelBasisTemplate", data = "ANY"),
          function(x, data, ...) {
            X <- as.matrix(data)
            L <- template_loadings(x)
            offset <- numeric(0)
            X_proj <- X
            if (isTRUE(x$center)) {
              offset <- colMeans(X)
              X_proj <- sweep(X, 2L, offset, "-")
            }
            proj <- X_proj %*% L
            coeff <- t(Matrix::solve(x$gram_factor, Matrix::t(proj)))
            list(
              coefficients = Matrix::Matrix(as.matrix(coeff), sparse = FALSE),
              offset = offset
            )
          })

#' @export
#' @rdname save_template
setMethod("save_template", signature(template = "ParcelBasisTemplate"),
          function(template, file, compress = "xz", ...) {
            saveRDS(template, file = file, compress = compress)
            invisible(normalizePath(file, winslash = "/", mustWork = FALSE))
          })

# --- Spec and encode_spec for shared parcel basis ----------------------------

#' Spatial parcel-basis spec (shared/template-based)
#'
#' Creates a spec for projecting data onto a pre-built
#' \code{\link{ParcelBasisTemplate}}. The loadings are fixed (shared across
#' subjects); only the temporal scores and per-subject offset vary.
#'
#' @param template A \code{\link{ParcelBasisTemplate}} object built by
#'   \code{\link{parcel_basis_template}()}.
#' @return A \code{spec_space_parcel} object for \code{\link{encode}()}.
#'
#' @examples
#' \dontrun{
#' tmpl <- parcel_basis_template(atlas, basis_slepian(k = 8))
#' lvec_s1 <- encode(bold_s1, spec_space_parcel(tmpl))
#' lvec_s2 <- encode(bold_s2, spec_space_parcel(tmpl))
#' # Same loadings; different basis (scores) and offset
#' }
#'
#' @seealso \code{\link{parcel_basis_template}}, \code{\link{encode}}
#' @export
spec_space_parcel <- function(template) {
  if (!inherits(template, "ParcelBasisTemplate")) {
    stop("template must be a ParcelBasisTemplate (from parcel_basis_template()).",
         call. = FALSE)
  }
  structure(list(template = template), class = "spec_space_parcel")
}

#' @exportS3Method
encode_spec.spec_space_parcel <- function(x, spec, mask, reduction, materialize, label, ...) {
  tmpl <- spec$template
  L <- template_loadings(tmpl)
  n_time <- nrow(x)
  n_vox <- ncol(x)
  template_mask <- template_mask(tmpl)
  template_mask_arr <- as.array(template_mask)
  supplied_mask_arr <- .mask_to_array(mask, "encode_spec.spec_space_parcel")

  if (n_vox != nrow(L)) {
    stop("Data has ", n_vox, " voxels but template loadings have ", nrow(L),
         " rows.", call. = FALSE)
  }
  if (!identical(supplied_mask_arr, template_mask_arr)) {
    stop("mask does not match the template mask. Shared parcel templates require identical voxel support and ordering.",
         call. = FALSE)
  }

  proj <- template_project(tmpl, x)
  basis <- proj$coefficients

  spc <- neuroim2::NeuroSpace(c(dim(template_mask_arr), n_time))
  meta_tmpl <- template_meta(tmpl)

  meta <- list(
    family = "parcel_basis",
    basis_family = meta_tmpl$family,
    k = meta_tmpl$k,
    center = tmpl$center,
    n_parcels = length(tmpl$reduction@cluster_ids),
    label_map = meta_tmpl$label_map,
    cluster_map = meta_tmpl$cluster_map
  )

  LatentNeuroVec(
    basis = basis,
    loadings = L,
    space = spc,
    mask = template_mask,
    offset = proj$offset,
    label = label,
    meta = meta
  )
}
