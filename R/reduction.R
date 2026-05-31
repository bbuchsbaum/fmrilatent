#' @include all_generic.R all_class.R
NULL

#' Graph reduction scaffolds (abstract)
#'
#' These classes describe how voxels are grouped or coarsened before a basis is
#' computed and lifted back to ambient space. Implementations of `lift()` for
#' specific combinations (e.g., supervoxel + Slepian, parcel PCA) live in
#' external packages or downstream code; fmrilatent ships only the contracts.
#'
#' @slot mask `LogicalNeuroVol` defining the ambient domain.
#' @slot info Optional list for implementation-specific metadata.
#' @export
setClass("GraphReduction",
  slots = c(
    mask = "LogicalNeuroVol",
    info = "list"
  )
)

#' Cluster-based reduction (e.g., supervoxels or atlas)
#'
#' @slot map Integer vector (voxel order) mapping each voxel to a cluster id.
#' @slot cluster_ids Unique cluster ids present in `map`.
#' @export
setClass("ClusterReduction",
  contains = "GraphReduction",
  slots = c(
    map = "integer",
    cluster_ids = "integer"
  )
)

#' Coarsened graph reduction (e.g., prolongation from coarse to fine)
#'
#' @slot P_matrix Sparse prolongation matrix (fine x coarse).
#' @slot coarse_adj Optional sparse adjacency on the coarse graph.
#' @export
setClass("CoarsenedReduction",
  contains = "GraphReduction",
  slots = c(
    P_matrix = "dgCMatrix",
    coarse_adj = "dgCMatrix"
  )
)

setValidity("CoarsenedReduction", function(object) {
  errors <- character()
  n_active <- sum(as.logical(as.array(object@mask)))
  if (nrow(object@P_matrix) != n_active) {
    errors <- c(
      errors,
      sprintf("P_matrix must have one row per active mask voxel (%d); got %d.",
              n_active, nrow(object@P_matrix))
    )
  }
  if (ncol(object@P_matrix) < 1L) {
    errors <- c(errors, "P_matrix must have at least one coarse column.")
  }
  n_coarse <- ncol(object@P_matrix)
  if (nrow(object@coarse_adj) != n_coarse || ncol(object@coarse_adj) != n_coarse) {
    errors <- c(
      errors,
      sprintf("coarse_adj must be square with dimensions matching ncol(P_matrix) (%d).",
              n_coarse)
    )
  }
  if (length(errors)) errors else TRUE
})

#' Create a coarsened graph reduction
#'
#' @param mask `LogicalNeuroVol` or array-like mask defining the fine domain.
#' @param P_matrix Fine-by-coarse sparse prolongation matrix.
#' @param coarse_adj Optional coarse-by-coarse sparse adjacency matrix.
#' @param info Optional metadata list.
#' @return A valid `CoarsenedReduction`.
#' @export
make_coarsened_reduction <- function(mask, P_matrix, coarse_adj = NULL, info = list()) {
  mask_arr <- .extract_mask_array(mask, "make_coarsened_reduction")
  mask_vol <- .mask_volume_from_array(mask, mask_arr, "make_coarsened_reduction")
  P <- Matrix::Matrix(P_matrix, sparse = TRUE)
  if (!inherits(P, "dgCMatrix")) {
    P <- methods::as(P, "dgCMatrix")
  }
  coarse <- if (is.null(coarse_adj)) {
    Matrix::sparseMatrix(
      i = integer(),
      j = integer(),
      x = numeric(),
      dims = c(ncol(P), ncol(P))
    )
  } else {
    Matrix::Matrix(coarse_adj, sparse = TRUE)
  }
  if (!inherits(coarse, "dgCMatrix")) {
    coarse <- methods::as(coarse, "dgCMatrix")
  }
  new("CoarsenedReduction",
      mask = mask_vol,
      info = info,
      P_matrix = P,
      coarse_adj = coarse)
}

#' Create a Slepian basis specification
#'
#' `basis_slepian()` creates a lightweight descriptor for graph/Slepian basis
#' construction during `lift()` or spatial encoding. It records the requested
#' component count and Slepian flavor; the actual basis is computed later by a
#' `lift()` method for the supplied reduction.
#'
#' @param k Positive integer number of Slepian components.
#' @param type Character scalar naming the Slepian basis flavor. The built-in
#'   spatial methods use `"laplacian"`.
#' @return A list with class `spec_slepian` containing `k` and `type`.
#' @export
basis_slepian <- function(k = 3, type = "laplacian") {
  k <- .validate_positive_count(k, "k")
  structure(list(k = k, type = type), class = "spec_slepian")
}

#' Create a PCA basis specification
#'
#' `basis_pca()` creates a lightweight descriptor for parcel- or
#' cluster-local PCA bases. The descriptor is consumed by
#' `lift(ClusterReduction, spec_pca, data = ...)`, where `data` supplies the
#' time-by-voxel matrix used to estimate the components.
#'
#' @param k Positive integer number of PCA components.
#' @param whiten Logical scalar. `basis_pca()` records this request, but
#'   `lift(ClusterReduction, spec_pca)` returns unwhitened loadings and emits a
#'   warning when `whiten = TRUE`; the higher-level `encode(spec_space_pca())`
#'   path applies whitening after projection.
#' @return A list with class `spec_pca` containing `k` and `whiten`.
#' @export
basis_pca <- function(k = 3, whiten = FALSE) {
  k <- .validate_positive_count(k, "k")
  whiten <- .validate_flag_scalar(whiten, "whiten")
  structure(list(k = k, whiten = whiten), class = "spec_pca")
}

#' Create a flat basis specification
#'
#' `basis_flat()` creates a descriptor for methods that should lift a reduction
#' without estimating local components, typically as a simple parcel/cluster
#' indicator basis.
#'
#' @return An empty list with class `spec_flat`.
#' @export
basis_flat <- function() {
  structure(list(), class = "spec_flat")
}

#' Lift reduced bases back to voxel space (abstract generic)
#'
#' @param reduction A `GraphReduction` subclass describing topology.
#' @param basis_spec A basis specification (e.g., `basis_slepian()`).
#' @param data Optional data for data-driven bases (e.g., PCA).
#' @param ... Additional arguments passed to methods (e.g., k_neighbors).
#' @return Typically a voxel x components `Matrix` (often sparse) or an
#'   implementation-defined object for implicit decoders.
#' @export
setGeneric("lift", function(reduction, basis_spec, data = NULL, ...) {
  standardGeneric("lift")
})

#' Default lift method (placeholder)
#'
#' This method exists to provide a clear error when no concrete lift is
#' registered. External packages should implement methods for specific
#' (reduction, basis_spec) signatures.
#'
#' @param reduction A `GraphReduction` subclass.
#' @param basis_spec A basis specification object.
#' @param data Optional data for data-driven bases.
#' @param ... Additional arguments (unused in default method).
#' @return This method does not return; it aborts with a classed unsupported-operation error.
#' @export
setMethod("lift", signature(reduction = "GraphReduction", basis_spec = "ANY"),
  function(reduction, basis_spec, data = NULL, ...) {
    .encoder_cli_abort(
      paste0("No lift() implementation for this reduction/basis combination. ",
             "Provide an external method (e.g., supervoxel Slepian or parcel PCA)."),
      class = "fmrilatent_error_unsupported_operation", call = rlang::caller_env())
  }
)

# --- ClusterReduction constructors -------------------------------------------

.abort_na_cluster_map <- function(map, context) {
  if (anyNA(map)) {
    .encoder_cli_abort(
      paste0(context, " must not contain NA cluster ids."),
      class = "fmrilatent_error_invalid_cluster_map",
      call = rlang::caller_env()
    )
  }
}

#' Create a ClusterReduction from a mask and voxel-to-cluster map
#'
#' @param mask A \code{LogicalNeuroVol} or logical 3D array defining the brain mask.
#' @param map Integer vector (mask order) mapping each voxel to a cluster id.
#' @return A \code{ClusterReduction} object.
#' @export
make_cluster_reduction <- function(mask, map) {
  mask_vol <- if (inherits(mask, "LogicalNeuroVol")) {
    mask
  } else {
    LogicalNeuroVol(mask, neuroim2::NeuroSpace(dim(mask)))
  }
  n_vox <- sum(as.array(mask_vol))
  if (length(map) != n_vox) {
    .encoder_cli_abort(
      paste0("map must have length ", n_vox, " (the number of voxels in the mask), got ",
             length(map), "."),
      class = "fmrilatent_error_dim", call = rlang::caller_env()
    )
  }
  map <- as.integer(map)
  .abort_na_cluster_map(map, "map")
  new("ClusterReduction",
      mask = mask_vol,
      map = map,
      cluster_ids = as.integer(sort(unique(map))),
      info = list())
}

#' Convert a ClusteredNeuroVol to a ClusterReduction
#'
#' Bridges the \code{neuroim2::ClusteredNeuroVol} parcellation representation
#' to fmrilatent's \code{ClusterReduction} class, preserving label metadata.
#'
#' @param cvol A \code{ClusteredNeuroVol} object (from \pkg{neuroim2}).
#' @return A \code{ClusterReduction} object with label metadata in \code{info}.
#' @export
as_cluster_reduction <- function(cvol) {
  stopifnot(inherits(cvol, "ClusteredNeuroVol"))
  mask_vol <- neuroim2::mask(cvol)
  clusters <- as.integer(cvol@clusters)
  .abort_na_cluster_map(clusters, "ClusteredNeuroVol clusters")
  info <- list()
  if (!is.null(cvol@label_map) && length(cvol@label_map) > 0L) {
    info$label_map <- cvol@label_map
  }
  if (!is.null(cvol@cluster_map) && length(cvol@cluster_map) > 0L) {
    info$cluster_map <- cvol@cluster_map
  }
  new("ClusterReduction",
      mask = mask_vol,
      map = clusters,
      cluster_ids = as.integer(sort(unique(clusters))),
      info = info)
}
