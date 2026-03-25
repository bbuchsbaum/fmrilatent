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

#' Basis specifications (lightweight descriptors)
#' @param k Number of components.
#' @param type Basis flavor (implementation-dependent).
#' @param whiten Whether to whiten PCA scores.
#' @export
basis_slepian <- function(k = 3, type = "laplacian") {
  structure(list(k = k, type = type), class = "spec_slepian")
}

#' @rdname basis_slepian
#' @export
basis_pca <- function(k = 3, whiten = FALSE) {
  structure(list(k = k, whiten = whiten), class = "spec_pca")
}

#' @rdname basis_slepian
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
#' @export
setMethod("lift", signature(reduction = "GraphReduction", basis_spec = "ANY"),
  function(reduction, basis_spec, data = NULL, ...) {
    stop("No lift() implementation for this reduction/basis combination. ",
         "Provide an external method (e.g., supervoxel Slepian or parcel PCA).")
  }
)

# --- ClusterReduction constructors -------------------------------------------

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
    stop(
      "map must have length ", n_vox, " (the number of voxels in the mask), got ",
      length(map), ".",
      call. = FALSE
    )
  }
  new("ClusterReduction",
      mask = mask_vol,
      map = as.integer(map),
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
  clusters <- cvol@clusters
  info <- list()
  if (!is.null(cvol@label_map) && length(cvol@label_map) > 0L) {
    info$label_map <- cvol@label_map
  }
  if (!is.null(cvol@cluster_map) && length(cvol@cluster_map) > 0L) {
    info$cluster_map <- cvol@cluster_map
  }
  new("ClusterReduction",
      mask = mask_vol,
      map = as.integer(clusters),
      cluster_ids = as.integer(sort(unique(clusters))),
      info = info)
}
