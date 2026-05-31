# Shared loadings handle for heat-wavelet spatial dictionaries

#' Construct a shared LoadingsHandle via heat-wavelet lifting
#'
#' Wraps a heat-wavelet `lift()` call so multiple `LatentNeuroVec`
#' instances can share the same spatial dictionary without embedding the
#' full matrix in each object.
#'
#' @param reduction  Graph/cluster reduction used by `lift()`.
#' @param basis_spec Basis specification; defaults to `basis_heat_wavelet()`.
#' @param data       Ignored by heat-wavelet lifting; accepted only to keep
#'   the lifted-handle constructor signature aligned across families.
#' @param k_neighbors Number of neighbors used for local graph construction
#'   when materializing the lifted basis.
#' @param id         Optional registry id; provide a stable string to reuse
#'   across sessions. If NULL, a deterministic id is derived from the spec and
#'   reduction.
#' @param label      Optional human-readable label.
#'
#' @return A \code{LoadingsHandle}.
#' @export
heat_wavelet_loadings_handle <- function(reduction,
                                         basis_spec = basis_heat_wavelet(),
                                         data  = NULL,
                                         k_neighbors = 6L,
                                         id    = NULL,
                                         label = "heat-wavelet") {
  spec_payload <- list(
    family     = "heat_wavelet",
    reduction  = reduction,
    basis_spec = basis_spec,
    data       = data,
    k_neighbors = k_neighbors
  )
  if (is.null(id)) {
    id <- .latent_handle_id("heat-wavelet", spec_payload)
  }

  cached <- .latent_get_matrix(id, type = "loadings")
  if (!is.null(cached)) {
    handle <- new("LoadingsHandle",
      id    = id,
      dim   = as.integer(dim(cached)),
      kind  = "lifted",
      spec  = spec_payload,
      label = label)
    .latent_get_matrix(
      id, type = "loadings",
      fingerprint = .latent_handle_fingerprint(handle)
    )
    return(handle)
  }

  L <- lift(reduction, basis_spec, data = data, k_neighbors = k_neighbors)
  handle <- new("LoadingsHandle",
      id    = id,
      dim   = as.integer(dim(L)),
      kind  = "lifted",
      spec  = spec_payload,
      label = label)
  .latent_register_matrix(
    id, L, type = "loadings", overwrite = FALSE,
    fingerprint = .latent_handle_fingerprint(handle)
  )
  handle
}
