# Shared loadings handle for heat-wavelet spatial dictionaries

#' Construct a shared LoadingsHandle via heat-wavelet lifting
#'
#' Wraps `lift(reduction, basis_spec, data)` so multiple `LatentNeuroVec`
#' instances can share the same spatial dictionary without embedding the
#' full matrix in each object.
#'
#' @param reduction  Graph/cluster reduction used by `lift()`.
#' @param basis_spec Basis specification, e.g., from `basis_heat_wavelet()`.
#' @param data       Optional data passed through to `lift()` (often NULL).
#' @param id         Optional registry id; provide a stable string to reuse
#'   across sessions. If NULL, a random id is generated.
#' @param label      Optional human-readable label.
#'
#' @return A \code{LoadingsHandle}.
#' @export
heat_wavelet_loadings_handle <- function(reduction,
                                         basis_spec,
                                         data  = NULL,
                                         id    = NULL,
                                         label = "heat-wavelet") {
  if (is.null(id)) {
    id <- paste0("heat-wavelet-", sprintf("%08x", as.integer(stats::runif(1, 0, 2^31))))
  }

  # Materialize once to capture dimensions and seed the registry
  L <- lift(reduction, basis_spec, data = data)
  .latent_register_matrix(id, L, type = "loadings", overwrite = FALSE)

  new("LoadingsHandle",
      id    = id,
      dim   = as.integer(dim(L)),
      kind  = "lifted",
      spec  = list(
        family     = "heat_wavelet",
        reduction  = reduction,
        basis_spec = basis_spec,
        data       = data
      ),
      label = label)
}
