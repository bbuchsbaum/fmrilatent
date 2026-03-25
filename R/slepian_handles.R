# Slepian / DPSS handles for temporal and spatial components

#' Create a BasisHandle for temporal Slepians (DPSS)
#'
#' @param n_time Integer number of time points.
#' @param tr Repetition time (seconds).
#' @param bandwidth Half-bandwidth in Hz (default 0.1).
#' @param k Optional number of tapers/components.
#' @param backend Backend passed to `dpss_time_basis` ("tridiag" or "dense").
#' @param id Optional registry key (generated if NULL).
#' @param label Optional human-readable label.
#'
#' @return A \code{BasisHandle}.
#' @export
slepian_temporal_handle <- function(n_time,
                                    tr,
                                    bandwidth = 0.1,
                                    k = NULL,
                                    backend = c("tridiag", "dense"),
                                    id = NULL,
                                    label = NULL) {
  backend <- match.arg(backend)
  n_time <- as.integer(n_time)
  if (is.null(k)) {
    NW <- n_time * bandwidth * tr
    k <- floor(2 * NW) - 1L
  }
  k <- as.integer(k)
  if (is.null(id)) {
    id <- sprintf("slepian-t-%d-%.4f-%d-%s", n_time, bandwidth, k, backend)
  }
  if (is.null(label)) {
    label <- sprintf("Slepian_t(n=%d,k=%d,W=%.4f)", n_time, k, bandwidth)
  }
  new("BasisHandle",
      id    = id,
      dim   = as.integer(c(n_time, k)),
      kind  = "slepian_temporal",
      spec  = list(
        n_time = n_time,
        tr = tr,
        bandwidth = bandwidth,
        k = k,
        backend = backend
      ),
      label = label)
}

#' Create a LoadingsHandle for spatial Slepians (graph Laplacian)
#'
#' @param reduction Graph reduction (e.g., ClusterReduction).
#' @param basis_spec Slepian basis spec (from `basis_slepian()`).
#' @param data Optional data passed to `lift()` (if needed).
#' @param k_neighbors Number of neighbors used for local graph construction
#'   when materializing the lifted basis.
#' @param id Optional registry id; generated if NULL.
#' @param label Optional label.
#'
#' @return A \code{LoadingsHandle}.
#' @export
slepian_spatial_loadings_handle <- function(reduction,
                                            basis_spec = basis_slepian(),
                                            data = NULL,
                                            k_neighbors = 6L,
                                            id = NULL,
                                            label = "slepian-spatial") {
  if (is.null(id)) {
    id <- paste0("slepian-spatial-",
                 sprintf("%08x", as.integer(stats::runif(1, 0, 2^31))))
  }
  L <- lift(reduction, basis_spec, data = data, k_neighbors = k_neighbors)
  .latent_register_matrix(id, L, type = "loadings", overwrite = FALSE)

  new("LoadingsHandle",
      id    = id,
      dim   = as.integer(dim(L)),
      kind  = "slepian_spatial",
      spec  = list(
        family     = "slepian_spatial",
        reduction  = reduction,
        basis_spec = basis_spec,
        data       = data,
        k_neighbors = k_neighbors
      ),
      label = label)
}
