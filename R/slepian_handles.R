# Slepian / DPSS handles for temporal and spatial components

#' Create a BasisHandle for temporal Slepians (DPSS)
#'
#' @param n_time Integer number of time points.
#' @param tr Repetition time (seconds).
#' @param bandwidth Half-bandwidth in Hz (default 0.1).
#' @param k Optional number of tapers/components.
#' @param backend Backend passed to `dpss_time_basis`; only "tridiag" is
#'   currently supported.
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
  if (identical(backend, "dense")) {
    .encoder_cli_abort(
      "backend = \"dense\" is disabled because it can return a different DPSS subspace under eigenvalue degeneracy; use backend = \"tridiag\".",
      class = "fmrilatent_error_unsupported_dpss_backend"
    )
  }
  if (missing(n_time)) {
    .encoder_cli_abort("n_time must be a positive integer.",
                       class = "fmrilatent_error_invalid_count")
  }
  if (missing(tr)) {
    .encoder_cli_abort("tr must be a positive finite number.",
                       class = "fmrilatent_error_invalid_scalar")
  }
  n_time <- .validate_positive_count(n_time, "n_time")
  tr <- .validate_positive_scalar(tr, "tr")
  bandwidth <- .validate_positive_scalar(bandwidth, "bandwidth")
  if (is.null(k)) {
    NW <- n_time * bandwidth * tr
    k <- floor(2 * NW) - 1L
  } else {
    k <- .validate_nonnegative_count(k, "k")
  }
  k <- max(1L, min(as.integer(k), n_time))
  if (is.null(id)) {
    id <- paste0(
      "slepian-t-",
      digest::digest(
        list(
          kind = "slepian_temporal",
          n_time = n_time,
          tr = tr,
          bandwidth = bandwidth,
          k = k,
          backend = backend
        ),
        algo = "xxhash64"
      )
    )
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
#' @details
#' This constructor lifts the spatial dictionary eagerly so the returned handle
#' records the realized dimensions and registers a fingerprinted cache entry.
#' Repeated constructor calls may therefore recompute the lift even when later
#' `loadings_mat()` calls can reuse the registry cache.
#'
#' @return A \code{LoadingsHandle}.
#' @export
slepian_spatial_loadings_handle <- function(reduction,
                                            basis_spec = basis_slepian(),
                                            data = NULL,
                                            k_neighbors = 6L,
                                            id = NULL,
                                            label = "slepian-spatial") {
  spec_payload <- list(
    family     = "slepian_spatial",
    reduction  = reduction,
    basis_spec = basis_spec,
    data       = data,
    k_neighbors = k_neighbors
  )
  if (is.null(id)) {
    id <- .latent_handle_id("slepian-spatial", spec_payload)
  }
  cached <- .latent_get_matrix(id, type = "loadings")
  if (!is.null(cached)) {
    handle <- new("LoadingsHandle",
      id    = id,
      dim   = as.integer(dim(cached)),
      kind  = "slepian_spatial",
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
      kind  = "slepian_spatial",
      spec  = spec_payload,
      label = label)
  .latent_register_matrix(
    id, L, type = "loadings", overwrite = FALSE,
    fingerprint = .latent_handle_fingerprint(handle)
  )
  handle
}
