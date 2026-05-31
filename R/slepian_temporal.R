#' DPSS temporal basis (Slepian sequences)
#'
#' Generates Discrete Prolate Spheroidal Sequences (DPSS) for a given series
#' length and bandwidth. Uses an internal RcppEigen solver over the prolate
#' matrix; suitable for moderate `n_time`. For very long series, a tridiagonal
#' solver can replace the backend with the same interface.
#'
#' @param n_time Integer length of the time dimension.
#' @param tr Repetition time in seconds.
#' @param bandwidth Half-bandwidth in Hz (typical BOLD range 0.008–0.1).
#' @param k Optional number of tapers; defaults to \code{floor(2 * NW) - 1}
#'   where \code{NW = n_time * bandwidth * tr}. Clamped to \code{[1, n_time]}.
#' @param backend DPSS backend. Only \code{"tridiag"} is currently supported.
#'
#' @return Numeric matrix (n_time x k) with orthonormal columns.
#' @export
dpss_time_basis <- function(n_time, tr, bandwidth, k = NULL,
                            backend = c("tridiag", "dense")) {
  if (missing(n_time)) {
    .encoder_cli_abort("n_time must be positive", class = "fmrilatent_error_invalid_argument")
  }
  if (missing(tr)) {
    .encoder_cli_abort("tr must be positive (seconds)", class = "fmrilatent_error_invalid_argument")
  }
  if (missing(bandwidth)) {
    .encoder_cli_abort("bandwidth must be positive (Hz)", class = "fmrilatent_error_invalid_argument")
  }
  n_time <- .validate_positive_count(n_time, "n_time")
  tr <- .validate_positive_scalar(tr, "tr")
  bandwidth <- .validate_positive_scalar(bandwidth, "bandwidth")
  W <- bandwidth * tr                    # normalized half-bandwidth (cycles/sample)
  NW <- n_time * W
  if (is.null(k)) {
    k <- floor(2 * NW) - 1
  } else {
    k <- .validate_nonnegative_count(k, "k")
  }
  k <- max(1L, min(as.integer(k), n_time))
  backend <- match.arg(backend)
  if (identical(backend, "dense")) {
    .encoder_cli_abort(
      "backend = \"dense\" is disabled because it can return a different DPSS subspace under eigenvalue degeneracy; use backend = \"tridiag\".",
      class = "fmrilatent_error_unsupported_dpss_backend"
    )
  }
  if (backend == "tridiag") {
    generate_dpss_tridiag_rcpp(as.integer(n_time), as.numeric(NW), as.integer(k))
  } else {
    generate_dpss_basis_rcpp(as.integer(n_time), as.numeric(NW), as.integer(k))
  }
}

#' LatentNeuroVec using a temporal DPSS basis
#'
#' @param X Numeric matrix (time x voxels within mask).
#' @param mask `LogicalNeuroVol`.
#' @param tr Repetition time in seconds.
#' @param bandwidth Half-bandwidth in Hz (default 0.1).
#' @param k Optional number of tapers; see `dpss_time_basis`.
#' @param denoise If TRUE, truncate to `floor(2 * NW) - 1` (Shannon number).
#' @param backend DPSS computation backend passed to `dpss_time_basis`; only
#'   \code{"tridiag"} is currently supported.
#' @param label Optional character label.
#' @return `LatentNeuroVec` with DPSS temporal basis and voxel loadings.
#' @export
slepian_temporal_latent <- function(X, mask, tr, bandwidth = 0.1, k = NULL,
                                    denoise = TRUE, backend = c("tridiag", "dense"),
                                    label = "") {
  n_time <- nrow(X)
  if (is.null(n_time) || n_time < 1) {
    .encoder_cli_abort("X must have time in rows", class = "fmrilatent_error_invalid_argument")
  }
  if (missing(tr)) {
    .encoder_cli_abort("tr must be positive (seconds)", class = "fmrilatent_error_invalid_argument")
  }
  tr <- .validate_positive_scalar(tr, "tr")
  bandwidth <- .validate_positive_scalar(bandwidth, "bandwidth")
  denoise <- .validate_flag_scalar(denoise, "denoise")

  W <- bandwidth * tr
  NW <- n_time * W
  shannon_k <- max(1L, min(as.integer(floor(2 * NW) - 1L), n_time))
  if (is.null(k)) {
    k_use <- shannon_k
  } else {
    k_use <- .validate_nonnegative_count(k, "k")
  }
  k_use <- max(1L, min(as.integer(k_use), n_time))
  if (denoise) {
    k_use <- min(k_use, shannon_k)
  }

  B <- dpss_time_basis(n_time, tr = tr, bandwidth = bandwidth, k = k_use,
                       backend = match.arg(backend))
  offset <- colMeans(X)
  X_centered <- sweep(as.matrix(X), 2L, offset, "-")
  L <- .temporal_loadings_from_basis(
    X_centered,
    B,
    context = "slepian_temporal_latent"
  )
  meta <- list(
    family = "slepian_temporal",
    k = k_use,
    tr = tr,
    bandwidth = bandwidth,
    nw = NW,
    denoise = denoise
  )
  .make_latent_neurovector(
    X, mask, L,
    basis = Matrix::Matrix(B),
    offset = offset,
    label = label,
    meta = meta,
    location = "slepian_temporal_latent"
  )
}
