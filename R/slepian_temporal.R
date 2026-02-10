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
#' @param backend Either \code{"tridiag"} (default, O(n^2)) or \code{"dense"} (O(n^3), debugging).
#'
#' @return Numeric matrix (n_time x k) with orthonormal columns.
#' @export
dpss_time_basis <- function(n_time, tr, bandwidth, k = NULL,
                            backend = c("tridiag", "dense")) {
  if (missing(n_time) || n_time < 1) stop("n_time must be positive")
  if (missing(tr) || tr <= 0) stop("tr must be positive (seconds)")
  if (missing(bandwidth) || bandwidth <= 0) stop("bandwidth must be positive (Hz)")
  W <- bandwidth * tr                    # normalized half-bandwidth (cycles/sample)
  NW <- n_time * W
  if (is.null(k)) k <- floor(2 * NW) - 1
  k <- max(1L, min(as.integer(k), n_time))
  backend <- match.arg(backend)
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
#' @param backend DPSS computation backend passed to `dpss_time_basis`.
#' @param label Optional character label.
#' @return `LatentNeuroVec` with DPSS temporal basis and voxel loadings.
#' @export
slepian_temporal_latent <- function(X, mask, tr, bandwidth = 0.1, k = NULL,
                                    denoise = TRUE, backend = c("tridiag", "dense"),
                                    label = "") {
  mask_arr <- extract_mask_array_logical(mask, "slepian_temporal_latent")
  n_time <- nrow(X)
  if (is.null(n_time) || n_time < 1) stop("X must have time in rows")

  W <- bandwidth * tr
  NW <- n_time * W
  if (is.null(k)) {
    k_use <- floor(2 * NW) - 1
  } else {
    k_use <- k
  }
  k_use <- max(1L, min(as.integer(k_use), n_time))

  B <- dpss_time_basis(n_time, tr = tr, bandwidth = bandwidth, k = k_use,
                       backend = match.arg(backend))
  L <- Matrix::crossprod(X, B)           # voxels x k
  spc <- neuroim2::NeuroSpace(c(dim(mask_arr), n_time))
  meta <- list(
    family = "slepian_temporal",
    tr = tr,
    bandwidth = bandwidth,
    nw = NW
  )
  LatentNeuroVec(
    basis = Matrix::Matrix(B),
    loadings = Matrix::Matrix(L),
    space = spc,
    mask = mask,
    offset = colMeans(X),
    label = label,
    meta = meta
  )
}
