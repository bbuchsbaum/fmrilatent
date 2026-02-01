# Diffusion wavelets (randomized, implicit)

#' @include reduction.R
NULL

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Diffusion wavelet basis specification
#'
#' @param target_rank Cap on retained components per scale (keeps runtime bounded).
#' @param oversample Oversampling for randomized range finder.
#' @param threshold Absolute value threshold to enforce sparsity in compressed ops.
#' @param max_scales Maximum diffusion scales to compute.
#' @param epsilon Optional precision (unused in capped-rank path; kept for API parity).
#' @export
basis_diffusion_wavelet <- function(target_rank = 2000L, oversample = 20L,
                                    threshold = 1e-5, max_scales = 1L,
                                    epsilon = NULL) {
  structure(
    list(
      target_rank = as.integer(target_rank),
      oversample = as.integer(oversample),
      threshold = threshold,
      max_scales = as.integer(max_scales),
      epsilon = epsilon
    ),
    class = "spec_diffusion_wavelet"
  )
}

as_linear_operator <- function(T_mat) {
  if (is.function(T_mat)) return(T_mat)
  function(X) T_mat %*% X
}

row_normalize_sparse <- function(adj) {
  rs <- Matrix::rowSums(adj)
  rs[rs == 0] <- 1
  inv_rs <- 1 / rs
  diag_scale <- Matrix::Diagonal(x = inv_rs)
  diag_scale %*% adj
}

randomized_diffusion_step <- function(T_op, n, target_rank, oversample, threshold) {
  k_work <- as.integer(target_rank + oversample)
  Omega <- matrix(stats::rnorm(n * k_work), nrow = n, ncol = k_work)

  # implicit T^2 @ Omega
  Y <- T_op(T_op(Omega))

  qr_y <- base::qr(Y)
  Q <- base::qr.Q(qr_y, complete = FALSE)

  # small operator: B^T B where B = T Q
  B <- T_op(Q)
  T_compressed <- crossprod(B)
  if (!is.null(threshold) && threshold > 0) {
    T_compressed[abs(T_compressed) < threshold] <- 0
  }

  list(Q = Q, T_compressed = T_compressed)
}

diffusion_wavelet_loadings <- function(T_mat, spec) {
  target_rank <- spec$target_rank %||% 2000L
  oversample <- spec$oversample %||% 20L
  threshold <- spec$threshold %||% 1e-5
  max_scales <- spec$max_scales %||% 1L

  T_op <- as_linear_operator(T_mat)
  if (is.null(dim(T_mat))) {
    stop("T_mat must have dimensions or provide a linear operator with known size", call. = FALSE)
  }
  current_dim <- nrow(T_mat)

  phi_global_prev <- NULL
  loadings_list <- list()
  current_op <- T_op

  for (scale in seq_len(max_scales)) {
    step <- randomized_diffusion_step(current_op, current_dim, target_rank, oversample, threshold)
    Q_local <- step$Q
    phi_global <- if (is.null(phi_global_prev)) {
      Q_local
    } else {
      phi_global_prev %*% Q_local
    }
    loadings_list[[length(loadings_list) + 1L]] <- phi_global

    current_dim <- ncol(step$T_compressed)
    if (current_dim < 2L) break
    current_op <- function(X) step$T_compressed %*% X
    phi_global_prev <- phi_global
  }

  Matrix::Matrix(do.call(cbind, loadings_list), sparse = FALSE)
}

build_cluster_graph <- function(reduction, k_neighbors = 6L) {
  if (!requireNamespace("rgsp", quietly = TRUE)) {
    stop("rgsp not installed; install bbuchsbaum/rgsp to use diffusion wavelets", call. = FALSE)
  }

  mask_arr <- .mask_to_array(reduction@mask, "build_cluster_graph")

  coords <- which(as.logical(mask_arr), arr.ind = TRUE)
  ids <- reduction@cluster_ids
  grp <- factor(reduction@map, levels = ids)

  centers <- vapply(seq_len(ncol(coords)), function(d) {
    tapply(coords[, d], grp, mean)
  }, numeric(length(ids)))
  centers <- Matrix::Matrix(centers, sparse = FALSE)

  g <- rgsp::graph_knn(centers, k = k_neighbors, weight = "distance", sym = "union")
  row_normalize_sparse(g$adjacency)
}

#' Lift diffusion wavelets for clustered reduction
#'
#' @param reduction ClusterReduction describing voxel-to-cluster map.
#' @param basis_spec Diffusion wavelet spec (basis_diffusion_wavelet()).
#' @param data Optional (unused) for API symmetry.
#' @param k_neighbors k for graph building (cluster centroids).
#' @param ... Additional arguments (unused).
#' @return Matrix of loadings (voxels x components) concatenating scaling bases across scales.
#' @export
setMethod("lift", signature(reduction = "ClusterReduction", basis_spec = "spec_diffusion_wavelet"),
  function(reduction, basis_spec, data = NULL, k_neighbors = 6L, ...) {
    T_sparse <- build_cluster_graph(reduction, k_neighbors = k_neighbors)
    loadings_coarse <- diffusion_wavelet_loadings(T_sparse, basis_spec)

    # expand cluster loadings back to voxels
    map <- reduction@map
    n_vox <- length(map)
    i <- seq_len(n_vox)
    j <- match(map, reduction@cluster_ids)
    P <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(n_vox, length(reduction@cluster_ids)))
    Matrix::Matrix(P %*% loadings_coarse, sparse = FALSE)
  }
)

#' Diffusion wavelet latent constructor (explicit basis)
#'
#' @param X Matrix time x voxels (mask order).
#' @param mask LogicalNeuroVol or 3D logical array.
#' @param reduction ClusterReduction; if NULL, defaults to one cluster per voxel.
#' @param spec diffusion wavelet spec (basis_diffusion_wavelet()).
#' @param k_neighbors k for graph building.
#' @param label Optional label.
#' @export
diffusion_wavelet_latent <- function(X, mask, reduction = NULL,
                                     spec = basis_diffusion_wavelet(),
                                     k_neighbors = 6L, label = "") {
  mask_arr <- .mask_to_array(mask, "diffusion_wavelet_latent")
  if (is.null(reduction)) {
    map <- seq_len(sum(mask_arr))
    reduction <- make_cluster_reduction(mask, map)
  }

  loadings <- lift(reduction, spec, k_neighbors = k_neighbors)
  basis <- as.matrix(X) %*% as.matrix(loadings)
  spc <- neuroim2::NeuroSpace(c(dim(mask_arr), nrow(X)))
  mask_vol <- if (inherits(mask, "LogicalNeuroVol")) mask else LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dim(mask_arr)))
  meta <- list(family = "diffusion_wavelet", spec = spec, k_neighbors = k_neighbors)
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask_vol, label = label, meta = meta)
}
