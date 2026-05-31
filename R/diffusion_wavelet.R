# Diffusion wavelets (randomized, implicit)

#' @include reduction.R
NULL

#' Diffusion wavelet basis specification
#'
#' @param target_rank Cap on retained components per scale (keeps runtime bounded).
#' @param oversample Oversampling for randomized range finder.
#' @param sparsify_eps Absolute value threshold to enforce sparsity in
#'   compressed operators. Stored as `threshold` for compatibility.
#' @param threshold Deprecated alias for `sparsify_eps`.
#' @param max_scales Maximum diffusion scales to compute.
#' @param epsilon Optional precision (unused in capped-rank path; kept for API parity).
#' @param seed Optional integer seed for deterministic randomized range finding.
#' @return A `spec_diffusion_wavelet` basis specification descriptor.
#' @export
basis_diffusion_wavelet <- function(target_rank = 2000L, oversample = 20L,
                                    threshold = NULL, max_scales = 1L,
                                    epsilon = NULL, seed = 1L,
                                    sparsify_eps = NULL) {
  if (!is.null(sparsify_eps) && !is.null(threshold) &&
      !identical(sparsify_eps, threshold)) {
    .encoder_cli_abort(
      "Specify only one of sparsify_eps or threshold.",
      class = "fmrilatent_error_invalid_argument"
    )
  }
  sparsify_eps <- sparsify_eps %||% threshold %||% 1e-5
  spec <- .validate_diffusion_wavelet_spec(list(
    target_rank = target_rank,
    oversample = oversample,
    threshold = sparsify_eps,
    sparsify_eps = sparsify_eps,
    max_scales = max_scales,
    epsilon = epsilon,
    seed = seed
  ))
  structure(spec, class = "spec_diffusion_wavelet")
}

.validate_diffusion_wavelet_spec <- function(spec) {
  if (!is.list(spec)) {
    .encoder_cli_abort(
      "basis_spec must be a diffusion-wavelet specification list.",
      class = "fmrilatent_error_invalid_diffusion_spec"
    )
  }
  target_rank <- .validate_positive_count(spec$target_rank %||% 2000L, "target_rank")
  oversample <- .validate_nonnegative_count(spec$oversample %||% 20L, "oversample")
  threshold <- .validate_nonnegative_scalar(spec$sparsify_eps %||% spec$threshold %||% 1e-5,
                                            "sparsify_eps")
  max_scales <- .validate_positive_count(spec$max_scales %||% 1L, "max_scales")
  epsilon <- spec$epsilon
  if (!is.null(epsilon)) {
    epsilon <- .validate_positive_scalar(epsilon, "epsilon")
  }
  seed <- spec$seed
  if (!is.null(seed)) {
    seed <- .validate_nonnegative_count(seed, "seed")
  }

  list(
    target_rank = target_rank,
    oversample = oversample,
    threshold = threshold,
    sparsify_eps = threshold,
    max_scales = max_scales,
    epsilon = epsilon,
    seed = seed
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

  # Galerkin restriction of T onto range(Q): Q^T T Q.
  B <- T_op(Q)
  T_compressed <- crossprod(Q, B)
  if (!is.null(threshold) && threshold > 0) {
    T_compressed[abs(T_compressed) < threshold] <- 0
  }

  list(Q = Q, T_compressed = T_compressed)
}

diffusion_wavelet_loadings <- function(T_mat, spec, n = NULL) {
  spec <- .validate_diffusion_wavelet_spec(spec)
  target_rank <- spec$target_rank %||% 2000L
  oversample <- spec$oversample %||% 20L
  threshold <- spec$threshold %||% 1e-5
  max_scales <- spec$max_scales %||% 1L
  seed <- spec$seed
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  T_op <- as_linear_operator(T_mat)
  T_dim <- dim(T_mat)
  if (is.null(T_dim)) {
    if (is.null(n)) {
      .encoder_cli_abort(
        "Function T_mat operators require explicit n.",
        class = "fmrilatent_error_invalid_diffusion_operator"
      )
    }
    current_dim <- .validate_positive_count(n, "n")
  } else {
    if (length(T_dim) != 2L || T_dim[1L] != T_dim[2L]) {
      .encoder_cli_abort(
        "T_mat must be a square matrix-like operator.",
        class = "fmrilatent_error_invalid_diffusion_operator"
      )
    }
    current_dim <- nrow(T_mat)
  }
  if (!is.function(T_op)) {
    .encoder_cli_abort(
      "T_mat must be matrix-like or a function.",
      class = "fmrilatent_error_invalid_diffusion_operator"
    )
  }

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
    current_op <- local({
      tc <- step$T_compressed
      function(X) tc %*% X
    })
    phi_global_prev <- phi_global
  }

  Matrix::Matrix(do.call(cbind, loadings_list), sparse = FALSE)
}

build_cluster_graph <- function(reduction, k_neighbors = 6L) {
  .build_cluster_centroid_graph(reduction, k_neighbors = k_neighbors)
}

#' Lift diffusion wavelets for clustered reduction
#'
#' @param reduction ClusterReduction describing voxel-to-cluster map.
#' @param basis_spec Diffusion wavelet spec (basis_diffusion_wavelet()).
#' @param data Ignored for this graph-only spatial dictionary; accepted only
#'   for the shared `lift()` method signature.
#' @param k_neighbors k for graph building (cluster centroids).
#' @param ... Additional arguments (unused).
#' @return Matrix of loadings (voxels x components) concatenating scaling bases across scales.
#' @export
setMethod("lift", signature(reduction = "ClusterReduction", basis_spec = "spec_diffusion_wavelet"),
  function(reduction, basis_spec, data = NULL, k_neighbors = 6L, ...) {
    T_sparse <- build_cluster_graph(reduction, k_neighbors = k_neighbors)
    loadings_coarse <- diffusion_wavelet_loadings(T_sparse, basis_spec)

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
#' @return A `LatentNeuroVec` object.
#' @export
diffusion_wavelet_latent <- function(X, mask, reduction = NULL,
                                     spec = basis_diffusion_wavelet(),
                                     k_neighbors = 6L, label = "") {
  mask_arr <- .extract_mask_array(mask, "diffusion_wavelet_latent")
  if (is.null(reduction)) {
    map <- seq_len(sum(mask_arr))
    reduction <- make_cluster_reduction(mask, map)
  }

  loadings <- lift(reduction, spec, k_neighbors = k_neighbors)
  meta <- list(family = "diffusion_wavelet", spec = spec, k_neighbors = k_neighbors)
  .make_latent_neurovector(X, mask, loadings, label = label, meta = meta,
                           location = "diffusion_wavelet_latent",
                           expect_dense = TRUE)
}
