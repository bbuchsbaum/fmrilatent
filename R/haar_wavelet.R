# Haar lifting wavelet utilities (mask-adaptive, Morton-ordered)
#
# Pure-R implementation with optional Rcpp acceleration. Exposes
# user-facing helpers to run forward/inverse transforms and a light
# latent wrapper.

# Internal null-coalescing helper (kept local to avoid namespace noise)
`%||%` <- function(x, y) if (is.null(x)) y else x

use_haar_rcpp <- function() {
  opt_new <- getOption("fmrilatent.haar.use_rcpp", NULL)
  opt_old <- getOption("fmrilatent.hwt.use_rcpp", NULL)
  flag <- opt_new %||% opt_old %||% TRUE
  isTRUE(flag) && exists("forward_lift_rcpp", mode = "function")
}

as_logical_mask <- function(arr, location = "haar_wavelet:mask") {
  if (inherits(arr, "LogicalNeuroVol")) {
    arr <- as.array(arr)
  }
  dims <- dim(arr)
  if (is.null(dims) || length(dims) != 3L) {
    stop(sprintf("mask must be a 3D array (%s)", location), call. = FALSE)
  }
  array(as.logical(arr), dim = dims)
}

# -- Morton ordering helpers -------------------------------------------------

#' Morton-ordered voxel indices from a 3D mask
#'
#' @keywords internal
get_morton_ordered_indices <- function(mask_3d_array, z_order_seed = 42L) {
  mask_logical <- as_logical_mask(mask_3d_array, "get_morton_ordered_indices")

  if (use_haar_rcpp() && exists("get_morton_ordered_indices_rcpp", mode = "function")) {
    return(get_morton_ordered_indices_rcpp(mask_logical))
  }

  vox_coords <- which(mask_logical, arr.ind = TRUE)
  if (nrow(vox_coords) == 0L) {
    return(integer(0))
  }

  dims <- dim(mask_logical)
  x <- vox_coords[, 1] - 1L
  y <- vox_coords[, 2] - 1L
  z <- vox_coords[, 3] - 1L

  max_dim <- max(dims)
  bits <- ceiling(log2(max_dim))
  codes <- integer(nrow(vox_coords))
  for (b in seq_len(bits)) {
    shift <- b - 1L
    codes <- bitwOr(codes,
                    bitwShiftL(bitwAnd(bitwShiftR(x, shift), 1L), 3L * shift))
    codes <- bitwOr(codes,
                    bitwShiftL(bitwAnd(bitwShiftR(y, shift), 1L), 3L * shift + 1L))
    codes <- bitwOr(codes,
                    bitwShiftL(bitwAnd(bitwShiftR(z, shift), 1L), 3L * shift + 2L))
  }

  ordering <- do.call(order, list(codes, x, y, z, seq_len(nrow(vox_coords))))
  ordered_coords <- vox_coords[ordering, , drop = FALSE]

  linear_idx <- (ordered_coords[, 1] - 1L) +
    (ordered_coords[, 2] - 1L) * dims[1L] +
    (ordered_coords[, 3] - 1L) * dims[1L] * dims[2L] + 1L

  as.integer(linear_idx)
}

#' Precompute scaling factors for Haar lifting
#'
#' @keywords internal
precompute_haar_scalings <- function(mask_3d_array, levels) {
  mask_logical <- as_logical_mask(mask_3d_array, "precompute_haar_scalings")
  levels <- as.integer(levels)
  if (length(levels) != 1L || is.na(levels) || levels < 1L) {
    stop("levels must be a positive integer", call. = FALSE)
  }

  if (use_haar_rcpp() && exists("precompute_haar_scalings_rcpp", mode = "function")) {
    return(precompute_haar_scalings_rcpp(mask_logical, levels))
  }

  scalings <- vector("list", levels)
  current <- mask_logical
  for (lvl in seq_len(levels)) {
    dims <- dim(current)
    x_seq <- seq(1L, dims[1], by = 2L)
    y_seq <- seq(1L, dims[2], by = 2L)
    z_seq <- seq(1L, dims[3], by = 2L)
    counts <- integer(length(x_seq) * length(y_seq) * length(z_seq))
    idx <- 1L
    for (x in x_seq) {
      for (y in y_seq) {
        for (z in z_seq) {
          block <- current[
            x:min(x + 1L, dims[1]),
            y:min(y + 1L, dims[2]),
            z:min(z + 1L, dims[3])
          ]
          counts[idx] <- sum(block)
          idx <- idx + 1L
        }
      }
    }
    scalings[[lvl]] <- list(
      sqrt_nvalid = sqrt(as.numeric(counts)),
      sqrt_nvalid_div_8 = sqrt(as.numeric(counts) / 8)
    )
    current <- array(counts > 0,
                     dim = c(length(x_seq), length(y_seq), length(z_seq)))
  }
  scalings
}

#' Morton codes of valid finest-level blocks
#'
#' @keywords internal
get_valid_finest_blocks <- function(mask_3d_array) {
  mask_logical <- as_logical_mask(mask_3d_array, "get_valid_finest_blocks")

  if (use_haar_rcpp() && exists("get_valid_finest_blocks_rcpp", mode = "function")) {
    return(get_valid_finest_blocks_rcpp(mask_logical))
  }

  dims <- dim(mask_logical)
  x_seq <- seq(1L, dims[1], by = 2L)
  y_seq <- seq(1L, dims[2], by = 2L)
  z_seq <- seq(1L, dims[3], by = 2L)
  block_mask <- array(FALSE,
                      dim = c(length(x_seq), length(y_seq), length(z_seq)))

  for (i in seq_along(x_seq)) {
    for (j in seq_along(y_seq)) {
      for (k in seq_along(z_seq)) {
        block <- mask_logical[
          x_seq[i]:min(x_seq[i] + 1L, dims[1]),
          y_seq[j]:min(y_seq[j] + 1L, dims[2]),
          z_seq[k]:min(z_seq[k] + 1L, dims[3])
        ]
        block_mask[i, j, k] <- any(block)
      }
    }
  }

  if (!any(block_mask)) return(integer(0))

  coords <- which(block_mask, arr.ind = TRUE)
  x <- coords[, 1] - 1L
  y <- coords[, 2] - 1L
  z <- coords[, 3] - 1L

  max_dim <- max(dim(block_mask))
  bits <- ceiling(log2(max_dim))
  codes <- integer(nrow(coords))
  for (b in seq_len(bits)) {
    shift <- b - 1L
    codes <- bitwOr(codes,
                    bitwShiftL(bitwAnd(bitwShiftR(x, shift), 1L), 3L * shift))
    codes <- bitwOr(codes,
                    bitwShiftL(bitwAnd(bitwShiftR(y, shift), 1L), 3L * shift + 1L))
    codes <- bitwOr(codes,
                    bitwShiftL(bitwAnd(bitwShiftR(z, shift), 1L), 3L * shift + 2L))
  }

  sort(as.integer(codes))
}

# -- Core lifting (pure R) ---------------------------------------------------

forward_lift_R <- function(data_masked_morton_ordered,
                           mask_flat_morton_ordered,
                           mask_dims,
                           levels,
                           scaling_factors_per_level) {
  current <- data_masked_morton_ordered
  detail_list <- vector("list", levels)
  for (lvl in seq_len(levels)) {
    sc <- scaling_factors_per_level[[lvl]]
    counts <- as.integer(round(sc$sqrt_nvalid^2))
    next_data <- numeric(sum(counts > 0L))
    sum_counts <- sum(counts)
    dvec <- numeric(sum_counts)
    idx_in <- 1L
    idx_out <- 1L
    idx_lp  <- 1L
    for (b in seq_along(counts)) {
      nv <- counts[b]
      if (nv > 0) {
        vals <- current[idx_in:(idx_in + nv - 1L)]
        avg <- mean(vals)
        next_data[idx_lp] <- avg * sc$sqrt_nvalid[b]
        dvec[idx_out:(idx_out + nv - 1L)] <- (vals - avg) * sc$sqrt_nvalid_div_8[b]
        idx_in <- idx_in + nv
        idx_out <- idx_out + nv
        idx_lp <- idx_lp + 1L
      }
    }
    detail_list[[lvl]] <- dvec
    current <- next_data
  }
  list(root_coeff = current, detail_coeffs_by_level = detail_list)
}

inverse_lift_R <- function(root_coeff,
                           detail_coeffs_by_level,
                           mask_flat_morton_ordered,
                           mask_dims,
                           levels,
                           scaling_factors_per_level) {
  current <- root_coeff
  for (lvl in rev(seq_len(levels))) {
    sc <- scaling_factors_per_level[[lvl]]
    counts <- as.integer(round(sc$sqrt_nvalid^2))
    total <- sum(counts)
    next_data <- numeric(total)
    idx_out <- 1L
    idx_det <- 1L
    idx_lp  <- 1L
    detail_vec <- detail_coeffs_by_level[[lvl]]
    if (length(detail_vec) != total) {
      stop(sprintf("Detail vector for level %d has incorrect length", lvl - 1L),
           call. = FALSE)
    }
    for (b in seq_along(counts)) {
      nv <- counts[b]
      if (nv > 0) {
        if (idx_lp > length(current)) {
          stop(sprintf("Lowpass underflow at level %d", lvl - 1L), call. = FALSE)
        }
        avg <- current[idx_lp] / sc$sqrt_nvalid[b]
        vals <- detail_vec[idx_det:(idx_det + nv - 1L)] /
          sc$sqrt_nvalid_div_8[b] + avg
        next_data[idx_out:(idx_out + nv - 1L)] <- vals
        idx_out <- idx_out + nv
        idx_det <- idx_det + nv
        idx_lp  <- idx_lp + 1L
      }
    }
    current <- next_data
  }
  current
}

#' Perform forward Haar lifting analysis
#'
#' @keywords internal
perform_haar_lift_analysis <- function(data_matrix_T_x_Nmask, mask_3d_array,
                                       levels, z_order_seed = 42L) {
  if (!is.matrix(data_matrix_T_x_Nmask)) {
    stop("data_matrix_T_x_Nmask must be a matrix", call. = FALSE)
  }

  mask_logical <- as_logical_mask(mask_3d_array, "perform_haar_lift_analysis")
  mask_dims <- dim(mask_logical)
  expected_cols <- sum(mask_logical)
  if (ncol(data_matrix_T_x_Nmask) != expected_cols) {
    if (nrow(data_matrix_T_x_Nmask) == expected_cols) {
      data_matrix_T_x_Nmask <- t(data_matrix_T_x_Nmask)
    } else {
      stop(sprintf("data matrix must have %d columns (one per masked voxel); got %d",
                   expected_cols, ncol(data_matrix_T_x_Nmask)), call. = FALSE)
    }
  }

  morton_idx <- get_morton_ordered_indices(mask_logical, z_order_seed)
  if (length(morton_idx) == 0L) {
    stop("mask must contain at least one voxel", call. = FALSE)
  }
  mask_linear <- which(mask_logical)
  perm <- match(morton_idx, mask_linear)
  data_morton <- data_matrix_T_x_Nmask[, perm, drop = FALSE]

  full_order <- get_morton_ordered_indices(array(TRUE, dim(mask_logical)),
                                           z_order_seed)
  mask_flat <- as.vector(mask_logical)
  mask_flat_morton <- mask_flat[full_order]

  scalings <- precompute_haar_scalings(mask_logical, levels)

  counts_by_level <- lapply(scalings, function(sc) as.integer(round(sc$sqrt_nvalid^2)))
  top_blocks <- length(counts_by_level[[levels]])
  if (top_blocks != 1L) {
    stop("Top-level lowpass must be a single block; increase 'levels'", call. = FALSE)
  }

  detail_lengths <- vapply(counts_by_level, sum, integer(1))

  n_time <- nrow(data_morton)
  root_coeffs <- matrix(0, nrow = n_time, ncol = top_blocks)
  detail_coeffs <- lapply(seq_len(levels), function(lv) {
    matrix(0, nrow = n_time, ncol = detail_lengths[lv])
  })

  use_rcpp <- use_haar_rcpp()

  for (tt in seq_len(n_time)) {
    if (use_rcpp) {
      res <- forward_lift_rcpp(
        data_morton[tt, ],
        mask_flat_morton,
        mask_dims,
        levels,
        scalings
      )
    } else {
      res <- forward_lift_R(
        data_morton[tt, ],
        mask_flat_morton,
        mask_dims,
        levels,
        scalings
      )
    }
    root_coeffs[tt, 1] <- res$root_coeff
    for (lvl in seq_len(levels)) {
      detail_vec <- res$detail_coeffs_by_level[[lvl]]
      if (length(detail_vec) != detail_lengths[lvl]) {
        stop(sprintf("Forward lift detail length mismatch at level %d", lvl - 1L),
             call. = FALSE)
      }
      detail_coeffs[[lvl]][tt, ] <- detail_vec
    }
  }

  list(root = root_coeffs, detail = detail_coeffs)
}

#' Perform inverse Haar lifting synthesis
#'
#' @keywords internal
perform_haar_lift_synthesis <- function(coeff_list, mask_3d_array, levels,
                                         z_order_seed = 42L) {
  if (is.null(coeff_list$root) || !is.matrix(coeff_list$root)) {
    stop("coeff_list$root must be a matrix", call. = FALSE)
  }
  if (!is.list(coeff_list$detail) || length(coeff_list$detail) < levels) {
    stop("coeff_list$detail must be a list of length 'levels'", call. = FALSE)
  }

  mask_logical <- as_logical_mask(mask_3d_array, "perform_haar_lift_synthesis")
  mask_dims <- dim(mask_logical)
  morton_idx <- get_morton_ordered_indices(mask_logical, z_order_seed)
  if (length(morton_idx) == 0L) {
    stop("mask must contain at least one voxel", call. = FALSE)
  }
  mask_linear <- which(mask_logical)
  perm <- match(morton_idx, mask_linear)
  inv_perm <- order(perm)

  full_order <- get_morton_ordered_indices(array(TRUE, dim(mask_logical)),
                                           z_order_seed)
  mask_flat <- as.vector(mask_logical)
  mask_flat_morton <- mask_flat[full_order]

  scalings <- precompute_haar_scalings(mask_logical, levels)
  counts_by_level <- lapply(scalings, function(sc) as.integer(round(sc$sqrt_nvalid^2)))
  top_blocks <- length(counts_by_level[[levels]])
  if (ncol(coeff_list$root) != top_blocks) {
    stop("Root coefficient matrix has wrong number of columns", call. = FALSE)
  }
  detail_expected <- vapply(counts_by_level, sum, integer(1))
  for (lv in seq_len(levels)) {
    if (ncol(coeff_list$detail[[lv]]) != detail_expected[lv]) {
      stop(sprintf("Detail matrix for level %d must have %d columns", lv - 1L,
                   detail_expected[lv]), call. = FALSE)
    }
  }

  n_time <- nrow(coeff_list$root)
  reco_morton <- matrix(0, nrow = n_time, ncol = length(morton_idx))

  use_rcpp <- use_haar_rcpp()
  for (tt in seq_len(n_time)) {
    detail_vecs <- lapply(seq_len(levels), function(lvl) coeff_list$detail[[lvl]][tt, ])
    if (use_rcpp) {
      reco_morton[tt, ] <- inverse_lift_rcpp(
        coeff_list$root[tt, 1],
        detail_vecs,
        mask_flat_morton,
        mask_dims,
        levels,
        scalings
      )
    } else {
      reco_morton[tt, ] <- inverse_lift_R(
        coeff_list$root[tt, 1],
        detail_vecs,
        mask_flat_morton,
        mask_dims,
        levels,
        scalings
      )
    }
  }

  reco_morton[, inv_perm, drop = FALSE]
}

# -- ROI mapping -------------------------------------------------------------

encode_morton3d <- function(x, y, z, bits) {
  code <- 0L
  for (b in seq_len(bits)) {
    shift <- b - 1L
    code <- bitwOr(code, bitwShiftL(bitwAnd(bitwShiftR(x, shift), 1L), 3L * shift))
    code <- bitwOr(code, bitwShiftL(bitwAnd(bitwShiftR(y, shift), 1L), 3L * shift + 1L))
    code <- bitwOr(code, bitwShiftL(bitwAnd(bitwShiftR(z, shift), 1L), 3L * shift + 2L))
  }
  as.integer(code)
}

compute_block_map <- function(mask_3d_array, levels) {
  mask_logical <- as_logical_mask(mask_3d_array, "compute_block_map")
  mapping <- vector("list", levels)
  current <- mask_logical
  for (lvl in seq_len(levels)) {
    dims <- dim(current)
    x_seq <- seq(1L, dims[1], by = 2L)
    y_seq <- seq(1L, dims[2], by = 2L)
    z_seq <- seq(1L, dims[3], by = 2L)
    nbx <- length(x_seq); nby <- length(y_seq); nbz <- length(z_seq)
    nblocks <- nbx * nby * nbz
    counts <- integer(nblocks)
    codes <- integer(nblocks)
    bits <- ceiling(log2(max(nbx, nby, nbz)))
    idx <- 1L
    for (i in seq_along(x_seq)) {
      for (j in seq_along(y_seq)) {
        for (k in seq_along(z_seq)) {
          block <- current[
            x_seq[i]:min(x_seq[i] + 1L, dims[1]),
            y_seq[j]:min(y_seq[j] + 1L, dims[2]),
            z_seq[k]:min(z_seq[k] + 1L, dims[3])
          ]
          counts[idx] <- sum(block)
          codes[idx] <- encode_morton3d(i - 1L, j - 1L, k - 1L, bits)
          idx <- idx + 1L
        }
      }
    }
    start <- c(0L, cumsum(counts)[-length(counts)])
    mapping[[lvl]] <- list(code = codes, count = counts, start = start)
    current <- array(counts > 0L, dim = c(nbx, nby, nbz))
  }
  mapping
}

get_roi_detail_indices <- function(roi_mask, mask_3d_array, levels) {
  if (!is.array(roi_mask) || !identical(dim(roi_mask), dim(mask_3d_array))) {
    stop("roi_mask must match mask dimensions", call. = FALSE)
  }

  mapping <- compute_block_map(mask_3d_array, levels)
  roi_codes_finest <- get_valid_finest_blocks(roi_mask)
  if (length(roi_codes_finest) == 0L) {
    return(lapply(seq_len(levels), function(x) integer()))
  }

  codes_by_level <- vector("list", levels)
  codes_by_level[[levels]] <- roi_codes_finest
  if (levels > 1) {
    for (lvl in seq_len(levels - 1)) {
      shift <- 3L * (levels - lvl)
      codes_by_level[[lvl]] <- unique(roi_codes_finest %/% (2L ^ shift))
    }
  }

  indices_by_level <- vector("list", levels)
  for (lvl in seq_len(levels)) {
    codes <- codes_by_level[[lvl]]
    map <- mapping[[lvl]]
    pos <- match(codes, map$code, nomatch = 0L)
    keep <- pos > 0L & map$count[pos] > 0L
    start <- map$start[pos[keep]]
    cnt <- map$count[pos[keep]]
    idx <- unlist(mapply(function(s, c) seq.int(s + 1L, s + c), start, cnt,
                         SIMPLIFY = FALSE), use.names = FALSE)
    indices_by_level[[lvl]] <- idx
  }
  indices_by_level
}

# -- Public API --------------------------------------------------------------

#' Forward Haar wavelet transform (mask-adaptive, Morton order)
#'
#' @param X Numeric matrix with time in rows and voxels (masked) in columns.
#' @param mask 3D logical array or LogicalNeuroVol defining the mask.
#' @param levels Number of decomposition levels; defaults to ceil(log2(max(dim(mask)))).
#' @param z_seed Integer seed for Morton ordering (tie-breaking).
#' @param threshold List with fields `type` ("none" | "absolute" | "relative_to_root_std") and `value`.
#' @return List with elements `coeff` (list(root, detail)), `meta` (counts, Morton hash, etc.).
#' @export
haar_wavelet_forward <- function(X, mask, levels = NULL, z_seed = 42L,
                                 threshold = list(type = "none", value = 0)) {
  mask_arr <- as_logical_mask(mask, "haar_wavelet_forward")
  if (is.null(levels)) {
    levels <- ceiling(log2(max(dim(mask_arr))))
  }
  levels <- as.integer(levels)
  coeffs <- perform_haar_lift_analysis(X, mask_arr, levels, z_seed)

  thresh_type <- threshold$type %||% "none"
  if (!identical(thresh_type, "none")) {
    thresh_val <- as.numeric(threshold$value %||% 0)
    if (is.na(thresh_val)) thresh_val <- 0
    actual_threshold <- switch(thresh_type,
      absolute = thresh_val,
      relative_to_root_std = thresh_val * stats::sd(c(as.vector(coeffs$root),
                                                      as.vector(coeffs$detail[[levels]]))),
      0
    )
    if (actual_threshold > 0) {
      for (lv in seq_len(levels)) {
        m <- coeffs$detail[[lv]]
        m[abs(m) < actual_threshold] <- 0
        coeffs$detail[[lv]] <- m
      }
    }
  }

  vox_coords <- which(mask_arr, arr.ind = TRUE)
  bbox <- c(min(vox_coords[, 1]) - 1L, max(vox_coords[, 1]) - 1L,
            min(vox_coords[, 2]) - 1L, max(vox_coords[, 2]) - 1L,
            min(vox_coords[, 3]) - 1L, max(vox_coords[, 3]) - 1L)
  morton_idx <- get_morton_ordered_indices(mask_arr, z_seed)
  morton_hash <- paste0("sha1:", digest::digest(morton_idx, algo = "sha1", serialize = TRUE))
  scalings <- precompute_haar_scalings(mask_arr, levels)
  counts_by_level <- lapply(scalings, function(sc) as.integer(round(sc$sqrt_nvalid^2)))
  lowpass_counts <- integer(levels + 1L)
  detail_counts <- integer(levels)
  for (lv in seq_len(levels)) {
    counts <- counts_by_level[[lv]]
    lowpass_counts[lv] <- length(counts)
    detail_counts[lv] <- sum(counts)
  }
  lowpass_counts[levels + 1L] <- length(counts_by_level[[levels]])

  valid_blocks <- get_valid_finest_blocks(mask_arr)

  meta <- list(
    levels = levels,
    z_seed = z_seed,
    mask_dims = dim(mask_arr),
    morton_hash_mask_indices = morton_hash,
    num_voxels_in_mask = length(morton_idx),
    octree_bounding_box_mask_space = bbox,
    num_coeffs_per_level = list(lowpass = lowpass_counts, detail = detail_counts),
    valid_finest_blocks = valid_blocks,
    morton_indices = morton_idx
  )

  list(coeff = coeffs, meta = meta)
}

#' Inverse Haar wavelet transform
#'
#' @param coeff Output from `haar_wavelet_forward()` (list with root/detail matrices).
#' @param mask 3D logical array or LogicalNeuroVol.
#' @param levels Integer levels (defaults to `meta$levels` if present).
#' @param z_seed Integer Morton seed (defaults to `meta$z_seed`).
#' @param roi_mask Optional ROI mask (same dims as mask) to subset voxels before returning.
#' @param time_idx Optional integer vector of time indices to subset rows.
#' @param levels_keep Optional integer vector of detail levels to include (1 = finest, `levels` = coarsest).
#'   Levels not listed are zeroed, enabling coarse/partial reconstruction without full detail.
#' @return Matrix (time x voxels) reconstructed (or subsetted if ROI/time specified).
#' @export
haar_wavelet_inverse <- function(coeff, mask, levels = NULL, z_seed = NULL,
                                 roi_mask = NULL, time_idx = NULL,
                                 levels_keep = NULL) {
  mask_arr <- as_logical_mask(mask, "haar_wavelet_inverse")
  meta <- coeff$meta %||% list()
  levels <- as.integer(levels %||% meta$levels %||% stop("levels must be provided", call. = FALSE))
  z_seed <- as.integer(z_seed %||% meta$z_seed %||% 42L)

  root <- coeff$coeff$root %||% coeff$root
  detail <- coeff$coeff$detail %||% coeff$detail
  coeff_list <- list(root = root, detail = detail)

  if (!is.null(time_idx)) {
    coeff_list$root <- coeff_list$root[time_idx, , drop = FALSE]
    coeff_list$detail <- lapply(coeff_list$detail, function(m) m[time_idx, , drop = FALSE])
  }

  if (!is.null(levels_keep)) {
    levels_keep <- intersect(as.integer(levels_keep), seq_len(levels))
    scalings <- precompute_haar_scalings(mask_arr, levels)
    detail_expected <- vapply(
      scalings,
      function(sc) sum(as.integer(round(sc$sqrt_nvalid^2))),
      integer(1)
    )
    coeff_list$detail <- lapply(seq_len(levels), function(lv) {
      m <- coeff_list$detail[[lv]]
      ncol_expected <- detail_expected[lv]
      n_time <- nrow(coeff_list$root)
      if (is.null(m) || !length(m)) {
        m <- matrix(0, nrow = n_time, ncol = ncol_expected)
      } else if (ncol(m) != ncol_expected) {
        # resize with zeros to expected width
        m0 <- matrix(0, nrow = n_time, ncol = ncol_expected)
        m0[, seq_len(min(ncol(m), ncol_expected))] <- m[, seq_len(min(ncol(m), ncol_expected)), drop = FALSE]
        m <- m0
      }
      if (lv %in% levels_keep) m else {
        m[] <- 0
        m
      }
    })
  }

  reco <- perform_haar_lift_synthesis(coeff_list, mask_arr, levels, z_seed)

  if (!is.null(roi_mask)) {
    # roi_mask is in full volume order; map to mask order indices
    global_idx <- which(as.logical(mask_arr))
    roi_global <- which(as.logical(roi_mask))
    col_keep <- which(global_idx %in% roi_global)
    reco <- reco[, col_keep, drop = FALSE]
  }
  reco
}

#' Build Haar latent representation
#'
#' @param X Numeric matrix (time x voxels within mask).
#' @param mask LogicalNeuroVol or 3D logical array.
#' @param levels Decomposition levels (optional).
#' @param z_seed Integer Morton seed (default 42).
#' @param threshold Threshold list passed to `haar_wavelet_forward`.
#' @return An object of class `HaarLatent` containing coefficients and metadata.
#' @export
haar_latent <- function(X, mask, levels = NULL, z_seed = 42L,
                        threshold = list(type = "none", value = 0)) {
  mask_arr <- tryCatch(extract_mask_array_logical(mask, "haar_latent"),
                       error = function(e) as_logical_mask(mask, "haar_latent"))
  fw <- haar_wavelet_forward(X, mask_arr, levels = levels, z_seed = z_seed,
                             threshold = threshold)
  meta <- fw$meta
  meta$family <- "haar"
  decoder <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL) {
    haar_wavelet_inverse(list(coeff = fw$coeff, meta = meta), mask_arr,
                         levels = meta$levels, z_seed = meta$z_seed,
                         roi_mask = roi_mask, time_idx = time_idx,
                         levels_keep = levels_keep)
  }
  obj <- implicit_latent(fw$coeff, decoder, meta, mask_arr)
  class(obj) <- c("HaarLatent", class(obj))
  obj
}

#' Convert HaarLatent to matrix
#' @param x HaarLatent object
#' @param ... Additional arguments (unused)
#' @return Matrix of reconstructed values
#' @method as.matrix HaarLatent
#' @export
as.matrix.HaarLatent <- function(x, ...) {
  haar_wavelet_inverse(x, x$mask, levels = x$meta$levels, z_seed = x$meta$z_seed)
}

#' Predict method for HaarLatent
#' @param object HaarLatent object
#' @param roi_mask Optional ROI mask
#' @param time_idx Optional time indices
#' @param levels_keep Optional levels to keep
#' @param ... Additional arguments passed to decoder
#' @return Matrix of predicted values
#' @export
predict.HaarLatent <- function(object, roi_mask = NULL, time_idx = NULL,
                               levels_keep = NULL, ...) {
  predict.ImplicitLatent(object, roi_mask = roi_mask, time_idx = time_idx,
                         levels_keep = levels_keep, ...)
}
