#' Hierarchical radial basis functions (HRBF) for latent fMRI
#'
#' Utilities to generate analytic HRBF bases, project data, reconstruct,
#' and build `LatentNeuroVec` objects. Parameters are kept simple and
#' in-R only (no descriptors or HDF5).
#'
#' @param params List with fields:
#'   - `sigma0` (numeric, default 6)
#'   - `levels` (integer, default 3)
#'   - `radius_factor` (numeric, default 2.5)
#'   - `num_extra_fine_levels` (integer, default 0)
#'   - `kernel_type` (\"gaussian\" or \"wendland_c6\", alias \"wendland_c4\")
#'   - `seed` (integer) for deterministic Poisson sampling
#' @param mask `LogicalNeuroVol` mask defining voxel locations.
#' @return For `hrbf_generate_basis`, a sparse matrix with one row per
#'   HRBF atom and columns matching mask voxels.
#' @export
hrbf_generate_basis <- function(params, mask) {
  params_clean <- params
  params_clean$centres_stored <- FALSE
  if (is.null(params_clean$seed)) params_clean$seed <- 1L
  hrbf_basis_from_params(params_clean, mask)
}

#' @rdname hrbf_generate_basis
#' @param X Numeric matrix with time in rows and voxels in columns.
#' @export
hrbf_project_matrix <- function(X, mask, params) {
  B <- hrbf_generate_basis(params, mask)
  X_dense <- as.matrix(X)
  as.matrix(Matrix::tcrossprod(X_dense, B))
}

#' @rdname hrbf_generate_basis
#' @param coeff Coefficient matrix with rows = time points.
#' @export
hrbf_reconstruct_matrix <- function(coeff, mask, params) {
  B <- hrbf_generate_basis(params, mask)
  coeff_dense <- as.matrix(coeff)
  as.matrix(coeff_dense %*% B)
}

#' Partially reconstruct selected voxels/timepoints
#'
#' Computes only the requested voxels (and optional timepoints) without
#' materializing the full dense volume. Useful for fast indexed access.
#'
#' @param coeff Coefficient matrix (time x atoms).
#' @param mask `LogicalNeuroVol`.
#' @param params HRBF parameter list.
#' @param voxel_idx Integer vector of voxel linear indices (within the mask grid).
#' @param time_idx Optional integer vector of time indices; defaults to all rows.
#' @return Matrix of shape length(time_idx) x length(voxel_idx).
#' @export
hrbf_reconstruct_partial <- function(coeff, mask, params,
                                     voxel_idx, time_idx = NULL) {
  if (missing(voxel_idx) || length(voxel_idx) == 0) {
    stop("voxel_idx must be provided")
  }
  mask_arr <- extract_mask_array_logical(mask, "hrbf_reconstruct_partial")
  n_vox_total <- length(mask_arr)
  if (any(voxel_idx < 1 | voxel_idx > n_vox_total)) {
    stop("voxel_idx out of bounds")
  }
  time_idx <- time_idx %||% seq_len(nrow(coeff))
  if (any(time_idx < 1 | time_idx > nrow(coeff))) {
    stop("time_idx out of bounds")
  }
  B <- hrbf_generate_basis(params, mask) # atoms x voxels
  B_sub <- B[, voxel_idx, drop = FALSE]
  coeff_sub <- coeff[time_idx, , drop = FALSE]
  out <- coeff_sub %*% B_sub
  as.matrix(out)
}

#' Build a LatentNeuroVec using an HRBF basis
#'
#' @param X Numeric matrix (time x voxels within mask).
#' @param mask `LogicalNeuroVol`.
#' @param params HRBF parameter list (see `hrbf_generate_basis`).
#' @param label Optional character label.
#' @return A `LatentNeuroVec` with `basis = coefficients`,
#'   `loadings = t(HRBF_basis)`.
#' @export
hrbf_latent <- function(X, mask, params = list(), label = "") {
  mask_arr <- extract_mask_array_logical(mask, "hrbf_latent")
  n_time <- nrow(X)
  if (is.null(n_time) || n_time < 1) stop("X must have time in rows")
  B_atoms_vox <- hrbf_generate_basis(params, mask)
  coeff <- hrbf_project_matrix(X, mask, params)
  loadings <- Matrix::t(B_atoms_vox) # voxels x atoms
  spc <- neuroim2::NeuroSpace(c(dim(mask_arr), n_time))
  meta <- list(family = "hrbf", params = params)
  LatentNeuroVec(basis = coeff, loadings = loadings, space = spc,
                 mask = mask, label = label, meta = meta)
}

#' Attach HRBF metadata to an existing LatentNeuroVec
#' @param lvec A LatentNeuroVec object
#' @param params HRBF parameters list
#' @param centres Optional matrix of HRBF centres
#' @param sigmas Optional vector of HRBF sigmas
#' @return The LatentNeuroVec with HRBF metadata attached
#' @export
as_hrbf_latent <- function(lvec, params, centres = NULL, sigmas = NULL) {
  stopifnot(inherits(lvec, "LatentNeuroVec"))
  meta <- list(family = "hrbf", params = params, centres = centres, sigmas = sigmas)
  lvec@meta <- meta
  lvec
}

#' Retrieve HRBF metadata if present
#' @param x A LatentNeuroVec object
#' @return HRBF metadata list or NULL if not HRBF-tagged
#' @export
hrbf_meta <- function(x) {
  if (!inherits(x, "LatentNeuroVec")) return(NULL)
  m <- x@meta %||% NULL
  if (is.null(m)) return(NULL)
  if (!identical(m$family, "hrbf")) return(NULL)
  m
}

#' Check if latent object carries HRBF metadata
#' @param x Object to test
#' @return Logical indicating if x has HRBF metadata
#' @export
is_hrbf_latent <- function(x) {
  !is.null(hrbf_meta(x))
}

# ---- Internal helpers ----------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

use_hrbf_rcpp <- function() {
  isTRUE(getOption("fmrilatent.hrbf.use_rcpp", FALSE))
}

extract_mask_array_logical <- function(mask_neurovol, location = "mask") {
  mask_arr <- tryCatch(as.array(mask_neurovol), error = function(e) NULL)
  if (is.null(mask_arr) || is.list(mask_arr)) {
    if (is.list(mask_neurovol) && !is.null(mask_neurovol$arr)) {
      mask_arr <- mask_neurovol$arr
    } else {
      stop(sprintf("Cannot extract array from mask (%s)", location))
    }
  }
  if (is.null(dim(mask_arr))) {
    stop(sprintf("mask must have dimensions (%s)", location))
  }
  array(as.logical(mask_arr), dim = dim(mask_arr))
}

extract_spacing_origin <- function(mask_neurovol, ndims, location) {
  spacing_default <- rep(1, ndims)
  origin_default <- rep(0, ndims)
  sp <- tryCatch(neuroim2::space(mask_neurovol), error = function(e) NULL)
  spacing_vec <- if (is.null(sp)) spacing_default else tryCatch(neuroim2::spacing(sp), error = function(e) spacing_default)
  origin_vec <- if (is.null(sp)) origin_default else tryCatch(neuroim2::origin(sp), error = function(e) origin_default)
  list(spacing = spacing_vec, origin = origin_vec)
}

label_components <- function(mask_arr_3d) {
  dims <- dim(mask_arr_3d)
  visited <- array(FALSE, dim = dims)
  labels <- array(0L, dim = dims)
  comp_id <- 0L
  neighbours <- matrix(c(1,0,0,-1,0,0,0,1,0,0,-1,0,0,0,1,0,0,-1),
                       ncol = 3, byrow = TRUE)
  for (i in seq_len(dims[1])) for (j in seq_len(dims[2])) for (k in seq_len(dims[3])) {
    if (isTRUE(mask_arr_3d[i,j,k]) && !visited[i,j,k]) {
      comp_id <- comp_id + 1L
      q <- list(c(i,j,k))
      while (length(q) > 0) {
        pt <- q[[1]]; q <- q[-1]
        ii <- pt[1]; jj <- pt[2]; kk <- pt[3]
        if (visited[ii,jj,kk]) next
        visited[ii,jj,kk] <- TRUE
        labels[ii,jj,kk] <- comp_id
        for (n in seq_len(nrow(neighbours))) {
          nn <- pt + neighbours[n, ]
          ni <- nn[1]; nj <- nn[2]; nk <- nn[3]
          if (ni >= 1 && ni <= dims[1] &&
              nj >= 1 && nj <= dims[2] &&
              nk >= 1 && nk <= dims[3] &&
              isTRUE(mask_arr_3d[ni,nj,nk]) && !visited[ni,nj,nk]) {
            q[[length(q) + 1L]] <- c(ni,nj,nk)
          }
        }
      }
    }
  }
  list(count = comp_id, labels = labels)
}

poisson_disk_sample_neuroim2 <- function(mask_neurovol, radius_mm, seed, component_id_for_seed_offset = 0) {
  mask_arr <- extract_mask_array_logical(mask_neurovol, "poisson_disk_sample_neuroim2")
  so <- extract_spacing_origin(mask_neurovol, length(dim(mask_arr)), "poisson_disk_sample_neuroim2")
  spacing_vec <- so$spacing
  comp_info <- label_components(mask_arr)
  radius_vox <- radius_mm / mean(spacing_vec)
  r2 <- radius_vox^2

  sample_component <- function(coords, base_seed) {
    set.seed(as.integer(base_seed))
    remaining <- coords
    selected <- matrix(numeric(0), ncol = 3)
    while (nrow(remaining) > 0) {
      cand <- remaining[1, , drop = FALSE]
      remaining <- remaining[-1, , drop = FALSE]
      if (nrow(selected) == 0) {
        selected <- rbind(selected, cand)
      } else {
        d2 <- rowSums((selected - matrix(cand, nrow = nrow(selected), ncol = 3, byrow = TRUE))^2)
        if (all(d2 >= r2)) selected <- rbind(selected, cand)
      }
    }
    if (nrow(coords) <= 8L && nrow(selected) == 0) {
      selected <- matrix(round(colMeans(coords)), nrow = 1)
    }
    selected
  }

  gather <- function(coords, seed_val) {
    out <- sample_component(coords, seed_val)
    if (nrow(coords) <= 8L) {
      centroid_pt <- round(colMeans(coords))
      out <- matrix(centroid_pt, nrow = 1)
    }
    out
  }

  if (comp_info$count > 1L && component_id_for_seed_offset == 0) {
    centres <- lapply(seq_len(comp_info$count), function(id) {
      coords <- which(comp_info$labels == id, arr.ind = TRUE)
      if (nrow(coords) == 0) return(matrix(numeric(0), ncol = 3))
      gather(coords, seed + id)
    })
    res <- do.call(rbind, centres)
    res <- if (nrow(res) > 0) matrix(as.numeric(res), ncol = 3) else matrix(numeric(0), ncol = 3)
    colnames(res) <- c("i","j","k")
    return(res)
  }

  vox_coords <- which(mask_arr, arr.ind = TRUE)
  if (nrow(vox_coords) == 0) {
    return(matrix(integer(0), ncol = 3, dimnames = list(NULL, c("i","j","k"))))
  }

  remaining <- vox_coords  # deterministic order
  selected <- matrix(numeric(0), ncol = 3)
  while (nrow(remaining) > 0) {
    cand <- remaining[1, , drop = FALSE]
    remaining <- remaining[-1, , drop = FALSE]
    if (nrow(selected) == 0) {
      selected <- rbind(selected, cand)
    } else {
      d2 <- rowSums((selected - matrix(cand, nrow = nrow(selected), ncol = 3, byrow = TRUE))^2)
      if (all(d2 >= r2)) {
        selected <- rbind(selected, cand)
      }
    }
  }
  selected <- if (nrow(selected) > 0) matrix(as.numeric(selected), ncol = 3) else selected
  colnames(selected) <- c("i","j","k")
  selected
}

generate_hrbf_atom <- function(mask_coords_world, mask_linear_indices,
                               centre_coord_world, sigma_mm,
                               current_level_j, total_levels, params,
                               normalize_over_mask = TRUE) {
  mask_coords_world <- as.matrix(mask_coords_world)
  centre_coord_world <- as.numeric(centre_coord_world)
  p_kernel_type <- params$kernel_type %||% "gaussian"
  if (identical(p_kernel_type, "wendland_c4")) p_kernel_type <- "wendland_c6"
  p_kernel_type_fine <- params$kernel_type_fine_levels %||% "wendland_c6"
  if (identical(p_kernel_type_fine, "wendland_c4")) p_kernel_type_fine <- "wendland_c6"
  num_alt <- params$num_fine_levels_alt_kernel %||% 0L
  total_levels_effective <- (params$levels %||% total_levels) + (params$num_extra_fine_levels %||% 0L)
  use_alt_kernel <- (num_alt > 0L) && (current_level_j >= (total_levels_effective - num_alt + 1L))
  eff_kernel <- if (use_alt_kernel) p_kernel_type_fine else p_kernel_type

  diffs <- sweep(mask_coords_world, 2, centre_coord_world, FUN = "-")
  dist_mm <- sqrt(rowSums(diffs^2))

  if (eff_kernel == "gaussian") {
    phi <- exp(-(dist_mm^2) / (2 * sigma_mm^2))
  } else {
    r <- dist_mm / sigma_mm
    base <- pmax(0, 1 - r)
    phi <- base^8 * (32 * r^3 + 25 * r^2 + 8 * r + 1)
    phi[r >= 1] <- 0
  }

  if (normalize_over_mask) {
    norm_val <- sqrt(sum(phi^2))
    if (norm_val > 0) phi <- phi / norm_val
  }

  list(values = phi, indices = mask_linear_indices)
}

hrbf_basis_from_params <- function(params, mask_neurovol,
                                   mask_world_coords = NULL,
                                   mask_arr = NULL,
                                   mask_linear_indices = NULL) {
  sigma0 <- params$sigma0 %||% 6
  levels <- params$levels %||% 3L
  radius_factor <- params$radius_factor %||% 2.5
  kernel_type <- params$kernel_type %||% "gaussian"
  if (identical(kernel_type, "wendland_c4")) kernel_type <- "wendland_c6"
  num_extra_fine_levels <- params$num_extra_fine_levels %||% 0L
  seed <- params$seed

  voxel_to_world <- function(vox_mat) {
    so <- extract_spacing_origin(mask_neurovol, ncol(vox_mat),
                                 "hrbf_basis_from_params:mask_space")
    sweep(vox_mat - 1, 2, so$spacing, `*`) +
      matrix(so$origin, nrow(vox_mat), ncol(vox_mat), byrow = TRUE)
  }

  if (is.null(mask_arr)) mask_arr <- extract_mask_array_logical(mask_neurovol, "hrbf_basis_from_params")
  if (is.null(mask_world_coords)) {
    mask_coords_vox_early <- which(mask_arr, arr.ind = TRUE)
    mask_world_coords <- voxel_to_world(mask_coords_vox_early)
  }
  if (is.null(mask_linear_indices)) {
    mask_linear_indices <- as.integer(which(mask_arr))
  }
  mask_coords_world <- mask_world_coords
  n_total_vox <- length(mask_linear_indices)

  if (n_total_vox == 0) {
    return(Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                                dims = c(0, 0)))
  }

  # Tiny-ROI fallback: identity-like basis for exact reconstruction
  if (!is.null(seed) && n_total_vox <= 64L) {
    centres_all <- mask_coords_world
    small_sigma <- max(1e-3, sigma0 / 100)
    use_rcpp <- use_hrbf_rcpp() && exists("hrbf_atoms_rcpp")
    if (use_rcpp) {
      B_try <- tryCatch(
        hrbf_atoms_rcpp(
          as.matrix(mask_world_coords),
          as.matrix(centres_all),
          rep_len(small_sigma, n_total_vox),
          kernel_type,
          value_threshold = 1e-12
        ),
        error = function(e) NULL
      )
      if (!is.null(B_try)) {
        norms <- sqrt(Matrix::rowSums(B_try^2))
        norms[norms == 0] <- 1
        B_try <- Matrix::Diagonal(x = 1 / norms) %*% B_try
        return(B_try)
      }
    }
    triplet_i_list <- vector("list", n_total_vox)
    triplet_j_list <- vector("list", n_total_vox)
    triplet_x_list <- vector("list", n_total_vox)
    for (kk in seq_len(n_total_vox)) {
      atom <- generate_hrbf_atom(mask_world_coords, mask_linear_indices,
                                 centres_all[kk, ], small_sigma,
                                 0L, 0L, params)
      nz <- which(abs(atom$values) > 1e-12)
      if (length(nz) > 0) {
        triplet_i_list[[kk]] <- rep.int(kk, length(nz))
        triplet_j_list[[kk]] <- mask_linear_indices[nz]
        triplet_x_list[[kk]] <- atom$values[nz]
      } else {
        triplet_i_list[[kk]] <- integer(); triplet_j_list[[kk]] <- integer(); triplet_x_list[[kk]] <- numeric()
      }
    }
    i_idx <- unlist(triplet_i_list, use.names = FALSE)
    j_idx <- unlist(triplet_j_list, use.names = FALSE)
    x_val <- unlist(triplet_x_list, use.names = FALSE)
    return(Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_val,
                                dims = c(n_total_vox, n_total_vox)))
  }

  if (is.null(seed)) stop("hrbf_basis_from_params requires a seed for centre generation")

  centres_list <- list(); sigs <- numeric(); level_vec <- integer()
  for (j in seq_len(levels + 1L) - 1L) {
    sigma_j <- sigma0 / (2^j)
    r_j <- radius_factor * sigma_j
    vox_centres <- poisson_disk_sample_neuroim2(mask_neurovol, r_j, seed + j,
                                                component_id_for_seed_offset = j)
    if (nrow(vox_centres) > 0) {
      centres_list[[length(centres_list) + 1L]] <- voxel_to_world(vox_centres)
      n_new <- nrow(vox_centres)
      sigs <- c(sigs, rep(sigma_j, n_new))
      level_vec <- c(level_vec, rep(j, n_new))
    }
  }
  if (num_extra_fine_levels > 0L) {
    for (j_extra in seq_len(num_extra_fine_levels)) {
      j_level <- levels + j_extra
      sigma_new <- sigma0 / (2^(levels + j_extra))
      r_new <- radius_factor * sigma_new
      vox_centres <- poisson_disk_sample_neuroim2(
        mask_neurovol, r_new, seed + levels + j_extra,
        component_id_for_seed_offset = j_level
      )
      if (nrow(vox_centres) > 0) {
        centres_list[[length(centres_list) + 1L]] <- voxel_to_world(vox_centres)
        n_new <- nrow(vox_centres)
        sigs <- c(sigs, rep(sigma_new, n_new))
        level_vec <- c(level_vec, rep(j_level, n_new))
      }
    }
  }
  C_total <- if (length(centres_list) > 0) do.call(rbind, centres_list)
             else matrix(numeric(0), ncol = 3)
  sigma_vec <- sigs
  k_actual <- nrow(C_total)

  if (k_actual == 0) {
    return(Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                                dims = c(0, n_total_vox)))
  }

  use_rcpp <- use_hrbf_rcpp() && exists("hrbf_atoms_rcpp")
  if (use_rcpp) {
    B_try <- tryCatch(
      hrbf_atoms_rcpp(
        as.matrix(mask_coords_world),
        as.matrix(C_total),
        as.numeric(sigma_vec),
        kernel_type,
        value_threshold = 1e-8
      ),
      error = function(e) NULL
    )
    if (!is.null(B_try)) {
      norms <- sqrt(Matrix::rowSums(B_try^2))
      norms[norms == 0] <- 1
      B_try <- Matrix::Diagonal(x = 1 / norms) %*% B_try
      return(B_try)
    }
  }

  triplet_i_list <- vector("list", k_actual)
  triplet_j_list <- vector("list", k_actual)
  triplet_x_list <- vector("list", k_actual)
  value_threshold <- 1e-8
  for (kk in seq_len(k_actual)) {
    atom <- generate_hrbf_atom(mask_coords_world, mask_linear_indices,
                               C_total[kk, ], sigma_vec[kk],
                               level_vec[kk], levels, params)
    keep <- abs(atom$values) >= value_threshold
    if (any(keep)) {
      triplet_i_list[[kk]] <- rep.int(kk, sum(keep))
      triplet_j_list[[kk]] <- atom$indices[keep]
      triplet_x_list[[kk]] <- atom$values[keep]
    } else {
      triplet_i_list[[kk]] <- integer()
      triplet_j_list[[kk]] <- integer()
      triplet_x_list[[kk]] <- numeric()
    }
  }
  i_idx <- unlist(triplet_i_list, use.names = FALSE)
  j_idx <- unlist(triplet_j_list, use.names = FALSE)
  x_val <- unlist(triplet_x_list, use.names = FALSE)
  Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_val,
                       dims = c(k_actual, n_total_vox))
}

compute_hrbf_coefficients <- function(X_time_vox, basis_atoms_vox) {
  stopifnot(ncol(X_time_vox) == ncol(basis_atoms_vox))
  coeff <- NULL
  # attempt square invert
  if (nrow(basis_atoms_vox) == ncol(basis_atoms_vox) && nrow(basis_atoms_vox) > 0) {
    inv_try <- tryCatch(Matrix::solve(basis_atoms_vox, Matrix::Diagonal(ncol(basis_atoms_vox))),
                        error = function(e) NULL)
    if (!is.null(inv_try)) {
      coeff_prod <- X_time_vox %*% inv_try
      coeff <- as.matrix(coeff_prod)
    }
  }
  if (is.null(coeff)) {
    proj <- Matrix::tcrossprod(X_time_vox, basis_atoms_vox)       # T x K
    gram <- tryCatch(Matrix::tcrossprod(basis_atoms_vox), error = function(e) NULL) # K x K
    if (!is.null(gram)) {
      d <- diag(gram)
      lambda <- if (length(d) > 0) 1e-6 * (mean(as.numeric(d)) + 1) else 1e-6
      gram_reg <- tryCatch(gram + Matrix::Diagonal(n = nrow(gram), x = lambda), error = function(e) gram)
      gram_inv <- tryCatch(Matrix::solve(gram_reg), error = function(e) NULL)
      if (!is.null(gram_inv)) {
        coeff_ls <- proj %*% gram_inv
        coeff <- as.matrix(coeff_ls)
      }
    }
    if (is.null(coeff)) coeff <- as.matrix(proj)
  }
  coeff
}
