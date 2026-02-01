# Hierarchical template and encoding helpers

#' Check whether an object is a HierarchicalBasisTemplate
#' @param x object to test
#' @export
is_hierarchical_template <- function(x) inherits(x, "HierarchicalBasisTemplate")

#' Build a hierarchical Laplacian template (offline)
#'
#' Constructs a multi-level spatial basis from nested parcellations using
#' graph Laplacian eigenvectors. The resulting template can be reused to
#' efficiently encode multiple datasets that share the same mask geometry.
#'
#' @param mask LogicalNeuroVol or logical array (3D) defining the domain.
#' @param parcellations List of integer vectors (one per level) of length = #voxels in mask.
#'   Levels must be nested: each child parcel maps to exactly one parent parcel in the
#'   previous level.
#' @param k_per_level Integer vector giving #modes per parcel at each level.
#' @param k_neighbors k for local graph construction inside parcels.
#' @param ridge Small diagonal ridge added to \eqn{G = t(B) \%*\% B} for stability.
#' @param solver Solver choice: "chol" (default) or "qr" fallback.
#' @param label Optional label stored in meta.
#' @return HierarchicalBasisTemplate (primal basis B + cached solver and metadata).
#' @export
build_hierarchical_template <- function(mask, parcellations, k_per_level,
                                        k_neighbors = 6L, ridge = 1e-8,
                                        solver = c("chol", "qr"),
                                        label = "hierarchical_laplacian") {
  solver <- match.arg(solver)

  mask_arr <- .mask_to_array(mask, "build_hierarchical_template")
  mask_vol <- if (inherits(mask, "LogicalNeuroVol")) mask else LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dim(mask_arr)))

  n_vox <- sum(mask_arr)

  if (!is.list(parcellations) || length(parcellations) == 0L) {
    stop("parcellations must be a non-empty list of integer vectors")
  }
  if (length(parcellations) != length(k_per_level)) {
    stop("length(k_per_level) must match number of parcellation levels")
  }

  # Normalize parcellations to integer vectors
  parcellations <- lapply(parcellations, function(v) as.integer(v))
  for (lvl in seq_along(parcellations)) {
    if (length(parcellations[[lvl]]) != n_vox) {
      stop("parcellation level ", lvl, " has length ", length(parcellations[[lvl]]), " but mask has ", n_vox, " voxels")
    }
  }

  # Check nesting and derive parent maps
  parents <- vector("list", length(parcellations))
  parents[[1L]] <- integer(0)
  if (length(parcellations) > 1L) {
    for (lvl in 2:length(parcellations)) {
      parents[[lvl]] <- .infer_parent_map(child = parcellations[[lvl]], parent = parcellations[[lvl - 1L]], lvl = lvl)
    }
  }

  # Build blocks per level
  block_list <- list()
  atom_rows <- list()
  col_offset <- 0L

  for (lvl in seq_along(parcellations)) {
    k_lvl <- as.integer(k_per_level[[lvl]])
    block_info <- .build_level_block(
      mask = mask_vol,
      labels = parcellations[[lvl]],
      k = k_lvl,
      k_neighbors = k_neighbors,
      level = lvl,
      parent_map = parents[[lvl]],
      col_offset = col_offset
    )

    block_list[[length(block_list) + 1L]] <- block_info$block
    atom_rows[[length(atom_rows) + 1L]] <- block_info$atoms
    col_offset <- block_info$next_col
  }

  if (col_offset == 0L) stop("No atoms were generated; check parcellations and k_per_level")

  B <- do.call(cbind, block_list)

  G <- Matrix::crossprod(B)
  if (ridge > 0) {
    G <- G + ridge * Matrix::Diagonal(n = ncol(G))
  }

  G_factor <- NULL
  if (solver == "chol") {
    G_factor <- try(Matrix::Cholesky(G, perm = TRUE, LDL = TRUE), silent = TRUE)
  }
  if (inherits(G_factor, "try-error") || is.null(G_factor)) {
    if (solver == "chol") message("Cholesky failed; falling back to sparse QR")
    G_factor <- qr(G)
    solver <- "qr"
  }

  atom_df <- do.call(rbind, atom_rows)

  spc <- neuroim2::NeuroSpace(c(dim(mask_arr), 1L))

  meta <- list(
    family = "hierarchical_laplacian",
    label = label,
    k_per_level = k_per_level,
    ridge = ridge,
    solver = solver
  )

  new("HierarchicalBasisTemplate",
    mask = mask_vol,
    space = spc,
    levels = parcellations,
    parents = parents,
    loadings = B,
    gram_factor = G_factor,
    atoms = atom_df,
    meta = meta
  )
}

#' Encode data using a hierarchical template
#' @param X matrix time x voxels (mask order) matching template mask
#' @param template HierarchicalBasisTemplate
#' @param label Optional label for the resulting LatentNeuroVec (defaults to template label)
#' @return LatentNeuroVec with basis = time x atoms coefficients, loadings = template loadings
#' @export
encode_hierarchical <- function(X, template, label = NULL) {
  if (!is_hierarchical_template(template)) stop("template must be a HierarchicalBasisTemplate")
  X_mat <- Matrix::Matrix(X, sparse = FALSE)
  n_time <- nrow(X_mat)
  n_vox <- ncol(X_mat)

  mask_arr <- as.array(template@mask)
  if (n_vox != sum(mask_arr)) {
    stop("X has ", n_vox, " voxels, but template mask has ", sum(mask_arr))
  }

  B <- template@loadings
  G_factor <- template@gram_factor

  proj <- Matrix::crossprod(B, Matrix::t(X_mat))  # atoms x time

  coeff <- if (inherits(G_factor, "sparseQR")) {
    Matrix::solve(G_factor, proj)
  } else {
    Matrix::solve(G_factor, proj)
  }

  coeff_t <- Matrix::t(coeff)  # time x atoms

  spc <- neuroim2::NeuroSpace(c(dim(mask_arr), n_time))
  lbl <- label %||% template@meta$label %||% "hierarchical_latent"
  meta <- list(family = "hierarchical_laplacian", template = template@meta)

  LatentNeuroVec(
    basis = coeff_t,
    loadings = template@loadings,
    space = spc,
    mask = template@mask,
    offset = numeric(0),
    label = lbl,
    meta = meta
  )
}

#' Project coefficients only (no LatentNeuroVec wrapper)
#' @param template HierarchicalBasisTemplate to use for projection
#' @param X matrix time x voxels (mask order) matching template mask
#' @return Matrix of coefficients (time x atoms)
#' @export
project_hierarchical <- function(template, X) {
  if (!is_hierarchical_template(template)) stop("template must be a HierarchicalBasisTemplate")
  X_mat <- Matrix::Matrix(X, sparse = FALSE)
  mask_arr <- as.array(template@mask)
  if (ncol(X_mat) != sum(mask_arr)) stop("X voxel dimension does not match template mask")

  B <- template@loadings
  G_factor <- template@gram_factor
  proj <- Matrix::crossprod(B, Matrix::t(X_mat))
  Matrix::t(Matrix::solve(G_factor, proj))  # time x atoms
}

#' Save a hierarchical template to disk
#' @param template HierarchicalBasisTemplate
#' @param file Path to .rds file
#' @param compress Compression passed to saveRDS (default "xz")
#' @export
save_hierarchical_template <- function(template, file, compress = "xz") {
  if (!is_hierarchical_template(template)) stop("template must be a HierarchicalBasisTemplate")
  saveRDS(template, file = file, compress = compress)
  invisible(normalizePath(file, winslash = "/", mustWork = FALSE))
}

#' Load a hierarchical template from disk
#' @param file Path to .rds produced by save_hierarchical_template
#' @return HierarchicalBasisTemplate
#' @export
load_hierarchical_template <- function(file) {
  obj <- readRDS(file)
  if (!is_hierarchical_template(obj)) stop("RDS does not contain a HierarchicalBasisTemplate")
  obj
}

# ---- helpers ---------------------------------------------------------------

.infer_parent_map <- function(child, parent, lvl) {
  child_ids <- sort(unique(child))
  parent_map <- integer(length(child_ids))
  names(parent_map) <- as.character(child_ids)
  for (cid in child_ids) {
    vox_idx <- which(child == cid)
    if (length(vox_idx) == 0L) next
    p_ids <- unique(parent[vox_idx])
    p_ids <- p_ids[!is.na(p_ids)]
    if (length(p_ids) != 1L) {
      stop("Parcellations are not nested at level ", lvl, ": child ", cid, " maps to parents ", paste(p_ids, collapse = ","))
    }
    parent_map[as.character(cid)] <- p_ids
  }
  parent_map
}

.build_level_block <- function(mask, labels, k, k_neighbors, level, parent_map, col_offset) {
  ids <- sort(unique(labels))
  n_vox <- length(labels)

  i_list <- list(); j_list <- list(); x_list <- list()
  atom_rows <- list()
  current_col <- col_offset

  for (cid in ids) {
    vox_idx <- which(labels == cid)
    n_loc <- length(vox_idx)
    if (n_loc == 0L) next

    k_use <- min(k, n_loc)
    vecs <- NULL
    if (n_loc == 1L) {
      vecs <- Matrix::sparseMatrix(i = 1L, j = 1L, x = 1, dims = c(1, 1))
    } else if (n_loc < 3L) {
      # RSpectra requires dimension >= 3; use base eigen() for small parcels
      g <- voxel_subset_to_gsp(mask, vox_idx, k_neighbors = k_neighbors)
      L <- as.matrix(g$laplacian)
      eig <- eigen(L, symmetric = TRUE)
      # Select k_use smallest eigenvalue eigenvectors
      ord <- order(eig$values)
      vecs <- Matrix::Matrix(eig$vectors[, ord[seq_len(k_use)], drop = FALSE], sparse = TRUE)
    } else {
      if (!requireNamespace("RSpectra", quietly = TRUE)) {
        stop("RSpectra is required to build hierarchical templates", call. = FALSE)
      }
      g <- voxel_subset_to_gsp(mask, vox_idx, k_neighbors = k_neighbors)
      L <- g$laplacian
      eig <- RSpectra::eigs(L, k = k_use, which = "SM")
      vecs <- Matrix::Matrix(eig$vectors, sparse = TRUE)
    }

    trip <- Matrix::summary(as(vecs, "TsparseMatrix"))
    if (nrow(trip) > 0) {
      i_list[[length(i_list) + 1L]] <- vox_idx[trip$i]
      j_list[[length(j_list) + 1L]] <- current_col + trip$j
      x_list[[length(x_list) + 1L]] <- trip$x
    }

    # atom metadata
    parent_id <- if (length(parent_map) > 0) parent_map[as.character(cid)] else NA_integer_
    atom_rows[[length(atom_rows) + 1L]] <- data.frame(
      col_id = seq.int(current_col + 1L, current_col + k_use),
      level = level,
      parcel_id = cid,
      parent_id = parent_id,
      mode = seq_len(k_use),
      label = paste0("L", level, ":P", cid)
    )

    current_col <- current_col + ncol(vecs)
  }

  if (current_col == col_offset) {
    block <- Matrix::Matrix(0, nrow = n_vox, ncol = 0, sparse = TRUE)
  } else {
    block <- Matrix::sparseMatrix(
      i = unlist(i_list, use.names = FALSE),
      j = unlist(j_list, use.names = FALSE),
      x = unlist(x_list, use.names = FALSE),
      dims = c(n_vox, current_col)
    )
    block <- block[, seq.int(col_offset + 1L, current_col), drop = FALSE]
  }

  list(block = block, atoms = do.call(rbind, atom_rows), next_col = current_col)
}
