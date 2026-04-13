# Indexing methods for LatentNeuroVec
# matricized_access, linear_access, [[, [

#' @include latent_neurovector.R
NULL

# --- matricized_access Methods ---

#' Internal Matrix-Based Access for LatentNeuroVec
#'
#' @description
#' Internal method providing efficient matrix-based access to elements.
#'
#' @param x A \code{LatentNeuroVec}.
#' @param i Index specification.
#' @return Computed values.
#'
#' @keywords internal
#' @noRd
setMethod(
  f = "matricized_access",
  signature = signature(x = "LatentNeuroVec", i = "matrix"),
  definition = function(x, i) {
    if (!is.numeric(i) || ncol(i) != 2L) {
      stop("`i` must be a numeric matrix with 2 columns (time, spatial-index)")
    }

    ## -- 1. split and sanity-check the two index columns
    t_idx <- as.integer(i[, 1L]) # time rows in @basis
    s_idx <- as.integer(i[, 2L]) # spatial indices 1..X*Y*Z

    nt <- .latent_basis_dim(x@basis)[1L]
    nxy <- prod(dim(x)[1:3])

    if (any(t_idx < 1L | t_idx > nt)) {
      stop("time index out of bounds")
    }
    if (any(s_idx < 1L | s_idx > nxy)) {
      stop("spatial index out of bounds")
    }

    ## -- 2. map spatial -> mask rows (0 means 'outside the mask')
    v_idx <- lookup(x@map, s_idx) # 0 / 1..nVox

    inside <- v_idx > 0L
    out <- numeric(length(t_idx)) # zeros by default

    if (any(inside)) {
      ## -- 3. gather the relevant rows
      b1 <- basis_mat(x, i = t_idx[inside])
      b2 <- loadings_mat(x, i = v_idx[inside])

      ## -- 4. pair-wise dot product + offset
      if (is.vector(b1) || (inherits(b1, "Matrix") && prod(dim(b1)) == length(b1))) {
        b1 <- as.matrix(b1)
        if (is.vector(b1)) b1 <- matrix(b1, nrow = 1)
      }
      if (is.vector(b2) || (inherits(b2, "Matrix") && prod(dim(b2)) == length(b2))) {
        b2 <- as.matrix(b2)
        if (is.vector(b2)) b2 <- matrix(b2, nrow = 1)
      }
      if (inherits(b1, "Matrix")) b1 <- as.matrix(b1)
      if (inherits(b2, "Matrix")) b2 <- as.matrix(b2)
      dot_products <- rowSums(b1 * b2)
      dot_products[is.na(dot_products)] <- 0

      if (length(x@offset) > 0) {
        out[inside] <- dot_products + x@offset[v_idx[inside]]
      } else {
        out[inside] <- dot_products
      }
    }

    out
  }
)

#' @keywords internal
#' @noRd
setMethod(
  f = "matricized_access",
  signature = signature(x = "LatentNeuroVec", i = "integer"),
  definition = function(x, i) {
    stopifnot(
      "[matricized_access,LatentNeuroVec,integer] Number of components mismatch." =
        .latent_basis_dim(x@basis)[2L] == .latent_loadings_dim(x@loadings)[2L]
    )

    if (any(i < 1) || any(i > .latent_loadings_dim(x@loadings)[1L])) {
      stop("Index out of bounds for 'loadings'")
    }
    b1 <- basis_mat(x)
    b2 <- loadings_mat(x, i = as.integer(i))

    if (inherits(b1, "Matrix")) b1 <- as.matrix(b1)
    if (inherits(b2, "Matrix")) b2 <- as.matrix(b2)

    out <- b1 %*% t(b2)

    if (length(x@offset) > 0) {
      as.matrix(sweep(out, 2, x@offset[as.integer(i)], "+"))
    } else {
      as.matrix(out)
    }
  }
)

#' @keywords internal
#' @noRd
setMethod(
  f = "matricized_access",
  signature = signature(x = "LatentNeuroVec", i = "numeric"),
  definition = function(x, i) {
    matricized_access(x, as.integer(i))
  }
)

# --- linear_access Methods ---

#' Linear access to LatentNeuroVec elements
#'
#' @description
#' Access elements of a \code{LatentNeuroVec} using linear (1D) indices
#' into the 4D array representation.
#'
#' @param x A \code{LatentNeuroVec} object
#' @param i Numeric index vector
#' @return The reconstructed values at the specified indices
#'
#' @name linear_access-methods
NULL

#' @export
#' @rdname linear_access-methods
setMethod(
  f = "linear_access",
  signature = signature(x = "LatentNeuroVec", i = "numeric"),
  definition = function(x, i) {
    linear_access(x, as.integer(i))
  }
)

#' @export
#' @rdname linear_access-methods
setMethod(
  f = "linear_access",
  signature = signature(x = "LatentNeuroVec", i = "integer"),
  definition = function(x, i) {
    dims_full <- dim(x)
    nels_4d <- prod(dims_full)
    nels_3d <- prod(dims_full[1:3])
    n_time <- dims_full[4]

    if (!is.numeric(i) || any(is.na(i))) {
      stop("[linear_access,LatentNeuroVec] Index `i` must be numeric without NA values")
    }
    if (any(i < 1) || any(i > nels_4d)) {
      stop(paste0("[linear_access,LatentNeuroVec] Index out of bounds for 4D volume [1..", nels_4d, "]"))
    }

    # Convert 4D linear indices to 3D spatial index + time index
    time_idx <- ceiling(i / nels_3d)
    spatial_idx_3d <- i %% nels_3d
    spatial_idx_3d[spatial_idx_3d == 0] <- nels_3d

    # Map 3D spatial indices to mask indices
    rowmap <- lookup(x@map, spatial_idx_3d)

    # Identify unique needed mask indices and time indices
    unique_valid_mask_idx <- unique(rowmap[rowmap > 0])
    unique_time_idx <- unique(time_idx)

    # If no requested indices fall within the mask, return zeros
    if (length(unique_valid_mask_idx) == 0) {
      return(numeric(length(i)))
    }

    # Pre-calculate
    loadings_subset <- loadings_mat(x, i = unique_valid_mask_idx)
    basis_subset <- basis_mat(x, i = unique_time_idx)

    # Shortcut for single time point query
    if (length(unique_time_idx) == 1L) {
      computed_vals <- drop(basis_subset %*% t(loadings_subset))
      if (length(x@offset) > 0) {
        computed_vals <- computed_vals + x@offset[unique_valid_mask_idx]
      }
      ovals <- numeric(length(i))
      original_indices_for_computed <- match(rowmap, unique_valid_mask_idx)
      valid_output_positions <- !is.na(original_indices_for_computed)
      ovals[valid_output_positions] <- computed_vals[original_indices_for_computed[valid_output_positions]]
      return(ovals)
    }

    # Calculate the required data block for multiple time points
    data_block <- basis_subset %*% t(loadings_subset)

    # Add offsets if they exist
    if (length(x@offset) > 0) {
      data_block <- sweep(data_block, 2, x@offset[unique_valid_mask_idx], "+")
    }

    # Create the output vector
    ovals <- numeric(length(i))

    # Map results back
    time_map <- match(time_idx, unique_time_idx)
    mask_map <- match(rowmap, unique_valid_mask_idx)

    in_mask_selector <- which(rowmap > 0)
    row_indices_in_block <- time_map[in_mask_selector]
    col_indices_in_block <- mask_map[in_mask_selector]
    linear_indices_in_block <- row_indices_in_block + (col_indices_in_block - 1) * nrow(data_block)
    ovals[in_mask_selector] <- data_block[linear_indices_in_block]

    ovals
  }
)

# --- Extract Methods ---

#' Extract Elements from LatentNeuroVec
#'
#' @description
#' Extract elements from a \code{LatentNeuroVec}. Use \code{[[} to extract a single
#' volume as a \code{SparseNeuroVol}, or \code{[} for 4D subsetting.
#'
#' @param x A \code{\link{LatentNeuroVec-class}} object.
#' @param i Numeric index for first dimension (x-axis) or volume index for \code{[[}.
#' @param j Numeric index for second dimension (y-axis).
#' @param k Numeric index for third dimension (z-axis).
#' @param l Numeric index for fourth dimension (time).
#' @param ... Additional arguments (unused).
#' @param drop Logical, whether to drop dimensions (default TRUE).
#'
#' @return For \code{[[}: A \code{\link[neuroim2]{SparseNeuroVol-class}} containing the computed volume.
#'   For \code{[}: An array of extracted values.
#'
#' @examples
#' \dontrun{
#' # Extract volumes
#' vol1 <- lvec[[1]]
#' vol_mid <- lvec[[dim(lvec)[4] / 2]]
#' vol_last <- lvec[[dim(lvec)[4]]]
#' }
#'
#' @importFrom neuroim2 SparseNeuroVol
#' @rdname extract-methods
#' @export
setMethod(
  "[[", signature(x = "LatentNeuroVec", i = "numeric"),
  function(x, i) {
    if (length(i) != 1L || is.na(i) || !is.finite(i)) stop("Index must be a single finite number")
    if (!isTRUE(all.equal(i, as.integer(i)))) stop("Index must be integer-valued")
    i <- as.integer(i)
    if (i < 1L || i > dim(x)[4]) stop("Index out of range")

    b1 <- basis_mat(x, i = i)
    b2 <- loadings_mat(x)
    dat <- as.numeric(b1 %*% t(b2))
    if (length(x@offset) > 0) {
      dat <- dat + x@offset
    }

    newdim <- dim(x)[1:3]
    bspace <- NeuroSpace(newdim,
      spacing = neuroim2::spacing(x),
      origin = neuroim2::origin(x),
      axes = neuroim2::axes(x@space),
      trans = neuroim2::trans(x)
    )

    SparseNeuroVol(dat, bspace, indices = neuroim2::indices(x))
  }
)

#' @rdname extract-methods
#' @export
setMethod(
  f = "[",
  signature = signature(x = "LatentNeuroVec", i = "numeric", j = "numeric", drop = "ANY"),
  definition = function(x, i, j, k, l, ..., drop = TRUE) {
    dims_full <- dim(x)
    if (missing(i)) i <- seq_len(dims_full[1])
    if (missing(j)) j <- seq_len(dims_full[2])
    if (missing(k)) k <- seq_len(dims_full[3])
    if (missing(l)) l <- seq_len(dims_full[4])

    # Convert logical indices to positions
    if (is.logical(i)) i <- which(i)
    if (is.logical(j)) j <- which(j)
    if (is.logical(k)) k <- which(k)
    if (is.logical(l)) l <- which(l)

    if (any(i < 1 | i > dims_full[1]) || any(j < 1 | j > dims_full[2]) ||
        any(k < 1 | k > dims_full[3]) || any(l < 1 | l > dims_full[4])) {
      stop("Subscript out of range for LatentNeuroVec.")
    }

    i <- as.integer(i)
    j <- as.integer(j)
    k <- as.integer(k)
    l <- as.integer(l)

    out_dim <- c(length(i), length(j), length(k), length(l))
    n_vox_req <- prod(out_dim[1:3])
    n_time_req <- out_dim[4]

    dims_space_3d <- dims_full[1:3]
    linear_idx_3d <- integer(n_vox_req)

    ij_grid <- outer(i, j, function(ii, jj) ii + (jj - 1L) * dims_space_3d[1])
    ij_step <- length(i) * length(j)

    for (kk in seq_along(k)) {
      start <- (kk - 1L) * ij_step + 1L
      end <- kk * ij_step
      linear_idx_3d[start:end] <-
        ij_grid + (k[kk] - 1L) * dims_space_3d[1] * dims_space_3d[2]
    }

    rowmap <- lookup(x@map, linear_idx_3d)
    valid_mask_indices <- rowmap[rowmap > 0]
    map_req_to_valid <- which(rowmap > 0)

    result <- array(0, dim = out_dim)

    if (length(valid_mask_indices) == 0) {
      if (drop) {
        return(drop(result))
      }
      return(result)
    }

    basis_subset <- basis_mat(x, i = l)
    b1 <- basis_subset
    b2 <- loadings_mat(x, i = valid_mask_indices)
    result_valid <- b1 %*% t(b2)

    if (length(x@offset) > 0) {
      result_valid <- sweep(result_valid, 2, x@offset[valid_mask_indices], "+")
    }

    n_vox_per_slice <- prod(out_dim[1:3])
    for (col_idx in seq_along(valid_mask_indices)) {
      voxel_linear_idx_in_slice <- map_req_to_valid[col_idx]
      linear_indices_in_result <- voxel_linear_idx_in_slice + (seq_len(n_time_req) - 1L) * n_vox_per_slice
      result[linear_indices_in_result] <- result_valid[, col_idx]
    }

    if (drop) result <- drop(result)
    result
  }
)

#' @rdname extract-methods
#' @export
setMethod(
  f = "[",
  signature = signature(x = "LatentNeuroVec", i = "ANY", j = "ANY", drop = "ANY"),
  definition = function(x, i, j, k, l, ..., drop = TRUE) {
    dims <- dim(x)

    normalize_subscript <- function(idx, extent, name) {
      if (is.logical(idx)) {
        return(which(idx))
      }
      if (!is.numeric(idx)) {
        stop("Subscript '", name, "' must be numeric or logical for LatentNeuroVec.",
             call. = FALSE)
      }
      idx_int <- as.integer(idx)
      if (anyNA(idx_int) || any(idx_int < 1L | idx_int > extent)) {
        stop("Subscript out of range for LatentNeuroVec.", call. = FALSE)
      }
      idx_int
    }

    # Handle matrix indexing: each row is a coordinate tuple
    if (!missing(i) && is.matrix(i) && missing(j)) {
      nc <- ncol(i)
      if (nc == 4L) {
        linear_idx <- i[, 1] + (i[, 2] - 1L) * dims[1] +
          (i[, 3] - 1L) * dims[1] * dims[2] +
          (i[, 4] - 1L) * dims[1] * dims[2] * dims[3]
        return(linear_access(x, as.integer(linear_idx)))
      } else if (nc == 2L) {
        return(matricized_access(x, i))
      } else {
        stop("Matrix index must have 2 or 4 columns for LatentNeuroVec")
      }
    }

    if (missing(i)) i <- seq_len(dims[1])
    if (missing(j)) j <- seq_len(dims[2])
    if (missing(k)) k <- seq_len(dims[3])
    if (missing(l)) l <- seq_len(dims[4])

    i <- normalize_subscript(i, dims[1], "i")
    j <- normalize_subscript(j, dims[2], "j")
    k <- normalize_subscript(k, dims[3], "k")
    l <- normalize_subscript(l, dims[4], "l")

    out_dims <- c(length(i), length(j), length(k), length(l))
    result <- array(0, dim = out_dims)

    spatial_indices <- array(0, dim = out_dims[1:3])
    idx <- 1
    for (kk in seq_along(k)) {
      for (jj in seq_along(j)) {
        for (ii in seq_along(i)) {
          lin_idx <- i[ii] + (j[jj] - 1) * dims[1] + (k[kk] - 1) * dims[1] * dims[2]
          spatial_indices[ii, jj, kk] <- lin_idx
        }
      }
    }

    spatial_vec <- as.vector(spatial_indices)
    mask_indices <- lookup(x@map, spatial_vec)
    valid_mask_idx <- unique(mask_indices[mask_indices > 0])

    if (length(valid_mask_idx) > 0) {
      basis_sub <- basis_mat(x, i = l)
      loadings_sub <- loadings_mat(x, i = valid_mask_idx)
      values <- basis_sub %*% t(loadings_sub)

      if (length(x@offset) > 0) {
        values <- sweep(values, 2, x@offset[valid_mask_idx], "+")
      }

      for (t_idx in seq_along(l)) {
        slice_3d <- array(0, dim = out_dims[1:3])
        for (v_idx in seq_along(valid_mask_idx)) {
          positions <- which(mask_indices == valid_mask_idx[v_idx])
          if (length(positions) > 0) {
            slice_3d[positions] <- values[t_idx, v_idx]
          }
        }
        result[, , , t_idx] <- slice_3d
      }
    }

    if (drop) {
      result <- drop(result)
    }

    result
  }
)
