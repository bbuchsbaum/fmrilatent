# Methods for LatentNeuroVec
# series, concat, show, as.matrix, as.array

#' @include latent_neurovector.R
NULL

# --- series Methods ---

#' Extract time series from LatentNeuroVec
#'
#' @description
#' Extract time series data from a \code{LatentNeuroVec} at specified
#' spatial coordinates or voxel indices.
#'
#' @param x A \code{LatentNeuroVec} object
#' @param i,j,k Spatial indices (x, y, z coordinates)
#' @param ... Additional arguments passed to methods
#' @param drop Logical; drop dimensions of length 1 (default TRUE)
#' @return A matrix or vector of time series values
#'
#' @name series-methods
NULL

#' @export
#' @rdname series-methods
setMethod(
  f = "series",
  signature = c(x = "LatentNeuroVec", i = "integer"),
  definition = function(x, i, j, k, ..., drop = TRUE) {
    if (missing(j)) j <- NULL
    if (missing(k)) k <- NULL

    dots <- list(...)
    if ("drop" %in% names(dots)) drop <- dots$drop

    has_j <- !is.null(j)
    has_k <- !is.null(k)

    nTime <- dim(x)[4]
    nels3d <- prod(dim(x)[1:3])

    # CASE A: user gave only i -> interpret as multiple 3D voxel indices
    if (!has_j && !has_k) {
      if (any(i < 1 | i > nels3d)) {
        stop("Some voxel index in 'i' is out of range [1..(X*Y*Z)].")
      }
      n_vox_req <- length(i)

      rowmap <- lookup(x@map, i)
      valid_mask_indices <- rowmap[rowmap > 0]
      map_req_to_valid <- which(rowmap > 0)

      if (length(valid_mask_indices) == 0) {
        if (n_vox_req == 1) {
          return(numeric(nTime))
        } else {
          return(matrix(0, nrow = nTime, ncol = n_vox_req))
        }
      }

      b1 <- basis_mat(x)
      b2 <- loadings_mat(x, i = valid_mask_indices)
      valid_vox_series <- b1 %*% t(b2)

      if (length(x@offset) > 0) {
        valid_vox_series <- sweep(valid_vox_series, 2, x@offset[valid_mask_indices], "+")
      }

      out_mat <- matrix(0, nrow = nTime, ncol = n_vox_req)
      out_mat[, map_req_to_valid] <- as.numeric(valid_vox_series)

      if (drop && n_vox_req == 1) {
        return(drop(out_mat))
      } else {
        return(out_mat)
      }
    } else {
      # CASE B: user gave i,j,k => single voxel
      if (!(length(i) == 1 && length(j) == 1 && length(k) == 1)) {
        stop("series(x, i,j,k): i,j,k must each be a single integer for one voxel.")
      }
      idx_1d <- i + (j - 1) * dim(x)[1] + (k - 1) * dim(x)[1] * dim(x)[2]
      if (idx_1d < 1 || idx_1d > nels3d) {
        stop("Voxel subscript (i,j,k) out of range for LatentNeuroVec.")
      }
      mr <- lookup(x@map, idx_1d)
      if (mr < 1) {
        return(numeric(nTime))
      }

      b1 <- basis_mat(x)
      b2 <- loadings_mat(x, i = mr)
      out_vec <- b1 %*% t(b2)
      if (length(x@offset) > 0) {
        out_vec <- out_vec + x@offset[mr]
      }

      as.vector(out_vec)
    }
  }
)

#' @export
#' @rdname series-methods
setMethod(
  f = "series",
  signature = signature(x = "LatentNeuroVec", i = "numeric"),
  definition = function(x, i, j, k, ..., drop = TRUE) {
    i <- as.integer(i)
    if (!missing(j)) {
      j <- as.integer(j)
      if (!missing(k)) {
        k <- as.integer(k)
        series(x, i, j, k, ..., drop = drop)
      } else {
        series(x, i, j, ..., drop = drop)
      }
    } else {
      series(x, i, ..., drop = drop)
    }
  }
)

#' @export
#' @rdname series-methods
setMethod(
  f = "series",
  signature = signature(x = "LatentNeuroVec"),
  definition = function(x, i, j, k, ..., drop = TRUE) {
    if (missing(i)) {
      stop("series requires at least one index argument 'i'")
    }
    if (is.numeric(i)) i <- as.integer(i)

    if (!missing(j)) {
      j <- as.integer(j)
      if (!missing(k)) {
        k <- as.integer(k)
        series(x, i, j, k, ..., drop = drop)
      } else {
        series(x, i, j, ..., drop = drop)
      }
    } else {
      series(x, i, ..., drop = drop)
    }
  }
)

# --- concat Method ---

#' Concatenate LatentNeuroVec Objects
#'
#' @description
#' Concatenates two or more \code{LatentNeuroVec} objects along the temporal dimension.
#'
#' @param x First \code{LatentNeuroVec}.
#' @param y Second \code{LatentNeuroVec}.
#' @param ... Additional \code{LatentNeuroVec} objects to concatenate.
#'
#' @return A new \code{LatentNeuroVec} if all objects are compatible, otherwise
#'   a \code{\link[neuroim2]{NeuroVecSeq-class}}.
#'
#' @examples
#' \dontrun{
#' combined <- concat(lvec1, lvec2, lvec3)
#' dim(combined)
#' }
#'
#' @importFrom neuroim2 NeuroVecSeq
#' @rdname concat-methods
#' @export
setMethod(
  f = "concat",
  signature = signature(x = "LatentNeuroVec", y = "LatentNeuroVec"),
  definition = function(x, y, ...) {
    additional <- list(...)
    all_objects <- c(list(x, y), additional)
    all_lvecs <- all(sapply(all_objects, is, "LatentNeuroVec"))

    if (!all_lvecs) {
      return(do.call(NeuroVecSeq, list(x, y, ...)))
    }

    # Check dimensions
    x_space <- space(x)
    x_dims_3d <- dim(x_space)[1:3]
    compatible_dims <- TRUE

    for (obj in all_objects[-1]) {
      obj_space <- space(obj)
      obj_dims_3d <- dim(obj_space)[1:3]
      validation_result <- validate_same_dims(x_dims_3d, obj_dims_3d, dims_to_compare = 1:3)
      if (!is.null(validation_result)) {
        compatible_dims <- FALSE
        break
      }
    }

    if (!compatible_dims) {
      return(do.call(NeuroVecSeq, list(x, y, ...)))
    }

    # Check masks
    x_mask_array <- as.array(mask(x))
    compatible_masks <- TRUE

    for (obj in all_objects[-1]) {
      obj_mask_array <- as.array(mask(obj))
      if (!identical(x_mask_array, obj_mask_array)) {
        compatible_masks <- FALSE
        break
      }
    }

    if (!compatible_masks) {
      return(do.call(NeuroVecSeq, list(x, y, ...)))
    }

    # Check k values
    x_k <- .latent_loadings_dim(x@loadings)[2L]
    compatible_k <- TRUE

    for (obj in all_objects[-1]) {
      obj_k <- .latent_loadings_dim(obj@loadings)[2L]
      if (x_k != obj_k) {
        compatible_k <- FALSE
        break
      }
    }

    if (!compatible_k) {
      return(do.call(NeuroVecSeq, list(x, y, ...)))
    }

    # Check loadings - must be identical
    compatible_loadings <- TRUE
    x_loadings <- loadings_mat(x)

    for (obj in all_objects[-1]) {
      if (!identical(as.matrix(x_loadings), as.matrix(loadings_mat(obj)))) {
        compatible_loadings <- FALSE
        break
      }
    }

    if (!compatible_loadings) {
      return(do.call(NeuroVecSeq, list(x, y, ...)))
    }

    # All compatible - create new LatentNeuroVec with concatenated basis
    all_basis <- lapply(all_objects, function(obj) basis_mat(obj))
    time_dims <- sapply(all_objects, function(obj) dim(obj@space)[4])
    total_time <- sum(time_dims)

    new_space_dims <- dim(x_space)
    new_space_dims[4] <- total_time
    new_space <- NeuroSpace(
      dim = new_space_dims,
      spacing = spacing(x_space),
      origin = origin(x_space),
      axes = axes(x_space),
      trans = trans(x_space)
    )

    if (all(sapply(all_basis, is, "sparseMatrix"))) {
      new_basis <- do.call(rbind, all_basis)
    } else {
      all_basis_matrix <- lapply(all_basis, as.matrix)
      new_basis_matrix <- do.call(rbind, all_basis_matrix)
      new_basis <- Matrix::Matrix(new_basis_matrix, sparse = (Matrix::nnzero(new_basis_matrix) / length(new_basis_matrix) < 0.5))
    }

    new_label <- paste0(
      x@label,
      ifelse(nchar(x@label) > 0, "_plus_", ""),
      length(all_objects) - 1,
      "_more"
    )

    LatentNeuroVec(
      basis = new_basis,
      loadings = x@loadings,
      space = new_space,
      mask = x@mask,
      offset = x@offset,
      label = new_label
    )
  }
)

# --- show Method ---

#' Display a LatentNeuroVec object
#'
#' @description
#' Print a formatted summary of a \code{LatentNeuroVec} object showing
#' dimensions, components, memory usage, and other relevant information.
#'
#' @param object A \code{LatentNeuroVec} object
#' @return Invisibly returns NULL; called for side effect of printing
#'
#' @name show-methods
NULL

#' @importFrom crayon bold green blue yellow silver
#' @importFrom utils object.size
#' @importMethodsFrom methods show
#' @rdname show-methods
#' @export
setMethod(
  f = "show",
  signature = "LatentNeuroVec",
  definition = function(object) {
    cat("\n", crayon::bold(crayon::blue("LatentNeuroVec Object")), "\n")
    cat(crayon::silver("======================\n"))

    dims <- dim(object)
    spatial_dims <- paste(dims[1:3], collapse = " x ")
    cat("\n", crayon::yellow("Dimensions:"), "\n")
    cat(" ", crayon::silver("*"), " Spatial: ", crayon::green(spatial_dims), "\n")
    cat(" ", crayon::silver("*"), " Temporal: ", crayon::green(dims[4]), "\n")

    basis_dim <- .latent_basis_dim(object@basis)
    loadings_dim <- .latent_loadings_dim(object@loadings)
    n_components <- basis_dim[2L]
    first_basis_coeffs <- format(basis_mat(object, i = 1:min(5, basis_dim[1L]), j = 1), digits = 3)
    cat("\n", crayon::yellow("Components:"), "\n")
    cat(" ", crayon::silver("*"), " Number: ", crayon::green(n_components), "\n")
    cat(
      " ", crayon::silver("*"), " First component, first 5 coeffs: ",
      crayon::green(paste(first_basis_coeffs, collapse = ", ")),
      if (length(first_basis_coeffs) < 5) "" else "...", "\n"
    )

    basis_size <- if (is(object@basis, "Matrix")) {
      format(object.size(object@basis), units = "auto")
    } else {
      "handle (lazy)"
    }
    loadings_size <- if (is(object@loadings, "Matrix")) {
      format(object.size(object@loadings), units = "auto")
    } else {
      "handle (lazy)"
    }
    total_size <- format(object.size(object), units = "auto")

    cat("\n", crayon::yellow("Memory Usage:"), "\n")
    cat(" ", crayon::silver("*"), " Basis: ", crayon::green(basis_size), "\n")
    cat(" ", crayon::silver("*"), " Loadings: ", crayon::green(loadings_size), "\n")
    cat(" ", crayon::silver("*"), " Total: ", crayon::green(total_size), "\n")

    n_nonzero <- sum(object@mask)
    sparsity <- round(100 * n_nonzero / prod(dims[1:3]), 2)
    cat("\n", crayon::yellow("Sparsity:"), "\n")
    cat(" ", crayon::silver("*"), " Non-zero voxels: ", crayon::green(n_nonzero), "\n")
    cat(" ", crayon::silver("*"), " Coverage: ", crayon::green(sparsity), "%\n")

    sp <- space(object)
    spacing_str <- paste(round(spacing(sp), 2), collapse = " x ")
    origin_str <- paste(round(origin(sp), 2), collapse = " x ")
    cat("\n", crayon::yellow("Space Information:"), "\n")
    cat(" ", crayon::silver("*"), " Spacing: ", crayon::green(spacing_str), "\n")
    cat(" ", crayon::silver("*"), " Origin:  ", crayon::green(origin_str), "\n")

    cat("\n", crayon::bold("Data Access:"), "\n")
    cat("\n", crayon::yellow("Reconstructed Space Access:"), "\n")
    cat(
      " ", crayon::silver("*"), " Extract volume: ",
      crayon::blue("object[[i]]"),
      crayon::silver("  # 3D volume at timepoint i\n")
    )
    cat(
      " ", crayon::silver("*"), " Get value: ",
      crayon::blue("object[i]"),
      crayon::silver("  # Value at linear index i\n")
    )
    cat(
      " ", crayon::silver("*"), " Subset: ",
      crayon::blue("object[mask]"),
      crayon::silver("  # Values at mask positions\n")
    )

    cat("\n", crayon::yellow("Latent Space Access:"), "\n")
    cat(
      " ", crayon::silver("*"), " Basis vectors: ",
      crayon::blue("basis(object)"),
      crayon::silver("  # nxk temporal basis\n")
    )
    cat(
      " ", crayon::silver("*"), " Loadings: ",
      crayon::blue("loadings(object)"),
      crayon::silver("  # pxk spatial loadings\n")
    )

    cat("\n", crayon::yellow("Conversions:"), "\n")
    cat(
      " ", crayon::silver("*"), " as.array(object): ",
      crayon::silver("4D reconstruction\n")
    )
    cat(
      " ", crayon::silver("*"), " as.matrix(object): ",
      crayon::silver("nxp matrix of reconstructed values\n")
    )

    cat(
      "\n", crayon::silver("Note: All access methods reconstruct data (X = B x L^T + offset)"),
      "\n", crayon::silver("unless you're explicitly accessing latent space."), "\n\n"
    )
  }
)

# --- Conversions ---

setMethod(
  f = "as.matrix",
  signature = signature(x = "LatentNeuroVec"),
  definition = function(x, ...) {
    B <- basis_mat(x)
    L <- loadings_mat(x)
    mat <- B %*% t(L)
    if (length(x@offset) > 0) {
      mat <- sweep(mat, 2, x@offset, "+")
    }
    as.matrix(mat)
  }
)

setMethod(
  f = "as.array",
  signature = signature(x = "LatentNeuroVec"),
  definition = function(x, ...) {
    dims <- dim(x)
    mat <- as.matrix(x)
    arr <- array(0, dim = dims)
    mask_idx <- which(as.logical(x@mask))
    for (t in seq_len(dims[4])) {
      slice <- numeric(prod(dims[1:3]))
      slice[mask_idx] <- mat[t, ]
      arr[, , , t] <- array(slice, dim = dims[1:3])
    }
    arr
  }
)
