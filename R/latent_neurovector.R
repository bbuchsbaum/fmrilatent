#' @include all_generic.R latent_neurovec_materialize.R
#' @importFrom neuroim2 matricized_access concat axes indices space origin spacing trans
#' @importFrom neuroim2 NeuroSpace LogicalNeuroVol SparseNeuroVol IndexLookupVol NeuroVecSeq
#' @importFrom neuroim2 drop_dim lookup series linear_access
#' @importFrom Matrix Matrix t nnzero
#' @importFrom methods setMethod signature new is as callNextMethod setValidity setGeneric isGeneric
#' @importFrom crayon bold blue green yellow silver
#' @importFrom utils object.size
NULL

# Ensure generics exist (defensive in case Collate order is bypassed)
if (!isGeneric("map")) {
  setGeneric("map", function(x, ...) standardGeneric("map"))
}

#' Helper function to check dimension consistency
#' @keywords internal
#' @noRd
check_same_dims <- function(x, y, dims_to_compare = NULL, msg = "Dimension mismatch") {
  x_dims <- if (is.numeric(x)) x else dim(x)
  y_dims <- if (is.numeric(y)) y else dim(y)

  if (!is.null(dims_to_compare)) {
    x_dims <- x_dims[dims_to_compare]
    y_dims <- y_dims[dims_to_compare]
  }

  if (!identical(x_dims, y_dims)) {
    stop(paste0(msg, ": ", paste(x_dims, collapse = "x"), " vs ", paste(y_dims, collapse = "x")))
  }
  invisible(TRUE)
}

#' Helper function to validate dimension consistency (returns message or NULL)
#' @keywords internal
#' @noRd
validate_same_dims <- function(x, y, dims_to_compare = NULL, msg = "Dimension mismatch") {
  x_dims <- if (is.numeric(x)) x else dim(x)
  y_dims <- if (is.numeric(y)) y else dim(y)

  if (!is.null(dims_to_compare)) {
    x_dims <- x_dims[dims_to_compare]
    y_dims <- y_dims[dims_to_compare]
  }

  if (!identical(x_dims, y_dims)) {
    return(paste0(msg, ": ", paste(x_dims, collapse = "x"), " vs ", paste(y_dims, collapse = "x")))
  }
  NULL
}

#' Construct a LatentNeuroVec Object
#'
#' @title Create a Latent Space Representation of Neuroimaging Data
#' @description
#' Constructs a \code{\link{LatentNeuroVec-class}} object, which provides a memory-efficient
#' representation of neuroimaging data using matrix factorization. This is particularly useful
#' for dimensionality reduction techniques (e.g., PCA or ICA).
#'
#' @param basis A numeric or \code{Matrix} object (\eqn{n \times k}) containing the temporal basis.
#'   Each row corresponds to a time point, each column to a component.
#' @param loadings A numeric or \code{Matrix} object (\eqn{p \times k}) containing spatial loadings.
#'   Each row corresponds to a voxel within the mask, each column to a component.
#' @param space A \code{\link[neuroim2]{NeuroSpace-class}} defining the spatial/temporal dimensions.
#'   Must be 4-dimensional.
#' @param mask A \code{\link[neuroim2]{LogicalNeuroVol-class}} defining the brain mask.
#'   The number of TRUE values must equal nrow(loadings).
#' @param offset Optional numeric vector of length \eqn{p} (voxel-wise offsets).
#'   If NULL, defaults to zero offset.
#' @param label Optional character label for the object.
#'
#' @return A new \code{\link{LatentNeuroVec-class}} instance.
#'
#' @details
#' The data is represented as the product:
#' \deqn{X = B \times L^T + c}
#' where:
#' \itemize{
#'   \item B is the basis matrix (\eqn{n \times k})
#'   \item L is the loadings matrix (\eqn{p \times k})
#'   \item c is an optional offset vector (length p)
#'   \item n is the number of time points
#'   \item p is the number of voxels in the mask
#'   \item k is the number of components
#' }
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' library(neuroim2)
#'
#' # Example data
#' n_timepoints <- 100
#' n_components <- 10
#' n_voxels <- 1000
#'
#' # Create basis & loadings
#' basis <- Matrix(rnorm(n_timepoints * n_components),
#'   nrow = n_timepoints,
#'   ncol = n_components
#' )
#' loadings <- Matrix(rnorm(n_voxels * n_components),
#'   nrow = n_voxels,
#'   ncol = n_components,
#'   sparse = TRUE
#' )
#'
#' # Create space (10x10x10 volume, 100 timepoints)
#' spc <- NeuroSpace(c(10, 10, 10, n_timepoints))
#'
#' # Create mask
#' mask_array <- array(TRUE, dim = c(10, 10, 10))
#' mask_vol <- LogicalNeuroVol(mask_array, NeuroSpace(c(10, 10, 10)))
#'
#' # Construct LatentNeuroVec
#' lvec <- LatentNeuroVec(
#'   basis = basis,
#'   loadings = loadings,
#'   space = spc,
#'   mask = mask_vol
#' )
#' }
#'
#' @param meta Optional list of metadata (e.g., HRBF params or centres).
#' @export
LatentNeuroVec <- function(basis, loadings, space, mask, offset = NULL, label = "", meta = list()) {
  # Validate 'space'
  if (!inherits(space, "NeuroSpace")) {
    stop("'space' must be a NeuroSpace object")
  }
  # Validate 'basis' / 'loadings' (accept Matrix or handle)
  if (!is.matrix(basis) && !inherits(basis, c("Matrix", "BasisHandle"))) {
    stop("'basis' must be a matrix, Matrix, or BasisHandle object")
  }
  if (!is.matrix(loadings) && !inherits(loadings, c("Matrix", "LoadingsHandle"))) {
    stop("'loadings' must be a matrix, Matrix, or LoadingsHandle object")
  }

  # Ensure we have a LogicalNeuroVol for mask
  # Initialize space_for_map using the main 3D space first
  space_3d <- drop_dim(space)
  space_for_map <- space_3d


  if (!inherits(mask, "LogicalNeuroVol")) {
    # Construct 3D space from the main 4D space
    mask <- LogicalNeuroVol(as.logical(mask), space_3d)
  } else {
    # If mask is already a LogicalNeuroVol, ensure its space matches the 3D part of the main space
    mask_space <- neuroim2::space(mask)
    # Use check_same_dims for the dimension check
    check_same_dims(mask_space, space_3d,
      msg = "Space dimensions of provided mask do not match dimensions derived from the main 4D space."
    )
    # Keep the check for space equality separate for now
    if (!isTRUE(all.equal(mask_space, space_3d))) {
      stop("Space object of provided mask does not match the space derived from the main 4D space. Cannot create IndexLookupVol.")
    }
  }

  cardinality <- sum(mask)

  # Handle offset: treat NULL or empty numeric offset as no offset
  if (is.null(offset) || length(offset) == 0) {
    # No offset provided or empty offset, use numeric(0)
    offset <- numeric(0)
  } else {
    # Check that provided offset matches loadings rows
    if (length(offset) != nrow(loadings)) {
      stop("'offset' length must match number of rows in 'loadings'")
    }
  }

  # Dimension checks - allow handles
  b_dim <- .latent_basis_dim(basis)
  l_dim <- .latent_loadings_dim(loadings)

  if (b_dim[2L] != l_dim[2L]) {
    stop("'basis' and 'loadings' must have the same number of columns (components)")
  }
  if (b_dim[1L] != dim(space)[4]) {
    stop("'basis' must have ", dim(space)[4], " rows (the 4th dimension of space)")
  }
  if (l_dim[1L] != cardinality) {
    stop(paste0("'loadings' must have ", cardinality, " rows (i.e. #non-zero in mask)"))
  }

  # Ensure all numeric inputs are finite
  if (is(basis, "Matrix") && !all(is.finite(basis))) {
    stop("'basis' must contain only finite values")
  }
  if (is(loadings, "Matrix") && !all(is.finite(loadings))) {
    stop("'loadings' must contain only finite values")
  }
  if (length(offset) > 0 && !all(is.finite(offset))) {
    stop("'offset' must contain only finite values")
  }

  # Convert basis/loadings to Matrix objects when provided as base matrices
  if (is.matrix(basis) && !is(basis, "Matrix")) {
    density_basis <- sum(basis != 0) / length(basis)
    if (density_basis > 0.5) {
      message("Input 'basis' is dense (", round(density_basis * 100), "% non-zero); storing as dense dgeMatrix.")
      basis <- Matrix::Matrix(basis, sparse = FALSE)
    } else {
      basis <- Matrix::Matrix(basis)
    }
  }

  if (is.matrix(loadings) && !is(loadings, "Matrix")) {
    density_loadings <- sum(loadings != 0) / length(loadings)
    if (density_loadings > 0.5) {
      message("Input 'loadings' is dense (", round(density_loadings * 100), "% non-zero); storing as dense dgeMatrix.")
      loadings <- Matrix::Matrix(loadings, sparse = FALSE)
    } else {
      loadings <- Matrix::Matrix(loadings)
    }
  }

  # Check component count
  k <- b_dim[2L]
  if (k < 1) {
    stop("Number of components (columns in basis and loadings) must be >= 1")
  }

  if (is.null(meta)) meta <- list()
  if (!is.list(meta)) stop("'meta' must be a list")

  # Create the object
  new("LatentNeuroVec",
    basis    = basis,
    loadings = loadings,
    space    = space,
    mask     = mask,
    map      = IndexLookupVol(space_for_map, as.integer(which(mask@.Data))),
    offset   = offset,
    label    = label,
    meta     = meta
  )
}

# --- LatentNeuroVec-specific materialization methods ---
# (base generics and Matrix/Handle methods are in latent_neurovec_materialize.R)

setMethod(
  "basis_mat", "LatentNeuroVec",
  function(x, i = NULL, j = NULL, ...) {
    basis_mat(x@basis, i = i, j = j, ...)
  }
)

setMethod(
  "loadings_mat", "LatentNeuroVec",
  function(x, i = NULL, j = NULL, ...) {
    loadings_mat(x@loadings, i = i, j = j, ...)
  }
)

# --- Accessor Methods ---

#' @export
#' @rdname basis-methods
setMethod("basis", "LatentNeuroVec", function(x) basis_mat(x))

#' @export
#' @rdname loadings-methods
setMethod("loadings", "LatentNeuroVec", function(x) loadings_mat(x))

#' @export
#' @rdname offset-methods
setMethod("offset", "LatentNeuroVec", function(object) object@offset)

#' @export
#' @rdname mask-methods
setMethod("mask", "LatentNeuroVec", function(x) x@mask)

#' @export
#' @rdname map-methods
setMethod("map", "LatentNeuroVec", function(x) x@map)

# Note: matricized_access, linear_access, [[, [ methods are in latent_indexing.R
# Note: series, concat, show, as.matrix, as.array methods are in latent_methods.R

# --- Validity Check ---

#' Validity Check for LatentNeuroVec Objects
#'
#' @param object A LatentNeuroVec object
#' @return TRUE if the object is valid, otherwise a character vector of error messages.
#' @keywords internal
#' @noRd
.validate_LatentNeuroVec <- function(object) {
  errors <- character()

  if (!inherits(object@basis, c("Matrix", "matrix", "BasisHandle"))) {
    errors <- c(errors, "Slot @basis must be Matrix, matrix, or BasisHandle.")
  }
  if (!inherits(object@loadings, c("Matrix", "matrix", "LoadingsHandle"))) {
    errors <- c(errors, "Slot @loadings must be Matrix, matrix, or LoadingsHandle.")
  }
  if (!is.numeric(object@offset)) {
    errors <- c(errors, "Slot @offset must be numeric.")
  }
  if (!inherits(object@mask, "LogicalNeuroVol")) {
    errors <- c(errors, "Slot @mask must be a LogicalNeuroVol object.")
  }
  if (!inherits(object@map, "IndexLookupVol")) {
    errors <- c(errors, "Slot @map must be an IndexLookupVol object.")
  }
  if (!is.character(object@label) || length(object@label) != 1) {
    errors <- c(errors, "Slot @label must be a single character string.")
  }
  if (!inherits(object@space, "NeuroSpace")) {
    errors <- c(errors, "Slot @space must be a NeuroSpace object.")
    return(errors)
  }

  if (length(errors) == 0) {
    s_dims <- dim(object@space)
    if (length(s_dims) != 4) {
      errors <- c(errors, "Slot @space must have 4 dimensions.")
    } else {
      b_dim <- .latent_basis_dim(object@basis)
      l_dim <- .latent_loadings_dim(object@loadings)

      if (b_dim[2L] != l_dim[2L]) {
        errors <- c(errors, paste0(
          "Component mismatch: k_b = ", b_dim[2L],
          " != k_l = ", l_dim[2L]
        ))
      }

      if (b_dim[1L] != s_dims[4]) {
        errors <- c(errors, paste0(
          "Time mismatch: n_time(basis) = ", b_dim[1L],
          " != dim(@space)[4] = ", s_dims[4]
        ))
      }

      dim_check_result <- validate_same_dims(
        object@mask,
        object@space,
        dims_to_compare = 1:3,
        msg = "[.validate_LatentNeuroVec] Mask/Space dim mismatch:"
      )
      if (!is.null(dim_check_result)) {
        errors <- c(errors, dim_check_result)
      }

      nVox_mask <- sum(object@mask)
      if (l_dim[1L] != nVox_mask) {
        errors <- c(errors, paste0(
          "Loadings rows (", l_dim[1L],
          ") mismatch non-zero count in mask (", nVox_mask, ")"
        ))
      }

      if (length(object@offset) > 0 && length(object@offset) != l_dim[1L]) {
        errors <- c(errors, paste0(
          "Offset length (", length(object@offset),
          ") mismatch number of rows in loadings (", l_dim[1L], ")"
        ))
      }

      if (length(object@map@indices) != nVox_mask) {
        errors <- c(errors, paste0(
          "Map indices length (", length(object@map@indices),
          ") mismatch non-zero count in mask (", nVox_mask, ")"
        ))
      }
    }
  }

  if (length(errors) == 0) TRUE else errors
}

#' @keywords internal
setValidity("LatentNeuroVec", .validate_LatentNeuroVec)
