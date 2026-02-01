#' @importFrom methods setGeneric
NULL

#' Get the basis matrix (temporal components)
#'
#' @description
#' Extract the basis matrix from a latent space representation. For
#' \code{LatentNeuroVec} objects, this returns the temporal basis matrix
#' with dimensions (nTime x k) where k is the number of components.
#'
#' @param x An object containing a basis matrix (e.g., \code{LatentNeuroVec})
#' @param ... Additional arguments (currently unused)
#' @return The basis matrix (typically time x components)
#'
#' @examples
#' \dontrun{
#' # For LatentNeuroVec:
#' b_matrix <- basis(lvec)
#' print(dim(b_matrix)) # nTime x nComponents
#' }
#'
#' @export
#' @rdname basis-methods
setGeneric("basis", function(x, ...) standardGeneric("basis"))

#' Get the loadings matrix (spatial components)
#'
#' @description
#' Extract the loadings matrix from a latent space representation. For
#' \code{LatentNeuroVec} objects, this returns the spatial loadings matrix
#' with dimensions (nVoxels x k) where k is the number of components and
#' nVoxels is the number of voxels within the mask.
#'
#' @param x An object containing loadings (e.g., \code{LatentNeuroVec})
#' @param ... Additional arguments (currently unused)
#' @return The loadings matrix (typically voxels x components)
#'
#' @examples
#' \dontrun{
#' # For LatentNeuroVec:
#' l_matrix <- loadings(lvec)
#' print(dim(l_matrix)) # nVoxels x nComponents
#' }
#'
#' @export
#' @rdname loadings-methods
setGeneric("loadings", function(x, ...) standardGeneric("loadings"))

#' Get the offset vector
#'
#' @description
#' Extract the offset vector from a latent space representation. For
#' \code{LatentNeuroVec} objects, this returns the voxel-wise offset
#' (mean or intercept) that is added after the basis x loadings reconstruction.
#'
#' @param object An object containing an offset (e.g., \code{LatentNeuroVec})
#' @return The offset vector (length = number of voxels in mask)
#'
#' @examples
#' \dontrun{
#' # For LatentNeuroVec:
#' off_vector <- offset(lvec)
#' print(length(off_vector)) # Should equal number of voxels in mask
#' }
#'
#' @export
#' @rdname offset-methods
#' @name offset
NULL

# Note: stats::offset exists with signature function(object).
# Calling setGeneric("offset", ...) will create an S4 generic using that signature
# (and will print a message); we then provide an S4 method for LatentNeuroVec.
if (!isGeneric("offset", where = topenv(parent.frame()))) {
  setGeneric("offset", function(object) standardGeneric("offset"))
}

#' Get the map object
#'
#' @description
#' Extract the index map from a latent space representation. For
#' \code{LatentNeuroVec} objects, this returns the \code{IndexLookupVol}
#' that maps 3D coordinates to linear mask indices.
#'
#' @param x An object containing a map (e.g., \code{LatentNeuroVec})
#' @return The index lookup object
#'
#' @export
#' @rdname map-methods
#' @name map
NULL

# Provide generics when not already available
if (!isGeneric("map", where = topenv(parent.frame()))) {
  setGeneric("map", function(x, ...) standardGeneric("map"))
}

#' Get the mask
#'
#' @description
#' Extract the brain mask from a latent space representation. For
#' \code{LatentNeuroVec} objects, this returns the \code{LogicalNeuroVol}
#' defining which voxels are included in the representation.
#'
#' @param x An object containing a mask (e.g., \code{LatentNeuroVec})
#' @return The logical mask volume
#'
#' @export
#' @importMethodsFrom neuroim2 mask
#' @rdname mask-methods
#' @name mask
NULL

# mask generic is imported from neuroim2

#' Linear access to LatentNeuroVec elements
#'
#' @description
#' Access elements of a \code{LatentNeuroVec} using linear (1D) indexing.
#' Reconstructs values on-the-fly from the basis and loadings factorization.
#'
#' @param x A \code{LatentNeuroVec} object
#' @param i Integer indices into the 4D volume (linear indexing)
#' @return Numeric vector of reconstructed values
#'
#' @name linear_access
#' @rdname linear_access-methods
NULL

#' Extract time series from LatentNeuroVec
#'
#' @description
#' Extract time series data from a \code{LatentNeuroVec}. Reconstructs
#' values on-the-fly from the basis and loadings factorization.
#'
#' @param x A \code{LatentNeuroVec} object
#' @param i,j,k Spatial indices (3D) or linear indices
#' @param ... Additional arguments passed to methods
#' @param drop Logical, whether to drop dimensions
#' @return Matrix or vector of time series values
#'
#' @name series
#' @rdname series-methods
NULL

#' Display a LatentNeuroVec object
#'
#' @description
#' Print a summary of a \code{LatentNeuroVec} object including dimensions,
#' component count, memory usage, and sparsity information.
#'
#' @param object A \code{LatentNeuroVec} object
#' @return Called for side effect (printing). Returns \code{NULL} invisibly.
#'
#' @name show
#' @rdname show-methods
NULL
