#' @include latent_handles.R
#' @importFrom neuroim2 NeuroVec IndexLookupVol LogicalNeuroVol NeuroSpace
#' @importFrom Matrix Matrix
#' @importClassesFrom Matrix Matrix dgCMatrix dgeMatrix
#' @importFrom methods setClass setOldClass
NULL

# Register S3 spec classes for S4 method dispatch
setOldClass("spec_slepian")
setOldClass("spec_heat_wavelet")
setOldClass("spec_diffusion_wavelet")
setOldClass("spec_pca")
setOldClass("spec_flat")
setOldClass("ImplicitLatent")
setOldClass("ParcelBasisTemplate")

#' LatentNeuroVec Class
#'
#' @description
#' A class that represents a 4-dimensional neuroimaging array using a latent space
#' decomposition. It stores the data as a set of basis functions (dictionary) and
#' a corresponding set of loadings (coefficients), enabling efficient representation
#' and manipulation of high-dimensional data.
#'
#' @slot basis A \code{Matrix} or \code{BasisHandle} where each column represents a basis
#'   vector in the latent space. Dimensions are (nTime x k) where k is the number of components.
#' @slot loadings A \code{Matrix} or \code{LoadingsHandle} (often sparse) containing the
#'   coefficients for each basis vector across the spatial dimensions. Dimensions are
#'   (nVoxels x k).
#' @slot offset A \code{numeric} vector representing a constant offset term for
#'   each voxel or spatial location. Length is nVoxels (within mask).
#' @slot map A \code{IndexLookupVol} object representing the mapping from 3D coordinates
#'   to linear indices within the mask.
#' @slot label A \code{character} string representing the label for the latent vector.
#' @slot meta A \code{list} storing lightweight metadata (e.g., HRBF params, centers).
#'
#' @details
#' \code{LatentNeuroVec} inherits from \code{\link[neuroim2]{NeuroVec-class}}
#' and \code{\link[neuroim2]{AbstractSparseNeuroVec-class}}. The original 4D data
#' can be reconstructed as:
#' \deqn{data[v,t] = \sum_k \bigl(basis[t,k] \times loadings[v,k]\bigr) + offset[v]}
#' where v indexes voxels within the mask and t indexes time points.
#'
#' This approach is especially useful for large datasets where storing the full
#' 4D array is expensive. Common use cases include:
#' \itemize{
#'   \item PCA-reduced fMRI data
#'   \item ICA decompositions
#'   \item Dictionary learning representations
#'   \item Any matrix factorization of neuroimaging data
#' }
#'
#' @section Inheritance:
#' \code{LatentNeuroVec} inherits from:
#' \itemize{
#'   \item \code{\link[neuroim2]{NeuroVec-class}} - base 4D neuroimaging class
#'   \item \code{\link[neuroim2]{AbstractSparseNeuroVec-class}} - sparse representation framework
#' }
#'
#' @seealso
#' \code{\link[neuroim2]{NeuroVec-class}},
#' \code{\link[neuroim2]{AbstractSparseNeuroVec-class}},
#' \code{\link{LatentNeuroVec}} constructor function.
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' library(neuroim2)
#'
#' # Example dimensions
#' n_timepoints <- 100
#' n_components <- 10
#' n_voxels <- 1000
#'
#' # Create basis (temporal) & loadings (spatial)
#' basis <- Matrix(rnorm(n_timepoints * n_components),
#'   nrow = n_timepoints, ncol = n_components
#' )
#' loadings <- Matrix(rnorm(n_voxels * n_components),
#'   nrow = n_voxels, ncol = n_components, sparse = TRUE
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
#'
#' print(lvec)
#' print(dim(lvec))
#' }
#'
#' @export
#' @rdname LatentNeuroVec-class
setClass("LatentNeuroVec",
  slots = c(
    basis = "MatrixOrBasisHandle",
    loadings = "MatrixOrLoadingsHandle",
    offset = "numeric",
    map = "IndexLookupVol",
    label = "character",
    meta = "list"
  ),
  contains = c("NeuroVec", "AbstractSparseNeuroVec")
)

#' HierarchicalBasisTemplate Class
#'
#' @description
#' A template-only, data-agnostic container for hierarchical Laplacian frames.
#' Stores the sparse spatial dictionary (primal basis), a cached solver for the
#' Gram matrix, and the parcellation hierarchy metadata. Intended to be built
#' offline for a fixed template (e.g., MNI) and reused for encoding fMRI data.
#'
#' @slot mask LogicalNeuroVol defining the domain (3D).
#' @slot space NeuroSpace (typically 4D with a singleton time dim) matching the mask.
#' @slot levels List of integer vectors (one per level) giving parcel ids per voxel (mask order).
#' @slot parents List mapping child parcel ids to parent ids for each level > 1.
#' @slot loadings Sparse Matrix (voxels x atoms) containing concatenated atoms B.
#' @slot gram_factor Cached factorization of \eqn{G = t(B) \%*\% B} (e.g., dCHMsimpl from Matrix::Cholesky).
#' @slot atoms data.frame describing each atom (col_id, level, parcel_id, parent_id, mode, label).
#' @slot meta List for auxiliary metadata (atlas names, k_per_level, ridge, version).
#' @export
setClass("HierarchicalBasisTemplate",
  slots = c(
    mask = "LogicalNeuroVol",
    space = "NeuroSpace",
    levels = "list",
    parents = "list",
    loadings = "Matrix",
    gram_factor = "ANY",
    atoms = "data.frame",
    meta = "list"
  )
)
