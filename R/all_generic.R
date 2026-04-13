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

#' Reconstruct a latent object as a matrix
#'
#' @description
#' Provides a common reconstruction interface for both explicit
#' \code{LatentNeuroVec} objects and implicit decoder-backed latent objects.
#'
#' @param x A latent object.
#' @param time_idx Optional integer time indices to keep.
#' @param roi_mask Optional logical ROI mask for spatial subsetting.
#' @param ... Additional arguments passed to methods.
#' @return Numeric matrix with rows = time and columns = voxels within the
#'   requested mask support.
#' @export
setGeneric("reconstruct_matrix", function(x, time_idx = NULL, roi_mask = NULL, ...) {
  standardGeneric("reconstruct_matrix")
})

#' Reconstruct a latent object as a 4D array
#'
#' @description
#' Reconstructs a latent object into a 4D array over its spatial support.
#'
#' @param x A latent object.
#' @param time_idx Optional integer time indices to keep.
#' @param roi_mask Optional logical ROI mask; voxels outside the ROI are zero.
#' @param ... Additional arguments passed to methods.
#' @return Numeric 4D array.
#' @export
setGeneric("reconstruct_array", function(x, time_idx = NULL, roi_mask = NULL, ...) {
  standardGeneric("reconstruct_array")
})

#' Get lightweight metadata from a latent object
#'
#' @param x A latent object.
#' @param ... Additional arguments (unused).
#' @return Metadata list.
#' @export
setGeneric("latent_meta", function(x, ...) standardGeneric("latent_meta"))

#' Test whether a latent object is explicit
#'
#' @param x A latent object.
#' @param ... Additional arguments (unused).
#' @return Logical scalar; \code{TRUE} when the object stores explicit basis and
#'   loadings matrices.
#' @export
setGeneric("is_explicit_latent", function(x, ...) standardGeneric("is_explicit_latent"))

#' Extract coefficient time series from a latent object
#'
#' @param x A latent object.
#' @param coordinates Coordinate system for the returned coefficients.
#' @param ... Additional arguments passed to methods.
#' @return Numeric matrix with rows = time and columns = latent coefficients.
#' @export
setGeneric("coef_time", function(x, coordinates = c("analysis", "raw"), ...) {
  standardGeneric("coef_time")
})

#' Extract coefficient-space metric from a latent object
#'
#' @param x A latent object.
#' @param coordinates Coordinate system for the requested metric.
#' @param ... Additional arguments passed to methods.
#' @return Matrix-like metric representation or \code{NULL}. For
#'   transport-backed latent objects, analysis coordinates are Euclidean by
#'   contract in v1. Raw-coordinate metrics are returned when the
#'   raw-to-analysis transform exposes a linear matrix representation.
#' @export
setGeneric("coef_metric", function(x, coordinates = c("raw", "analysis"), ...) {
  standardGeneric("coef_metric")
})

#' Describe the transform from raw to analysis coordinates
#'
#' @param x A latent object.
#' @param ... Additional arguments passed to methods.
#' @return Transform descriptor or \code{NULL}. Downstream model-fitting code
#'   should ordinarily consume \code{coef_time(x, "analysis")} rather than
#'   reasoning about raw coordinates directly.
#' @export
setGeneric("analysis_transform", function(x, ...) standardGeneric("analysis_transform"))

#' Extract the basis asset behind a latent object
#'
#' @param x A latent object.
#' @param ... Additional arguments passed to methods.
#' @return Basis asset object or \code{NULL}.
#' @export
setGeneric("basis_asset", function(x, ...) standardGeneric("basis_asset"))

#' Get a decoder view for a latent object
#'
#' @param x A latent object.
#' @param space Output space to decode into.
#' @param coordinates Coordinate system consumed by the decoder.
#' @param ... Additional arguments passed to methods.
#' @return Decoder view object.
#' @export
setGeneric("decoder", function(x, space = c("native", "template"),
                                 coordinates = c("analysis", "raw"), ...) {
  standardGeneric("decoder")
})

#' Decode coefficient-space vectors into an output space
#'
#' @param x A latent object.
#' @param gamma Coefficient-space vector or matrix.
#' @param space Output space to decode into.
#' @param coordinates Coordinate system used by \code{gamma}.
#' @param wrap If \code{"auto"}, wrap the decoded values with
#'   \code{\link{wrap_decoded}()} so the result is domain-aware (for
#'   example a \code{NeuroVol} for volumetric targets or a surface field
#'   for surface targets). Defaults to \code{"none"}, which returns the
#'   raw numeric vector or matrix.
#' @param ... Additional arguments passed to methods.
#' @return Numeric vector or matrix in the requested output space, or a
#'   domain-aware wrapper when \code{wrap = "auto"}.
#' @export
setGeneric("decode_coefficients", function(x, gamma, space = c("native", "template"),
                                            coordinates = c("analysis", "raw"),
                                            wrap = c("none", "auto"), ...) {
  standardGeneric("decode_coefficients")
})

#' Push coefficient covariance into an output space
#'
#' @param x A latent object.
#' @param Sigma Coefficient covariance matrix.
#' @param space Output space to decode into.
#' @param coordinates Coordinate system used by \code{Sigma}.
#' @param diag_only Logical; if \code{TRUE}, return only the diagonal.
#' @param ... Additional arguments passed to methods.
#' @return Numeric vector or matrix in the requested output space.
#' @export
setGeneric("decode_covariance", function(x, Sigma, space = c("native", "template"),
                                          coordinates = c("analysis", "raw"),
                                          diag_only = TRUE, ...) {
  standardGeneric("decode_covariance")
})

#' Compatibility wrapper for decoder-based coefficient projection
#'
#' @param x A latent object.
#' @param gamma Coefficient-space vector or matrix.
#' @param space Output space to decode into.
#' @param coordinates Coordinate system used by \code{gamma}.
#' @param ... Additional arguments passed to methods.
#' @return Numeric vector or matrix in the requested output space.
#' @export
setGeneric("project_effect", function(x, gamma, space = c("native", "template"),
                                       coordinates = c("analysis", "raw"), ...) {
  standardGeneric("project_effect")
})

#' Compatibility wrapper for decoder-based covariance pushforward
#'
#' @param x A latent object.
#' @param Sigma Coefficient covariance matrix.
#' @param space Output space to decode into.
#' @param coordinates Coordinate system used by \code{Sigma}.
#' @param diag_only Logical; if \code{TRUE}, return only the diagonal.
#' @param ... Additional arguments passed to methods.
#' @return Numeric vector or matrix in the requested output space.
#' @export
setGeneric("project_vcov", function(x, Sigma, space = c("native", "template"),
                                     coordinates = c("analysis", "raw"),
                                     diag_only = TRUE, ...) {
  standardGeneric("project_vcov")
})

#' Extract template loadings
#'
#' @param x A template object.
#' @param ... Additional arguments (unused).
#' @return Spatial loadings matrix stored by the template.
#' @export
setGeneric("template_loadings", function(x, ...) standardGeneric("template_loadings"))

#' Extract template mask
#'
#' @param x A template object.
#' @param ... Additional arguments (unused).
#' @return Mask associated with the template.
#' @export
setGeneric("template_mask", function(x, ...) standardGeneric("template_mask"))

#' Extract template metadata
#'
#' @param x A template object.
#' @param ... Additional arguments (unused).
#' @return Metadata list.
#' @export
setGeneric("template_meta", function(x, ...) standardGeneric("template_meta"))

#' Extract a basis decoder from a template asset
#'
#' @param template A template object.
#' @param ... Additional arguments passed to methods.
#' @return Basis decoder view object.
#' @export
setGeneric("basis_decoder", function(template, ...) standardGeneric("basis_decoder"))

#' Query the rank of a template basis
#'
#' @param template A template object.
#' @param ... Additional arguments passed to methods.
#' @return Integer scalar rank.
#' @export
setGeneric("template_rank", function(template, ...) standardGeneric("template_rank"))

#' Extract the domain associated with a template
#'
#' @param template A template object.
#' @param ... Additional arguments passed to methods.
#' @return Domain object or identifier.
#' @export
setGeneric("template_domain", function(template, ...) standardGeneric("template_domain"))

#' Extract the support associated with a template
#'
#' @param template A template object.
#' @param ... Additional arguments passed to methods.
#' @return Support object describing the decoded support of the template.
#' @export
setGeneric("template_support", function(template, ...) standardGeneric("template_support"))

#' Extract optional measure or mass information for a template
#'
#' @param template A template object.
#' @param ... Additional arguments passed to methods.
#' @return Optional measure or mass object, or \code{NULL}.
#' @export
setGeneric("template_measure", function(template, ...) standardGeneric("template_measure"))

#' Extract the spatial roughness operator for a template asset
#'
#' @param template A shared basis asset.
#' @param coordinates Coordinate system for the returned roughness operator.
#' @param ... Additional arguments passed to methods.
#' @return Matrix-like roughness operator or \code{NULL}.
#' @export
setGeneric("template_roughness", function(template, coordinates = c("analysis", "raw"), ...) {
  standardGeneric("template_roughness")
})

#' Project data onto a template
#'
#' @param x A template object.
#' @param data Numeric matrix (time x voxels in mask order).
#' @param ... Additional arguments passed to methods.
#' @return A list with \code{coefficients} (time x atoms/components) and
#'   \code{offset} (voxel means or \code{numeric(0)}).
#' @export
setGeneric("template_project", function(x, data, ...) standardGeneric("template_project"))

#' Save a template object to disk
#'
#' @param template A template object.
#' @param file Output path.
#' @param compress Compression passed to \code{saveRDS()}.
#' @param ... Additional arguments passed to methods.
#' @return Normalized output path, invisibly.
#' @export
setGeneric("save_template", function(template, file, compress = "xz", ...) {
  standardGeneric("save_template")
})

#' Extract the decoded domain associated with a latent object
#'
#' @param x A latent object.
#' @param ... Additional arguments passed to methods.
#' @return Domain object or identifier.
#' @export
setGeneric("latent_domain", function(x, ...) standardGeneric("latent_domain"))

#' Extract the decoded support associated with a latent object
#'
#' @param x A latent object.
#' @param ... Additional arguments passed to methods.
#' @return Support object for the latent representation.
#' @export
setGeneric("latent_support", function(x, ...) standardGeneric("latent_support"))

#' Wrap flat decoded outputs into a domain-native representation
#'
#' @param x A latent object.
#' @param values Flat decoded values, typically as a vector or matrix.
#' @param time_idx Optional integer time indices associated with \code{values}.
#' @param space Output space to wrap into.
#' @param ... Additional arguments passed to methods.
#' @return Domain-native wrapped representation.
#' @export
setGeneric("wrap_decoded", function(x, values, time_idx = NULL,
                                     space = c("native", "template"), ...) {
  standardGeneric("wrap_decoded")
})
