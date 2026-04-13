# Surface basis templates and domain-aware surface helpers

#' @include all_generic.R
#' @importFrom Matrix Cholesky Matrix crossprod solve
#' @importFrom methods setMethod
NULL

.require_neurosurf <- function(context = "surface support") {
  if (!requireNamespace("neurosurf", quietly = TRUE)) {
    stop(context, " requires the 'neurosurf' package.", call. = FALSE)
  }
}

#' Build a shared surface basis template
#'
#' @param geometry A \code{neurosurf::SurfaceGeometry} or \code{neurosurf::SurfaceSet}.
#' @param loadings A matrix-like decoder basis with rows aligned to \code{support}
#'   and columns aligned to latent coefficients.
#' @param support Surface support as vertex indices or a logical vector over all
#'   surface nodes. Defaults to the full surface.
#' @param roughness Optional coefficient-space roughness matrix.
#' @param measure Optional support-aligned weighting or mass information.
#' @param ridge Small diagonal ridge added to the Gram matrix if needed.
#' @param label Optional label stored in metadata.
#' @param meta Optional additional metadata.
#' @return A \code{SurfaceBasisTemplate} object.
#' @export
surface_basis_template <- function(geometry, loadings, support = NULL, roughness = NULL,
                                   measure = NULL, ridge = 1e-8,
                                   label = "surface_basis", meta = list()) {
  .require_neurosurf("surface_basis_template")
  if (!(methods::is(geometry, "SurfaceGeometry") || methods::is(geometry, "SurfaceSet"))) {
    stop("geometry must be a neurosurf::SurfaceGeometry or neurosurf::SurfaceSet.",
         call. = FALSE)
  }

  n_nodes <- length(neurosurf::nodes(geometry))
  if (is.null(support)) {
    support <- seq_len(n_nodes)
  }
  support <- .normalize_surface_support(
    support,
    domain = geometry,
    context = "surface_basis_template support"
  )

  loadings <- Matrix::Matrix(as.matrix(loadings), sparse = FALSE)
  if (nrow(loadings) != length(support)) {
    stop("loadings must have ", length(support),
         " rows to match the surface support cardinality.", call. = FALSE)
  }

  if (!is.null(measure)) {
    if (is.atomic(measure) && is.null(dim(measure))) {
      if (length(measure) != length(support)) {
        stop("measure must have length ", length(support),
             " when supplied as a vector.", call. = FALSE)
      }
    } else {
      measure <- as.matrix(measure)
      if (!identical(dim(measure), c(length(support), length(support)))) {
        stop("measure must be either a support-length vector or a ",
             length(support), "x", length(support), " matrix.", call. = FALSE)
      }
    }
  }

  if (!is.null(roughness)) {
    roughness <- Matrix::Matrix(as.matrix(roughness), sparse = FALSE)
    if (!identical(dim(roughness), c(ncol(loadings), ncol(loadings)))) {
      stop("roughness must be a ", ncol(loadings), "x", ncol(loadings),
           " matrix.", call. = FALSE)
    }
  }

  G <- Matrix::crossprod(loadings)
  gram_factor <- tryCatch(
    suppressWarnings(Matrix::Cholesky(G, perm = TRUE)),
    error = function(e) {
      suppressWarnings(Matrix::Cholesky(G + ridge * Matrix::Diagonal(n = ncol(G)), perm = TRUE))
    }
  )

  structure(
    list(
      geometry = geometry,
      support = support,
      loadings = loadings,
      roughness = roughness,
      measure = measure,
      gram_factor = gram_factor,
      meta = utils::modifyList(list(
        family = "surface_basis",
        label = label,
        k = ncol(loadings),
        ridge = ridge
      ), meta)
    ),
    class = "SurfaceBasisTemplate"
  )
}

#' Print method for SurfaceBasisTemplate
#'
#' @param x A \code{SurfaceBasisTemplate}.
#' @param ... Ignored.
#' @export
print.SurfaceBasisTemplate <- function(x, ...) {
  cat("SurfaceBasisTemplate\n")
  cat("  Basis family:", x$meta$family %||% "surface_basis", "\n")
  cat("  Atoms:", ncol(x$loadings), "\n")
  cat("  Support vertices:", length(x$support), "\n")
  cat("  Domain class:", class(x$geometry)[1], "\n")
  invisible(x)
}

#' Test whether an object is a surface basis template
#'
#' @param x Object to test.
#' @return Logical scalar.
#' @export
is_surface_template <- function(x) inherits(x, "SurfaceBasisTemplate")

#' @export
#' @rdname template_loadings
setMethod("template_loadings", "SurfaceBasisTemplate", function(x, ...) x$loadings)

#' @export
#' @rdname template_mask
setMethod("template_mask", "SurfaceBasisTemplate",
          function(x, ...) {
            stop("template_mask() is not defined for SurfaceBasisTemplate.", call. = FALSE)
          })

#' @export
#' @rdname template_meta
setMethod("template_meta", "SurfaceBasisTemplate", function(x, ...) x$meta %||% list())

#' @export
#' @rdname basis_decoder
setMethod("basis_decoder", "SurfaceBasisTemplate",
          function(template, ...) {
            decoder_map <- .template_custom_field(template, "decoder_map") %||%
              template_meta(template)$decoder_map %||% NULL
            if (!is.null(decoder_map)) {
              return(.normalize_linear_map(decoder_map, context = "basis decoder"))
            }
            payload <- .template_coordinate_payload(
              raw_loadings = template_loadings(template),
              measure = template_measure(template),
              analysis_transform = .template_custom_field(template, "analysis_transform") %||%
                template_meta(template)$analysis_transform %||% NULL,
              default_measure = "null"
            )
            B <- payload$analysis_loadings
            .linear_map_from_matrix(
              B,
              source_domain_id = paste0("latent:", template_meta(template)$family %||% "surface"),
              target_domain_id = digest::digest(list(
                support = template_support(template),
                meta = template_meta(template)
              )),
              provenance = list(
                basis_asset_class = "SurfaceBasisTemplate",
                basis_family = template_meta(template)$family %||% "surface_basis",
                basis_id = digest::digest(B),
                target_support = template_support(template)
              )
            )
          })

#' @export
#' @rdname template_rank
setMethod("template_rank", "SurfaceBasisTemplate",
          function(template, ...) ncol(template_loadings(template)))

#' @export
#' @rdname template_domain
setMethod("template_domain", "SurfaceBasisTemplate",
          function(template, ...) template$geometry)

#' @export
#' @rdname template_support
setMethod("template_support", "SurfaceBasisTemplate",
          function(template, ...) template$support)

#' @export
#' @rdname template_measure
setMethod("template_measure", "SurfaceBasisTemplate",
          function(template, ...) template$measure %||% NULL)

#' @export
#' @rdname template_roughness
setMethod("template_roughness", "SurfaceBasisTemplate",
          function(template, coordinates = c("analysis", "raw"), ...) {
            coordinates <- match.arg(coordinates)
            .transform_quadratic_form(
              template$roughness %||% NULL,
              .template_coordinate_payload(
                raw_loadings = template_loadings(template),
                measure = template_measure(template),
                analysis_transform = .template_custom_field(template, "analysis_transform") %||%
                  template_meta(template)$analysis_transform %||% NULL,
                default_measure = "null"
              )$analysis_transform,
              coordinates = coordinates
            )
          })

#' @export
#' @rdname template_project
setMethod("template_project", signature(x = "SurfaceBasisTemplate", data = "ANY"),
          function(x, data, ...) {
            X <- as.matrix(data)
            if (ncol(X) != length(template_support(x))) {
              stop("data must have ", length(template_support(x)),
                   " columns to match the surface support cardinality.", call. = FALSE)
            }
            .template_projection_payload(
              data = X,
              raw_loadings = template_loadings(x),
              measure = template_measure(x),
              center = FALSE,
              analysis_transform = .template_custom_field(x, "analysis_transform") %||%
                template_meta(x)$analysis_transform %||% NULL,
              default_measure = "null"
            )
          })

#' @export
#' @rdname save_template
setMethod("save_template", signature(template = "SurfaceBasisTemplate"),
          function(template, file, compress = TRUE, ...) {
            saveRDS(template, file = file, compress = compress)
            invisible(file)
          })
