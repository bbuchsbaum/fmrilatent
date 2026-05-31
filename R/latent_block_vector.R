#' @include all_generic.R transport_latent.R latent_shared_validation.R
#' @importFrom methods setMethod setValidity new
#' @importFrom Matrix Matrix
NULL

.normalize_block_latent_blocks <- function(blocks, context = "BlockLatentNeuroVector") {
  if (!is.list(blocks) || length(blocks) == 0L) {
    .encoder_cli_abort(paste0(context, " requires a non-empty list of latent blocks."),
                       class = "fmrilatent_error_value", call = rlang::caller_env())
  }
  nm <- names(blocks)
  if (is.null(nm) || any(!nzchar(nm))) {
    names(blocks) <- paste0("block", seq_along(blocks))
  }
  blocks
}

.validate_block_latent_blocks <- function(blocks, context = "BlockLatentNeuroVector",
                                          compare_values = TRUE) {
  blocks <- .normalize_block_latent_blocks(blocks, context = context)
  errs <- character()

  # First pass: ensure each block is an explicit latent. Constructors compare
  # materialized basis values once; S4 validity only compares dimensions so
  # validObject() does not trigger expensive handle materialization.
  bases <- list()
  basis_dims <- list()
  labels <- character()
  for (nm in names(blocks)) {
    block <- blocks[[nm]]
    ok_explicit <- FALSE
    explicit_try <- tryCatch(
      is_explicit_latent(block),
      error = function(e) e
    )
    if (!inherits(explicit_try, "error")) {
      ok_explicit <- isTRUE(explicit_try)
    }
    if (!ok_explicit) {
      dispatch_msg <- if (inherits(explicit_try, "error")) {
        paste0(" is_explicit_latent() dispatch failed: ", conditionMessage(explicit_try))
      } else {
        ""
      }
      errs <- c(errs, paste0("Block '", nm, "' must be an explicit latent object.", dispatch_msg))
      next
    }

    if (isTRUE(compare_values)) {
      block_basis <- try(as.matrix(basis(block)), silent = TRUE)
      if (inherits(block_basis, "try-error")) {
        errs <- c(errs, paste0("Block '", nm, "' does not expose a basis() matrix."))
        next
      }
      bases[[length(bases) + 1L]] <- block_basis
    } else {
      block_dim <- try(.explicit_latent_basis_dim(block), silent = TRUE)
      if (inherits(block_dim, "try-error")) {
        errs <- c(errs, paste0("Block '", nm, "' does not expose basis dimensions."))
        next
      }
      basis_dims[[length(basis_dims) + 1L]] <- block_dim
    }
    labels <- c(labels, nm)
  }

  basis_errs <- if (isTRUE(compare_values)) {
    .validate_shared_basis(
      bases,
      labels = labels,
      tolerance = 1e-8,
      dim_msg = "Block '%s' basis dimensions do not match the reference block.",
      value_msg = "Block '%s' basis matrix must match the shared block basis."
    )
  } else {
    .validate_shared_basis_dims(
      basis_dims,
      labels = labels,
      dim_msg = "Block '%s' basis dimensions do not match the reference block."
    )
  }
  if (!isTRUE(basis_errs)) {
    errs <- c(errs, basis_errs)
  }

  if (length(errs) == 0L) TRUE else errs
}

.block_latent_sizes <- function(x) {
  vapply(x@blocks, function(block) nrow(as.matrix(loadings(block))), integer(1))
}

.block_latent_split_values <- function(x, values, context = "wrap_decoded.BlockLatentNeuroVector") {
  sizes <- .block_latent_sizes(x)
  total <- sum(sizes)
  block_names <- names(x@blocks)

  if (is.atomic(values) && is.null(dim(values))) {
    if (length(values) != total) {
      .encoder_cli_abort(
        paste0(context, " vector length does not match total block support cardinality."),
        class = "fmrilatent_error_dim", call = rlang::caller_env())
    }
    starts <- cumsum(c(1L, utils::head(sizes, -1L)))
    ends <- cumsum(sizes)
    out <- stats::setNames(vector("list", length(sizes)), block_names)
    for (i in seq_along(sizes)) {
      out[[i]] <- values[seq.int(starts[i], ends[i])]
    }
    return(out)
  }

  values <- as.matrix(values)
  if (ncol(values) != total) {
    .encoder_cli_abort(
      paste0(context, " matrix column count does not match total block support cardinality."),
      class = "fmrilatent_error_dim", call = rlang::caller_env())
  }
  starts <- cumsum(c(1L, utils::head(sizes, -1L)))
  ends <- cumsum(sizes)
  out <- stats::setNames(vector("list", length(sizes)), block_names)
  for (i in seq_along(sizes)) {
    out[[i]] <- values[, seq.int(starts[i], ends[i]), drop = FALSE]
  }
  out
}

#' Construct a block-domain latent object
#'
#' @param blocks Named list of explicit latent objects sharing a common basis.
#' @param label Optional label.
#' @param meta Optional metadata list.
#' @return A \code{BlockLatentNeuroVector}.
#' @export
BlockLatentNeuroVector <- function(blocks, label = "", meta = list()) {
  blocks <- .normalize_block_latent_blocks(blocks)
  validity <- .validate_block_latent_blocks(blocks)
  if (!isTRUE(validity)) {
    .encoder_cli_abort(paste(validity, collapse = "\n"),
                       class = "fmrilatent_error_value", call = rlang::caller_env())
  }
  if (!is.list(meta)) {
    .encoder_cli_abort("'meta' must be a list.",
                       class = "fmrilatent_error_type", call = rlang::caller_env())
  }
  new("BlockLatentNeuroVector",
      blocks = blocks,
      label = label,
      meta = meta)
}

#' @export
#' @rdname latent_meta
setMethod("latent_meta", "BlockLatentNeuroVector",
          function(x, ...) {
            utils::modifyList(
              list(
                family = "block_explicit",
                blocks = lapply(x@blocks, latent_meta)
              ),
              x@meta %||% list()
            )
          })

#' @export
#' @rdname latent_domain
setMethod("latent_domain", "BlockLatentNeuroVector",
          function(x, ...) {
            structure(
              lapply(x@blocks, latent_domain),
              class = "BlockLatentDomain"
            )
          })

#' @export
#' @rdname latent_support
setMethod("latent_support", "BlockLatentNeuroVector",
          function(x, ...) lapply(x@blocks, latent_support))

#' @export
#' @rdname basis-methods
setMethod("basis", "BlockLatentNeuroVector", function(x) basis(x@blocks[[1L]]))

#' @export
#' @rdname loadings-methods
setMethod("loadings", "BlockLatentNeuroVector",
          function(x) {
            do.call(rbind, lapply(x@blocks, loadings))
          })

#' @export
#' @rdname offset-methods
setMethod("offset", "BlockLatentNeuroVector",
          function(object, ...) unlist(lapply(object@blocks, offset), use.names = FALSE))

#' @export
#' @rdname reconstruct_matrix
setMethod("reconstruct_matrix", "BlockLatentNeuroVector",
          function(x, time_idx = NULL, roi_mask = NULL, ...) {
            block_names <- names(x@blocks)
            if (is.null(roi_mask)) {
              roi_list <- stats::setNames(vector("list", length(x@blocks)), block_names)
            } else if (is.list(roi_mask)) {
              roi_list <- roi_mask
              roi_names <- names(roi_list)
              if (is.null(roi_names) || any(!nzchar(roi_names))) {
                if (length(roi_list) != length(x@blocks)) {
                  .encoder_cli_abort(
                    "roi_mask list must be named by block or have one entry per block.",
                    class = "fmrilatent_error_invalid_argument"
                  )
                }
                names(roi_list) <- block_names
              }
              unknown <- setdiff(names(roi_list), block_names)
              if (length(unknown) > 0L) {
                .encoder_cli_abort(
                  paste0("roi_mask contains unknown block names: ",
                         paste(unknown, collapse = ", ")),
                  class = "fmrilatent_error_invalid_argument"
                )
              }
            } else if (length(x@blocks) == 1L) {
              roi_list <- stats::setNames(list(roi_mask), block_names)
            } else {
              .encoder_cli_abort(
                "roi_mask must be a named list for multi-block latent vectors.",
                class = "fmrilatent_error_invalid_argument"
              )
            }
            pieces <- Map(
              function(block, nm) {
                reconstruct_matrix(block,
                                   time_idx = time_idx,
                                   roi_mask = roi_list[[nm]] %||% NULL,
                                   ...)
              },
              x@blocks,
              names(x@blocks)
            )
            do.call(cbind, pieces)
          })

#' @export
#' @rdname reconstruct_array
setMethod("reconstruct_array", "BlockLatentNeuroVector",
          function(x, time_idx = NULL, roi_mask = NULL, ...) {
            .encoder_cli_abort(
              paste0("reconstruct_array() is not defined for block-domain latent objects. ",
                     "Use reconstruct_matrix() plus wrap_decoded()."),
              class = "fmrilatent_error_unsupported_operation", call = rlang::caller_env())
          })

#' @export
#' @rdname wrap_decoded
setMethod("wrap_decoded", "BlockLatentNeuroVector",
          function(x, values, time_idx = NULL, space = c("native", "template"), ...) {
            space <- match.arg(space)
            if (space != "native") {
              .encoder_cli_abort(
                "wrap_decoded() for BlockLatentNeuroVector currently supports only native-space wrapping.",
                class = "fmrilatent_error_unsupported_operation", call = rlang::caller_env())
            }
            split_vals <- .block_latent_split_values(x, values)
            out <- stats::setNames(vector("list", length(x@blocks)), names(x@blocks))
            for (nm in names(x@blocks)) {
              out[[nm]] <- wrap_decoded(x@blocks[[nm]], split_vals[[nm]], space = "native", ...)
            }
            class(out) <- c("BlockDecodedLatent", "list")
            out
          })

#' @export
#' @rdname as.matrix-LatentNeuroVec-method
setMethod("as.matrix", "BlockLatentNeuroVector",
          function(x, ...) reconstruct_matrix(x, ...))

#' @export
#' @rdname show-methods
setMethod("show", "BlockLatentNeuroVector",
          function(object) {
            cat("BlockLatentNeuroVector\n")
            cat("  Blocks:", paste(names(object@blocks), collapse = ", "), "\n")
            cat("  Time points:", nrow(basis(object)), "\n")
            cat("  Components:", ncol(basis(object)), "\n")
            sizes <- .block_latent_sizes(object)
            for (i in seq_along(sizes)) {
              cat("  ", names(sizes)[i], " support:", sizes[[i]], "\n", sep = "")
            }
            invisible(object)
          })

#' @export
#' @rdname coef_time
setMethod("coef_time", "BlockLatentNeuroVector",
          function(x, coordinates = c("analysis", "raw"), ...) as.matrix(basis(x)))

#' @export
#' @rdname coef_metric
setMethod("coef_metric", "BlockLatentNeuroVector",
          function(x, coordinates = c("analysis", "raw"), ...) diag(ncol(basis(x))))

#' @export
#' @rdname analysis_transform
setMethod("analysis_transform", "BlockLatentNeuroVector",
          function(x, ...) .transport_identity_transform(ncol(basis(x))))

#' @export
#' @rdname basis_asset
setMethod("basis_asset", "BlockLatentNeuroVector", function(x, ...) NULL)

#' @export
#' @rdname decoder
setMethod("decoder", "BlockLatentNeuroVector",
          function(x, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), ...) {
            space <- match.arg(space)
            coordinates <- match.arg(coordinates)
            if (space == "template") {
              .encoder_cli_warn(
                "BlockLatentNeuroVector has no separate template domain; returning the stacked block decoder.",
                class = "fmrilatent_warning_no_template_domain", call = rlang::caller_env())
            }
            .latent_loadings_map(x)
          })

#' @export
#' @rdname decode_coefficients
setMethod("decode_coefficients", "BlockLatentNeuroVector",
          function(x, gamma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"),
                   wrap = c("none", "auto"), ...) {
            .decode_coefficients_via_decoder(x, gamma, space = space,
                                             coordinates = coordinates,
                                             wrap = wrap, ...)
          })

#' @export
#' @rdname decode_covariance
setMethod("decode_covariance", "BlockLatentNeuroVector",
          function(x, Sigma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), diag_only = TRUE, ...) {
            map <- decoder(x, space = space, coordinates = coordinates, ...)
            Sigma <- .as_square_matrix(Sigma, map$n_source, context = "Sigma")
            if (isTRUE(diag_only)) {
              .project_covariance_diag(map, Sigma)
            } else {
              D <- .materialize_linear_map(map)
              D %*% Sigma %*% t(D)
            }
          })

#' @export
#' @rdname project_effect
setMethod("project_effect", "BlockLatentNeuroVector",
          function(x, gamma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), ...) {
            decode_coefficients(x, gamma, space = space, coordinates = coordinates, ...)
          })

#' @export
#' @rdname project_vcov
setMethod("project_vcov", "BlockLatentNeuroVector",
          function(x, Sigma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), diag_only = TRUE, ...) {
            decode_covariance(x, Sigma, space = space, coordinates = coordinates,
                              diag_only = diag_only, ...)
          })

.validate_BlockLatentNeuroVector <- function(object) {
  .validate_block_latent_blocks(object@blocks, compare_values = FALSE)
}

setValidity("BlockLatentNeuroVector", .validate_BlockLatentNeuroVector)
