# Transport-backed implicit latent helpers and protocol methods

#' @include implicit_latent.R
#' @importFrom methods setMethod
NULL

#' Portable linear-map contract
#'
#' The portable linear map is the seam `fmrilatent` uses to consume subject
#' field operators from neurofunctor (or any other producer) without pulling
#' in file-format-specific warp code. A portable linear map represents a
#' single linear operator \eqn{A : \mathbb{R}^{n_\text{source}} \to
#' \mathbb{R}^{n_\text{target}}} together with enough metadata to compose,
#' adjoint, and (optionally) materialize it.
#'
#' An object satisfies the contract if it is either
#' \enumerate{
#'   \item a base \code{matrix} or \code{Matrix::Matrix}, or
#'   \item a \code{list} with the following fields:
#'   \describe{
#'     \item{\code{n_source}, \code{n_target}}{Positive integer scalar
#'       dimensions of the source and target sample vectors.}
#'     \item{\code{forward(x, ...)}}{Function applying the operator to
#'       source-space data. Accepts a vector of length \code{n_source} or a
#'       matrix with \code{n_source} rows; returns a vector/matrix with
#'       \code{n_target} rows.}
#'     \item{\code{adjoint_apply(y, ...)} (alias \code{adjoint})}{Function
#'       applying the discrete adjoint (or a declared alternative) to
#'       target-space data.}
#'     \item{\code{source_domain_id}, \code{target_domain_id}}{Stable
#'       character identifiers for the source and target domains. Used for
#'       composition safety and provenance digests. \code{.compose_linear_maps()}
#'       requires matching domain ids across the junction unless at least one
#'       side is an empty string ("unspecified").}
#'     \item{\code{source_support}, \code{target_support}}{Optional
#'       descriptors of the sample layout on each side (for example a
#'       \code{LogicalNeuroVol}, a 3D logical mask, surface vertex indices, or
#'       a vector of sample ids). When present these let \code{fmrilatent}
#'       reconstruct domain-aware outputs; when absent it falls back to
#'       \code{provenance$target_support} / \code{provenance$target_mask}.}
#'     \item{\code{adjoint_convention}}{Character tag identifying the
#'       adjoint convention. Defaults to \code{"euclidean_discrete"} (the
#'       discrete transpose under the standard Euclidean inner product).
#'       Covariance pushforward paths such as \code{decode_covariance()} and
#'       \code{.project_covariance_diag()} currently require this value;
#'       operators using another convention must be re-normalized before
#'       being handed to the pushforward code.}
#'     \item{\code{provenance}}{Optional named list of provenance metadata.
#'       \code{fmrilatent} interprets a small set of reserved keys:
#'       \describe{
#'         \item{\code{target_mask}, \code{target_support}, \code{target_domain}}{
#'           Legacy fallback locations read by
#'           \code{.resolve_transport_target_support()} when the top-level
#'           fields are absent. New producers should set the top-level
#'           fields instead.}
#'         \item{\code{source_support}}{Legacy fallback for the top-level
#'           field of the same name.}
#'         \item{\code{coordinates}}{Appended with value \code{"raw"} when
#'           the map has been wrapped with
#'           \code{.transform_linear_map_coordinates()}.}
#'       }
#'       All other keys are advisory and carried through composition
#'       verbatim. Producers should prefer top-level fields over provenance
#'       entries for anything fmrilatent actually reads.}
#'     \item{\code{materialize(...)} (optional)}{Returns the dense
#'       \code{n_target x n_source} matrix representation. Used only when a
#'       caller explicitly asks for it; the main encode and decode paths are
#'       matrix-free and only call \code{forward}/\code{adjoint_apply}.}
#'     \item{\code{contract_version}}{Integer scalar identifying the
#'       contract revision. Auto-filled by \code{as_portable_linear_map()}
#'       to \code{1L} in this release. Future revisions will bump this tag
#'       so producers and consumers can negotiate compatibility.}
#'   }
#' }
#'
#' \strong{Producer-form vs canonical form on \code{TransportLatent}.}
#' \code{transport_latent()} records the raw producer object at
#' \code{x$field_operator} (and its legacy alias \code{x$observation_operator})
#' as a back-reference for debugging and round-trip. The canonical
#' portable-linear-map representation used internally by all encode/decode
#' math is at \code{x$transport$field_operator}. Consumers that need to
#' introspect the operator should always read the canonical form.
#'
#' Producers (such as \code{neurofunctor::compile_observation_operator}) can
#' return either a bare callback list or a wrapped object. Call
#' \code{as_portable_linear_map()} to coerce any accepted input to the
#' canonical contract and \code{validate_portable_linear_map()} to verify
#' that an input already satisfies it without coercing.
#'
#' @name portable_linear_map
NULL

#' Coerce an operator to the portable linear-map contract
#'
#' Normalizes an incoming object to the canonical
#' \code{\link{portable_linear_map}} representation used internally by
#' \code{\link{encode_operator}()}, \code{\link{decoder}()}, and related
#' routines. Matrices and \code{Matrix::Matrix} objects are wrapped with
#' closure-based \code{forward}/\code{adjoint_apply} callbacks; list inputs
#' are validated and re-normalized with top-level \code{source_support} /
#' \code{target_support} preferred over \code{provenance$*} fallbacks.
#'
#' @param x An input operator. Either a base \code{matrix},
#'   \code{Matrix::Matrix}, or a list satisfying the
#'   \code{\link{portable_linear_map}} contract.
#' @param source_domain_id,target_domain_id Optional character tags applied
#'   only when \code{x} is a matrix or \code{Matrix::Matrix}. Ignored for
#'   list inputs (which carry their own).
#' @param source_support,target_support Optional support descriptors applied
#'   only when \code{x} is a matrix or \code{Matrix::Matrix}. Typical values
#'   are a \code{LogicalNeuroVol}, a 3D logical mask, integer vertex
#'   indices, or \code{NULL}. Ignored for list inputs.
#' @param provenance Optional provenance list applied only when \code{x} is
#'   a matrix or \code{Matrix::Matrix}. Ignored for list inputs.
#' @param context Optional label used in error messages.
#' @return A list in the canonical portable-linear-map form.
#' @seealso \code{\link{portable_linear_map}},
#'   \code{\link{validate_portable_linear_map}}.
#' @examples
#' m <- matrix(1:6, nrow = 2)
#' op <- as_portable_linear_map(m, target_domain_id = "native")
#' op$forward(c(1, 1, 1))
#' @export
as_portable_linear_map <- function(x,
                                   source_domain_id = "",
                                   target_domain_id = "",
                                   source_support = NULL,
                                   target_support = NULL,
                                   provenance = list(),
                                   context = "portable linear map") {
  if (is.matrix(x) || inherits(x, "Matrix")) {
    built <- .linear_map_from_matrix(
      as.matrix(x),
      source_domain_id = source_domain_id,
      target_domain_id = target_domain_id,
      source_support = source_support,
      target_support = target_support,
      provenance = provenance
    )
    return(.validate_linear_map_contract(built, context = context))
  }
  .normalize_linear_map(x, context = context)
}

#' Validate an object against the portable linear-map contract
#'
#' Returns the normalized canonical form invisibly when \code{x} satisfies
#' the \code{\link{portable_linear_map}} contract. The canonical return
#' value lets callers use the result directly rather than paying
#' normalization a second time via \code{\link{as_portable_linear_map}()}.
#' Set \code{error = FALSE} to receive \code{FALSE} on failure instead of
#' a stop.
#'
#' @param x Object to validate.
#' @param context Optional label used in error messages.
#' @param error If \code{TRUE} (the default), signal an error on failure.
#'   If \code{FALSE}, return \code{FALSE} silently on failure.
#' @return Invisibly: the normalized portable-linear-map list on success,
#'   or \code{FALSE} when \code{error = FALSE} and validation fails. On
#'   failure with \code{error = TRUE} (the default), an error is raised.
#' @seealso \code{\link{portable_linear_map}},
#'   \code{\link{as_portable_linear_map}}.
#' @examples
#' op <- validate_portable_linear_map(matrix(1:4, 2, 2))
#' op$n_source
#' op$n_target
#' validate_portable_linear_map(list(forward = function(x) x), error = FALSE)
#' @export
validate_portable_linear_map <- function(x, context = "portable linear map",
                                         error = TRUE) {
  res <- tryCatch(.normalize_linear_map(x, context = context),
                  error = function(e) e)
  if (inherits(res, "error")) {
    if (isTRUE(error)) stop(conditionMessage(res), call. = FALSE)
    return(invisible(FALSE))
  }
  invisible(res)
}

.transport_identity_transform <- function(k) {
  list(
    type = "identity",
    dim = as.integer(k),
    to_analysis = function(x) x,
    to_raw = function(x) x,
    matrix = diag(k)
  )
}

.transport_raw_metric <- function(transform, k) {
  mat <- transform$matrix %||% NULL
  if (is.null(mat)) {
    return(NULL)
  }
  mat <- as.matrix(mat)
  if (!identical(dim(mat), c(k, k))) {
    stop("analysis transform matrix must have dimensions ", k, "x", k, ".",
         call. = FALSE)
  }
  crossprod(mat)
}

.is_transport_latent <- function(x) {
  inherits(x, "TransportLatent") ||
    identical(x$meta$family %||% "", "transport") ||
    grepl("^transport", x$meta$family %||% "")
}

.transport_vector_or_matrix_input <- function(data, expected_n, context) {
  input_was_vector <- is.atomic(data) && is.null(dim(data))
  if (input_was_vector) {
    data <- matrix(data, ncol = 1L)
  } else {
    data <- as.matrix(data)
  }
  if (nrow(data) != expected_n) {
    stop(context, " expects ", expected_n, " rows but got ", nrow(data), ".", call. = FALSE)
  }
  list(data = data, input_was_vector = input_was_vector)
}

.transport_finalize_output <- function(result, input_was_vector) {
  if (input_was_vector && is.matrix(result) && ncol(result) == 1L) {
    as.vector(result)
  } else {
    result
  }
}

# Current version of the portable linear-map contract.
.PORTABLE_LINEAR_MAP_VERSION <- 1L

.linear_map_from_matrix <- function(mat, source_domain_id = "", target_domain_id = "",
                                    source_support = NULL, target_support = NULL,
                                    provenance = list(), adjoint_convention = "euclidean_discrete") {
  mat <- as.matrix(mat)
  list(
    mode = "matrix",
    contract_version = .PORTABLE_LINEAR_MAP_VERSION,
    n_source = ncol(mat),
    n_target = nrow(mat),
    source_domain_id = source_domain_id,
    target_domain_id = target_domain_id,
    source_support = source_support,
    target_support = target_support,
    adjoint_convention = adjoint_convention,
    provenance = provenance,
    forward = function(data, ...) {
      prep <- .transport_vector_or_matrix_input(data, ncol(mat), "linear map")
      result <- mat %*% prep$data
      .transport_finalize_output(result, prep$input_was_vector)
    },
    adjoint_apply = function(data, ...) {
      prep <- .transport_vector_or_matrix_input(data, nrow(mat), "linear map adjoint")
      result <- t(mat) %*% prep$data
      .transport_finalize_output(result, prep$input_was_vector)
    },
    materialize = function(...) mat
  )
}

.require_euclidean_adjoint <- function(map, context = "adjoint convention") {
  conv <- map$adjoint_convention %||% "euclidean_discrete"
  if (!identical(conv, "euclidean_discrete")) {
    stop(context, " requires an operator with adjoint_convention = ",
         "'euclidean_discrete' (the Euclidean discrete transpose). ",
         "Got '", conv, "'. The covariance pushforward math depends on ",
         "this convention; supply an operator that normalizes its ",
         "adjoint to the Euclidean inner product, or override ",
         "adjoint_convention explicitly if the current value is ",
         "incorrectly set.",
         call. = FALSE)
  }
  invisible(map)
}

.validate_linear_map_contract <- function(map, context = "linear map") {
  for (nm in c("n_source", "n_target")) {
    value <- map[[nm]]
    if (length(value) != 1L || !is.finite(value) || value <= 0) {
      stop(context, " must define a positive scalar ", nm, ".", call. = FALSE)
    }
    map[[nm]] <- as.integer(value)
  }
  if (!is.function(map$forward)) {
    stop(context, " must define forward().", call. = FALSE)
  }
  if (!is.function(map$adjoint_apply)) {
    stop(context, " must define adjoint_apply().", call. = FALSE)
  }
  map
}

.build_normalized_linear_map <- function(map, forward_fn, adjoint_fn) {
  prov <- map$provenance %||% list()
  list(
    mode = map$mode %||% "callbacks",
    contract_version = map$contract_version %||% .PORTABLE_LINEAR_MAP_VERSION,
    n_source = as.integer(map$n_source),
    n_target = as.integer(map$n_target),
    source_domain_id = map$source_domain_id %||% "",
    target_domain_id = map$target_domain_id %||% "",
    source_support = map$source_support %||% prov$source_support %||% NULL,
    target_support = map$target_support %||% prov$target_support %||% NULL,
    adjoint_convention = map$adjoint_convention %||% "euclidean_discrete",
    provenance = prov,
    forward = forward_fn,
    adjoint_apply = adjoint_fn,
    materialize = map$materialize %||% NULL
  )
}

.normalize_linear_map <- function(map, context = "linear map") {
  if (is.matrix(map) || inherits(map, "Matrix")) {
    return(.validate_linear_map_contract(.linear_map_from_matrix(as.matrix(map)), context = context))
  }

  if (is.list(map)) {
    # Accept `adjoint` as a legacy alias for `adjoint_apply`: promote it
    # before the presence check so the alias path remains reachable.
    if (is.null(map$adjoint_apply) && !is.null(map$adjoint)) {
      map$adjoint_apply <- map$adjoint
    }
    required <- c("forward", "adjoint_apply", "n_source", "n_target")
    present <- vapply(required, function(nm) !is.null(map[[nm]]), logical(1))
    missing <- required[!present]
    if (any(present) && length(missing) > 0L) {
      stop(context, " is missing required fields: ", paste(missing, collapse = ", "), ".",
           call. = FALSE)
    }
    if (length(missing) == 0L) {
      return(.validate_linear_map_contract(
        .build_normalized_linear_map(map,
                                     forward_fn = map$forward,
                                     adjoint_fn = map$adjoint_apply),
        context = context
      ))
    }
  }

  stop("Unsupported ", context, " contract.", call. = FALSE)
}

.materialize_linear_map <- function(map) {
  map <- .normalize_linear_map(map)
  if (is.function(map$materialize)) {
    return(as.matrix(map$materialize()))
  }
  eye <- diag(map$n_source)
  as.matrix(map$forward(eye))
}

.compose_linear_maps <- function(first, second, context = "composed linear map",
                                  strict = getOption("fmrilatent.strict_compose", FALSE)) {
  first <- .normalize_linear_map(first, context = paste(context, "first map"))
  second <- .normalize_linear_map(second, context = paste(context, "second map"))
  if (first$n_target != second$n_source) {
    stop(context, " dimensions do not align: ", first$n_target, " vs ", second$n_source, ".",
         call. = FALSE)
  }

  # Domain-id compatibility: when both sides carry a non-empty id and the
  # caller opts in to strict mode, warn on mismatches. The check is opt-in
  # because fmrilatent's own basis_decoder methods use content digests for
  # target_domain_id while external field operators typically use readable
  # names; those two conventions live in different namespaces and are not
  # meaningfully comparable. Users who control both sides of a composition
  # (e.g. bench tests) can pass `strict = TRUE` or set
  # `options(fmrilatent.strict_compose = TRUE)` to enable the check.
  if (isTRUE(strict)) {
    first_out <- first$target_domain_id %||% ""
    second_in <- second$source_domain_id %||% ""
    if (nzchar(first_out) && nzchar(second_in) && !identical(first_out, second_in)) {
      warning(context, " domain ids do not match: first$target_domain_id = '",
              first_out, "' but second$source_domain_id = '", second_in,
              "'. Silencing this warning requires either matching ids ",
              "or setting one side to \"\" (unspecified).",
              call. = FALSE)
    }
  }

  # Adjoint conventions: warn on mismatched non-euclidean conventions so the
  # user has a signal at composition time. The actual math gate lives in
  # .project_covariance_diag(), which errors at the point of use when the
  # composed map is not euclidean_discrete.
  first_conv <- first$adjoint_convention %||% "euclidean_discrete"
  second_conv <- second$adjoint_convention %||% "euclidean_discrete"
  if (!identical(first_conv, second_conv) &&
      !identical(first_conv, "euclidean_discrete") &&
      !identical(second_conv, "euclidean_discrete")) {
    warning(context, " has incompatible adjoint_conventions: '", first_conv,
            "' and '", second_conv,
            "'. decode_covariance() will reject this map at use time. ",
            "Re-normalize one side to euclidean_discrete before composing ",
            "if you intend to push a covariance through it.",
            call. = FALSE)
  }
  composed_conv <- if (identical(first_conv, "euclidean_discrete")) {
    second_conv
  } else {
    first_conv
  }

  # Composition exposes only the outer source/target supports and domain
  # ids; intermediate supports are intentionally dropped (they are not
  # observable from the composed operator).
  .build_normalized_linear_map(
    list(
      mode = "callbacks",
      contract_version = max(
        first$contract_version %||% .PORTABLE_LINEAR_MAP_VERSION,
        second$contract_version %||% .PORTABLE_LINEAR_MAP_VERSION
      ),
      n_source = first$n_source,
      n_target = second$n_target,
      source_domain_id = first$source_domain_id %||% "",
      target_domain_id = second$target_domain_id %||% "",
      source_support = first$source_support %||% NULL,
      target_support = second$target_support %||% NULL,
      adjoint_convention = composed_conv,
      provenance = c(first$provenance %||% list(), second$provenance %||% list()),
      materialize = function(...) {
        .materialize_linear_map(second) %*% .materialize_linear_map(first)
      }
    ),
    forward_fn = function(data, ...) second$forward(first$forward(data, ...), ...),
    adjoint_fn = function(data, ...) first$adjoint_apply(second$adjoint_apply(data, ...), ...)
  )
}

.normalize_basis_decoder_map <- function(template) {
  .normalize_linear_map(basis_decoder(template), context = "basis decoder")
}

.resolve_field_operator <- function(field_operator = NULL, observation_operator = NULL,
                                    context = "subject field operator") {
  if (is.null(field_operator) && is.null(observation_operator)) {
    stop(context, " is required.", call. = FALSE)
  }
  if (!is.null(field_operator) && !is.null(observation_operator) &&
      !isTRUE(identical(field_operator, observation_operator))) {
    stop("Provide only one of field_operator or observation_operator.", call. = FALSE)
  }
  field_operator %||% observation_operator
}

.normalize_field_operator_map <- function(operator, context = "subject field operator") {
  .normalize_linear_map(operator, context = context)
}

.normalize_observation_operator_map <- function(operator, context = "subject observation operator") {
  .normalize_field_operator_map(operator, context = context)
}

.mask_to_neurovol <- function(mask, location = "unknown function") {
  mask_arr <- .mask_to_array(mask, location)
  if (inherits(mask, "LogicalNeuroVol")) {
    return(list(mask = mask, mask_arr = mask_arr))
  }
  list(
    mask = neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dim(mask_arr))),
    mask_arr = mask_arr
  )
}

.support_cardinality <- function(support, domain = NULL, context = "support") {
  if (inherits(support, "LogicalNeuroVol")) {
    return(sum(as.array(support)))
  }
  if (is.logical(support) && length(dim(support)) == 3L) {
    return(sum(support))
  }
  if (.is_surface_domain(domain)) {
    support_norm <- .normalize_surface_support(support, domain = domain, context = context)
    return(length(support_norm))
  }
  if (is.logical(support) && is.null(dim(support))) {
    return(sum(support))
  }
  if (is.atomic(support) && is.null(dim(support))) {
    return(length(support))
  }
  stop(context, " cardinality is not defined for this support object.", call. = FALSE)
}

.subset_reconstruction_by_support <- function(rec_mat, support, roi_mask = NULL,
                                              domain = NULL, context = "support subset") {
  if (is.null(roi_mask)) {
    return(rec_mat)
  }
  if (inherits(support, "LogicalNeuroVol") || (is.logical(support) && length(dim(support)) == 3L)) {
    mask_arr <- if (inherits(support, "LogicalNeuroVol")) as.array(support) else support
    return(roi_subset_columns(rec_mat, mask_arr, roi_mask))
  }
  support_idx <- if (.is_surface_domain(domain)) {
    .normalize_surface_support(support, domain = domain, context = context)
  } else {
    as.integer(support)
  }
  if (is.logical(roi_mask) && is.null(dim(roi_mask))) {
    keep <- if (length(roi_mask) == length(support_idx)) {
      roi_mask
    } else if (.is_surface_domain(domain) && length(roi_mask) == length(neurosurf::nodes(domain))) {
      roi_mask[support_idx]
    } else {
      stop(context, " logical roi_mask has incompatible length.", call. = FALSE)
    }
  } else {
    roi_idx <- as.integer(roi_mask)
    keep <- support_idx %in% roi_idx
  }
  rec_mat[, keep, drop = FALSE]
}

.resolve_transport_target_support <- function(mask = NULL, domain = NULL, support = NULL,
                                              field_operator,
                                              location = "transport-backed encoding") {
  op_map <- .normalize_field_operator_map(
    field_operator,
    context = paste(location, "field operator")
  )

  if (!is.null(mask)) {
    out <- .mask_to_neurovol(mask, location)
    support_out <- out$mask
    domain_out <- neuroim2::space(out$mask)
    source <- "explicit"
  } else {
    support_use <- support %||% op_map$target_support %||%
      op_map$provenance$target_support %||% NULL
    domain_use <- domain %||% op_map$provenance$target_domain %||% NULL
    target_mask <- op_map$provenance$target_mask %||% NULL

    if (!is.null(support_use)) {
      if (.is_surface_domain(domain_use)) {
        support_out <- .normalize_surface_support(
          support_use,
          domain = domain_use,
          context = paste(location, "target support")
        )
      } else {
        support_out <- support_use
      }
      domain_out <- domain_use
      source <- if (!is.null(support)) "explicit_support" else "operator_provenance_support"
    } else if (!is.null(target_mask)) {
      out <- .mask_to_neurovol(target_mask, location)
      support_out <- out$mask
      domain_out <- neuroim2::space(out$mask)
      source <- "operator_provenance"
    } else {
      stop(location,
           " requires an explicit target mask or field_operator$provenance$target_mask. ",
           "Alternatively, provide an explicit target support/domain.",
           call. = FALSE)
    }
  }

  target_card <- .support_cardinality(support_out, domain = domain_out, context = location)
  support_label <- if (inherits(support_out, "LogicalNeuroVol") ||
    (is.logical(support_out) && length(dim(support_out)) == 3L)) {
    "target mask cardinality "
  } else {
    "target support cardinality "
  }
  if (target_card != op_map$n_target) {
    stop(
      location,
      " ",
      support_label,
      target_card,
      " does not match field operator target dimension ",
      op_map$n_target,
      ".",
      call. = FALSE
    )
  }

  list(
    mask = if (inherits(support_out, "LogicalNeuroVol") ||
      (is.logical(support_out) && length(dim(support_out)) == 3L)) {
      if (inherits(support_out, "LogicalNeuroVol")) support_out else .mask_to_neurovol(support_out, location)$mask
    } else NULL,
    support = support_out,
    domain = domain_out,
    source = source
  )
}

.transform_linear_map_coordinates <- function(map, transform, coordinates = c("analysis", "raw")) {
  coordinates <- match.arg(coordinates)
  map <- .normalize_linear_map(map)
  transform <- transform %||% .transport_identity_transform(map$n_source)

  if (coordinates == "analysis") {
    return(map)
  }

  if (!is.function(transform$to_analysis) || !is.function(transform$to_raw)) {
    stop("analysis transform must provide to_analysis() and to_raw() callables.", call. = FALSE)
  }

  .build_normalized_linear_map(
    list(
      mode = "callbacks",
      contract_version = map$contract_version %||% .PORTABLE_LINEAR_MAP_VERSION,
      n_source = map$n_source,
      n_target = map$n_target,
      source_domain_id = map$source_domain_id %||% "",
      target_domain_id = map$target_domain_id %||% "",
      source_support = map$source_support %||% NULL,
      target_support = map$target_support %||% NULL,
      adjoint_convention = map$adjoint_convention %||% "euclidean_discrete",
      provenance = c(map$provenance %||% list(), list(coordinates = "raw")),
      materialize = function(...) {
        base <- .materialize_linear_map(map)
        analysis_matrix <- transform$matrix %||% diag(map$n_source)
        base %*% analysis_matrix
      }
    ),
    forward_fn = function(data, ...) map$forward(transform$to_analysis(data), ...),
    adjoint_fn = function(data, ...) transform$to_raw(map$adjoint_apply(data, ...))
  )
}

.as_coeff_matrix <- function(x, k, context = "coefficients") {
  if (is.atomic(x) && is.null(dim(x))) {
    x <- matrix(x, ncol = 1L)
  } else {
    x <- as.matrix(x)
  }
  if (nrow(x) != k) {
    stop(context, " must have ", k, " rows but has ", nrow(x), ".", call. = FALSE)
  }
  x
}

.as_square_matrix <- function(x, k, context = "covariance matrix") {
  x <- as.matrix(x)
  if (!identical(dim(x), c(k, k))) {
    stop(context, " must have dimensions ", k, "x", k, ".", call. = FALSE)
  }
  x
}

.covariance_factor <- function(Sigma, tol = 1e-10) {
  Sigma <- 0.5 * (Sigma + t(Sigma))
  eig <- eigen(Sigma, symmetric = TRUE)
  keep <- Re(eig$values) > tol
  if (!any(keep)) {
    return(matrix(0, nrow = nrow(Sigma), ncol = 0L))
  }
  eig$vectors[, keep, drop = FALSE] %*% diag(sqrt(Re(eig$values[keep])), nrow = sum(keep))
}

.project_covariance_diag <- function(map, Sigma, block_size = 32L) {
  map <- .normalize_linear_map(map)
  .require_euclidean_adjoint(map, context = "covariance diagonal pushforward")
  Sigma <- .as_square_matrix(Sigma, map$n_source, context = "Sigma")
  factor <- .covariance_factor(Sigma)
  if (ncol(factor) == 0L) {
    return(rep(0, map$n_target))
  }

  out <- numeric(map$n_target)
  block_size <- max(1L, as.integer(block_size))
  for (start in seq.int(1L, ncol(factor), by = block_size)) {
    stop_idx <- min(start + block_size - 1L, ncol(factor))
    block <- factor[, start:stop_idx, drop = FALSE]
    prep <- .transport_vector_or_matrix_input(block, map$n_source, "covariance diagonal pushforward")
    decoded <- as.matrix(map$forward(prep$data))
    if (nrow(decoded) != map$n_target) {
      stop("covariance diagonal pushforward returned ", nrow(decoded),
           " rows; expected ", map$n_target, ".", call. = FALSE)
    }
    out <- out + rowSums(decoded^2)
  }
  out
}

.latent_meta_field <- function(x, field, default = NULL) {
  meta <- tryCatch(latent_meta(x), error = function(e) list())
  meta[[field]] %||% default
}

.explicit_latent_analysis_transform <- function(x, k) {
  .latent_meta_field(x, "analysis_transform", .transport_identity_transform(k))
}

.explicit_latent_basis_asset <- function(x) {
  .latent_meta_field(x, "basis_asset", NULL)
}

.explicit_coef_time <- function(x, basis_matrix, coordinates = c("analysis", "raw")) {
  coordinates <- match.arg(coordinates)
  if (coordinates == "analysis") {
    return(basis_matrix)
  }
  transform <- .explicit_latent_analysis_transform(x, ncol(basis_matrix))
  t(transform$to_raw(t(basis_matrix)))
}

.explicit_coef_metric <- function(x, k, coordinates = c("raw", "analysis")) {
  coordinates <- match.arg(coordinates)
  transform <- .explicit_latent_analysis_transform(x, k)
  if (coordinates == "analysis") {
    diag(k)
  } else {
    .transport_raw_metric(transform, k)
  }
}

# Shared helpers for decode/project methods across LatentNeuroVec,
# LatentNeuroSurfaceVector, and ImplicitLatent.

.decode_coefficients_via_decoder <- function(x, gamma, space, coordinates,
                                              wrap = c("none", "auto"), ...) {
  wrap <- match.arg(wrap)
  map <- decoder(x, space = space, coordinates = coordinates)
  gamma <- .as_coeff_matrix(gamma, map$n_source, context = "gamma")
  result <- map$forward(gamma)
  values <- if (ncol(result) == 1L) as.vector(result) else result
  if (wrap == "none") {
    return(values)
  }
  if (identical(space, "template")) {
    stop("decode_coefficients(wrap = \"auto\") is only implemented for ",
         "space = \"native\". No wrap_decoded() method supports template ",
         "space yet; call with wrap = \"none\" for raw template-space values.",
         call. = FALSE)
  }
  # decoder forward() returns target-samples x ncol(gamma); wrap_decoded()
  # methods expect time-first (rows = frames), so transpose when gamma has
  # more than one column.
  wrap_values <- if (is.matrix(values)) t(values) else values
  wrap_decoded(x, wrap_values, space = space)
}

.decode_covariance_via_decoder <- function(x, Sigma, space, coordinates, diag_only, ...) {
  map <- decoder(x, space = space, coordinates = coordinates)
  Sigma <- .as_square_matrix(Sigma, map$n_source, context = "Sigma")
  if (isTRUE(diag_only)) {
    .project_covariance_diag(map, Sigma)
  } else {
    D <- .materialize_linear_map(map)
    D %*% Sigma %*% t(D)
  }
}

.latent_loadings_map <- function(x) {
  L <- as.matrix(loadings(x))
  .linear_map_from_matrix(
    L,
    source_domain_id = paste0("latent:", latent_meta(x)$family %||% "explicit"),
    target_domain_id = "stored_spatial_domain",
    provenance = list(
      basis_asset_class = "LatentNeuroVec",
      basis_family = latent_meta(x)$family %||% "explicit",
      basis_id = digest::digest(L)
    )
  )
}

#' Construct a transport-backed implicit latent object
#'
#' @param coeff_raw Raw coefficient matrix (time x k).
#' @param basis_asset Shared basis asset.
#' @param field_operator Subject field operator satisfying the
#'   contract documented in \code{\link{encode_operator}()}.
#' @param observation_operator Legacy alias for \code{field_operator}.
#' @param mask Target-domain mask for volumetric reconstruction.
#' @param domain Optional target domain for non-volumetric reconstruction.
#' @param support Optional target support for non-volumetric reconstruction.
#' @param coeff_analysis Optional analysis-coordinate coefficient matrix (time x k).
#' @param analysis_transform Optional transform descriptor from raw to analysis coordinates.
#' @param offset Optional voxel offset vector in target space.
#' @param run_info Optional run metadata.
#' @param meta Optional metadata list.
#' @details
#' `transport_latent()` stores the `fmrilatent` side of the operator-backed
#' workflow: shared basis asset, subject field operator, raw and
#' analysis-space coefficients, and decoder views for native or template
#' projection. It does not fit GLMs or carry statistical summaries. Downstream
#' model-fitting code should ordinarily consume
#' \code{coef_time(x, coordinates = "analysis")}, then return coefficient-space
#' effects and covariance matrices for projection with
#' \code{decode_coefficients()} and \code{decode_covariance()}.
#' @return A \code{TransportLatent} object inheriting from \code{ImplicitLatent}.
#' @export
transport_latent <- function(coeff_raw, basis_asset, field_operator = NULL, observation_operator = NULL, mask = NULL,
                             domain = NULL, support = NULL,
                             coeff_analysis = NULL, analysis_transform = NULL,
                             offset = numeric(0), run_info = NULL, meta = list()) {
  coeff_raw <- as.matrix(coeff_raw)
  coeff_analysis <- as.matrix(coeff_analysis %||% coeff_raw)
  if (!identical(dim(coeff_raw), dim(coeff_analysis))) {
    stop("coeff_raw and coeff_analysis must have identical dimensions.", call. = FALSE)
  }

  field_operator <- .resolve_field_operator(
    field_operator = field_operator,
    observation_operator = observation_operator,
    context = "transport_latent() field_operator"
  )

  basis_map <- .normalize_basis_decoder_map(basis_asset)
  obs_map <- .normalize_field_operator_map(field_operator)
  if (basis_map$n_target != obs_map$n_source) {
    stop("basis decoder output dimension ", basis_map$n_target,
         " does not match field operator source dimension ", obs_map$n_source, ".",
         call. = FALSE)
  }

  analysis_transform <- analysis_transform %||% .transport_identity_transform(ncol(coeff_analysis))
  # Store unwrapped decoders; coordinate transforms applied on demand in decoder()
  native_decoder <- .compose_linear_maps(basis_map, obs_map, context = "subject decoder")
  template_decoder <- .normalize_linear_map(basis_map, context = "template decoder")
  native_decoder_analysis <- .transform_linear_map_coordinates(
    native_decoder,
    analysis_transform,
    coordinates = "analysis"
  )

  support_use <- support %||% mask
  target_card <- .support_cardinality(support_use, domain = domain, context = "transport_latent")
  if (native_decoder$n_target != target_card) {
    stop("Target support cardinality ", target_card,
         " does not match decoder target dimension ", native_decoder$n_target, ".",
         call. = FALSE)
  }
  if (length(offset) > 0L && length(offset) != native_decoder$n_target) {
    stop("offset length must match decoder target dimension.", call. = FALSE)
  }

  decoder_fn <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL, ...) {
    coeff_mat <- coeff_analysis
    if (!is.null(time_idx)) {
      coeff_mat <- coeff_mat[as.integer(time_idx), , drop = FALSE]
    }
    rec <- t(native_decoder_analysis$forward(t(coeff_mat)))
    if (length(offset) > 0L) {
      rec <- sweep(rec, 2L, offset, "+")
    }
    .subset_reconstruction_by_support(
      rec,
      support = support_use,
      roi_mask = roi_mask,
      domain = domain,
      context = "transport_latent decoder"
    )
  }

  meta_defaults <- list(
    family = "transport",
    coordinate_mode = "analysis_euclidean",
    basis_asset_class = class(basis_asset)[1],
    target_mask_source = "explicit",
    native_domain_id = obs_map$target_domain_id %||% "",
    template_domain_id = basis_map$target_domain_id %||% "",
    operator_adjoint_convention = obs_map$adjoint_convention %||% "euclidean_discrete",
    run_info = run_info %||% list(n_runs = 1L, run_lengths = nrow(coeff_analysis))
  )
  meta <- utils::modifyList(meta_defaults, meta)

  obj <- implicit_latent(
    coeff = list(
      raw = coeff_raw,
      analysis = coeff_analysis,
      offset = offset
    ),
    decoder = decoder_fn,
    meta = meta,
    mask = mask,
    domain = domain,
    support = support_use
  )
  obj$basis_asset <- basis_asset
  obj$field_operator <- field_operator
  obj$observation_operator <- field_operator
  obj$transport <- list(
    basis_decoder = basis_map,
    field_operator = obs_map,
    observation_operator = obs_map,
    native_decoder = native_decoder,
    template_decoder = template_decoder,
    analysis_transform = analysis_transform,
    run_info = run_info %||% list(n_runs = 1L, run_lengths = nrow(coeff_analysis))
  )
  class(obj) <- c("TransportLatent", class(obj))
  obj
}

#' Test whether an object is a transport-backed latent object
#'
#' @param x Object to test.
#' @return Logical scalar.
#' @export
is_transport_latent <- function(x) .is_transport_latent(x)

#' @export
#' @rdname coef_time
setMethod("coef_time", "ImplicitLatent",
          function(x, coordinates = c("analysis", "raw"), ...) {
            if (!.is_transport_latent(x)) {
              stop("coef_time() is only defined for transport-backed ImplicitLatent objects.",
                   call. = FALSE)
            }
            coordinates <- match.arg(coordinates)
            x$coeff[[coordinates]]
          })

#' @export
#' @rdname coef_metric
setMethod("coef_metric", "ImplicitLatent",
          function(x, coordinates = c("raw", "analysis"), ...) {
            if (!.is_transport_latent(x)) {
              return(NULL)
            }
            coordinates <- match.arg(coordinates)
            k <- ncol(x$coeff$analysis)
            if (coordinates == "analysis") {
              diag(k)
            } else {
              .transport_raw_metric(x$transport$analysis_transform, k)
            }
          })

#' @export
#' @rdname analysis_transform
setMethod("analysis_transform", "ImplicitLatent",
          function(x, ...) {
            if (!.is_transport_latent(x)) return(NULL)
            x$transport$analysis_transform
          })

#' @export
#' @rdname basis_asset
setMethod("basis_asset", "ImplicitLatent",
          function(x, ...) {
            if (!.is_transport_latent(x)) return(NULL)
            x$basis_asset %||% NULL
          })

#' @export
#' @rdname decoder
setMethod("decoder", "ImplicitLatent",
          function(x, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), ...) {
            if (!.is_transport_latent(x)) {
              stop("decoder() is only defined for transport-backed ImplicitLatent objects.",
                   call. = FALSE)
            }
            space <- match.arg(space)
            coordinates <- match.arg(coordinates)
            base_map <- switch(
              space,
              native = x$transport$native_decoder,
              template = x$transport$template_decoder
            )
            .transform_linear_map_coordinates(base_map, x$transport$analysis_transform,
                                              coordinates = coordinates)
          })

#' @export
#' @rdname decode_coefficients
setMethod("decode_coefficients", "ImplicitLatent",
          function(x, gamma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"),
                   wrap = c("none", "auto"), ...) {
            if (!.is_transport_latent(x)) {
              stop("decode_coefficients() is only defined for transport-backed ImplicitLatent objects.",
                   call. = FALSE)
            }
            .decode_coefficients_via_decoder(x, gamma, space = space,
                                             coordinates = coordinates,
                                             wrap = wrap, ...)
          })

#' @export
#' @rdname decode_covariance
setMethod("decode_covariance", "ImplicitLatent",
          function(x, Sigma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), diag_only = TRUE, ...) {
            if (!.is_transport_latent(x)) {
              stop("decode_covariance() is only defined for transport-backed ImplicitLatent objects.",
                   call. = FALSE)
            }
            .decode_covariance_via_decoder(x, Sigma, space = space, coordinates = coordinates,
                                           diag_only = diag_only, ...)
          })

#' @export
#' @rdname project_effect
setMethod("project_effect", "ImplicitLatent",
          function(x, gamma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), ...) {
            decode_coefficients(x, gamma, space = space, coordinates = coordinates, ...)
          })

#' @export
#' @rdname project_vcov
setMethod("project_vcov", "ImplicitLatent",
          function(x, Sigma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), diag_only = TRUE, ...) {
            decode_covariance(x, Sigma, space = space, coordinates = coordinates,
                              diag_only = diag_only, ...)
          })

#' @export
#' @rdname coef_time
setMethod("coef_time", "LatentNeuroVec",
          function(x, coordinates = c("analysis", "raw"), ...) {
            .explicit_coef_time(x, as.matrix(basis(x)), coordinates = coordinates)
          })

#' @export
#' @rdname coef_metric
setMethod("coef_metric", "LatentNeuroVec",
          function(x, coordinates = c("raw", "analysis"), ...) {
            .explicit_coef_metric(x, ncol(basis(x)), coordinates = coordinates)
          })

#' @export
#' @rdname analysis_transform
setMethod("analysis_transform", "LatentNeuroVec",
          function(x, ...) .explicit_latent_analysis_transform(x, ncol(basis(x))))

#' @export
#' @rdname basis_asset
setMethod("basis_asset", "LatentNeuroVec", function(x, ...) .explicit_latent_basis_asset(x))

#' @export
#' @rdname decoder
setMethod("decoder", "LatentNeuroSurfaceVector",
          function(x, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), ...) {
            space <- match.arg(space)
            coordinates <- match.arg(coordinates)
            if (space == "template") {
              warning("LatentNeuroSurfaceVector has no separate template domain; returning the stored surface decoder.",
                      call. = FALSE)
            }
            .transform_linear_map_coordinates(
              .latent_loadings_map(x),
              analysis_transform(x),
              coordinates = coordinates
            )
          })

#' @export
#' @rdname decode_coefficients
setMethod("decode_coefficients", "LatentNeuroSurfaceVector",
          function(x, gamma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"),
                   wrap = c("none", "auto"), ...) {
            .decode_coefficients_via_decoder(x, gamma, space = space,
                                             coordinates = coordinates,
                                             wrap = wrap, ...)
          })

#' @export
#' @rdname decode_covariance
setMethod("decode_covariance", "LatentNeuroSurfaceVector",
          function(x, Sigma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), diag_only = TRUE, ...) {
            .decode_covariance_via_decoder(x, Sigma, space = space, coordinates = coordinates,
                                           diag_only = diag_only, ...)
          })

#' @export
#' @rdname project_effect
setMethod("project_effect", "LatentNeuroSurfaceVector",
          function(x, gamma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), ...) {
            decode_coefficients(x, gamma, space = space, coordinates = coordinates, ...)
          })

#' @export
#' @rdname project_vcov
setMethod("project_vcov", "LatentNeuroSurfaceVector",
          function(x, Sigma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), diag_only = TRUE, ...) {
            decode_covariance(x, Sigma, space = space, coordinates = coordinates,
                              diag_only = diag_only, ...)
          })

#' @export
#' @rdname decoder
setMethod("decoder", "LatentNeuroVec",
          function(x, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), ...) {
            space <- match.arg(space)
            coordinates <- match.arg(coordinates)
            if (space == "template") {
              warning("LatentNeuroVec has no separate template domain; returning the stored spatial decoder.",
                      call. = FALSE)
            }
            .transform_linear_map_coordinates(
              .latent_loadings_map(x),
              analysis_transform(x),
              coordinates = coordinates
            )
          })

#' @export
#' @rdname decode_coefficients
setMethod("decode_coefficients", "LatentNeuroVec",
          function(x, gamma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"),
                   wrap = c("none", "auto"), ...) {
            .decode_coefficients_via_decoder(x, gamma, space = space,
                                             coordinates = coordinates,
                                             wrap = wrap, ...)
          })

#' @export
#' @rdname decode_covariance
setMethod("decode_covariance", "LatentNeuroVec",
          function(x, Sigma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), diag_only = TRUE, ...) {
            .decode_covariance_via_decoder(x, Sigma, space = space, coordinates = coordinates,
                                           diag_only = diag_only, ...)
          })

#' @export
#' @rdname project_effect
setMethod("project_effect", "LatentNeuroVec",
          function(x, gamma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), ...) {
            decode_coefficients(x, gamma, space = space, coordinates = coordinates, ...)
          })

#' @export
#' @rdname project_vcov
setMethod("project_vcov", "LatentNeuroVec",
          function(x, Sigma, space = c("native", "template"),
                   coordinates = c("analysis", "raw"), diag_only = TRUE, ...) {
            decode_covariance(x, Sigma, space = space, coordinates = coordinates,
                              diag_only = diag_only, ...)
          })
