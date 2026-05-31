# Shared-structure contracts for reusable latent representations

.shared_check_list <- function(x, name) {
  if (!is.list(x)) {
    .encoder_cli_abort(paste0(name, " must be a list."),
                       class = "fmrilatent_error_invalid_type")
  }
  x
}

.shared_check_id <- function(id, name = "id") {
  if (length(id) != 1L || is.na(id) || !nzchar(id)) {
    .encoder_cli_abort(paste0(name, " must be a non-empty character scalar."),
                       class = "fmrilatent_error_invalid_id")
  }
  as.character(id)
}

.shared_check_positive_integer <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x != as.integer(x) || x < 1L) {
    .encoder_cli_abort(paste0(name, " must be a positive integer scalar."),
                       class = "fmrilatent_error_invalid_argument")
  }
  as.integer(x)
}

.shared_matrix_like <- function(x, name) {
  if (inherits(x, "SharedReference")) {
    x <- resolve_shared_reference(x)
  }
  if (inherits(x, "GroupDeltaLoadings")) {
    x <- materialize_group_delta_loadings(x)
  }
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    .encoder_cli_abort(paste0(name, " must be a matrix-like object."),
                       class = "fmrilatent_error_invalid_type")
  }
  x <- Matrix::Matrix(x, sparse = inherits(x, "sparseMatrix"))
  if (!.matrix_like_all_finite(x)) {
    .encoder_cli_abort(paste0(name, " must contain only finite values."),
                       class = "fmrilatent_error_invalid_value")
  }
  x
}

.shared_support_cardinality <- function(support, domain = NULL, context = "support") {
  if (is.null(support)) {
    return(NULL)
  }
  tryCatch(
    .support_cardinality(support, domain = domain, context = context),
    error = function(e) {
      if (is.atomic(support) && is.null(dim(support))) {
        length(support)
      } else if (is.logical(support)) {
        sum(support)
      } else {
        stop(e)  # re-raise to preserve the classed inner condition
      }
    }
  )
}

.shared_reference_env <- local({
  env <- new.env(parent = emptyenv())
  function(reset = FALSE) {
    if (isTRUE(reset)) {
      rm(list = ls(env, all.names = TRUE), envir = env)
    }
    env
  }
})

#' Create an in-session shared object reference
#'
#' @description
#' `shared_reference()` records a reusable object behind a stable in-session
#' identifier. It deliberately does not create files, archive locators, HDF5
#' datasets, checksums, or persistent registries; those responsibilities belong
#' to archive packages such as neuroarchive.
#'
#' @param value Optional object to cache immediately.
#' @param id Optional reference id. If omitted, a deterministic digest is used
#'   when `value` is present; otherwise a unique session id is generated.
#' @param kind Human-readable object kind.
#' @param materialize Optional zero-argument function used to create the value
#'   on first resolution.
#' @param meta Optional advisory metadata.
#' @param register If `TRUE`, cache `value` in the session registry.
#' @return A `SharedReference` descriptor.
#' @export
shared_reference <- function(value = NULL, id = NULL, kind = "object",
                             materialize = NULL, meta = list(),
                             register = TRUE) {
  if (!is.null(materialize) && !is.function(materialize)) {
    .encoder_cli_abort("materialize must be NULL or a function.",
                       class = "fmrilatent_error_invalid_argument")
  }
  meta <- .shared_check_list(meta, "meta")
  kind <- .shared_check_id(kind, "kind")
  if (is.null(id)) {
    id <- if (!is.null(value)) {
      paste0("shared:", digest::digest(list(kind = kind, value = value)))
    } else {
      paste0("shared:", digest::digest(list(kind = kind, time = Sys.time(), nonce = stats::runif(1))))
    }
  }
  id <- .shared_check_id(id)

  if (!is.null(value) && isTRUE(register)) {
    assign(id, value, envir = .shared_reference_env())
  }

  structure(
    list(
      id = id,
      kind = kind,
      has_cached_value = exists(id, .shared_reference_env(), inherits = FALSE),
      materialize = materialize,
      meta = meta,
      persistence = "session"
    ),
    class = "SharedReference"
  )
}

#' Test whether an object is a shared reference
#'
#' @param x Object to test.
#' @return Logical scalar.
#' @export
is_shared_reference <- function(x) inherits(x, "SharedReference")

#' Resolve an in-session shared reference
#'
#' @param x A `SharedReference`, or any object returned unchanged.
#' @param cache Cache materialized values in the session registry.
#' @return The referenced object.
#' @export
resolve_shared_reference <- function(x, cache = TRUE) {
  if (!inherits(x, "SharedReference")) {
    return(x)
  }
  env <- .shared_reference_env()
  if (exists(x$id, env, inherits = FALSE)) {
    return(get(x$id, envir = env, inherits = FALSE))
  }
  if (!is.function(x$materialize)) {
    .encoder_cli_abort(
      paste0("Shared reference '", x$id, "' has no cached value or materialize() function."),
      class = "fmrilatent_error_unresolved_reference"
    )
  }
  value <- x$materialize()
  if (isTRUE(cache)) {
    assign(x$id, value, envir = env)
  }
  value
}

#' Summarize a shared reference without resolving it
#'
#' @param x A `SharedReference`.
#' @return A list with id, kind, cache state, persistence, and metadata.
#' @export
shared_reference_info <- function(x) {
  if (!inherits(x, "SharedReference")) {
    .encoder_cli_abort("x must be a SharedReference.",
                       class = "fmrilatent_error_invalid_type")
  }
  list(
    id = x$id,
    kind = x$kind,
    has_cached_value = exists(x$id, .shared_reference_env(), inherits = FALSE),
    persistence = x$persistence,
    meta = x$meta
  )
}

#' Reset the in-session shared reference cache
#'
#' @return Invisibly, `TRUE`.
#' @export
shared_reference_clear <- function() {
  .shared_reference_env(reset = TRUE)
  invisible(TRUE)
}

#' Build a method-neutral shared component contract
#'
#' @description
#' A shared component contract describes a reusable spatial component matrix
#' without committing to a particular decomposition method. The contract records
#' dimensions, domain/support metadata, a digest, and optional advisory fields.
#'
#' @param loadings Matrix-like component loadings (`features x components`), a
#'   `SharedReference` to such a matrix, or a `GroupDeltaLoadings` object.
#' @param id Optional component id.
#' @param family Method-neutral or method-family label.
#' @param domain_id Stable domain identifier.
#' @param support Optional support descriptor whose cardinality must match
#'   `nrow(loadings)` when supplied.
#' @param measure Optional support measure. Vectors must have one value per
#'   feature; matrices must be square over features.
#' @param meta Optional advisory metadata.
#' @return A `SharedComponentContract` object.
#' @export
shared_component_contract <- function(loadings, id = NULL, family = "shared_component",
                                      domain_id = "", support = NULL,
                                      measure = NULL, meta = list()) {
  L <- .shared_matrix_like(loadings, "loadings")
  family <- .shared_check_id(family, "family")
  domain_id <- as.character(domain_id %||% "")
  if (length(domain_id) != 1L || is.na(domain_id)) {
    .encoder_cli_abort("domain_id must be a character scalar.",
                       class = "fmrilatent_error_invalid_argument")
  }
  meta <- .shared_check_list(meta, "meta")
  support_n <- .shared_support_cardinality(support, context = "shared component support")
  if (!is.null(support_n) && support_n != nrow(L)) {
    .encoder_cli_abort(
      paste0("support cardinality ", support_n, " does not match loadings rows ",
             nrow(L), "."),
      class = "fmrilatent_error_dimension_mismatch"
    )
  }
  measure_use <- .template_measure_resolve(
    measure,
    n_support = nrow(L),
    default = "null",
    context = "shared component measure"
  )
  if (is.null(id)) {
    id <- paste0(
      "component:",
      digest::digest(list(
        family = family,
        domain_id = domain_id,
        dim = dim(L),
        support = support,
        loadings = L
      ))
    )
  }
  id <- .shared_check_id(id)

  out <- list(
    id = id,
    contract_version = 1L,
    family = family,
    n_features = as.integer(nrow(L)),
    n_components = as.integer(ncol(L)),
    domain_id = domain_id,
    support = support,
    measure = measure_use,
    loadings_ref = if (inherits(loadings, "SharedReference")) loadings else NULL,
    digest = digest::digest(list(dim = dim(L), loadings = L, measure = measure_use)),
    meta = meta
  )
  structure(out, class = "SharedComponentContract")
}

#' Validate a shared component contract
#'
#' @param x Object to validate.
#' @param error If `TRUE`, throw on failure; otherwise return `FALSE`.
#' @return Invisibly, `x` on success, or `FALSE` when `error = FALSE`.
#' @export
validate_shared_component_contract <- function(x, error = TRUE) {
  res <- tryCatch({
    if (!inherits(x, "SharedComponentContract")) {
      .encoder_cli_abort("x must be a SharedComponentContract.",
                         class = "fmrilatent_error_invalid_type")
    }
    .shared_check_id(x$id)
    if (!identical(as.integer(x$contract_version), 1L)) {
      .encoder_cli_abort("unsupported shared component contract_version.",
                         class = "fmrilatent_error_unsupported_version")
    }
    .shared_check_positive_integer(x$n_features, "n_features")
    .shared_check_positive_integer(x$n_components, "n_components")
    .shared_check_id(x$family, "family")
    support_n <- .shared_support_cardinality(x$support, context = "shared component support")
    if (!is.null(support_n) && support_n != x$n_features) {
      .encoder_cli_abort("support cardinality does not match n_features.",
                         class = "fmrilatent_error_dimension_mismatch")
    }
    x
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (isTRUE(error)) stop(res)  # re-raise to preserve the classed inner condition
    return(invisible(FALSE))
  }
  invisible(res)
}

#' Validate the reusable template protocol
#'
#' @description
#' Checks that a template asset exposes the method-neutral surface expected by
#' shared-structure consumers: loadings, rank, support/domain, measure,
#' roughness, decoder, and projection behavior.
#'
#' @param template Template object.
#' @param error If `TRUE`, throw on failure; otherwise return `FALSE`.
#' @return Invisibly, a manifest describing the verified protocol.
#' @export
validate_template_protocol <- function(template, error = TRUE) {
  res <- tryCatch({
    L <- .shared_matrix_like(template_loadings(template), "template_loadings(template)")
    rank <- template_rank(template)
    rank <- .shared_check_positive_integer(rank, "template_rank(template)")
    if (rank != ncol(L)) {
      .encoder_cli_abort(
        "template_rank(template) does not match ncol(template_loadings(template)).",
        class = "fmrilatent_error_dimension_mismatch"
      )
    }
    domain <- template_domain(template)
    support <- template_support(template)
    support_n <- .shared_support_cardinality(
      support,
      domain = domain,
      context = "template_support(template)"
    )
    if (!is.null(support_n) && support_n != nrow(L)) {
      .encoder_cli_abort(
        "template_support(template) cardinality does not match loadings rows.",
        class = "fmrilatent_error_dimension_mismatch"
      )
    }
    measure <- template_measure(template)
    .template_measure_resolve(
      measure,
      n_support = nrow(L),
      default = "null",
      context = "template_measure(template)"
    )
    roughness <- template_roughness(template, coordinates = "analysis")
    if (!is.null(roughness)) {
      roughness <- as.matrix(roughness)
      if (!identical(dim(roughness), c(rank, rank))) {
        .encoder_cli_abort(
          "template_roughness(template) must have dimensions rank x rank.",
          class = "fmrilatent_error_dimension_mismatch"
        )
      }
    }
    decoder_map <- .normalize_linear_map(basis_decoder(template), context = "basis decoder")
    if (decoder_map$n_source != rank || decoder_map$n_target != nrow(L)) {
      .encoder_cli_abort(
        "basis_decoder(template) dimensions do not match template loadings.",
        class = "fmrilatent_error_dimension_mismatch"
      )
    }
    list(
      class = class(template),
      rank = rank,
      n_features = as.integer(nrow(L)),
      domain = domain,
      support_cardinality = support_n,
      decoder_contract_version = decoder_map$contract_version %||% 1L,
      has_measure = !is.null(measure),
      has_roughness = !is.null(roughness)
    )
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (isTRUE(error)) stop(res)  # re-raise to preserve the classed inner condition
    return(invisible(FALSE))
  }
  invisible(res)
}

#' Represent subject loadings as group shared loadings plus a delta
#'
#' @param group Matrix-like shared/group loadings or a `SharedReference`.
#' @param delta Optional subject-specific loading delta. If `NULL`, a zero
#'   delta is used.
#' @param scale Numeric multiplier applied to `delta`.
#' @param id Optional id.
#' @param meta Optional advisory metadata.
#' @return A `GroupDeltaLoadings` object.
#' @export
group_delta_loadings <- function(group, delta = NULL, scale = 1, id = NULL,
                                 meta = list()) {
  group_mat <- .shared_matrix_like(group, "group")
  if (is.null(delta)) {
    delta_mat <- Matrix::Matrix(0, nrow = nrow(group_mat), ncol = ncol(group_mat),
                                sparse = TRUE)
  } else {
    delta_mat <- .shared_matrix_like(delta, "delta")
  }
  if (!identical(dim(group_mat), dim(delta_mat))) {
    .encoder_cli_abort("group and delta must have identical dimensions.",
                       class = "fmrilatent_error_dimension_mismatch")
  }
  if (length(scale) != 1L || !is.finite(scale)) {
    .encoder_cli_abort("scale must be a finite numeric scalar.",
                       class = "fmrilatent_error_invalid_argument")
  }
  meta <- .shared_check_list(meta, "meta")
  if (is.null(id)) {
    id <- paste0("group-delta:", digest::digest(list(group = group_mat, delta = delta_mat, scale = scale)))
  }
  structure(
    list(
      id = .shared_check_id(id),
      group = if (inherits(group, "SharedReference")) group else shared_reference(group_mat, kind = "loadings"),
      delta = if (inherits(delta, "SharedReference")) delta else shared_reference(delta_mat, kind = "loadings_delta"),
      scale = as.numeric(scale),
      dim = as.integer(dim(group_mat)),
      meta = meta
    ),
    class = "GroupDeltaLoadings"
  )
}

#' Materialize group-plus-delta loadings
#'
#' @param x A `GroupDeltaLoadings` object.
#' @return Matrix-like materialized loadings.
#' @export
materialize_group_delta_loadings <- function(x) {
  if (!inherits(x, "GroupDeltaLoadings")) {
    .encoder_cli_abort("x must be a GroupDeltaLoadings object.",
                       class = "fmrilatent_error_invalid_type")
  }
  group <- .shared_matrix_like(resolve_shared_reference(x$group), "group")
  delta <- .shared_matrix_like(resolve_shared_reference(x$delta), "delta")
  if (!identical(dim(group), dim(delta))) {
    .encoder_cli_abort("stored group and delta dimensions no longer match.",
                       class = "fmrilatent_error_dimension_mismatch")
  }
  group + x$scale * delta
}

#' Coerce group-plus-delta loadings to a dense matrix
#'
#' @param x A `GroupDeltaLoadings` object.
#' @param ... Ignored.
#' @return Dense matrix of materialized loadings.
#' @method as.matrix GroupDeltaLoadings
#' @export
as.matrix.GroupDeltaLoadings <- function(x, ...) {
  as.matrix(materialize_group_delta_loadings(x))
}

#' Construct a shared temporal basis descriptor
#'
#' @param kind Temporal basis kind: `"dct"`, `"bspline"`, `"slepian"`, or
#'   `"custom"`.
#' @param n_time Number of time points. Inferred from `basis` for custom specs.
#' @param rank Number of basis columns. Inferred from `basis` when omitted.
#' @param basis Matrix-like basis for `kind = "custom"`.
#' @param params Basis-specific parameters.
#' @param id Optional id.
#' @param meta Optional advisory metadata.
#' @return A `SharedTemporalSpec` descriptor.
#' @export
shared_temporal_spec <- function(kind = c("dct", "bspline", "slepian", "custom"),
                                 n_time = NULL, rank = NULL, basis = NULL,
                                 params = list(), id = NULL, meta = list()) {
  kind <- match.arg(kind)
  params <- .shared_check_list(params, "params")
  meta <- .shared_check_list(meta, "meta")
  basis_ref <- NULL
  if (!is.null(basis)) {
    B <- .shared_matrix_like(basis, "basis")
    n_time <- n_time %||% nrow(B)
    rank <- rank %||% ncol(B)
    if (n_time != nrow(B) || rank != ncol(B)) {
      .encoder_cli_abort("basis dimensions must match n_time and rank.",
                         class = "fmrilatent_error_dimension_mismatch")
    }
    basis_ref <- if (inherits(basis, "SharedReference")) {
      basis
    } else {
      shared_reference(B, kind = "temporal_basis")
    }
  }
  n_time <- .shared_check_positive_integer(n_time, "n_time")
  rank <- .shared_check_positive_integer(rank, "rank")
  if (rank > n_time) {
    .encoder_cli_abort("rank cannot exceed n_time.",
                       class = "fmrilatent_error_invalid_argument")
  }
  if (is.null(id)) {
    id <- paste0("temporal:", digest::digest(list(kind = kind, n_time = n_time, rank = rank, params = params, basis = basis_ref)))
  }
  structure(
    list(
      id = .shared_check_id(id),
      contract_version = 1L,
      kind = kind,
      n_time = n_time,
      rank = rank,
      basis = basis_ref,
      params = params,
      meta = meta
    ),
    class = "SharedTemporalSpec"
  )
}

#' Materialize a shared temporal specification
#'
#' @param spec A `SharedTemporalSpec`.
#' @return Dense temporal basis matrix.
#' @export
materialize_shared_temporal_spec <- function(spec) {
  if (!inherits(spec, "SharedTemporalSpec")) {
    .encoder_cli_abort("spec must be a SharedTemporalSpec.",
                       class = "fmrilatent_error_invalid_type")
  }
  if (!is.null(spec$basis)) {
    return(as.matrix(resolve_shared_reference(spec$basis)))
  }
  if (identical(spec$kind, "dct")) {
    return(as.matrix(build_dct_basis(
      spec$n_time,
      spec$rank,
      norm = spec$params$norm %||% "ortho"
    )))
  }
  if (identical(spec$kind, "bspline")) {
    return(as.matrix(build_bspline_basis(
      n_time = spec$n_time,
      k = spec$rank,
      degree = spec$params$degree %||% 3L,
      knots = spec$params$knots %||% NULL,
      boundary_knots = spec$params$boundary_knots %||% NULL,
      include_intercept = spec$params$include_intercept %||% FALSE,
      orthonormalize = spec$params$orthonormalize %||% TRUE
    )))
  }
  if (identical(spec$kind, "slepian")) {
    return(as.matrix(dpss_time_basis(
      n_time = spec$n_time,
      tr = spec$params$tr %||% 1,
      bandwidth = spec$params$bandwidth %||% 0.1,
      k = spec$rank,
      backend = spec$params$backend %||% "tridiag"
    )))
  }
  .encoder_cli_abort(
    paste0("Shared temporal spec kind '", spec$kind,
           "' is a descriptor only and cannot be materialized without a custom basis."),
    class = "fmrilatent_error_unsupported_operation"
  )
}

#' Construct a reusable event shape dictionary
#'
#' @param shapes Named list of numeric shape vectors.
#' @param n_time Optional temporal length used for validation.
#' @param meta Optional advisory metadata.
#' @return A `SharedEventDictionary`.
#' @export
shared_event_dictionary <- function(shapes = list(impulse = 1), n_time = NULL,
                                    meta = list()) {
  shapes <- .shared_check_list(shapes, "shapes")
  if (length(shapes) == 0L || is.null(names(shapes)) || any(!nzchar(names(shapes)))) {
    .encoder_cli_abort("shapes must be a non-empty named list.",
                       class = "fmrilatent_error_invalid_argument")
  }
  shapes <- lapply(shapes, function(shape) {
    shape <- as.numeric(shape)
    if (length(shape) < 1L || any(!is.finite(shape))) {
      .encoder_cli_abort("each event shape must be a non-empty finite numeric vector.",
                         class = "fmrilatent_error_invalid_value")
    }
    shape
  })
  if (!is.null(n_time)) {
    n_time <- .shared_check_positive_integer(n_time, "n_time")
  }
  structure(
    list(
      id = paste0("events:", digest::digest(list(shapes = shapes, n_time = n_time))),
      contract_version = 1L,
      shapes = shapes,
      n_time = n_time,
      meta = .shared_check_list(meta, "meta")
    ),
    class = "SharedEventDictionary"
  )
}

#' Render sparse events with a shared event dictionary
#'
#' @param dictionary A `SharedEventDictionary`.
#' @param events Data frame with columns `atom`, `time`, `amplitude`, and
#'   optional `shape_id`.
#' @param n_atoms Number of event rows/atoms.
#' @param n_time Number of time points. Defaults to the dictionary value.
#' @return Sparse `n_atoms x n_time` event matrix.
#' @export
render_shared_events <- function(dictionary, events, n_atoms, n_time = NULL) {
  if (!inherits(dictionary, "SharedEventDictionary")) {
    .encoder_cli_abort("dictionary must be a SharedEventDictionary.",
                       class = "fmrilatent_error_invalid_type")
  }
  n_atoms <- .shared_check_positive_integer(n_atoms, "n_atoms")
  n_time <- n_time %||% dictionary$n_time
  n_time <- .shared_check_positive_integer(n_time, "n_time")
  if (nrow(events) == 0L) {
    return(Matrix::Matrix(0, nrow = n_atoms, ncol = n_time, sparse = TRUE))
  }
  required <- c("atom", "time", "amplitude")
  missing <- setdiff(required, names(events))
  if (length(missing) > 0L) {
    # `missing` is a subset of the fixed `required` names, so no glue hazard.
    .encoder_cli_abort(
      paste0("events is missing required columns: ", paste(missing, collapse = ", "), "."),
      class = "fmrilatent_error_invalid_events"
    )
  }
  shape_id <- events$shape_id %||% rep.int(names(dictionary$shapes)[[1L]], nrow(events))
  n_events <- nrow(events)
  row_chunks <- vector("list", n_events)
  col_chunks <- vector("list", n_events)
  value_chunks <- vector("list", n_events)
  for (idx in seq_len(n_events)) {
    shape <- dictionary$shapes[[as.character(shape_id[[idx]])]]
    if (is.null(shape)) {
      # shape_id is user-supplied free-form text; escape glue braces so a stray
      # "{" cannot break cli_abort's message interpolation.
      sid <- gsub("}", "}}", gsub("{", "{{", as.character(shape_id[[idx]]), fixed = TRUE), fixed = TRUE)
      .encoder_cli_abort(paste0("unknown event shape_id '", sid, "'."),
                         class = "fmrilatent_error_invalid_events")
    }
    atom <- as.integer(events$atom[[idx]])
    start <- as.integer(events$time[[idx]])
    amp <- as.numeric(events$amplitude[[idx]])
    if (atom < 1L || atom > n_atoms || start < 1L || start > n_time || !is.finite(amp)) {
      .encoder_cli_abort(
        "events contains an out-of-range atom/time or non-finite amplitude.",
        class = "fmrilatent_error_invalid_events"
      )
    }
    cols_i <- seq.int(start, length.out = length(shape))
    keep <- cols_i <= n_time
    n_keep <- sum(keep)
    row_chunks[[idx]] <- rep.int(atom, n_keep)
    col_chunks[[idx]] <- cols_i[keep]
    value_chunks[[idx]] <- amp * shape[keep]
  }
  rows <- unlist(row_chunks, use.names = FALSE)
  cols <- unlist(col_chunks, use.names = FALSE)
  values <- unlist(value_chunks, use.names = FALSE)
  Matrix::sparseMatrix(i = rows, j = cols, x = values,
                       dims = c(n_atoms, n_time))
}

#' Create the fmrilatent-to-neuroarchive handoff contract
#'
#' @description
#' The handoff contract records the boundary between in-memory representation
#' responsibilities owned by fmrilatent and persistent archive responsibilities
#' owned by neuroarchive. It is intentionally a manifest, not a lazy archive
#' reader or file locator.
#'
#' @param representation Optional fmrilatent representation object.
#' @param components Optional list of shared component contracts.
#' @param templates Optional list of reusable template assets.
#' @param references Optional list of in-session shared references.
#' @param meta Optional advisory metadata.
#' @return A `NeuroarchiveHandoffContract`.
#' @export
neuroarchive_handoff_contract <- function(representation = NULL,
                                          components = list(),
                                          templates = list(),
                                          references = list(),
                                          meta = list()) {
  components <- .shared_check_list(components, "components")
  templates <- .shared_check_list(templates, "templates")
  references <- .shared_check_list(references, "references")
  components <- lapply(components, validate_shared_component_contract)
  template_manifests <- lapply(templates, validate_template_protocol)
  if (!all(vapply(references, is_shared_reference, logical(1)))) {
    .encoder_cli_abort("references must contain only SharedReference objects.",
                       class = "fmrilatent_error_invalid_type")
  }
  structure(
    list(
      contract_version = 1L,
      representation_class = if (is.null(representation)) character(0) else class(representation),
      fmrilatent_owns = c(
        "method-neutral representation contracts",
        "reusable templates",
        "projection and reconstruction behavior",
        "in-session shared references",
        "group-plus-delta abstractions",
        "temporal specs and event dictionaries"
      ),
      neuroarchive_owns = c(
        "persistent shared stores",
        "HDF5 artifacts",
        "persistent registries",
        "checksums",
        "locators",
        "lazy archive readers"
      ),
      components = components,
      templates = templates,
      template_manifests = template_manifests,
      references = lapply(references, shared_reference_info),
      persistent = FALSE,
      archive_locator = NULL,
      checksum = NULL,
      meta = .shared_check_list(meta, "meta")
    ),
    class = "NeuroarchiveHandoffContract"
  )
}

#' Validate the fmrilatent-to-neuroarchive handoff contract
#'
#' @param x Object to validate.
#' @param error If `TRUE`, throw on failure; otherwise return `FALSE`.
#' @return Invisibly, `x` on success, or `FALSE` when `error = FALSE`.
#' @export
validate_neuroarchive_handoff_contract <- function(x, error = TRUE) {
  res <- tryCatch({
    if (!inherits(x, "NeuroarchiveHandoffContract")) {
      .encoder_cli_abort("x must be a NeuroarchiveHandoffContract.",
                         class = "fmrilatent_error_invalid_type")
    }
    if (!identical(as.integer(x$contract_version), 1L)) {
      .encoder_cli_abort("unsupported neuroarchive handoff contract_version.",
                         class = "fmrilatent_error_unsupported_version")
    }
    if (!identical(x$persistent, FALSE) || !is.null(x$archive_locator) || !is.null(x$checksum)) {
      .encoder_cli_abort(
        "fmrilatent handoff contracts must not contain persistent archive locators or checksums.",
        class = "fmrilatent_error_invalid_handoff"
      )
    }
    x
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (isTRUE(error)) stop(res)  # re-raise to preserve the classed inner condition
    return(invisible(FALSE))
  }
  invisible(res)
}
