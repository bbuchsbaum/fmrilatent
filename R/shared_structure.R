# Shared-structure contracts for reusable latent representations

.shared_check_list <- function(x, name) {
  if (!is.list(x)) {
    stop(name, " must be a list.", call. = FALSE)
  }
  x
}

.shared_check_id <- function(id, name = "id") {
  if (length(id) != 1L || is.na(id) || !nzchar(id)) {
    stop(name, " must be a non-empty character scalar.", call. = FALSE)
  }
  as.character(id)
}

.shared_check_positive_integer <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x != as.integer(x) || x < 1L) {
    stop(name, " must be a positive integer scalar.", call. = FALSE)
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
    stop(name, " must be a matrix-like object.", call. = FALSE)
  }
  x <- Matrix::Matrix(x, sparse = inherits(x, "sparseMatrix"))
  if (any(!is.finite(as.matrix(x)))) {
    stop(name, " must contain only finite values.", call. = FALSE)
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
        stop(conditionMessage(e), call. = FALSE)
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
    stop("materialize must be NULL or a function.", call. = FALSE)
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
    stop("Shared reference '", x$id, "' has no cached value or materialize() function.",
         call. = FALSE)
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
    stop("x must be a SharedReference.", call. = FALSE)
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
    stop("domain_id must be a character scalar.", call. = FALSE)
  }
  meta <- .shared_check_list(meta, "meta")
  support_n <- .shared_support_cardinality(support, context = "shared component support")
  if (!is.null(support_n) && support_n != nrow(L)) {
    stop("support cardinality ", support_n, " does not match loadings rows ",
         nrow(L), ".", call. = FALSE)
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
      stop("x must be a SharedComponentContract.", call. = FALSE)
    }
    .shared_check_id(x$id)
    if (!identical(as.integer(x$contract_version), 1L)) {
      stop("unsupported shared component contract_version.", call. = FALSE)
    }
    .shared_check_positive_integer(x$n_features, "n_features")
    .shared_check_positive_integer(x$n_components, "n_components")
    .shared_check_id(x$family, "family")
    support_n <- .shared_support_cardinality(x$support, context = "shared component support")
    if (!is.null(support_n) && support_n != x$n_features) {
      stop("support cardinality does not match n_features.", call. = FALSE)
    }
    x
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (isTRUE(error)) stop(conditionMessage(res), call. = FALSE)
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
      stop("template_rank(template) does not match ncol(template_loadings(template)).",
           call. = FALSE)
    }
    domain <- template_domain(template)
    support <- template_support(template)
    support_n <- .shared_support_cardinality(
      support,
      domain = domain,
      context = "template_support(template)"
    )
    if (!is.null(support_n) && support_n != nrow(L)) {
      stop("template_support(template) cardinality does not match loadings rows.",
           call. = FALSE)
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
        stop("template_roughness(template) must have dimensions rank x rank.",
             call. = FALSE)
      }
    }
    decoder_map <- .normalize_linear_map(basis_decoder(template), context = "basis decoder")
    if (decoder_map$n_source != rank || decoder_map$n_target != nrow(L)) {
      stop("basis_decoder(template) dimensions do not match template loadings.",
           call. = FALSE)
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
    if (isTRUE(error)) stop(conditionMessage(res), call. = FALSE)
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
    stop("group and delta must have identical dimensions.", call. = FALSE)
  }
  if (length(scale) != 1L || !is.finite(scale)) {
    stop("scale must be a finite numeric scalar.", call. = FALSE)
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
    stop("x must be a GroupDeltaLoadings object.", call. = FALSE)
  }
  group <- .shared_matrix_like(resolve_shared_reference(x$group), "group")
  delta <- .shared_matrix_like(resolve_shared_reference(x$delta), "delta")
  if (!identical(dim(group), dim(delta))) {
    stop("stored group and delta dimensions no longer match.", call. = FALSE)
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
      stop("basis dimensions must match n_time and rank.", call. = FALSE)
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
    stop("rank cannot exceed n_time.", call. = FALSE)
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
    stop("spec must be a SharedTemporalSpec.", call. = FALSE)
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
  stop("Shared temporal spec kind '", spec$kind,
       "' is a descriptor only and cannot be materialized without a custom basis.",
       call. = FALSE)
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
    stop("shapes must be a non-empty named list.", call. = FALSE)
  }
  shapes <- lapply(shapes, function(shape) {
    shape <- as.numeric(shape)
    if (length(shape) < 1L || any(!is.finite(shape))) {
      stop("each event shape must be a non-empty finite numeric vector.",
           call. = FALSE)
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
    stop("dictionary must be a SharedEventDictionary.", call. = FALSE)
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
    stop("events is missing required columns: ", paste(missing, collapse = ", "),
         ".", call. = FALSE)
  }
  shape_id <- events$shape_id %||% rep.int(names(dictionary$shapes)[[1L]], nrow(events))
  rows <- integer()
  cols <- integer()
  values <- numeric()
  for (idx in seq_len(nrow(events))) {
    shape <- dictionary$shapes[[as.character(shape_id[[idx]])]]
    if (is.null(shape)) {
      stop("unknown event shape_id '", shape_id[[idx]], "'.", call. = FALSE)
    }
    atom <- as.integer(events$atom[[idx]])
    start <- as.integer(events$time[[idx]])
    amp <- as.numeric(events$amplitude[[idx]])
    if (atom < 1L || atom > n_atoms || start < 1L || start > n_time || !is.finite(amp)) {
      stop("events contains an out-of-range atom/time or non-finite amplitude.",
           call. = FALSE)
    }
    cols_i <- seq.int(start, length.out = length(shape))
    keep <- cols_i <= n_time
    rows <- c(rows, rep.int(atom, sum(keep)))
    cols <- c(cols, cols_i[keep])
    values <- c(values, amp * shape[keep])
  }
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
  invisible(lapply(components, validate_shared_component_contract))
  invisible(lapply(templates, validate_template_protocol))
  if (!all(vapply(references, is_shared_reference, logical(1)))) {
    stop("references must contain only SharedReference objects.", call. = FALSE)
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
      templates = lapply(templates, function(template) validate_template_protocol(template)),
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
      stop("x must be a NeuroarchiveHandoffContract.", call. = FALSE)
    }
    if (!identical(as.integer(x$contract_version), 1L)) {
      stop("unsupported neuroarchive handoff contract_version.", call. = FALSE)
    }
    if (!identical(x$persistent, FALSE) || !is.null(x$archive_locator) || !is.null(x$checksum)) {
      stop("fmrilatent handoff contracts must not contain persistent archive locators or checksums.",
           call. = FALSE)
    }
    x
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (isTRUE(error)) stop(conditionMessage(res), call. = FALSE)
    return(invisible(FALSE))
  }
  invisible(res)
}
