#' @include latent_utils.R
NULL

# Handle classes and registry for shared/implicit dictionaries

#' @keywords internal
setClass(
  "BasisHandle",
  slots = c(
    id    = "character",
    dim   = "integer",
    kind  = "character",
    spec  = "list",
    label = "character"
  )
)

.handle_kind_registries <- local({
  basis <- new.env(parent = emptyenv())
  loadings <- new.env(parent = emptyenv())
  function(type = c("basis", "loadings")) {
    type <- match.arg(type)
    if (identical(type, "basis")) basis else loadings
  }
})

.handle_kind_names <- function(type = c("basis", "loadings")) {
  sort(ls(.handle_kind_registries(type)))
}

.handle_kind_materializer <- function(type = c("basis", "loadings"), kind) {
  type <- match.arg(type)
  env <- .handle_kind_registries(type)
  if (!exists(kind, envir = env, inherits = FALSE)) {
    .encoder_cli_abort(
      paste0("Unknown ", type, " handle kind: ", kind),
      class = "fmrilatent_error_value",
      call = rlang::caller_env()
    )
  }
  get(kind, envir = env, inherits = FALSE)
}

#' Register a lazy handle materializer
#'
#' @param kind Character scalar handle kind.
#' @param materializer Function that accepts a handle and returns a matrix-like object.
#' @param type Registry to update: \code{"basis"} or \code{"loadings"}.
#' @return Invisibly, the registered kind.
#' @export
register_handle_kind <- function(kind, materializer, type = c("basis", "loadings")) {
  type <- match.arg(type)
  if (length(kind) != 1L || is.na(kind) || !nzchar(kind)) {
    .encoder_cli_abort("kind must be a non-empty character scalar.",
                       class = "fmrilatent_error_invalid_argument")
  }
  if (!is.function(materializer)) {
    .encoder_cli_abort("materializer must be a function.",
                       class = "fmrilatent_error_invalid_argument")
  }
  assign(kind, materializer, envir = .handle_kind_registries(type))
  invisible(kind)
}

.validate_handle_dim <- function(dim, class_name) {
  if (length(dim) != 2L) {
    return(paste0(class_name, "@dim must have length 2."))
  }
  if (any(is.na(dim)) || any(dim < 1L)) {
    return(paste0(class_name, "@dim must contain positive integers."))
  }
  TRUE
}

.validate_explicit_handle_spec <- function(spec, class_name) {
  if (is.null(spec$matrix)) {
    return(paste0(class_name, "(kind = 'explicit') requires spec$matrix."))
  }
  TRUE
}

setValidity("BasisHandle", function(object) {
  if (length(object@id) != 1L || is.na(object@id) || !nzchar(object@id)) {
    return("BasisHandle@id must be a non-empty character scalar.")
  }
  dim_valid <- .validate_handle_dim(object@dim, "BasisHandle")
  if (!isTRUE(dim_valid)) return(dim_valid)
  basis_kinds <- .handle_kind_names("basis")
  if (length(object@kind) != 1L || is.na(object@kind) || !(object@kind %in% basis_kinds)) {
    return(paste0(
      "BasisHandle@kind must be a registered basis handle kind. Registered kinds: ",
      paste(basis_kinds, collapse = ", "),
      "."
    ))
  }
  if (length(object@label) != 1L || is.na(object@label)) {
    return("BasisHandle@label must be a character scalar.")
  }
  if (!is.list(object@spec)) {
    return("BasisHandle@spec must be a list.")
  }
  if (identical(object@kind, "explicit")) {
    spec_valid <- .validate_explicit_handle_spec(object@spec, "BasisHandle")
    if (!isTRUE(spec_valid)) return(spec_valid)
  }
  TRUE
})

#' @keywords internal
setClass(
  "LoadingsHandle",
  slots = c(
    id    = "character",
    dim   = "integer",
    kind  = "character",
    spec  = "list",
    label = "character"
  )
)

setValidity("LoadingsHandle", function(object) {
  if (length(object@id) != 1L || is.na(object@id) || !nzchar(object@id)) {
    return("LoadingsHandle@id must be a non-empty character scalar.")
  }
  dim_valid <- .validate_handle_dim(object@dim, "LoadingsHandle")
  if (!isTRUE(dim_valid)) return(dim_valid)
  loadings_kinds <- .handle_kind_names("loadings")
  if (length(object@kind) != 1L || is.na(object@kind) || !(object@kind %in% loadings_kinds)) {
    return(paste0(
      "LoadingsHandle@kind must be a registered loadings handle kind. Registered kinds: ",
      paste(loadings_kinds, collapse = ", "),
      "."
    ))
  }
  if (length(object@label) != 1L || is.na(object@label)) {
    return("LoadingsHandle@label must be a character scalar.")
  }
  if (!is.list(object@spec)) {
    return("LoadingsHandle@spec must be a list.")
  }
  if (identical(object@kind, "explicit")) {
    spec_valid <- .validate_explicit_handle_spec(object@spec, "LoadingsHandle")
    if (!isTRUE(spec_valid)) return(spec_valid)
  }
  TRUE
})

# Unions so LatentNeuroVec can accept either Matrix or handle
setClassUnion("MatrixOrBasisHandle",    c("Matrix", "matrix", "BasisHandle"))
setClassUnion("MatrixOrLoadingsHandle", c("Matrix", "matrix", "LoadingsHandle"))

# Encapsulated cache for materialized dictionaries.
# NOTE: This is intentionally not a top-level mutable registry object. The cache
# is only accessible through helper functions and can be disabled via option.
.fmrilatent_cache_env <- local({
  basis_env <- NULL
  loadings_env <- NULL

  new_registry_env <- function() {
    env <- new.env(parent = emptyenv())
    attr(env, "fmrilatent.lru_order") <- character(0)
    env
  }

  init <- function() {
    basis_env <<- new_registry_env()
    loadings_env <<- new_registry_env()
    invisible(TRUE)
  }

  function(type = c("basis", "loadings"), reset = FALSE) {
    type <- match.arg(type)
    if (reset || is.null(basis_env) || is.null(loadings_env)) {
      init()
    }
    if (type == "basis") basis_env else loadings_env
  }
})

.latent_registry_order <- function(env) {
  attr(env, "fmrilatent.lru_order", exact = TRUE) %||% character(0)
}

.latent_registry_set_order <- function(env, order) {
  attr(env, "fmrilatent.lru_order") <- as.character(order)
  invisible(NULL)
}

.latent_registry_touch <- function(env, id) {
  .latent_registry_set_order(env, c(setdiff(.latent_registry_order(env), id), id))
}

.latent_registry_drop <- function(env, ids) {
  .latent_registry_set_order(env, setdiff(.latent_registry_order(env), ids))
}

.latent_registry_reset <- function(env) {
  .latent_registry_set_order(env, character(0))
}

.latent_registry_enabled <- function() {
  isTRUE(getOption("fmrilatent.registry.enabled", TRUE))
}

# Per-type maximum number of cached matrices before LRU eviction kicks in.
# Controlled by option `fmrilatent.registry.max_entries` (default 256). Set the
# option to `Inf` (or a non-positive value) for an unbounded cache.
.fmrilatent_registry_default_max <- 256L

.latent_registry_max_entries <- function() {
  val <- getOption("fmrilatent.registry.max_entries", .fmrilatent_registry_default_max)
  if (length(val) != 1L || is.na(val)) return(.fmrilatent_registry_default_max)
  if (is.infinite(val) || val <= 0) return(Inf)
  as.integer(val)
}

# Evict least-recently-used entries from a registry env until it fits the cap.
.latent_enforce_cap <- function(type) {
  cap <- .latent_registry_max_entries()
  if (is.infinite(cap)) return(invisible(NULL))
  env     <- .latent_get_registry_env(type)
  order   <- .latent_registry_order(env)
  n_over  <- length(order) - cap
  if (n_over > 0) {
    evict <- order[seq_len(n_over)]
    rm(list = evict, envir = env)
    .latent_registry_drop(env, evict)
  }
  invisible(NULL)
}

.latent_get_registry_env <- function(type = c("basis", "loadings")) {
  type <- match.arg(type)
  .fmrilatent_cache_env(type)
}

.latent_register_matrix <- function(id, value, type = c("basis", "loadings"),
                                    overwrite = FALSE, fingerprint = NULL) {
  stopifnot(is.character(id), length(id) == 1L)
  type <- match.arg(type)
  if (!.latent_registry_enabled()) {
    return(invisible(FALSE))
  }
  env  <- .latent_get_registry_env(type)
  if (!overwrite && exists(id, env, inherits = FALSE)) {
    existing <- get(id, envir = env, inherits = FALSE)
    existing_fingerprint <- attr(existing, "fmrilatent.handle_fingerprint", exact = TRUE)
    if (!is.null(fingerprint)) {
      # A fingerprinted registration must be provably the same entry. Refuse to
      # silently reuse an entry whose fingerprint differs OR is absent (an older
      # entry registered before fingerprinting cannot be verified to match).
      if (is.null(existing_fingerprint)) {
        .encoder_cli_abort(
          paste0("Object with id '", id, "' is already registered in ", type,
                 " registry without a fingerprint; cannot confirm it matches the ",
                 "requested handle. Clear the registry or register with overwrite = TRUE."),
          class = "fmrilatent_error_invalid_id", call = rlang::caller_env())
      }
      if (!identical(fingerprint, existing_fingerprint)) {
        .encoder_cli_abort(
          paste0("Object with id '", id, "' is already registered in ", type,
                 " registry with a different fingerprint."),
          class = "fmrilatent_error_invalid_id", call = rlang::caller_env())
      }
      # Same id, same fingerprint: reuse the cached entry and mark it as fresh.
      .latent_registry_touch(env, id)
      return(invisible(FALSE))
    }
    .encoder_cli_warn(
      paste0("Object with id '", id, "' already registered in ", type,
             " registry; set overwrite = TRUE to replace."),
      class = "fmrilatent_warning_encoder", call = rlang::caller_env())
    return(invisible(FALSE))
  }
  if (!is.null(fingerprint)) {
    attr(value, "fmrilatent.handle_fingerprint") <- fingerprint
  }
  assign(id, value, envir = env)
  .latent_registry_touch(env, id)
  .latent_enforce_cap(type)
  invisible(TRUE)
}

.latent_get_matrix <- function(id, type = c("basis", "loadings"), fingerprint = NULL) {
  stopifnot(is.character(id), length(id) == 1L)
  type <- match.arg(type)
  if (!.latent_registry_enabled()) {
    return(NULL)
  }
  env  <- .latent_get_registry_env(type)
  if (exists(id, env, inherits = FALSE)) {
    value <- get(id, envir = env, inherits = FALSE)
    existing_fingerprint <- attr(value, "fmrilatent.handle_fingerprint", exact = TRUE)
    if (!is.null(fingerprint) && !is.null(existing_fingerprint) &&
        !identical(fingerprint, existing_fingerprint)) {
      .encoder_cli_abort(
        paste0("Registry entry '", id, "' in ", type,
               " registry does not match the requested handle fingerprint."),
        class = "fmrilatent_error_invalid_id", call = rlang::caller_env())
    }
    .latent_registry_touch(env, id)
    value
  } else {
    NULL
  }
}

.latent_has_matrix <- function(id, type = c("basis", "loadings")) {
  !is.null(.latent_get_matrix(id, type = type))
}

.latent_handle_fingerprint <- function(handle) {
  digest::digest(
    list(
      class = class(handle),
      kind = handle@kind,
      dim = as.integer(handle@dim),
      spec = handle@spec
    ),
    algo = "xxhash64"
  )
}

.latent_handle_id <- function(prefix, payload) {
  paste0(prefix, "-", digest::digest(payload, algo = "xxhash64"))
}

# --- Public registry lifecycle API ---

#' Clear the fmrilatent handle registry
#'
#' Removes cached materialized matrices from the internal registry.
#' Use this to free memory or force re-materialization of handles.
#'
#' @param type Which registry to clear: "basis", "loadings", or "all" (default).
#' @return Invisibly, the number of entries removed.
#' @export
#' @examples
#' # Clear all cached matrices
#' fmrilatent_registry_clear()
#'
#' # Clear only basis matrices
#' fmrilatent_registry_clear("basis")
fmrilatent_registry_clear <- function(type = c("all", "basis", "loadings")) {

  type <- match.arg(type)
  removed <- 0L


  if (type %in% c("all", "basis")) {
    env <- .latent_get_registry_env("basis")
    removed <- removed + length(ls(env))
    rm(list = ls(env), envir = env)
    .latent_registry_reset(env)
  }

  if (type %in% c("all", "loadings")) {
    env <- .latent_get_registry_env("loadings")
    removed <- removed + length(ls(env))
    rm(list = ls(env), envir = env)
    .latent_registry_reset(env)
  }

  invisible(removed)
}

#' List entries in the fmrilatent handle registry
#'
#' Returns the IDs of all cached matrices in the registry.
#'
#' @param type Which registry to list: "basis", "loadings", or "all" (default).
#' @return Character vector of registered IDs.
#' @export
#' @examples
#' # See what's cached
#' fmrilatent_registry_list()
fmrilatent_registry_list <- function(type = c("all", "basis", "loadings")) {
  type <- match.arg(type)

  ids <- character(0)


  if (type %in% c("all", "basis")) {
    env <- .latent_get_registry_env("basis")
    basis_ids <- ls(env)
    if (length(basis_ids) > 0) {
      names(basis_ids) <- rep("basis", length(basis_ids))
    }
    ids <- c(ids, basis_ids)
  }

  if (type %in% c("all", "loadings")) {
    env <- .latent_get_registry_env("loadings")
    loadings_ids <- ls(env)
    if (length(loadings_ids) > 0) {
      names(loadings_ids) <- rep("loadings", length(loadings_ids))
    }
    ids <- c(ids, loadings_ids)
  }

  ids
}

#' Get registry statistics
#'
#' Returns count and approximate memory usage of cached matrices.
#'
#' @param type Which registry: "basis", "loadings", or "all" (default).
#' @return A list with count and bytes for each registry type.
#' @export
#' @examples
#' fmrilatent_registry_stats()
fmrilatent_registry_stats <- function(type = c("all", "basis", "loadings")) {
  type <- match.arg(type)

  get_env_stats <- function(env_type) {
    env <- .latent_get_registry_env(env_type)
    ids <- ls(env)
    count <- length(ids)
    bytes <- 0
    if (count > 0) {
      bytes <- sum(vapply(ids, function(id) {
        obj <- get(id, envir = env, inherits = FALSE)
        as.numeric(object.size(obj))
      }, numeric(1)))
    }
    list(count = count, bytes = bytes)
  }

  result <- list()

  if (type %in% c("all", "basis")) {
    result$basis <- get_env_stats("basis")
  }

  if (type %in% c("all", "loadings")) {
    result$loadings <- get_env_stats("loadings")
  }

  if (type == "all") {
    result$total <- list(
      count = (result$basis$count %||% 0L) + (result$loadings$count %||% 0L),
      bytes = (result$basis$bytes %||% 0) + (result$loadings$bytes %||% 0)
    )
  }

  result
}

#' Enable or disable the fmrilatent handle registry
#'
#' The handle registry caches materialized matrices for \code{BasisHandle} and
#' \code{LoadingsHandle} objects. This can improve performance and reduce memory
#' duplication when multiple \code{LatentNeuroVec} objects share the same handle
#' IDs.
#'
#' Set \code{fmrilatent_registry_disable()} to turn off caching (useful for
#' deterministic benchmarking or to avoid retaining large matrices).
#'
#' @section Capacity:
#' The cache is bounded. Each registry (basis and loadings) holds at most
#' \code{getOption("fmrilatent.registry.max_entries", 256)} materialized
#' matrices; once full, the least-recently-used entry is evicted on the next
#' registration. Both registration and retrieval count as a use. Set the option
#' to \code{Inf} (or a non-positive value) to restore an unbounded cache, or call
#' \code{\link{fmrilatent_registry_clear}()} to drop everything immediately.
#'
#' @return Invisibly, \code{TRUE}.
#' @export
#' @examples
#' fmrilatent_registry_disable()
#' fmrilatent_registry_enable()
fmrilatent_registry_enable <- function() {
  options(fmrilatent.registry.enabled = TRUE)
  invisible(TRUE)
}

#' @rdname fmrilatent_registry_enable
#' @export
fmrilatent_registry_disable <- function() {
  options(fmrilatent.registry.enabled = FALSE)
  invisible(TRUE)
}

#' @rdname fmrilatent_registry_enable
#' @export
fmrilatent_registry_enabled <- function() {
  .latent_registry_enabled()
}

# Dimension helpers that avoid materializing handles
.latent_basis_dim <- function(b) {
  if (is(b, "Matrix")) {
    dim(b)
  } else if (is.matrix(b)) {
    dim(b)
  } else if (is.array(b) && length(dim(b)) == 2L) {
    dim(b)
  } else if (is(b, "BasisHandle")) {
    b@dim
  } else {
    .encoder_cli_abort(
      paste0("Unsupported basis slot type: ", paste(class(b), collapse = ", ")),
      class = "fmrilatent_error_type", call = rlang::caller_env())
  }
}

.latent_loadings_dim <- function(L) {
  if (is(L, "Matrix")) {
    dim(L)
  } else if (is.matrix(L)) {
    dim(L)
  } else if (is.array(L) && length(dim(L)) == 2L) {
    dim(L)
  } else if (is(L, "LoadingsHandle")) {
    L@dim
  } else {
    .encoder_cli_abort(
      paste0("Unsupported loadings slot type: ", paste(class(L), collapse = ", ")),
      class = "fmrilatent_error_type", call = rlang::caller_env())
  }
}
