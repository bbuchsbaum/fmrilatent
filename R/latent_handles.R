# Handle classes and registry for shared/implicit dictionaries

#' @keywords internal
setClass(
  "BasisHandle",
  slots = c(
    id    = "character",  # registry key
    dim   = "integer",    # c(n_time, k)
    kind  = "character",  # e.g., "dct", "lifted", "explicit"
    spec  = "list",       # parameters for reconstruction
    label = "character"   # optional human-readable label
  )
)

#' @keywords internal
setClass(
  "LoadingsHandle",
  slots = c(
    id    = "character",  # registry key
    dim   = "integer",    # c(n_vox, k)
    kind  = "character",  # e.g., "lifted", "explicit"
    spec  = "list",       # parameters for reconstruction
    label = "character"   # optional human-readable label
  )
)

# Unions so LatentNeuroVec can accept either Matrix or handle
setClassUnion("MatrixOrBasisHandle",    c("Matrix", "matrix", "BasisHandle"))
setClassUnion("MatrixOrLoadingsHandle", c("Matrix", "matrix", "LoadingsHandle"))

# Encapsulated cache for materialized dictionaries.
# NOTE: This is intentionally not a top-level mutable registry object. The cache
# is only accessible through helper functions and can be disabled via option.
.fmrilatent_cache_env <- local({
  basis_env <- NULL
  loadings_env <- NULL

  init <- function() {
    basis_env <<- new.env(parent = emptyenv())
    loadings_env <<- new.env(parent = emptyenv())
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

.latent_registry_enabled <- function() {
  isTRUE(getOption("fmrilatent.registry.enabled", TRUE))
}

.latent_get_registry_env <- function(type = c("basis", "loadings")) {
  type <- match.arg(type)
  .fmrilatent_cache_env(type)
}

.latent_register_matrix <- function(id, value, type = c("basis", "loadings"),
                                    overwrite = FALSE) {
  stopifnot(is.character(id), length(id) == 1L)
  type <- match.arg(type)
  if (!.latent_registry_enabled()) {
    return(invisible(FALSE))
  }
  env  <- .latent_get_registry_env(type)
  if (!overwrite && exists(id, env, inherits = FALSE)) {
    warning("Object with id '", id, "' already registered in ", type,
            " registry; set overwrite = TRUE to replace.", call. = FALSE)
    return(invisible(FALSE))
  }
  assign(id, value, envir = env)
  invisible(TRUE)
}

.latent_get_matrix <- function(id, type = c("basis", "loadings")) {
  stopifnot(is.character(id), length(id) == 1L)
  type <- match.arg(type)
  if (!.latent_registry_enabled()) {
    return(NULL)
  }
  env  <- .latent_get_registry_env(type)
  if (exists(id, env, inherits = FALSE)) {
    get(id, envir = env, inherits = FALSE)
  } else {
    NULL
  }
}

.latent_has_matrix <- function(id, type = c("basis", "loadings")) {
  !is.null(.latent_get_matrix(id, type = type))
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
  }

  if (type %in% c("all", "loadings")) {
    env <- .latent_get_registry_env("loadings")
    removed <- removed + length(ls(env))
    rm(list = ls(env), envir = env)
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

# Utility: safe NULL coalesce
`%||%` <- function(x, y) if (is.null(x)) y else x

# Helper to convert mask to array with informative error messages
#' @keywords internal
.mask_to_array <- function(mask, location = "unknown function") {
  result <- tryCatch(
    as.array(mask),
    error = function(e) {
      stop(
        sprintf("In %s: mask must be array-like or LogicalNeuroVol. ", location),
        "Underlying error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )
  if (is.null(result)) {
    stop(
      sprintf("In %s: mask must be array-like or LogicalNeuroVol ", location),
      "(conversion returned NULL).",
      call. = FALSE
    )
  }
  result
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
    stop("Unsupported basis slot type: ", paste(class(b), collapse = ", "))
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
    stop("Unsupported loadings slot type: ", paste(class(L), collapse = ", "))
  }
}
